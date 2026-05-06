#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import yaml
from scipy import sparse

from src.stages.storage import DEFAULT_LOW_RES_GROUP, raw_dir


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Prepare FFPE Mouse Brain Sagittal Visium + annotated snRNA reference as a Stage1 real sample."
    )
    p.add_argument("--project_root", default=".")
    p.add_argument("--source_dir", default="data/raw/low_resolution/FFPE Mouse Brain Sagittal")
    p.add_argument("--target_sample", default="ffpe_mouse_brain_sagittal_real")
    p.add_argument(
        "--annotation_mode",
        choices=("broad", "refined"),
        default="broad",
        help="Use broad deterministic classes or original refined annotation_1 labels.",
    )
    p.add_argument(
        "--include_cell_types",
        nargs="*",
        default=None,
        help="Optional allow-list after annotation mapping. Use annotation_1 labels with --annotation_mode=refined.",
    )
    p.add_argument("--max_cells_total", type=int, default=12000)
    p.add_argument("--min_cells_per_type", type=int, default=20)
    p.add_argument(
        "--max_cells_per_type",
        type=int,
        default=0,
        help="Optional per-type cap before total downsampling; useful for balanced simulations.",
    )
    p.add_argument("--export_genes", type=int, default=0, help="Keep top variable common genes; 0 keeps all.")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def _decode_arr(arr: np.ndarray) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else str(x) for x in arr.tolist()]


def _broad_cell_type(annotation: str) -> str | None:
    if annotation.startswith("Ext_"):
        return "Excitatory neuron"
    if annotation.startswith("Inh_"):
        return "Inhibitory neuron"
    if annotation.startswith("Astro_"):
        return "Astrocyte"
    if annotation.startswith("Oligo_"):
        return "Oligodendrocyte"
    if annotation.startswith("OPC_"):
        return "Oligodendrocyte precursor cell"
    if annotation == "Micro":
        return "Microglia"
    if annotation == "Endo":
        return "Endothelial cell"
    if annotation.startswith("Nb_"):
        return "Neuroblast"
    return None


def _make_unique_gene_names(feature_ids: list[str], gene_symbols: list[str]) -> tuple[list[str], dict[str, int]]:
    seen: dict[str, int] = {}
    duplicate_counts: dict[str, int] = {}
    out: list[str] = []
    for feature_id, symbol in zip(feature_ids, gene_symbols, strict=True):
        n = seen.get(symbol, 0) + 1
        seen[symbol] = n
        if n == 1:
            out.append(symbol)
        else:
            duplicate_counts[symbol] = n
            out.append(f"{symbol}-dup{n}-{feature_id}")
    return out, duplicate_counts


def _read_visium_h5(path: Path) -> tuple[list[str], list[str], list[str], sparse.csc_matrix]:
    with h5py.File(path, "r") as f:
        g = f["matrix"]
        shape = tuple(int(x) for x in g["shape"][()])
        data = g["data"][()]
        indices = g["indices"][()]
        indptr = g["indptr"][()]
        barcodes = _decode_arr(g["barcodes"][()])
        feature_ids = _decode_arr(g["features"]["id"][()])
        gene_symbols = _decode_arr(g["features"]["name"][()])
    mat = sparse.csc_matrix((data, indices, indptr), shape=shape)
    return feature_ids, gene_symbols, barcodes, mat


def _sample_cells_by_type(
    rng: np.random.Generator,
    cell_types: np.ndarray,
    max_cells_total: int,
    min_cells_per_type: int,
) -> np.ndarray:
    type_to_idx = {t: np.where(cell_types == t)[0] for t in sorted(pd.unique(cell_types).tolist())}
    total = sum(len(v) for v in type_to_idx.values())
    if max_cells_total <= 0 or max_cells_total >= total:
        return np.arange(total, dtype=np.int32)

    base = {t: min(min_cells_per_type, len(idx)) for t, idx in type_to_idx.items()}
    remain = max_cells_total - sum(base.values())
    if remain < 0:
        raise ValueError("max_cells_total is smaller than required per-type minimum allocation.")
    caps = {t: len(idx) - base[t] for t, idx in type_to_idx.items()}
    cap_sum = sum(max(0, v) for v in caps.values())
    extra = {t: 0 for t in type_to_idx}
    if remain > 0 and cap_sum > 0:
        raw = {t: remain * max(0, caps[t]) / cap_sum for t in type_to_idx}
        extra = {t: int(np.floor(raw[t])) for t in type_to_idx}
        left = remain - sum(extra.values())
        order = sorted(type_to_idx, key=lambda t: raw[t] - extra[t], reverse=True)
        i = 0
        while left > 0 and order:
            t = order[i % len(order)]
            if extra[t] < caps[t]:
                extra[t] += 1
                left -= 1
            i += 1

    keep: list[np.ndarray] = []
    for t, idx in type_to_idx.items():
        k = base[t] + extra[t]
        keep.append(idx if k >= len(idx) else rng.choice(idx, size=k, replace=False))
    merged = np.concatenate(keep)
    merged.sort()
    return merged.astype(np.int32, copy=False)


def _write_expr_tsv(path: Path, genes: list[str], sample_ids: list[str], mat_g_by_s: sparse.spmatrix) -> None:
    mat_csr = mat_g_by_s.tocsr()
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["Gene", *sample_ids])
        for i, gene in enumerate(genes):
            row = np.asarray(mat_csr.getrow(i).toarray()).ravel()
            writer.writerow([gene, *row.tolist()])
            if (i + 1) % 1000 == 0 or i + 1 == len(genes):
                print(f"[WRITE] {path.name}: {i + 1}/{len(genes)}")


def _write_dataset_yaml(project_root: Path, sample: str) -> None:
    cfg = {
        "storage": {
            "group": DEFAULT_LOW_RES_GROUP,
        },
        "paths": {
            "sc_expr": "brca_scRNA_GEP.txt",
            "sc_meta": "brca_scRNA_celllabels.txt",
            "st_expr": "brca_STdata_GEP.txt",
            "st_meta": "brca_STdata_coordinates.txt",
            "svg_marker_whitelist": None,
        },
        "qc": {
            "sc_min_genes": 100,
            "sc_max_genes": 10000,
            "sc_max_mt": 15,
            "st_min_genes": 50,
            "st_max_genes": "Inf",
            "st_max_mt": 30,
            "hvg_nfeatures": 2000,
            "mt_pattern": "^(mt-|MT-)",
        },
        "gene_filter": {"min_cells_sc": 0, "min_cells_st": 0},
        "stage3": {
            "strong_th": 0.7,
            "weak_th": 0.4,
            "st_cluster_k": 30,
            "unknown_floor": 0.3,
            "min_cells_rare_type": 20,
            "eps": 1.0e-8,
            "plugin_genes_path": (
                f"data/processed/{DEFAULT_LOW_RES_GROUP}/{sample}/stage1_preprocess/hvg_genes.txt"
            ),
            "gene_weights_path": None,
            "auto_missing_detection": {
                "enable": True,
                "method": "adaptive_low_support",
                "min_cells": 50,
                "robust_z_th": -2.0,
                "soft_z_th": -0.9,
                "require_masked_for_soft": False,
                "require_masked_for_hard": True,
                "require_confirmation": True,
                "confirmation_max_support_score": None,
                "confirmation_use_masked_missing": True,
                "confirmation_use_marker_identity": True,
                "confirmation_marker_identity_z_th": -1.5,
                "confirmation_marker_support_score_th": 0.5,
                "max_fraction_types": 0.5,
                "max_types": 4,
                "action": "mark_unknown",
            },
            "masked_missing_detection": {
                "enable": True,
                "apply_to_auto_missing": True,
                "neighbor_cosine_th": 0.9,
                "neighbor_cell_ratio_min": 1.0,
                "marker_top_n": 20,
                "min_identity_markers": 3,
                "min_marker_specificity": 1.2,
                "min_marker_type_mean": 0.001,
                "min_marker_st_detect_frac": 0.005,
                "st_presence_quantile": 0.9,
                "identity_z_th": -0.8,
                "pressure_z_th": 1.0,
                "min_support_score_for_apply": 0.8,
                "max_types": 4,
            },
            "marker_identity_diagnostics": {
                "enable": True,
                "marker_top_n": 80,
                "min_identity_markers": 5,
                "min_all_specificity": 1.3,
                "min_neighbor_specificity": 1.1,
                "min_marker_type_mean": 0.001,
                "min_marker_st_detect_frac": 0.005,
                "st_presence_quantile": 0.9,
                "depleted_z_th": -0.9,
            },
        },
    }
    out = project_root / "configs" / "datasets" / f"{sample}.yaml"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True), encoding="utf-8")
    print(f"[CFG] wrote: {out}")


def main() -> int:
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    project_root = Path(args.project_root).resolve()
    source_dir = (project_root / args.source_dir).resolve()
    target_dir = raw_dir(project_root, args.target_sample, {"storage": {"group": DEFAULT_LOW_RES_GROUP}})
    if target_dir.exists():
        if not args.overwrite:
            raise FileExistsError(f"Target exists: {target_dir} (use --overwrite)")
        shutil.rmtree(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)

    sc_h5ad = source_dir / "all_cells_20200625.h5ad"
    ann_csv = source_dir / "snRNA_annotation_astro_subtypes_refined59_20200823.csv"
    st_h5 = source_dir / "CytAssist_FFPE_Sagittal_Mouse_Brain_filtered_feature_bc_matrix.h5"
    st_pos_csv = source_dir / "spatial" / "tissue_positions.csv"
    for path in (sc_h5ad, ann_csv, st_h5, st_pos_csv):
        if not path.exists():
            raise FileNotFoundError(path)

    print("[STEP] load annotated snRNA reference")
    adata = ad.read_h5ad(sc_h5ad)
    ann = pd.read_csv(ann_csv, index_col=0)
    ann.index = ann.index.astype(str)
    obs = adata.obs.copy()
    obs["cell_id"] = obs.index.astype(str)
    obs = obs.join(ann[["annotation_1"]], how="left")
    obs["annotation_1"] = obs["annotation_1"].astype(str)
    if args.annotation_mode == "broad":
        obs["cell_type"] = obs["annotation_1"].map(_broad_cell_type)
    else:
        obs["cell_type"] = obs["annotation_1"]
    obs = obs.loc[obs["cell_type"].notna()].copy()
    if args.include_cell_types:
        include = set(args.include_cell_types)
        missing = sorted(include - set(obs["cell_type"].unique().tolist()))
        if missing:
            raise ValueError(f"include_cell_types not found after mapping: {missing}")
        obs = obs.loc[obs["cell_type"].isin(include)].copy()
    if args.max_cells_per_type and args.max_cells_per_type > 0:
        capped_parts = []
        for _, part in obs.groupby("cell_type", sort=True):
            if len(part) > args.max_cells_per_type:
                capped_parts.append(part.sample(n=args.max_cells_per_type, random_state=args.seed))
            else:
                capped_parts.append(part)
        obs = pd.concat(capped_parts, axis=0).sort_index()
    if obs.empty:
        raise ValueError("No annotated cells after cell-type mapping/filtering.")

    keep_idx_in_obs = _sample_cells_by_type(
        rng,
        obs["cell_type"].to_numpy(dtype=object),
        args.max_cells_total,
        args.min_cells_per_type,
    )
    obs = obs.iloc[keep_idx_in_obs].copy()
    keep_cell_ids = obs["cell_id"].to_numpy(dtype=object)
    keep_cell_types = obs["cell_type"].to_numpy(dtype=object)
    sc_sel = adata[keep_cell_ids, :].copy()
    sc_feature_ids = sc_sel.var_names.astype(str).tolist()
    sc_symbols = sc_sel.var["SYMBOL"].astype(str).tolist()
    sc_genes, sc_dups = _make_unique_gene_names(sc_feature_ids, sc_symbols)
    if sc_dups:
        print(f"[INFO] renamed duplicate snRNA gene symbols: {sc_dups}")
    print(f"[SC] selected cells={len(keep_cell_ids)} types={len(pd.unique(keep_cell_types))}")

    print("[STEP] load Visium h5 + coordinates")
    st_feature_ids, st_symbols, st_spots, st_mat = _read_visium_h5(st_h5)
    st_genes, st_dups = _make_unique_gene_names(st_feature_ids, st_symbols)
    if st_dups:
        print(f"[INFO] renamed duplicate ST gene symbols: {st_dups}")
    pos = pd.read_csv(st_pos_csv)
    required_cols = {"barcode", "pxl_row_in_fullres", "pxl_col_in_fullres"}
    if not required_cols.issubset(pos.columns):
        raise ValueError(f"Invalid tissue positions columns: {pos.columns.tolist()}")
    coord_df = pos.set_index("barcode").reindex(st_spots)
    if coord_df["pxl_row_in_fullres"].isna().any():
        raise ValueError("Some ST barcodes are missing from tissue positions.")

    print("[STEP] align genes by Ensembl feature id")
    st_feature_to_idx = {feature_id: i for i, feature_id in enumerate(st_feature_ids)}
    sc_feature_to_idx = {feature_id: i for i, feature_id in enumerate(sc_feature_ids)}
    common_feature_ids = [feature_id for feature_id in sc_feature_ids if feature_id in st_feature_to_idx]
    common_genes = [sc_genes[sc_feature_to_idx[feature_id]] for feature_id in common_feature_ids]
    if len(common_genes) != len(set(common_genes)):
        raise ValueError("Gene names must be unique after duplicate-symbol renaming.")
    if len(common_genes) < 500:
        raise ValueError(f"Too few common genes: {len(common_genes)}")
    sc_idx = np.array([sc_feature_to_idx[feature_id] for feature_id in common_feature_ids], dtype=np.int32)
    st_idx = np.array([st_feature_to_idx[feature_id] for feature_id in common_feature_ids], dtype=np.int32)

    sc_x = sc_sel.X
    if not sparse.issparse(sc_x):
        sc_x = sparse.csr_matrix(sc_x)
    sc_common = sc_x[:, sc_idx].tocsr()
    st_common = st_mat[st_idx, :].tocsc()
    if args.export_genes and args.export_genes > 0 and args.export_genes < len(common_genes):
        mean = np.asarray(sc_common.mean(axis=0)).ravel()
        mean_sq = np.asarray(sc_common.multiply(sc_common).mean(axis=0)).ravel()
        var = mean_sq - mean * mean
        top = np.argsort(var)[-int(args.export_genes) :]
        top.sort()
        common_genes = [common_genes[i] for i in top.tolist()]
        sc_common = sc_common[:, top].tocsr()
        st_common = st_common[top, :].tocsc()
    print(f"[INFO] common_genes={len(common_genes)} cells={sc_common.shape[0]} spots={st_common.shape[1]}")

    print("[STEP] write Stage1 inputs")
    _write_expr_tsv(target_dir / "brca_scRNA_GEP.txt", common_genes, keep_cell_ids.tolist(), sc_common.T)
    pd.DataFrame({"cell_id": keep_cell_ids, "cell_type": keep_cell_types}).to_csv(
        target_dir / "brca_scRNA_celllabels.txt", sep="\t", index=False
    )
    _write_expr_tsv(target_dir / "brca_STdata_GEP.txt", common_genes, st_spots, st_common)
    pd.DataFrame(
        {
            "spot_id": st_spots,
            "row": coord_df["pxl_row_in_fullres"].to_numpy(dtype=int),
            "col": coord_df["pxl_col_in_fullres"].to_numpy(dtype=int),
        }
    ).to_csv(target_dir / "brca_STdata_coordinates.txt", sep="\t", index=False)

    audit = (
        obs.groupby(["annotation_1", "cell_type"], as_index=False)
        .size()
        .rename(columns={"size": "n_cells"})
        .sort_values(["cell_type", "annotation_1"])
    )
    audit.to_csv(target_dir / "cell_type_annotation_audit.csv", index=False)
    info = {
        "sample": args.target_sample,
        "source_dir": str(source_dir),
        "n_cells": int(sc_common.shape[0]),
        "n_spots": int(st_common.shape[1]),
        "n_genes": int(len(common_genes)),
        "annotation_mode": args.annotation_mode,
        "include_cell_types": args.include_cell_types,
        "annotation_source": "cell2location mouse brain snRNA annotation_1; broad classes are derived deterministically from annotation prefixes.",
        "dropped_annotation_prefixes": ["LowQ", "Unk"],
        "cell_type_counts": pd.Series(keep_cell_types).value_counts().to_dict(),
        "params": vars(args),
    }
    (target_dir / "real_input_info.json").write_text(json.dumps(info, indent=2, ensure_ascii=False), encoding="utf-8")
    _write_dataset_yaml(project_root, args.target_sample)
    print(f"[DONE] prepared: {target_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
