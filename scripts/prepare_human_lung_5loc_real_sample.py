#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import json
import shutil
import tempfile
import zipfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import yaml
from scipy import sparse


LUNG_BROAD_V1 = {
    "AT1": "Alveolar epithelial",
    "AT2": "Alveolar epithelial",
    "Ciliated": "Airway epithelial",
    "Epi_Basal": "Airway epithelial",
    "Epi_secretory": "Airway epithelial",
    "Submucosal_Glands": "Airway epithelial",
    "Endothelia_lymphatic": "Endothelial",
    "Endothelia_vascular": "Endothelial",
    "Fibroblast": "Mesenchymal",
    "Muscle": "Mesenchymal",
    "Chondrocyte": "Chondrocyte",
    "B_cell": "B cell",
    "CD4": "T/NK/ILC",
    "CD8": "T/NK/ILC",
    "ILC": "T/NK/ILC",
    "NK": "T/NK/ILC",
    "TNK": "T/NK/ILC",
    "DC": "Myeloid",
    "Macrophage_alveolar": "Myeloid",
    "Macrophage_other": "Myeloid",
    "Monocyte": "Myeloid",
    "Mast_cell": "Mast cell",
    "Erythrocyte": "Erythrocyte",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Prepare the Human Lung 5 Locations sc/snRNA + one Visium sample as standard raw TSV inputs."
    )
    p.add_argument("--project_root", default=".")
    p.add_argument("--source_dir", default="data/raw/Human Lung 5 Locations")
    p.add_argument("--target_sample", default="human_lung_5loc_broad_real")
    p.add_argument("--visium_sample", default="WSA_LngSP10193345")
    p.add_argument("--cell_type_col", default="Celltypes_master_higher_immune")
    p.add_argument(
        "--type_granularity",
        choices=("broad", "fine"),
        default="broad",
        help="Use curated broad lung groups or the original fine labels from --cell_type_col.",
    )
    p.add_argument(
        "--include_cell_types",
        nargs="*",
        default=None,
        help="Optional allow-list after type mapping. Use fine labels when --type_granularity=fine.",
    )
    p.add_argument("--max_cells_total", type=int, default=15000)
    p.add_argument("--min_cells_per_type", type=int, default=250)
    p.add_argument(
        "--max_cells_per_type",
        type=int,
        default=0,
        help="Optional per-type cap before total downsampling; useful for balanced simulations.",
    )
    p.add_argument("--export_genes", type=int, default=3000)
    p.add_argument("--confirmation_marker_identity_z_th", type=float, default=-1.2)
    p.add_argument("--confirmation_marker_max_support_score", type=float, default=0.75)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def _decode_arr(arr: np.ndarray) -> list[str]:
    out: list[str] = []
    for x in arr.tolist():
        out.append(x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else str(x))
    return out


def _read_visium_h5(path: Path) -> tuple[list[str], list[str], sparse.csc_matrix]:
    with h5py.File(path, "r") as f:
        g = f["matrix"]
        shape = tuple(int(x) for x in g["shape"][()])
        data = g["data"][()]
        indices = g["indices"][()]
        indptr = g["indptr"][()]
        barcodes = _decode_arr(g["barcodes"][()])
        genes = _decode_arr(g["features"]["name"][()])
    mat = sparse.csc_matrix((data, indices, indptr), shape=shape)
    return genes, barcodes, mat


def _write_expr_tsv(path: Path, genes: list[str], sample_ids: list[str], mat_g_by_s: np.ndarray) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Gene", *sample_ids])
        for i, gene in enumerate(genes):
            writer.writerow([gene, *mat_g_by_s[i, :].tolist()])
            if (i + 1) % 500 == 0 or i + 1 == len(genes):
                print(f"[WRITE] {path.name}: {i + 1}/{len(genes)}")


def _sample_cells_by_type(
    rng: np.random.Generator,
    cell_ids: np.ndarray,
    cell_types: np.ndarray,
    max_cells_total: int,
    min_cells_per_type: int,
) -> np.ndarray:
    type_to_idx = {t: np.where(cell_types == t)[0] for t in sorted(pd.unique(cell_types).tolist())}
    total = sum(len(v) for v in type_to_idx.values())
    if max_cells_total <= 0 or max_cells_total >= total:
        return np.arange(total, dtype=np.int32)

    n_types = len(type_to_idx)
    if max_cells_total < n_types:
        raise ValueError(f"max_cells_total({max_cells_total}) < n_types({n_types})")

    counts = {t: len(idx) for t, idx in type_to_idx.items()}
    base = {t: min(min_cells_per_type, counts[t]) for t in type_to_idx}
    remain = max_cells_total - sum(base.values())
    if remain < 0:
        target = {t: max(1, int(round(counts[t] / total * max_cells_total))) for t in type_to_idx}
        while sum(target.values()) > max_cells_total:
            key = max(target, key=target.get)
            if target[key] > 1:
                target[key] -= 1
            else:
                break
    else:
        caps = {t: counts[t] - base[t] for t in type_to_idx}
        cap_sum = sum(max(0, v) for v in caps.values())
        extra = {t: 0 for t in type_to_idx}
        if remain > 0 and cap_sum > 0:
            alloc = {t: max(0, caps[t]) / cap_sum * remain for t in type_to_idx}
            extra = {t: int(np.floor(v)) for t, v in alloc.items()}
            left = remain - sum(extra.values())
            order = sorted(type_to_idx, key=lambda t: alloc[t] - extra[t], reverse=True)
            i = 0
            while left > 0 and order:
                t = order[i % len(order)]
                if extra[t] < caps[t]:
                    extra[t] += 1
                    left -= 1
                i += 1
        target = {t: base[t] + extra[t] for t in type_to_idx}

    keep: list[np.ndarray] = []
    for t, idx in type_to_idx.items():
        k = min(target[t], len(idx))
        keep.append(idx if k >= len(idx) else rng.choice(idx, size=k, replace=False))
    merged = np.concatenate(keep)
    merged.sort()
    return merged.astype(np.int32, copy=False)


def _log_norm_dense(x: np.ndarray, scale: float = 10000.0) -> np.ndarray:
    x = x.astype(np.float32, copy=True)
    lib = x.sum(axis=1, keepdims=True)
    lib[lib <= 0] = 1.0
    x = x / lib * scale
    np.log1p(x, out=x)
    return x


def _write_dataset_yaml(
    project_root: Path,
    sample: str,
    confirmation_marker_identity_z_th: float,
    confirmation_marker_max_support_score: float,
) -> None:
    cfg = {
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
            "mt_pattern": "^(MT-|mt-)",
        },
        "gene_filter": {"min_cells_sc": 0, "min_cells_st": 0},
        "stage3": {
            "strong_th": 0.7,
            "weak_th": 0.4,
            "st_cluster_k": 30,
            "unknown_floor": 0.3,
            "min_cells_rare_type": 20,
            "eps": 1.0e-8,
            "plugin_genes_path": f"data/processed/low_resolution_experiments/{sample}/stage1_preprocess/hvg_genes.txt",
            "gene_weights_path": None,
            "auto_missing_detection": {
                "enable": True,
                "method": "adaptive_low_support",
                "min_cells": 50,
                "robust_z_th": -2.0,
                "soft_z_th": -0.9,
                "require_masked_for_soft": True,
                "require_masked_for_hard": True,
                "require_confirmation": True,
                "confirmation_use_masked_missing": True,
                "confirmation_use_marker_identity": True,
                "confirmation_marker_identity_z_th": float(confirmation_marker_identity_z_th),
                "confirmation_marker_support_score_th": 0.65,
                "confirmation_marker_max_support_score": float(confirmation_marker_max_support_score),
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
        },
        "storage": {"group": "low_resolution_experiments"},
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
    target_dir = project_root / "data" / "raw" / "low_resolution_experiments" / args.target_sample
    if target_dir.exists():
        if not args.overwrite:
            raise FileExistsError(f"Target exists: {target_dir} (use --overwrite)")
        shutil.rmtree(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)

    sc_h5ad = source_dir / "sc_reference" / "lung_5loc_sc_sn_raw_counts_cellxgene.h5ad"
    visium_zip = source_dir / "visium" / f"{args.visium_sample}_spaceranger_counts.zip"
    if not sc_h5ad.exists():
        raise FileNotFoundError(sc_h5ad)
    if not visium_zip.exists():
        raise FileNotFoundError(visium_zip)

    print("[STEP] load sc/snRNA h5ad metadata")
    adata = ad.read_h5ad(sc_h5ad, backed="r")
    if args.cell_type_col not in adata.obs.columns:
        raise ValueError(f"cell_type_col not found: {args.cell_type_col}")
    obs = adata.obs[[args.cell_type_col]].copy()
    obs["cell_id"] = obs.index.astype(str)
    obs["source_type"] = obs[args.cell_type_col].astype(str).str.strip()
    if args.type_granularity == "broad":
        obs["cell_type"] = obs["source_type"].map(LUNG_BROAD_V1)
    else:
        obs["cell_type"] = obs["source_type"]
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
        raise ValueError("No cells after type mapping/filtering.")
    print(f"[SC] {args.type_granularity} type source counts:")
    print(obs["cell_type"].value_counts().to_string())

    cell_ids_all = obs["cell_id"].to_numpy(dtype=object)
    cell_types_all = obs["cell_type"].to_numpy(dtype=object)
    keep_idx = _sample_cells_by_type(
        rng,
        cell_ids_all,
        cell_types_all,
        args.max_cells_total,
        args.min_cells_per_type,
    )
    keep_cell_ids = cell_ids_all[keep_idx]
    keep_cell_types = cell_types_all[keep_idx]
    sc_sel = adata[keep_cell_ids, :]
    sc_genes = sc_sel.var_names.astype(str).tolist()
    print(f"[SC] selected cells={len(keep_cell_ids)} types={len(pd.unique(keep_cell_types))}")
    print(pd.Series(keep_cell_types).value_counts().to_string())

    print("[STEP] extract Visium matrix + coordinates")
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        with zipfile.ZipFile(visium_zip) as zf:
            h5_member = f"{args.visium_sample}/filtered_feature_bc_matrix.h5"
            pos_member = f"{args.visium_sample}/spatial/tissue_positions_list.csv"
            h5_path = tmp / "filtered_feature_bc_matrix.h5"
            pos_path = tmp / "tissue_positions_list.csv"
            h5_path.write_bytes(zf.read(h5_member))
            pos_path.write_bytes(zf.read(pos_member))
        st_genes, st_spots, st_mat = _read_visium_h5(h5_path)
        pos = pd.read_csv(pos_path, header=None)

    if pos.shape[1] < 6:
        raise ValueError("Invalid tissue_positions_list.csv")
    pos.columns = ["spot_id", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    pos["spot_id"] = pos["spot_id"].astype(str)
    coords = pos.set_index("spot_id").reindex(st_spots)
    if coords["pxl_row"].isna().any():
        raise ValueError("Some ST barcodes are missing in tissue positions.")

    print("[STEP] align genes and select compact export panel")
    st_gene_to_idx = {g: i for i, g in enumerate(st_genes)}
    common_genes = [g for g in sc_genes if g in st_gene_to_idx]
    if len(common_genes) < args.export_genes:
        raise ValueError(f"Too few common genes: {len(common_genes)}")
    sc_gene_to_idx = {g: i for i, g in enumerate(sc_genes)}
    sc_idx = np.array([sc_gene_to_idx[g] for g in common_genes], dtype=np.int32)
    st_idx = np.array([st_gene_to_idx[g] for g in common_genes], dtype=np.int32)
    sc_x = sc_sel.X
    if sparse.issparse(sc_x):
        sc_common = sc_x[:, sc_idx].tocsr()
        mean = np.asarray(sc_common.mean(axis=0)).ravel()
        mean_sq = np.asarray(sc_common.multiply(sc_common).mean(axis=0)).ravel()
        var = mean_sq - mean * mean
    else:
        sc_common = np.asarray(sc_x[:, sc_idx], dtype=np.float32)
        var = sc_common.var(axis=0)
    top = np.argsort(var)[-int(args.export_genes) :]
    top.sort()
    export_genes = [common_genes[i] for i in top.tolist()]
    if sparse.issparse(sc_common):
        sc_counts = sc_common[:, top].T.toarray().astype(np.float32, copy=False)
    else:
        sc_counts = sc_common[:, top].T.astype(np.float32, copy=False)
    st_counts = st_mat[st_idx[top], :].toarray().astype(np.float32, copy=False)
    print(
        f"[INFO] common_genes={len(common_genes)} export_genes={len(export_genes)} "
        f"cells={sc_counts.shape[1]} spots={st_counts.shape[1]}"
    )

    print("[STEP] write standard raw files")
    _write_expr_tsv(target_dir / "brca_scRNA_GEP.txt", export_genes, keep_cell_ids.tolist(), sc_counts)
    pd.DataFrame({"cell_id": keep_cell_ids, "cell_type": keep_cell_types}).to_csv(
        target_dir / "brca_scRNA_celllabels.txt", sep="\t", index=False
    )
    _write_expr_tsv(target_dir / "brca_STdata_GEP.txt", export_genes, st_spots, st_counts)
    pd.DataFrame(
        {
            "spot_id": st_spots,
            "row": coords["pxl_row"].to_numpy(),
            "col": coords["pxl_col"].to_numpy(),
        }
    ).to_csv(target_dir / "brca_STdata_coordinates.txt", sep="\t", index=False)

    info = {
        "sample": args.target_sample,
        "source_dir": str(source_dir),
        "visium_sample": args.visium_sample,
        "cell_type_col": args.cell_type_col,
        "type_granularity": args.type_granularity,
        "include_cell_types": args.include_cell_types,
        "broad_mapping": LUNG_BROAD_V1,
        "n_cells": int(sc_counts.shape[1]),
        "n_spots": int(st_counts.shape[1]),
        "n_genes": int(len(export_genes)),
        "cell_type_counts": pd.Series(keep_cell_types).value_counts().to_dict(),
        "params": vars(args),
    }
    (target_dir / "real_input_info.json").write_text(json.dumps(info, indent=2, ensure_ascii=False), encoding="utf-8")
    _write_dataset_yaml(
        project_root,
        args.target_sample,
        args.confirmation_marker_identity_z_th,
        args.confirmation_marker_max_support_score,
    )
    print(f"[DONE] human lung sample prepared: {target_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
