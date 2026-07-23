#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import yaml
from scipy import sparse


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Convert one low_resolution h5ad + Visium dataset into Stage1 TSV inputs.")
    p.add_argument("--project_root", default=".")
    p.add_argument("--source_dir", required=True)
    p.add_argument("--target_sample", required=True)
    p.add_argument("--sc_h5ad", required=True)
    p.add_argument("--st_h5", required=True)
    p.add_argument("--st_pos_csv", default="spatial/tissue_positions_list.csv")
    p.add_argument("--cell_type_col", default="Celltype")
    p.add_argument("--drop_unknown", action="store_true")
    p.add_argument("--min_source_cells_per_type", type=int, default=200)
    p.add_argument("--max_types", type=int, default=10)
    p.add_argument("--max_cells_total", type=int, default=12000)
    p.add_argument("--min_cells_per_type", type=int, default=50)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def _decode_arr(arr: np.ndarray) -> list[str]:
    out: list[str] = []
    for x in arr.tolist():
        if isinstance(x, (bytes, bytearray)):
            out.append(x.decode("utf-8"))
        else:
            out.append(str(x))
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
    type_to_idx = {
        t: np.where(cell_types == t)[0]
        for t in sorted(pd.unique(cell_types).tolist())
    }
    total = sum(len(v) for v in type_to_idx.values())
    if max_cells_total <= 0 or max_cells_total >= total:
        return np.arange(total, dtype=np.int32)

    n_types = len(type_to_idx)
    if max_cells_total < n_types:
        raise ValueError(f"max_cells_total({max_cells_total}) < n_types({n_types})")

    counts = {t: len(idx) for t, idx in type_to_idx.items()}
    base = {t: min(min_cells_per_type, counts[t]) for t in type_to_idx}
    base_sum = sum(base.values())
    if base_sum > max_cells_total:
        target = {t: max(1, int(round(counts[t] / total * max_cells_total))) for t in type_to_idx}
        while sum(target.values()) > max_cells_total:
            key = max(target, key=target.get)
            if target[key] > 1:
                target[key] -= 1
            else:
                break
    else:
        remain = max_cells_total - base_sum
        caps = {t: counts[t] - base[t] for t in type_to_idx}
        cap_sum = sum(max(0, x) for x in caps.values())
        extra = {t: 0 for t in type_to_idx}
        if remain > 0 and cap_sum > 0:
            alloc = {t: (max(0, caps[t]) / cap_sum) * remain for t in type_to_idx}
            extra = {t: int(np.floor(v)) for t, v in alloc.items()}
            left = remain - sum(extra.values())
            order = sorted(type_to_idx.keys(), key=lambda t: alloc[t] - extra[t], reverse=True)
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
        k = target[t]
        keep.append(idx if k >= len(idx) else rng.choice(idx, size=k, replace=False))
    merged = np.concatenate(keep)
    merged.sort()
    return merged.astype(np.int32, copy=False)


def _write_dataset_yaml(project_root: Path, sample: str) -> None:
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
            "plugin_genes_path": f"data/processed/{sample}/stage1_preprocess/hvg_genes.txt",
            "gene_weights_path": None,
            "auto_missing_detection": {
                "enable": True,
                "method": "adaptive_low_support",
                "min_cells": 50,
                "robust_z_th": -2.0,
                "soft_z_th": -0.9,
                "require_masked_for_soft": True,
                "require_masked_for_hard": True,
                "max_fraction_types": 0.3,
                "max_types": 2,
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
                "max_types": 2,
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
    target_dir = project_root / "data" / "raw" / args.target_sample
    if target_dir.exists():
        if not args.overwrite:
            raise FileExistsError(f"Target exists: {target_dir} (use --overwrite)")
        shutil.rmtree(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)

    sc_h5ad = source_dir / args.sc_h5ad
    st_h5 = source_dir / args.st_h5
    st_pos_csv = source_dir / args.st_pos_csv
    for path in (sc_h5ad, st_h5, st_pos_csv):
        if not path.exists():
            raise FileNotFoundError(path)

    print("[STEP] load scRNA h5ad")
    adata = ad.read_h5ad(sc_h5ad)
    if args.cell_type_col not in adata.obs.columns:
        raise ValueError(f"cell_type_col not found: {args.cell_type_col}")
    obs = adata.obs.copy()
    obs["cell_id"] = obs.index.astype(str)
    obs["cell_type"] = obs[args.cell_type_col].astype(str).str.strip()
    if args.drop_unknown:
        bad = {"", "nan", "na", "none", "unknown"}
        obs = obs.loc[~obs["cell_type"].str.lower().isin(bad)].copy()
    if args.min_source_cells_per_type > 0:
        counts = obs["cell_type"].value_counts()
        keep_types = counts[counts >= args.min_source_cells_per_type].index.tolist()
        obs = obs.loc[obs["cell_type"].isin(keep_types)].copy()
        print(f"[SC] min_source_cells_per_type={args.min_source_cells_per_type}; types={len(keep_types)} cells={len(obs)}")
    if args.max_types > 0:
        counts = obs["cell_type"].value_counts()
        keep_types = counts.head(args.max_types).index.tolist()
        obs = obs.loc[obs["cell_type"].isin(keep_types)].copy()
        print(f"[SC] max_types={args.max_types}; types={len(keep_types)} cells={len(obs)}")
    if obs.empty:
        raise ValueError("No SC cells after filtering.")

    cell_ids_all = obs["cell_id"].to_numpy(dtype=object)
    cell_types_all = obs["cell_type"].to_numpy(dtype=object)
    keep_idx = _sample_cells_by_type(rng, cell_ids_all, cell_types_all, args.max_cells_total, args.min_cells_per_type)
    keep_cell_ids = cell_ids_all[keep_idx]
    keep_cell_types = cell_types_all[keep_idx]
    sc_sel = adata[keep_cell_ids, :].copy()
    sc_genes = sc_sel.var_names.astype(str).tolist()
    print(f"[SC] selected cells={len(keep_cell_ids)} types={len(pd.unique(keep_cell_types))}")

    print("[STEP] load Visium h5 + coordinates")
    st_genes, st_spots, st_mat = _read_visium_h5(st_h5)
    pos = pd.read_csv(st_pos_csv, header=None)
    if pos.shape[1] < 6:
        raise ValueError(f"Invalid tissue positions file: {st_pos_csv}")
    pos.columns = ["spot_id", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    pos["spot_id"] = pos["spot_id"].astype(str)
    coord_df = pos.set_index("spot_id").reindex(st_spots)
    if coord_df["pxl_row"].isna().any():
        raise ValueError("Some ST barcodes are missing from tissue positions.")

    print("[STEP] align genes")
    st_gene_to_idx = {g: i for i, g in enumerate(st_genes)}
    common_genes = [g for g in sc_genes if g in st_gene_to_idx]
    if len(common_genes) < 500:
        raise ValueError(f"Too few common genes: {len(common_genes)}")
    sc_idx = np.array([sc_genes.index(g) for g in common_genes], dtype=np.int32)
    st_idx = np.array([st_gene_to_idx[g] for g in common_genes], dtype=np.int32)

    sc_x = sc_sel.X
    if not sparse.issparse(sc_x):
        sc_x = sparse.csr_matrix(sc_x)
    sc_common = sc_x[:, sc_idx].tocsr()
    st_common = st_mat[st_idx, :].tocsc()
    print(f"[INFO] common_genes={len(common_genes)} cells={sc_common.shape[0]} spots={st_common.shape[1]}")

    sc_expr_path = target_dir / "brca_scRNA_GEP.txt"
    sc_meta_path = target_dir / "brca_scRNA_celllabels.txt"
    st_expr_path = target_dir / "brca_STdata_GEP.txt"
    st_meta_path = target_dir / "brca_STdata_coordinates.txt"

    _write_expr_tsv(sc_expr_path, common_genes, keep_cell_ids.tolist(), sc_common.T.toarray().astype(np.float32))
    pd.DataFrame({"cell_id": keep_cell_ids, "cell_type": keep_cell_types}).to_csv(sc_meta_path, sep="\t", index=False)
    _write_expr_tsv(st_expr_path, common_genes, st_spots, st_common.toarray().astype(np.float32))
    pd.DataFrame({"spot_id": st_spots, "row": coord_df["pxl_row"].to_numpy(), "col": coord_df["pxl_col"].to_numpy()}).to_csv(
        st_meta_path, sep="\t", index=False
    )

    info = {
        "sample": args.target_sample,
        "source_dir": str(source_dir),
        "n_cells": int(sc_common.shape[0]),
        "n_spots": int(st_common.shape[1]),
        "n_genes": int(len(common_genes)),
        "cell_type_col": args.cell_type_col,
        "cell_type_counts": pd.Series(keep_cell_types).value_counts().to_dict(),
        "params": vars(args),
    }
    (target_dir / "real_input_info.json").write_text(json.dumps(info, indent=2, ensure_ascii=False), encoding="utf-8")
    _write_dataset_yaml(project_root, args.target_sample)
    print(f"[DONE] real sample prepared: {target_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
