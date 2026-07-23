#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse

from src.stages.storage import DEFAULT_LOW_RES_GROUP, raw_dir


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Prepare Spatial Skin Atlas h5ad files as a Stage1 real sample.")
    p.add_argument("--project_root", default=".")
    p.add_argument("--source_dir", default="data/raw/low_resolution/Human Skin Atlas")
    p.add_argument("--target_sample", default="human_skin_atlas_real")
    p.add_argument("--cell_type_col", default="04_celltypes")
    p.add_argument("--max_types", type=int, default=12)
    p.add_argument("--min_source_cells_per_type", type=int, default=500)
    p.add_argument("--max_cells_total", type=int, default=12000)
    p.add_argument("--min_cells_per_type", type=int, default=50)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


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


def _raw_or_x(adata: ad.AnnData):
    if adata.raw is not None:
        return adata.raw.X, adata.raw.var_names.astype(str).tolist(), "raw.X"
    return adata.X, adata.var_names.astype(str).tolist(), "X"


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

    sc_h5ad = source_dir / "bcc_and_normal-CG_portal_fat.h5ad"
    st_h5ad = source_dir / "WSSKNKCLsp10767968.h5ad"
    for path in (sc_h5ad, st_h5ad):
        if not path.exists():
            raise FileNotFoundError(path)

    print("[STEP] load scRNA reference metadata")
    sc = ad.read_h5ad(sc_h5ad, backed="r")
    if args.cell_type_col not in sc.obs.columns:
        raise ValueError(f"cell_type_col not found: {args.cell_type_col}")
    obs = sc.obs.copy()
    obs["cell_id"] = obs.index.astype(str)
    obs["cell_type"] = obs[args.cell_type_col].astype(str).str.strip()
    counts_all = obs["cell_type"].value_counts()
    keep_types = counts_all[counts_all >= args.min_source_cells_per_type].head(args.max_types).index.tolist()
    obs = obs.loc[obs["cell_type"].isin(keep_types)].copy()
    keep_rows = _sample_cells_by_type(
        rng,
        obs["cell_type"].to_numpy(dtype=object),
        args.max_cells_total,
        args.min_cells_per_type,
    )
    obs = obs.iloc[keep_rows].copy()
    keep_cell_ids = obs["cell_id"].to_numpy(dtype=object)
    keep_cell_types = obs["cell_type"].to_numpy(dtype=object)
    print(f"[SC] selected cells={len(keep_cell_ids)} types={len(pd.unique(keep_cell_types))}")

    print("[STEP] load selected scRNA expression")
    sc_sel = sc[keep_cell_ids, :].to_memory()
    sc.file.close()
    st = ad.read_h5ad(st_h5ad)
    sc_x, sc_genes, sc_source = _raw_or_x(sc_sel)
    st_x, st_genes, st_source = _raw_or_x(st)
    print(f"[INFO] expression source: sc={sc_source} st={st_source}")
    if len(sc_genes) != len(set(sc_genes)) or len(st_genes) != len(set(st_genes)):
        raise ValueError("Gene names must be unique in both h5ad files.")

    print("[STEP] align genes")
    st_gene_to_idx = {g: i for i, g in enumerate(st_genes)}
    common_genes = [g for g in sc_genes if g in st_gene_to_idx]
    if len(common_genes) < 500:
        raise ValueError(f"Too few common genes: {len(common_genes)}")
    sc_gene_to_idx = {g: i for i, g in enumerate(sc_genes)}
    sc_idx = np.array([sc_gene_to_idx[g] for g in common_genes], dtype=np.int32)
    st_idx = np.array([st_gene_to_idx[g] for g in common_genes], dtype=np.int32)
    if not sparse.issparse(sc_x):
        sc_x = sparse.csr_matrix(sc_x)
    if not sparse.issparse(st_x):
        st_x = sparse.csr_matrix(st_x)
    sc_common = sc_x[:, sc_idx].tocsr()
    st_common = st_x[:, st_idx].tocsr()
    print(f"[INFO] common_genes={len(common_genes)} cells={sc_common.shape[0]} spots={st_common.shape[0]}")

    if "X_spatial" not in st.obsm.keys():
        raise ValueError("ST h5ad is missing obsm['X_spatial'].")
    coords = np.asarray(st.obsm["X_spatial"])
    st_spots = st.obs_names.astype(str).tolist()

    print("[STEP] write Stage1 inputs")
    _write_expr_tsv(target_dir / "brca_scRNA_GEP.txt", common_genes, keep_cell_ids.tolist(), sc_common.T)
    pd.DataFrame({"cell_id": keep_cell_ids, "cell_type": keep_cell_types}).to_csv(
        target_dir / "brca_scRNA_celllabels.txt", sep="\t", index=False
    )
    _write_expr_tsv(target_dir / "brca_STdata_GEP.txt", common_genes, st_spots, st_common.T)
    pd.DataFrame({"spot_id": st_spots, "row": coords[:, 1], "col": coords[:, 0]}).to_csv(
        target_dir / "brca_STdata_coordinates.txt", sep="\t", index=False
    )

    audit = (
        obs.groupby("cell_type", as_index=False)
        .size()
        .rename(columns={"size": "n_cells"})
        .sort_values("n_cells", ascending=False)
    )
    audit.to_csv(target_dir / "cell_type_annotation_audit.csv", index=False)
    info = {
        "sample": args.target_sample,
        "source_dir": str(source_dir),
        "st_h5ad": st_h5ad.name,
        "sc_h5ad": sc_h5ad.name,
        "cell_type_col": args.cell_type_col,
        "n_cells": int(sc_common.shape[0]),
        "n_spots": int(st_common.shape[0]),
        "n_genes": int(len(common_genes)),
        "cell_type_counts": pd.Series(keep_cell_types).value_counts().to_dict(),
        "dropped_cell_types": [t for t in counts_all.index.tolist() if t not in keep_types],
        "params": vars(args),
    }
    (target_dir / "real_input_info.json").write_text(json.dumps(info, indent=2, ensure_ascii=False), encoding="utf-8")
    _write_dataset_yaml(project_root, args.target_sample)
    print(f"[DONE] prepared: {target_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
