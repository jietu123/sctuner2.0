#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import yaml
from scipy import sparse

_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.stages.storage import DEFAULT_LOW_RES_GROUP, raw_dir


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Prepare full cell-type Human Intestine Cancer h5ad + Visium data as a standard raw sample."
    )
    p.add_argument("--project_root", default=".")
    p.add_argument("--source_dir", default="data/raw/low_resolution/Human Intestine Cancer")
    p.add_argument("--target_sample", default="human_intestine_cancer_full_real")
    p.add_argument("--cell_type_col", default="cell_type")
    p.add_argument("--max_cells_total", type=int, default=15000)
    p.add_argument("--min_cells_per_type", type=int, default=80)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def _decode_arr(arr: np.ndarray) -> list[str]:
    out: list[str] = []
    for x in arr.tolist():
        out.append(x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else str(x))
    return out


def _make_unique_gene_names(feature_ids: list[str], gene_names: list[str]) -> tuple[list[str], dict[str, int]]:
    seen: dict[str, int] = {}
    duplicate_counts: dict[str, int] = {}
    unique: list[str] = []
    for feature_id, gene_name in zip(feature_ids, gene_names, strict=True):
        n = seen.get(gene_name, 0) + 1
        seen[gene_name] = n
        if n == 1:
            unique.append(gene_name)
        else:
            duplicate_counts[gene_name] = n
            unique.append(f"{gene_name}__dup{n}_{feature_id}")
    return unique, duplicate_counts


def _read_10x_h5(path: Path) -> tuple[list[str], list[str], list[str], sparse.csc_matrix]:
    with h5py.File(path, "r") as f:
        g = f["matrix"]
        shape = tuple(int(x) for x in g["shape"][()])
        data = g["data"][()]
        indices = g["indices"][()]
        indptr = g["indptr"][()]
        barcodes = _decode_arr(g["barcodes"][()])
        feature_ids = _decode_arr(g["features"]["id"][()])
        genes_raw = _decode_arr(g["features"]["name"][()])
        genes, duplicate_counts = _make_unique_gene_names(feature_ids, genes_raw)
    if duplicate_counts:
        print(f"[INFO] renamed duplicate gene symbols in {path.name}: {duplicate_counts}")
    mat = sparse.csc_matrix((data, indices, indptr), shape=shape)
    return feature_ids, genes, barcodes, mat


def _sample_cells_by_type(
    rng: np.random.Generator,
    obs: pd.DataFrame,
    max_cells_total: int,
    min_cells_per_type: int,
) -> pd.DataFrame:
    groups = {t: idx.to_numpy() for t, idx in obs.groupby("cell_type", observed=True).groups.items()}
    total = len(obs)
    if max_cells_total <= 0 or max_cells_total >= total:
        return obs.copy()

    base: dict[str, int] = {}
    for t, idx in groups.items():
        base[t] = min(len(idx), min_cells_per_type)
    base_sum = sum(base.values())
    if base_sum > max_cells_total:
        raise ValueError("max_cells_total is too small to retain requested per-type minimum.")

    remain = max_cells_total - base_sum
    caps = {t: max(0, len(idx) - base[t]) for t, idx in groups.items()}
    cap_sum = sum(caps.values())
    extra = {t: 0 for t in groups}
    if remain > 0 and cap_sum > 0:
        raw = {t: remain * caps[t] / cap_sum for t in groups}
        extra = {t: int(np.floor(raw[t])) for t in groups}
        left = remain - sum(extra.values())
        order = sorted(groups, key=lambda t: raw[t] - extra[t], reverse=True)
        i = 0
        while left > 0 and order:
            t = order[i % len(order)]
            if extra[t] < caps[t]:
                extra[t] += 1
                left -= 1
            i += 1

    keep: list[np.ndarray] = []
    for t, idx in groups.items():
        k = min(len(idx), base[t] + extra[t])
        keep.append(idx if k >= len(idx) else rng.choice(idx, size=k, replace=False))
    keep_idx = np.concatenate(keep)
    keep_idx.sort()
    return obs.iloc[keep_idx].copy()


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
        "storage": {"group": DEFAULT_LOW_RES_GROUP},
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
            "sc_max_mt": 20,
            "st_min_genes": 50,
            "st_max_genes": "Inf",
            "st_max_mt": 40,
            "hvg_nfeatures": 2000,
            "mt_pattern": "^(MT-|mt-)",
        },
        "gene_filter": {"min_cells_sc": 0, "min_cells_st": 0},
        "stage3": {
            "strong_th": 0.7,
            "weak_th": 0.4,
            "st_cluster_k": 30,
            "unknown_floor": 0.3,
            "min_cells_rare_type": 10,
            "eps": 1.0e-8,
            "plugin_genes_path": (
                f"data/processed/{DEFAULT_LOW_RES_GROUP}/{sample}/stage1_preprocess/hvg_genes.txt"
            ),
            "gene_weights_path": None,
            "auto_missing_detection": {
                "enable": True,
                "method": "adaptive_low_support",
                "min_cells": 10,
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


def _read_positions(path: Path, spot_ids: list[str]) -> pd.DataFrame:
    positions = pd.read_csv(path, header=None)
    if positions.shape[1] < 6:
        raise ValueError(f"invalid tissue positions file: {path}")
    positions = positions.iloc[:, :6].copy()
    positions.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    positions["barcode"] = positions["barcode"].astype(str)
    coord_df = positions.set_index("barcode").reindex(spot_ids)
    if coord_df[["pxl_row", "pxl_col"]].isna().any().any():
        missing = int(coord_df[["pxl_row", "pxl_col"]].isna().any(axis=1).sum())
        raise ValueError(f"{missing} ST spots are missing tissue positions.")
    return pd.DataFrame(
        {
            "spot_id": spot_ids,
            "row": coord_df["pxl_row"].to_numpy(dtype=int),
            "col": coord_df["pxl_col"].to_numpy(dtype=int),
        }
    )


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

    sc_h5ad = source_dir / "sc_human_intestine_heca.h5ad"
    st_h5 = source_dir / "Visium_FFPE_Human_Intestinal_Cancer_filtered_feature_bc_matrix.h5"
    positions_csv = source_dir / "spatial" / "tissue_positions_list.csv"
    for path in (sc_h5ad, st_h5, positions_csv):
        if not path.exists():
            raise FileNotFoundError(path)

    print("[STEP] load scRNA metadata")
    sc = ad.read_h5ad(sc_h5ad, backed="r")
    if args.cell_type_col not in sc.obs.columns:
        raise ValueError(f"cell_type_col not found: {args.cell_type_col}")
    obs = pd.DataFrame(
        {
            "cell_id": sc.obs_names.astype(str),
            "cell_type": sc.obs[args.cell_type_col].astype(str).str.strip().to_numpy(),
        }
    )
    obs = obs.loc[(obs["cell_type"] != "") & (obs["cell_type"].str.lower() != "nan")].copy()
    obs = _sample_cells_by_type(rng, obs, args.max_cells_total, args.min_cells_per_type)
    obs = obs.reset_index(drop=True)
    cell_ids = obs["cell_id"].astype(str).tolist()
    cell_types = obs["cell_type"].astype(str).tolist()
    print(f"[SC] selected cells={len(cell_ids)} types={obs['cell_type'].nunique()}")
    print(obs["cell_type"].value_counts().to_string())

    print("[STEP] load selected scRNA expression")
    sc_sel = sc[cell_ids, :].to_memory()
    sc.file.close()
    sc_genes = sc_sel.var_names.astype(str).tolist()
    sc_mat = sc_sel.X
    if not sparse.issparse(sc_mat):
        sc_mat = sparse.csr_matrix(sc_mat)
    sc_mat = sc_mat.T.tocsr()

    print("[STEP] load Visium matrix")
    _st_feature_ids, st_genes, st_spots, st_mat = _read_10x_h5(st_h5)

    print("[STEP] align genes")
    sc_gene_to_idx = {g: i for i, g in enumerate(sc_genes)}
    st_gene_to_idx = {g: i for i, g in enumerate(st_genes)}
    common_genes = [g for g in st_genes if g in sc_gene_to_idx]
    if len(common_genes) < 1000:
        raise ValueError(f"Too few common genes: {len(common_genes)}")
    sc_idx = np.array([sc_gene_to_idx[g] for g in common_genes], dtype=np.int32)
    st_idx = np.array([st_gene_to_idx[g] for g in common_genes], dtype=np.int32)
    sc_common = sc_mat[sc_idx, :]
    st_common = st_mat[st_idx, :]
    print(f"[INFO] common_genes={len(common_genes)} cells={sc_common.shape[1]} spots={st_common.shape[1]}")

    print("[STEP] write Stage1 raw inputs")
    _write_expr_tsv(target_dir / "brca_scRNA_GEP.txt", common_genes, cell_ids, sc_common)
    pd.DataFrame({"cell_id": cell_ids, "cell_type": cell_types}).to_csv(
        target_dir / "brca_scRNA_celllabels.txt", sep="\t", index=False
    )
    _write_expr_tsv(target_dir / "brca_STdata_GEP.txt", common_genes, st_spots, st_common)
    _read_positions(positions_csv, st_spots).to_csv(
        target_dir / "brca_STdata_coordinates.txt", sep="\t", index=False
    )

    info = {
        "sample": args.target_sample,
        "source_dir": str(source_dir),
        "n_cells": int(sc_common.shape[1]),
        "n_spots": int(st_common.shape[1]),
        "n_genes": int(len(common_genes)),
        "n_types": int(obs["cell_type"].nunique()),
        "cell_type_counts": obs["cell_type"].value_counts().to_dict(),
        "params": vars(args),
    }
    (target_dir / "real_input_info.json").write_text(
        json.dumps(info, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    _write_dataset_yaml(project_root, args.target_sample)
    print(f"[DONE] prepared: {target_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
