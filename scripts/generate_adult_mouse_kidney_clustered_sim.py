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
    p = argparse.ArgumentParser(
        description=(
            "Generate one spatially-clustered simulation dataset from "
            "scRNA h5ad + Visium filtered matrix."
        )
    )
    p.add_argument("--project_root", default=".", help="Project root.")
    p.add_argument(
        "--source_dir",
        default="data/raw/low_resolution/Adult Mouse Kidney",
        help="Source dataset directory containing adata.h5ad + Visium h5 + spatial positions.",
    )
    p.add_argument("--sc_h5ad", default="adata.h5ad", help="scRNA h5ad file name under source_dir.")
    p.add_argument(
        "--st_h5",
        default="Visium_FFPE_Mouse_Kidney_filtered_feature_bc_matrix.h5",
        help="Visium filtered_feature_bc_matrix.h5 file name under source_dir.",
    )
    p.add_argument(
        "--st_pos_csv",
        default="spatial/tissue_positions_list.csv",
        help="Visium tissue positions csv under source_dir.",
    )
    p.add_argument(
        "--cell_type_col",
        default="Celltype",
        help="Cell type column in adata.obs.",
    )
    p.add_argument(
        "--drop_unknown",
        action="store_true",
        help="Drop cells with unknown-like labels.",
    )
    p.add_argument(
        "--min_source_cells_per_type",
        type=int,
        default=200,
        help="Drop cell types with source SC cell count below this threshold before simulation.",
    )
    p.add_argument(
        "--max_types",
        type=int,
        default=0,
        help="Keep only top-K cell types by source SC abundance after filtering (0 means no cap).",
    )
    p.add_argument(
        "--exclude_cell_types",
        nargs="*",
        default=[],
        help="Cell types to exclude before top-K filtering, useful for removing highly similar subtype pairs.",
    )
    p.add_argument(
        "--target_sample",
        default="adult_mouse_kidney_clustered_sim",
        help="Target sample id (output: data/sim/<sim_group>/<target_sample>).",
    )
    p.add_argument(
        "--sim_group",
        default="adult_mouse_kidney",
        help="Simulation group under data/sim/<sim_group>/.",
    )
    p.add_argument("--seed", type=int, default=42, help="Random seed.")
    p.add_argument(
        "--source_sample_label",
        default="Adult Mouse Kidney",
        help="Source sample label written into sim_info.json.",
    )
    p.add_argument(
        "--simulation_type_label",
        default="spatial_clustered_generic",
        help="Simulation type label written into sim_info.json.",
    )
    p.add_argument("--max_cells_total", type=int, default=12000, help="Max SC cells kept for simulation.")
    p.add_argument("--min_cells_per_type", type=int, default=50, help="Lower bound for each type during sampling.")
    p.add_argument("--centers_per_type", type=int, default=3, help="Number of spatial centers for each type.")
    p.add_argument("--distance_scale", type=float, default=0.65, help="Distance decay scale.")
    p.add_argument("--temperature", type=float, default=0.38, help="Softmax temperature.")
    p.add_argument("--noise_sigma", type=float, default=0.08, help="Noise sigma on type scores.")
    p.add_argument("--mix_alpha", type=float, default=0.85, help="Blend weight for assignment mixture.")
    p.add_argument("--depth_scale", type=float, default=1.0, help="Scale factor for ST depth.")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing target sample dir.")
    return p.parse_args()


def _softmax(x: np.ndarray, axis: int = -1) -> np.ndarray:
    x = x - np.max(x, axis=axis, keepdims=True)
    e = np.exp(x)
    return e / np.clip(np.sum(e, axis=axis, keepdims=True), 1e-12, None)


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
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Gene", *sample_ids])
        n_rows = mat_g_by_s.shape[0]
        for i, g in enumerate(genes):
            w.writerow([g, *mat_g_by_s[i, :].tolist()])
            if (i + 1) % 500 == 0 or i + 1 == n_rows:
                print(f"[WRITE] {path.name}: {i + 1}/{n_rows}")


def _sample_cells_by_type(
    rng: np.random.Generator,
    cell_ids: np.ndarray,
    cell_types: np.ndarray,
    max_cells_total: int,
    min_cells_per_type: int,
) -> np.ndarray:
    type_to_idx: dict[str, np.ndarray] = {}
    for t in sorted(pd.unique(cell_types).tolist()):
        idx = np.where(cell_types == t)[0]
        if len(idx) > 0:
            type_to_idx[t] = idx

    n_types = len(type_to_idx)
    total = sum(len(v) for v in type_to_idx.values())
    if max_cells_total <= 0 or max_cells_total >= total:
        return np.arange(total, dtype=np.int32)
    if max_cells_total < n_types:
        raise ValueError(f"max_cells_total({max_cells_total}) < n_types({n_types})")

    vc = {t: len(idx) for t, idx in type_to_idx.items()}
    base = {t: min(min_cells_per_type, vc[t]) for t in type_to_idx}
    base_sum = sum(base.values())

    if base_sum > max_cells_total:
        # fallback: proportional allocation with >=1 per type
        raw = {t: max(1, int(round(vc[t] / total * max_cells_total))) for t in type_to_idx}
        while sum(raw.values()) > max_cells_total:
            t = max(raw, key=raw.get)
            if raw[t] > 1:
                raw[t] -= 1
            else:
                break
        while sum(raw.values()) < max_cells_total:
            t = max(vc, key=vc.get)
            raw[t] += 1
        target = raw
    else:
        remain = max_cells_total - base_sum
        caps = {t: vc[t] - base[t] for t in type_to_idx}
        cap_sum = sum(max(0, x) for x in caps.values())
        extra = {t: 0 for t in type_to_idx}
        if remain > 0 and cap_sum > 0:
            alloc_float = {t: (max(0, caps[t]) / cap_sum) * remain for t in type_to_idx}
            extra = {t: int(np.floor(v)) for t, v in alloc_float.items()}
            used = sum(extra.values())
            left = remain - used
            if left > 0:
                order = sorted(type_to_idx.keys(), key=lambda k: alloc_float[k] - extra[k], reverse=True)
                i = 0
                while left > 0 and i < len(order):
                    t = order[i]
                    if extra[t] < caps[t]:
                        extra[t] += 1
                        left -= 1
                    i = (i + 1) % len(order)
        target = {t: base[t] + extra[t] for t in type_to_idx}

    keep_idx: list[np.ndarray] = []
    for t, idx in type_to_idx.items():
        k = target[t]
        if k <= 0:
            continue
        if k >= len(idx):
            keep_idx.append(idx)
        else:
            keep_idx.append(rng.choice(idx, size=k, replace=False))
    merged = np.concatenate(keep_idx)
    merged.sort()
    return merged.astype(np.int32, copy=False)


def _write_dataset_yaml(project_root: Path, sample: str, sim_group: str) -> None:
    cfg = {
        "storage": {
            "group": f"simulation_experiments/{sim_group}",
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
                f"data/processed/simulation_experiments/{sim_group}/"
                f"{sample}/stage1_preprocess/hvg_genes.txt"
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


def main() -> int:
    args = parse_args()
    rng = np.random.default_rng(args.seed)

    project_root = Path(args.project_root).resolve()
    source_dir = (project_root / args.source_dir).resolve()
    sc_h5ad = source_dir / args.sc_h5ad
    st_h5 = source_dir / args.st_h5
    st_pos_csv = source_dir / args.st_pos_csv

    for p in (sc_h5ad, st_h5, st_pos_csv):
        if not p.exists():
            raise FileNotFoundError(f"Missing input file: {p}")

    target_dir = project_root / "data" / "sim" / args.sim_group / args.target_sample
    if target_dir.exists():
        if not args.overwrite:
            raise FileExistsError(f"Target exists: {target_dir} (use --overwrite)")
        shutil.rmtree(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)

    print("[STEP] load scRNA h5ad")
    adata = ad.read_h5ad(sc_h5ad)
    if args.cell_type_col not in adata.obs.columns:
        raise ValueError(f"cell_type_col not found in adata.obs: {args.cell_type_col}")

    obs = adata.obs.copy()
    obs["cell_id"] = obs.index.astype(str)
    obs["cell_type"] = obs[args.cell_type_col].astype(str).str.strip()
    if args.drop_unknown:
        bad = {"", "nan", "na", "none", "unknown"}
        mask = ~obs["cell_type"].str.lower().isin(bad)
        obs = obs.loc[mask].copy()

    exclude_types = {str(x).strip() for x in args.exclude_cell_types if str(x).strip()}
    if exclude_types:
        before_exclude_n = len(obs)
        before_exclude_types = int(obs["cell_type"].nunique())
        obs = obs.loc[~obs["cell_type"].isin(exclude_types)].copy()
        print(
            "[SC] excluded cell types: "
            f"{sorted(exclude_types)}; types {before_exclude_types} -> {int(obs['cell_type'].nunique())}, "
            f"cells kept={len(obs)} / {before_exclude_n}"
        )

    if int(args.min_source_cells_per_type) > 0:
        vc_src = obs["cell_type"].value_counts()
        keep_types = vc_src[vc_src >= int(args.min_source_cells_per_type)].index.tolist()
        before_n_types = int(vc_src.shape[0])
        obs = obs.loc[obs["cell_type"].isin(keep_types)].copy()
        after_n_types = int(len(keep_types))
        print(
            "[SC] type prefilter: "
            f"min_source_cells_per_type={int(args.min_source_cells_per_type)} "
            f"types {before_n_types} -> {after_n_types}, cells kept={len(obs)}"
        )
    if int(args.max_types) > 0:
        vc_top = obs["cell_type"].value_counts()
        keep_top = vc_top.head(int(args.max_types)).index.tolist()
        before_top_types = int(vc_top.shape[0])
        obs = obs.loc[obs["cell_type"].isin(keep_top)].copy()
        print(
            "[SC] top-k type filter: "
            f"max_types={int(args.max_types)} types {before_top_types} -> {len(keep_top)}, "
            f"cells kept={len(obs)}"
        )
    if obs.empty:
        raise ValueError("No valid SC cells after filtering.")

    cell_ids_all = obs["cell_id"].to_numpy(dtype=object)
    cell_types_all = obs["cell_type"].to_numpy(dtype=object)

    keep_idx = _sample_cells_by_type(
        rng=rng,
        cell_ids=cell_ids_all,
        cell_types=cell_types_all,
        max_cells_total=int(args.max_cells_total),
        min_cells_per_type=int(args.min_cells_per_type),
    )
    keep_cell_ids = cell_ids_all[keep_idx]
    keep_cell_types = cell_types_all[keep_idx]
    print(
        f"[SC] selected cells: {len(keep_cell_ids)} / {len(cell_ids_all)}; "
        f"types={len(pd.unique(keep_cell_types))}"
    )

    sc_sel = adata[keep_cell_ids, :].copy()
    sc_genes = sc_sel.var_names.astype(str).tolist()

    print("[STEP] load Visium h5 + coordinates")
    st_genes, st_spots, st_mat = _read_visium_h5(st_h5)  # genes x spots
    pos = pd.read_csv(st_pos_csv, header=None)
    if pos.shape[1] < 6:
        raise ValueError(f"Invalid tissue positions file: {st_pos_csv}")
    pos.columns = ["spot_id", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    pos["spot_id"] = pos["spot_id"].astype(str)
    pos = pos.set_index("spot_id")
    coord_df = pos.reindex(st_spots)
    if coord_df["pxl_row"].isna().any():
        raise ValueError("Some ST barcodes are missing in tissue_positions_list.csv")
    coords = coord_df[["pxl_row", "pxl_col"]].to_numpy(dtype=np.float32)

    print("[STEP] align genes")
    st_gene_to_idx = {g: i for i, g in enumerate(st_genes)}
    common_genes = [g for g in sc_genes if g in st_gene_to_idx]
    if len(common_genes) < 500:
        raise ValueError(f"Too few common genes: {len(common_genes)}")
    sc_gene_to_idx = {g: i for i, g in enumerate(sc_genes)}
    sc_idx = np.array([sc_gene_to_idx[g] for g in common_genes], dtype=np.int32)
    st_idx = np.array([st_gene_to_idx[g] for g in common_genes], dtype=np.int32)

    sc_x = sc_sel.X
    if not sparse.issparse(sc_x):
        sc_x = sparse.csr_matrix(sc_x)
    sc_common_csr = sc_x[:, sc_idx].tocsr()  # cells x genes
    st_common = st_mat[st_idx, :].tocsc()  # genes x spots

    n_cells = sc_common_csr.shape[0]
    n_genes = sc_common_csr.shape[1]
    n_spots = st_common.shape[1]
    print(f"[INFO] common_genes={n_genes}, cells={n_cells}, spots={n_spots}")

    type_names = sorted(pd.unique(keep_cell_types).tolist())
    type_to_idx = {t: i for i, t in enumerate(type_names)}
    cell_type_idx = np.array([type_to_idx[t] for t in keep_cell_types], dtype=np.int32)
    n_types = len(type_names)

    print("[STEP] build clustered spot-type field")
    xy = coords.copy()
    xy[:, 0] = (xy[:, 0] - xy[:, 0].mean()) / (xy[:, 0].std() + 1e-6)
    xy[:, 1] = (xy[:, 1] - xy[:, 1].mean()) / (xy[:, 1].std() + 1e-6)
    diff = xy[:, None, :] - xy[None, :, :]
    dist = np.sqrt((diff * diff).sum(axis=2)).astype(np.float32, copy=False)
    del diff

    type_counts = pd.Series(keep_cell_types).value_counts().reindex(type_names).fillna(0).to_numpy(dtype=np.float32)
    type_prior = np.clip(type_counts / np.clip(type_counts.sum(), 1.0, None), 1e-8, None)

    scores = np.zeros((n_spots, n_types), dtype=np.float32)
    n_centers = max(1, int(args.centers_per_type))
    for t_idx in range(n_types):
        center_idx = rng.choice(n_spots, size=min(n_centers, n_spots), replace=False)
        min_dist = dist[:, center_idx].min(axis=1)
        scores[:, t_idx] = (
            -(min_dist / max(float(args.distance_scale), 1e-6))
            + np.log(type_prior[t_idx] + 1e-8)
            + rng.normal(0.0, float(args.noise_sigma), size=n_spots).astype(np.float32)
        )

    spot_type_probs = _softmax(scores / max(float(args.temperature), 1e-6), axis=1).astype(np.float32)
    del scores
    del dist

    print("[STEP] assign SC cells to spots")
    assigned_spot_idx = np.empty(n_cells, dtype=np.int32)
    for t_idx, t_name in enumerate(type_names):
        idx = np.where(cell_type_idx == t_idx)[0]
        p = spot_type_probs[:, t_idx].astype(np.float64)
        p = p / np.clip(p.sum(), 1e-12, None)
        assigned_spot_idx[idx] = rng.choice(n_spots, size=len(idx), replace=True, p=p)
        print(f"[ASSIGN] {t_name}: {len(idx)}")
    assigned_spot_ids = np.array(st_spots, dtype=object)[assigned_spot_idx]

    print("[STEP] derive truth tables")
    query_id = f"{args.target_sample}_q0"
    truth_query = pd.DataFrame(
        {
            "query_id": query_id,
            "cell_id": keep_cell_ids,
            "true_spot_id": assigned_spot_ids,
            "cell_type": keep_cell_types,
        }
    )
    count_tab = pd.crosstab(
        pd.Categorical(assigned_spot_ids, categories=st_spots),
        pd.Categorical(keep_cell_types, categories=type_names),
    )
    count_tab = count_tab.reindex(index=st_spots, fill_value=0).reindex(columns=type_names, fill_value=0)
    frac_tab = count_tab.div(np.clip(count_tab.sum(axis=1), 1, None), axis=0).astype(np.float32)

    print("[STEP] build type profiles + simulate ST counts")
    type_profile = np.zeros((n_genes, n_types), dtype=np.float32)
    for t_idx, _t_name in enumerate(type_names):
        idx = np.where(cell_type_idx == t_idx)[0]
        sub = sc_common_csr[idx, :]
        mean_vec = np.asarray(sub.mean(axis=0)).ravel().astype(np.float32)
        type_profile[:, t_idx] = mean_vec + 1e-3
    type_profile = type_profile / np.clip(type_profile.sum(axis=0, keepdims=True), 1e-8, None)

    weights = frac_tab.to_numpy(dtype=np.float32, copy=True)
    row_sum = weights.sum(axis=1)
    zero_mask = row_sum <= 0
    nz_mask = ~zero_mask
    alpha = float(np.clip(args.mix_alpha, 0.0, 1.0))
    if nz_mask.any():
        weights[nz_mask] = alpha * weights[nz_mask] + (1.0 - alpha) * spot_type_probs[nz_mask]
    if zero_mask.any():
        weights[zero_mask] = spot_type_probs[zero_mask]
    weights = weights / np.clip(weights.sum(axis=1, keepdims=True), 1e-8, None)

    lib_sizes = np.asarray(st_common.sum(axis=0)).ravel().astype(np.float64)
    lib_sizes = np.clip(lib_sizes, 1.0, None) * float(args.depth_scale)
    expected = type_profile @ weights.T
    expected = expected * lib_sizes[np.newaxis, :]
    sim_counts = rng.poisson(np.clip(expected, 0.0, None)).astype(np.int32)

    print("[STEP] write simulation files")
    sc_expr_path = target_dir / "brca_scRNA_GEP.txt"
    sc_meta_path = target_dir / "brca_scRNA_celllabels.txt"
    st_expr_path = target_dir / "brca_STdata_GEP.txt"
    st_meta_path = target_dir / "brca_STdata_coordinates.txt"
    truth_query_path = target_dir / "sim_truth_query_cell_spot.csv"
    truth_spot_soft_path = target_dir / "sim_truth_spot_type_fraction.csv"
    truth_spot_cell_path = target_dir / "sim_truth_spot_type_fraction_from_cells.csv"
    dominant_path = target_dir / "sim_truth_spot_dominant_type.csv"
    sim_info_path = target_dir / "sim_info.json"

    # SC expression: genes x cells
    sc_gxc = sc_common_csr.T.toarray().astype(np.int32, copy=False)
    _write_expr_tsv(sc_expr_path, common_genes, keep_cell_ids.tolist(), sc_gxc)

    sc_meta_df = pd.DataFrame({"cell_id": keep_cell_ids, "cell_type": keep_cell_types})
    sc_meta_df.to_csv(sc_meta_path, sep="\t", index=False)

    _write_expr_tsv(st_expr_path, common_genes, st_spots, sim_counts)

    st_meta_df = pd.DataFrame(
        {
            "spot_id": st_spots,
            "row": coords[:, 0],
            "col": coords[:, 1],
        }
    )
    st_meta_df.to_csv(st_meta_path, sep="\t", index=False)

    truth_query.to_csv(truth_query_path, index=False, encoding="utf-8")
    truth_spot_soft = pd.DataFrame(weights, index=st_spots, columns=type_names)
    truth_spot_soft.insert(0, "spot_id", truth_spot_soft.index)
    truth_spot_soft.to_csv(truth_spot_soft_path, index=False, encoding="utf-8")
    truth_spot_cell = frac_tab.copy()
    truth_spot_cell.insert(0, "spot_id", truth_spot_cell.index)
    truth_spot_cell.to_csv(truth_spot_cell_path, index=False, encoding="utf-8")
    dominant = pd.DataFrame(
        {
            "spot_id": st_spots,
            "dominant_type": np.array(type_names, dtype=object)[
                np.argmax(truth_spot_soft.drop(columns=["spot_id"]).to_numpy(dtype=np.float32), axis=1)
            ],
            "assigned_cells": count_tab.sum(axis=1).to_numpy(dtype=np.int32),
        }
    )
    dominant.to_csv(dominant_path, index=False, encoding="utf-8")

    sim_info = {
        "sample": args.target_sample,
        "source_sample": str(args.source_sample_label),
        "simulation_type": str(args.simulation_type_label),
        "seed": int(args.seed),
        "query_id": query_id,
        "missing_type": None,
        "replacement_type": None,
        "cells_per_spot": float(count_tab.sum(axis=1).mean()),
        "n_cells": int(n_cells),
        "n_spots": int(n_spots),
        "n_types": int(n_types),
        "n_genes": int(n_genes),
        "cell_type_col": args.cell_type_col,
        "params": {
            "drop_unknown": bool(args.drop_unknown),
            "exclude_cell_types": sorted(exclude_types),
            "min_source_cells_per_type": int(args.min_source_cells_per_type),
            "max_types": int(args.max_types),
            "max_cells_total": int(args.max_cells_total),
            "min_cells_per_type": int(args.min_cells_per_type),
            "centers_per_type": int(args.centers_per_type),
            "distance_scale": float(args.distance_scale),
            "temperature": float(args.temperature),
            "noise_sigma": float(args.noise_sigma),
            "mix_alpha": float(args.mix_alpha),
            "depth_scale": float(args.depth_scale),
        },
        "files": {
            "sc_expr": str(sc_expr_path),
            "sc_meta": str(sc_meta_path),
            "st_expr": str(st_expr_path),
            "st_meta": str(st_meta_path),
            "sim_truth_query_cell_spot": str(truth_query_path),
            "sim_truth_spot_type_fraction": str(truth_spot_soft_path),
            "sim_truth_spot_type_fraction_from_cells": str(truth_spot_cell_path),
            "sim_truth_spot_dominant_type": str(dominant_path),
        },
    }
    sim_info_path.write_text(json.dumps(sim_info, indent=2, ensure_ascii=False), encoding="utf-8")

    _write_dataset_yaml(project_root, args.target_sample, args.sim_group)

    print("[DONE] clustered simulation generated.")
    print(f"[OUT] {target_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
