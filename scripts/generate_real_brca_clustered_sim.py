#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.utils.sample_paths import resolve_sample_dir
from src.stages.storage import raw_dir, read_dataset_config


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate spatially-clustered simulation data from real_brca raw files."
    )
    p.add_argument("--project_root", default=".", help="Project root path.")
    p.add_argument("--source_sample", default="real_brca", help="Source sample id.")
    p.add_argument(
        "--target_sample",
        default="real_brca_clustered_sim",
        help="Target sample name under data/sim/<sim_group>/.",
    )
    p.add_argument(
        "--sim_group",
        default="real_brca",
        help="simulation group under data/sim/<sim_group>/ (default: real_brca)",
    )
    p.add_argument("--seed", type=int, default=42, help="Random seed.")
    p.add_argument(
        "--centers_per_type",
        type=int,
        default=3,
        help="Number of spatial centers per cell type (larger => more fragmented regions).",
    )
    p.add_argument(
        "--distance_scale",
        type=float,
        default=0.65,
        help="Distance decay scale in normalized coordinate space.",
    )
    p.add_argument(
        "--temperature",
        type=float,
        default=0.38,
        help="Softmax temperature for spot-type probabilities (smaller => stronger dominance).",
    )
    p.add_argument(
        "--noise_sigma",
        type=float,
        default=0.08,
        help="Noise level added to spatial type score.",
    )
    p.add_argument(
        "--mix_alpha",
        type=float,
        default=0.85,
        help="Blend weight for assignment-derived fractions when generating ST expression.",
    )
    p.add_argument(
        "--merge_cell_type",
        action="append",
        default=[],
        metavar="FROM=TO",
        help=(
            "Merge one source cell type into a coarse target label. "
            "Can be provided multiple times, e.g. --merge_cell_type 'B cell=B/Plasma B cell'."
        ),
    )
    p.add_argument(
        "--exclude_cell_types",
        nargs="*",
        default=[],
        help="Cell types to exclude after optional merging.",
    )
    p.add_argument(
        "--simulation_type_label",
        default="spatial_clustered_standard_raw",
        help="Simulation type label written to sim_info.json.",
    )
    p.add_argument(
        "--depth_scale",
        type=float,
        default=1.0,
        help="Scale factor for simulated ST library size.",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite target directory if it exists.",
    )
    return p.parse_args()


def _softmax(x: np.ndarray, axis: int = -1) -> np.ndarray:
    x = x - np.max(x, axis=axis, keepdims=True)
    e = np.exp(x)
    return e / np.clip(np.sum(e, axis=axis, keepdims=True), 1e-12, None)


def read_expr_tsv(path: Path) -> tuple[np.ndarray, list[str], np.ndarray]:
    print(f"[READ] {path}")
    df = pd.read_csv(path, sep="\t", low_memory=False)
    genes = df.iloc[:, 0].astype(str).to_numpy()
    sample_ids = [str(c) for c in df.columns[1:]]
    mat = df.iloc[:, 1:].to_numpy(dtype=np.float32, copy=False)
    return genes, sample_ids, mat


def write_expr_tsv(path: Path, genes: np.ndarray, spot_ids: list[str], mat: np.ndarray) -> None:
    print(f"[WRITE] {path}")
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Gene", *spot_ids])
        n_rows = mat.shape[0]
        for i, g in enumerate(genes):
            w.writerow([g, *mat[i, :].tolist()])
            if (i + 1) % 500 == 0 or (i + 1) == n_rows:
                print(f"[WRITE] rows: {i + 1}/{n_rows}")


def parse_merge_map(items: list[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"Invalid --merge_cell_type value, expected FROM=TO: {item}")
        src, dst = item.split("=", 1)
        src = src.strip()
        dst = dst.strip()
        if not src or not dst:
            raise ValueError(f"Invalid --merge_cell_type value, expected non-empty FROM=TO: {item}")
        out[src] = dst
    return out


def write_sim_dataset_yaml(project_root: Path, target_sample: str, sim_group: str, source_cfg: dict) -> None:
    cfg = dict(source_cfg or {})
    cfg.setdefault("paths", {
        "sc_expr": "brca_scRNA_GEP.txt",
        "sc_meta": "brca_scRNA_celllabels.txt",
        "st_expr": "brca_STdata_GEP.txt",
        "st_meta": "brca_STdata_coordinates.txt",
        "svg_marker_whitelist": None,
    })
    cfg.setdefault("qc", {
        "sc_min_genes": 100,
        "sc_max_genes": 10000,
        "sc_max_mt": 20,
        "st_min_genes": 50,
        "st_max_genes": "Inf",
        "st_max_mt": 40,
        "hvg_nfeatures": 2000,
        "mt_pattern": "^(mt-|MT-)",
    })
    cfg.setdefault("gene_filter", {"min_cells_sc": 0, "min_cells_st": 0})
    cfg.setdefault("stage3", {})
    cfg["storage"] = {"group": f"simulation_experiments/{sim_group}"}
    stage3 = cfg["stage3"]
    stage3["plugin_genes_path"] = (
        f"data/processed/simulation_experiments/{sim_group}/{target_sample}/stage1_preprocess/hvg_genes.txt"
    )
    stage3.setdefault("strong_th", 0.7)
    stage3.setdefault("weak_th", 0.4)
    stage3.setdefault("st_cluster_k", 30)
    stage3.setdefault("unknown_floor", 0.3)
    stage3.setdefault("min_cells_rare_type", 20)
    stage3.setdefault("eps", 1.0e-8)
    auto = stage3.setdefault("auto_missing_detection", {})
    auto.setdefault("enable", True)
    auto.setdefault("method", "adaptive_low_support")
    auto.setdefault("min_cells", 50)
    auto.setdefault("robust_z_th", -2.0)
    auto.setdefault("soft_z_th", -0.9)
    auto.setdefault("require_masked_for_soft", False)
    auto.setdefault("require_masked_for_hard", True)
    auto.setdefault("require_confirmation", True)
    auto.setdefault("confirmation_use_masked_missing", True)
    auto.setdefault("confirmation_use_marker_identity", True)
    auto.setdefault("confirmation_marker_identity_z_th", -1.5)
    auto.setdefault("confirmation_marker_support_score_th", 0.5)
    auto.setdefault("confirmation_marker_max_support_score", None)
    auto.setdefault("max_fraction_types", 0.3)
    auto.setdefault("max_types", 2)
    auto.setdefault("action", "mark_unknown")
    masked = stage3.setdefault("masked_missing_detection", {})
    masked.setdefault("enable", True)
    masked.setdefault("apply_to_auto_missing", True)
    masked.setdefault("neighbor_cosine_th", 0.9)
    masked.setdefault("neighbor_cell_ratio_min", 1.0)
    masked.setdefault("marker_top_n", 20)
    masked.setdefault("min_identity_markers", 3)
    masked.setdefault("min_marker_specificity", 1.2)
    masked.setdefault("min_marker_type_mean", 0.001)
    masked.setdefault("min_marker_st_detect_frac", 0.005)
    masked.setdefault("st_presence_quantile", 0.9)
    masked.setdefault("identity_z_th", -0.8)
    masked.setdefault("pressure_z_th", 1.0)
    masked.setdefault("min_support_score_for_apply", 0.8)
    masked.setdefault("max_types", 2)
    diag = stage3.setdefault("marker_identity_diagnostics", {})
    diag.setdefault("enable", True)
    diag.setdefault("marker_top_n", 80)
    diag.setdefault("min_identity_markers", 5)
    diag.setdefault("min_all_specificity", 1.3)
    diag.setdefault("min_neighbor_specificity", 1.1)
    diag.setdefault("min_marker_type_mean", 0.001)
    diag.setdefault("min_marker_st_detect_frac", 0.005)
    diag.setdefault("st_presence_quantile", 0.9)
    diag.setdefault("depleted_z_th", -0.9)

    out = project_root / "configs" / "datasets" / f"{target_sample}.yaml"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True), encoding="utf-8")
    print(f"[CFG] wrote: {out}")


def main() -> int:
    args = parse_args()
    rng = np.random.default_rng(args.seed)

    project_root = Path(args.project_root).resolve()
    source_cfg = read_dataset_config(project_root, args.source_sample)
    raw_source_dir = raw_dir(project_root, args.source_sample, source_cfg)
    src_dir = raw_source_dir if raw_source_dir.exists() else resolve_sample_dir(
        project_root,
        args.source_sample,
        sim_group=args.sim_group,
        must_exist=True,
    )
    dst_dir = project_root / "data" / "sim" / args.sim_group / args.target_sample

    if not src_dir.exists():
        raise FileNotFoundError(f"Source sample dir not found: {src_dir}")
    if dst_dir.exists():
        if not args.overwrite:
            raise FileExistsError(f"Target dir exists: {dst_dir} (use --overwrite)")
        shutil.rmtree(dst_dir)
    dst_dir.mkdir(parents=True, exist_ok=True)

    sc_expr_src = src_dir / "brca_scRNA_GEP.txt"
    sc_meta_src = src_dir / "brca_scRNA_celllabels.txt"
    st_expr_src = src_dir / "brca_STdata_GEP.txt"
    st_meta_src = src_dir / "brca_STdata_coordinates.txt"
    for p in [sc_expr_src, sc_meta_src, st_expr_src, st_meta_src]:
        if not p.exists():
            raise FileNotFoundError(f"Missing source file: {p}")

    print("[STEP] Load SC metadata")
    sc_meta = pd.read_csv(sc_meta_src, sep="\t")
    if sc_meta.shape[1] < 2:
        raise ValueError("sc metadata must have >=2 columns")
    sc_meta = sc_meta.rename(columns={sc_meta.columns[0]: "cell_id", sc_meta.columns[1]: "cell_type"})
    sc_meta = sc_meta[["cell_id", "cell_type"]].drop_duplicates("cell_id")
    sc_meta["cell_id"] = sc_meta["cell_id"].astype(str)
    sc_meta["cell_type"] = sc_meta["cell_type"].astype(str)
    merge_map = parse_merge_map(args.merge_cell_type)
    if merge_map:
        before_types = int(sc_meta["cell_type"].nunique())
        sc_meta["cell_type"] = sc_meta["cell_type"].replace(merge_map)
        print(
            "[SC] merged cell types: "
            f"{merge_map}; types {before_types} -> {int(sc_meta['cell_type'].nunique())}"
        )
    exclude_types = {str(x).strip() for x in args.exclude_cell_types if str(x).strip()}
    if exclude_types:
        before_n = len(sc_meta)
        before_types = int(sc_meta["cell_type"].nunique())
        sc_meta = sc_meta.loc[~sc_meta["cell_type"].isin(exclude_types)].copy()
        print(
            "[SC] excluded cell types: "
            f"{sorted(exclude_types)}; types {before_types} -> {int(sc_meta['cell_type'].nunique())}, "
            f"cells kept={len(sc_meta)} / {before_n}"
        )

    print("[STEP] Load SC expression")
    sc_genes, sc_cells, sc_mat = read_expr_tsv(sc_expr_src)
    sc_cells_arr = np.array(sc_cells, dtype=object)
    meta_map = dict(zip(sc_meta["cell_id"], sc_meta["cell_type"]))
    keep_mask = np.array([c in meta_map for c in sc_cells_arr], dtype=bool)
    if not keep_mask.any():
        raise ValueError("No overlapping cells between SC expression and metadata.")
    sc_cells_arr = sc_cells_arr[keep_mask]
    sc_mat = sc_mat[:, keep_mask]
    cell_types_arr = np.array([meta_map[c] for c in sc_cells_arr], dtype=object)

    print("[STEP] Load ST expression")
    st_genes, st_spots, st_mat = read_expr_tsv(st_expr_src)
    st_spots_arr = np.array(st_spots, dtype=object)

    print("[STEP] Load ST coordinates")
    st_meta = pd.read_csv(st_meta_src, sep="\t")
    if st_meta.shape[1] < 3:
        raise ValueError("st coordinates must have >=3 columns")
    st_meta = st_meta.rename(
        columns={st_meta.columns[0]: "spot_id", st_meta.columns[1]: "row", st_meta.columns[2]: "col"}
    )
    st_meta["spot_id"] = st_meta["spot_id"].astype(str)
    coord_df = st_meta.set_index("spot_id")
    missing_spots = [s for s in st_spots_arr if s not in coord_df.index]
    if missing_spots:
        raise ValueError(f"{len(missing_spots)} spots missing coordinates. example={missing_spots[:3]}")
    coords = coord_df.loc[st_spots_arr, ["row", "col"]].to_numpy(dtype=np.float32)

    print("[STEP] Align genes (SC ∩ ST)")
    sc_gene_to_idx = {g: i for i, g in enumerate(sc_genes)}
    common_genes = [g for g in st_genes if g in sc_gene_to_idx]
    if len(common_genes) < 1000:
        raise ValueError(f"Too few common genes: {len(common_genes)}")
    st_gene_to_idx = {g: i for i, g in enumerate(st_genes)}
    sc_idx = np.array([sc_gene_to_idx[g] for g in common_genes], dtype=np.int32)
    st_idx = np.array([st_gene_to_idx[g] for g in common_genes], dtype=np.int32)
    sc_common = sc_mat[sc_idx, :]
    st_common = st_mat[st_idx, :]
    common_genes_arr = np.array(common_genes, dtype=object)
    del sc_mat
    del st_mat

    type_names = sorted(pd.unique(cell_types_arr).tolist())
    type_to_idx = {t: i for i, t in enumerate(type_names)}
    cell_type_idx = np.array([type_to_idx[t] for t in cell_types_arr], dtype=np.int32)
    n_cells = len(sc_cells_arr)
    n_spots = len(st_spots_arr)
    n_types = len(type_names)

    print(
        f"[INFO] cells={n_cells}, spots={n_spots}, types={n_types}, genes_common={len(common_genes_arr)}"
    )

    print("[STEP] Build spatially clustered type field")
    xy = coords.copy()
    xy[:, 0] = (xy[:, 0] - xy[:, 0].mean()) / (xy[:, 0].std() + 1e-6)
    xy[:, 1] = (xy[:, 1] - xy[:, 1].mean()) / (xy[:, 1].std() + 1e-6)
    diff = xy[:, None, :] - xy[None, :, :]
    dist = np.sqrt((diff * diff).sum(axis=2)).astype(np.float32, copy=False)
    del diff

    type_counts = pd.Series(cell_types_arr).value_counts().reindex(type_names).fillna(0).to_numpy(dtype=np.float32)
    type_prior = np.clip(type_counts / np.clip(type_counts.sum(), 1.0, None), 1e-8, None)

    scores = np.zeros((n_spots, n_types), dtype=np.float32)
    n_centers = max(1, int(args.centers_per_type))
    for t_idx in range(n_types):
        center_idx = rng.choice(n_spots, size=min(n_centers, n_spots), replace=False)
        min_dist = dist[:, center_idx].min(axis=1)
        scores[:, t_idx] = (
            -(min_dist / max(args.distance_scale, 1e-6))
            + np.log(type_prior[t_idx] + 1e-8)
            + rng.normal(0.0, args.noise_sigma, size=n_spots).astype(np.float32)
        )

    spot_type_probs = _softmax(scores / max(args.temperature, 1e-6), axis=1).astype(np.float32)
    del scores
    del dist

    print("[STEP] Assign each SC cell to a spatial spot (cluster-aware)")
    assigned_spot_idx = np.empty(n_cells, dtype=np.int32)
    for t_idx, t_name in enumerate(type_names):
        mask = cell_type_idx == t_idx
        idx = np.where(mask)[0]
        p = spot_type_probs[:, t_idx].astype(np.float64)
        p_sum = p.sum()
        if p_sum <= 0:
            p = np.ones_like(p) / len(p)
        else:
            p = p / p_sum
        assigned_spot_idx[idx] = rng.choice(n_spots, size=len(idx), replace=True, p=p)
        print(f"[ASSIGN] {t_name}: n_cells={len(idx)}")

    assigned_spot_ids = st_spots_arr[assigned_spot_idx]

    print("[STEP] Build truth tables")
    query_id = f"{args.target_sample}_q0"
    truth_query = pd.DataFrame(
        {
            "query_id": query_id,
            "cell_id": sc_cells_arr,
            "true_spot_id": assigned_spot_ids,
            "cell_type": cell_types_arr,
        }
    )

    count_tab = pd.crosstab(
        pd.Categorical(assigned_spot_ids, categories=st_spots_arr),
        pd.Categorical(cell_types_arr, categories=type_names),
    )
    count_tab.index = count_tab.index.astype(str)
    count_tab.columns = count_tab.columns.astype(str)
    count_tab = count_tab.reindex(index=st_spots_arr, fill_value=0).reindex(columns=type_names, fill_value=0)
    frac_tab = count_tab.div(np.clip(count_tab.sum(axis=1), 1, None), axis=0).astype(np.float32)

    print("[STEP] Build cell-type expression profiles")
    one_hot = np.zeros((n_cells, n_types), dtype=np.float32)
    one_hot[np.arange(n_cells), cell_type_idx] = 1.0
    type_sum = sc_common @ one_hot
    type_profile = type_sum + 1e-3
    type_profile = type_profile / np.clip(type_profile.sum(axis=0, keepdims=True), 1e-8, None)

    print("[STEP] Simulate ST expression from clustered type mixtures")
    weights = frac_tab.to_numpy(dtype=np.float32, copy=True)
    row_sum = weights.sum(axis=1)
    zero_mask = row_sum <= 0
    nonzero_mask = ~zero_mask
    alpha = float(np.clip(args.mix_alpha, 0.0, 1.0))
    if nonzero_mask.any():
        weights[nonzero_mask] = alpha * weights[nonzero_mask] + (1.0 - alpha) * spot_type_probs[nonzero_mask]
    if zero_mask.any():
        weights[zero_mask] = spot_type_probs[zero_mask]
    weights = weights / np.clip(weights.sum(axis=1, keepdims=True), 1e-8, None)

    lib_sizes = st_common.sum(axis=0, dtype=np.float64)
    lib_sizes = np.clip(lib_sizes, 1.0, None) * float(args.depth_scale)

    expected = type_profile @ weights.T
    expected = expected * lib_sizes[np.newaxis, :]
    sim_counts = rng.poisson(np.clip(expected, 0.0, None)).astype(np.int32)

    print("[STEP] Write target sample files")
    sc_expr_dst = dst_dir / "brca_scRNA_GEP.txt"
    sc_meta_dst = dst_dir / "brca_scRNA_celllabels.txt"
    st_expr_dst = dst_dir / "brca_STdata_GEP.txt"
    st_meta_dst = dst_dir / "brca_STdata_coordinates.txt"
    truth_query_dst = dst_dir / "sim_truth_query_cell_spot.csv"
    truth_spot_dst = dst_dir / "sim_truth_spot_type_fraction.csv"
    truth_spot_from_cells_dst = dst_dir / "sim_truth_spot_type_fraction_from_cells.csv"
    sim_info_dst = dst_dir / "sim_info.json"
    dominant_spot_dst = dst_dir / "sim_truth_spot_dominant_type.csv"

    shutil.copy2(sc_expr_src, sc_expr_dst)
    pd.DataFrame({"cell_id": sc_cells_arr, "cell_type": cell_types_arr}).to_csv(
        sc_meta_dst,
        sep="\t",
        index=False,
        encoding="utf-8",
    )
    shutil.copy2(st_meta_src, st_meta_dst)
    write_expr_tsv(st_expr_dst, common_genes_arr, st_spots_arr.tolist(), sim_counts)

    truth_query.to_csv(truth_query_dst, index=False, encoding="utf-8")

    truth_spot_soft = pd.DataFrame(weights, index=st_spots_arr, columns=type_names)
    truth_spot_soft.insert(0, "spot_id", truth_spot_soft.index)
    truth_spot_soft.to_csv(truth_spot_dst, index=False, encoding="utf-8")

    truth_spot_from_cells = frac_tab.copy()
    truth_spot_from_cells.insert(0, "spot_id", truth_spot_from_cells.index)
    truth_spot_from_cells.to_csv(truth_spot_from_cells_dst, index=False, encoding="utf-8")

    dominant = pd.DataFrame(
        {
            "spot_id": st_spots_arr,
            "dominant_type": np.array(type_names, dtype=object)[
                np.argmax(truth_spot_soft.drop(columns=["spot_id"]).to_numpy(), axis=1)
            ],
            "assigned_cells": count_tab.sum(axis=1).to_numpy(dtype=np.int32),
        }
    )
    dominant.to_csv(dominant_spot_dst, index=False, encoding="utf-8")

    sim_info = {
        "sample": args.target_sample,
        "source_sample": args.source_sample,
        "simulation_type": args.simulation_type_label,
        "seed": args.seed,
        "query_id": query_id,
        "missing_type": None,
        "cells_per_spot": float(count_tab.sum(axis=1).mean()),
        "n_cells": int(n_cells),
        "n_spots": int(n_spots),
        "n_types": int(n_types),
        "n_genes": int(len(common_genes_arr)),
        "params": {
            "centers_per_type": int(args.centers_per_type),
            "distance_scale": float(args.distance_scale),
            "temperature": float(args.temperature),
            "noise_sigma": float(args.noise_sigma),
            "mix_alpha": float(args.mix_alpha),
            "depth_scale": float(args.depth_scale),
            "merge_cell_type": merge_map,
            "exclude_cell_types": sorted(exclude_types),
        },
        "files": {
            "sc_expr": str(sc_expr_dst),
            "sc_meta": str(sc_meta_dst),
            "st_expr": str(st_expr_dst),
            "st_meta": str(st_meta_dst),
            "sim_truth_query_cell_spot": str(truth_query_dst),
            "sim_truth_spot_type_fraction": str(truth_spot_dst),
            "sim_truth_spot_type_fraction_from_cells": str(truth_spot_from_cells_dst),
            "sim_truth_spot_dominant_type": str(dominant_spot_dst),
        },
    }
    sim_info_dst.write_text(json.dumps(sim_info, indent=2, ensure_ascii=False), encoding="utf-8")
    write_sim_dataset_yaml(project_root, args.target_sample, args.sim_group, source_cfg)

    print("[DONE] clustered simulation generated.")
    print(f"[OUT] {dst_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
