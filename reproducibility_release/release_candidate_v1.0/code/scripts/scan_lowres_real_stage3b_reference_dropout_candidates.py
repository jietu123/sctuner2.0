#!/usr/bin/env python
from __future__ import annotations

import argparse
from collections import deque
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_SAMPLES = [
    "adult_mouse_kidney_real",
    "ffpe_mouse_brain_sagittal_real",
    "human_breast_cancer_real",
    "human_breast_cancer_visium_ff_wta_real",
    "human_breast_cancer_wta_120_real",
    "human_cervical_cancer_real",
    "human_heart_ff_real",
    "human_intestine_cancer_real",
    "human_lymph_node_real",
    "mouse_embryo_real",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Scan low-resolution real samples for Stage3B reference-dropout candidates."
    )
    p.add_argument("--project_root", default=".")
    p.add_argument("--storage_group", default="low_resolution_experiments")
    p.add_argument("--out_dir", default="visualizations/stage3b_realdata_candidate_scan")
    p.add_argument("--min_sc_cells", type=int, default=20)
    p.add_argument("--n_marker_genes", type=int, default=30)
    p.add_argument("--samples", nargs="*", default=DEFAULT_SAMPLES)
    return p.parse_args()


def read_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str)
    return df.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def coord_columns(coords: pd.DataFrame) -> tuple[str, str]:
    cols = {c.lower(): c for c in coords.columns}
    if "row" in cols and "col" in cols:
        return cols["row"], cols["col"]
    if "y" in cols and "x" in cols:
        return cols["y"], cols["x"]
    numeric = [c for c in coords.columns if pd.api.types.is_numeric_dtype(coords[c])]
    if len(numeric) >= 2:
        return numeric[0], numeric[1]
    raise ValueError("Cannot infer coordinate columns")


def knn_indices(xy: np.ndarray, k: int = 6) -> np.ndarray:
    try:
        from scipy.spatial import cKDTree

        _, idx = cKDTree(xy).query(xy, k=min(k + 1, len(xy)))
        return idx[:, 1:]
    except Exception:
        dist = ((xy[:, None, :] - xy[None, :, :]) ** 2).sum(axis=2)
        order = np.argsort(dist, axis=1)
        return order[:, 1 : k + 1]


def moran_knn(score: np.ndarray, neigh: np.ndarray) -> float:
    z = score.astype(float) - float(np.mean(score))
    denom = float(np.sum(z * z))
    if denom <= 0 or neigh.size == 0:
        return 0.0
    num = float(np.sum(z[:, None] * z[neigh]))
    n = len(z)
    w = neigh.size
    return (n / w) * num / denom


def largest_component_fraction(mask: np.ndarray, neigh: np.ndarray) -> float:
    nodes = set(np.flatnonzero(mask).tolist())
    if not nodes:
        return 0.0
    seen: set[int] = set()
    largest = 0
    for start in list(nodes):
        if start in seen:
            continue
        q: deque[int] = deque([start])
        seen.add(start)
        size = 0
        while q:
            cur = q.popleft()
            size += 1
            for nb in neigh[cur]:
                nb = int(nb)
                if nb in nodes and nb not in seen:
                    seen.add(nb)
                    q.append(nb)
        largest = max(largest, size)
    return largest / len(nodes)


def top_region_metrics(score: pd.Series, neigh: np.ndarray, quantile: float) -> tuple[int, float, float, float]:
    threshold = float(score.quantile(quantile))
    mask = score.to_numpy() >= threshold
    n_top = int(mask.sum())
    if n_top <= 0:
        return 0, 0.0, 0.0, 0.0
    purity = float(mask[neigh].mean(axis=1)[mask].mean()) if neigh.size else 0.0
    frac = n_top / len(mask)
    enrichment = purity / frac if frac > 0 else 0.0
    largest = largest_component_fraction(mask, neigh)
    return n_top, purity, enrichment, largest


def marker_scores(st_expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in st_expr.columns]
    if not present:
        raise ValueError("No marker genes present in ST expression")
    return st_expr[present].mean(axis=1).astype(float)


def scan_sample(root: Path, storage_group: str, sample: str, min_sc_cells: int, n_marker_genes: int) -> list[dict]:
    export = root / "data" / "processed" / storage_group / sample / "stage1_preprocess" / "exported"
    if not export.exists():
        raise FileNotFoundError(f"Missing Stage1 export: {export}")

    meta = pd.read_csv(export / "sc_metadata.csv", index_col=0)
    meta.index = meta.index.astype(str)
    if "cell_type" not in meta.columns:
        raise ValueError(f"sc_metadata.csv lacks cell_type: {sample}")
    meta["cell_type"] = meta["cell_type"].astype(str)

    sc_expr = read_matrix(export / "sc_expression_normalized.csv")
    st_expr = read_matrix(export / "st_expression_normalized.csv")
    common_cells = sc_expr.index.intersection(meta.index)
    sc_expr = sc_expr.loc[common_cells]
    meta = meta.loc[common_cells]

    common_genes = [g for g in sc_expr.columns if g in st_expr.columns]
    if not common_genes:
        raise ValueError(f"No common SC/ST genes: {sample}")
    sc_expr = sc_expr[common_genes]
    st_expr = st_expr[common_genes]

    counts = meta["cell_type"].value_counts()
    keep_types = counts[counts >= min_sc_cells].index.tolist()
    if len(keep_types) < 2:
        return []

    grouped = sc_expr.groupby(meta["cell_type"]).mean()
    global_sum = sc_expr.sum(axis=0)
    n_total = len(sc_expr)

    coords = pd.read_csv(export / "st_coordinates.csv", index_col=0)
    coords.index = coords.index.astype(str)
    common_spots = st_expr.index.intersection(coords.index)
    st_expr = st_expr.loc[common_spots]
    coords = coords.loc[common_spots]
    y_col, x_col = coord_columns(coords)
    xy = coords[[x_col, y_col]].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy()
    neigh = knn_indices(xy, k=6)

    rows: list[dict] = []
    eps = 1e-6
    for cell_type in keep_types:
        n_type = int(counts[cell_type])
        type_mean = grouped.loc[cell_type]
        other_n = max(n_total - n_type, 1)
        other_mean = (global_sum - type_mean * n_type) / other_n
        specificity = (type_mean + eps) / (other_mean + eps)
        expressed = type_mean[type_mean > 0]
        if expressed.empty:
            continue
        marker_table = pd.DataFrame(
            {
                "gene": expressed.index,
                "type_mean": type_mean.loc[expressed.index],
                "specificity": specificity.loc[expressed.index],
            }
        )
        marker_table = marker_table[marker_table["gene"].isin(st_expr.columns)]
        marker_table = marker_table.sort_values(["specificity", "type_mean"], ascending=False)
        markers = marker_table["gene"].head(n_marker_genes).tolist()
        if len(markers) < 5:
            continue

        score = marker_scores(st_expr, markers)
        top15_n, top15_purity, top15_enrich, top15_largest = top_region_metrics(score, neigh, 0.85)
        top20_n, top20_purity, top20_enrich, top20_largest = top_region_metrics(score, neigh, 0.80)
        moran = moran_knn(score.to_numpy(), neigh)
        q50 = float(score.quantile(0.50))
        q85 = float(score.quantile(0.85))
        q95 = float(score.quantile(0.95))
        dynamic = q95 - q50
        specificity_median = float(marker_table["specificity"].head(n_marker_genes).median())
        candidate_score = (
            2.0 * max(moran, 0.0)
            + 0.35 * top15_enrich
            + 0.8 * top15_largest
            + 0.15 * np.log1p(max(dynamic, 0.0))
        )
        rows.append(
            {
                "pair_id": sample,
                "source_sample": sample,
                "storage_group": storage_group,
                "cell_type": cell_type,
                "sc_cells": n_type,
                "st_spots": int(len(st_expr)),
                "marker_genes": int(len(markers)),
                "marker_specificity_median": specificity_median,
                "marker_score_q50": q50,
                "marker_score_q85": q85,
                "marker_score_q95": q95,
                "marker_dynamic_q95_q50": dynamic,
                "moran_i": float(moran),
                "top15_spots": top15_n,
                "top15_neighbor_purity": top15_purity,
                "top15_neighbor_enrichment": top15_enrich,
                "top15_largest_component_fraction": top15_largest,
                "top20_spots": top20_n,
                "top20_neighbor_purity": top20_purity,
                "top20_neighbor_enrichment": top20_enrich,
                "top20_largest_component_fraction": top20_largest,
                "marker_genes_list": ";".join(markers),
                "candidate_score": float(candidate_score),
            }
        )
    return rows


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    all_rows: list[dict] = []
    for sample in args.samples:
        print(f"[scan] {sample}", flush=True)
        rows = scan_sample(root, args.storage_group, sample, args.min_sc_cells, args.n_marker_genes)
        print(f"[scan] {sample}: candidates={len(rows)}", flush=True)
        all_rows.extend(rows)

    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    out = pd.DataFrame(all_rows)
    out = out.sort_values(["pair_id", "candidate_score"], ascending=[True, False])
    out_path = out_dir / "lowres_masked_source_st_target_candidate_scan.csv"
    out.to_csv(out_path, index=False)
    print(f"[done] {out_path} rows={len(out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
