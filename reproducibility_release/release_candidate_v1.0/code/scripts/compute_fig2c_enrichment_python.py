from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io, sparse


EPS = 1.0e-12


def _read_assigned_expression(cytospace_dir: Path) -> tuple[sparse.csc_matrix, list[str], list[str]]:
    expr_dir = cytospace_dir / "assigned_expression"
    mat = io.mmread(expr_dir / "matrix.mtx").tocsc().astype(np.float32)
    genes = pd.read_csv(expr_dir / "genes.tsv", sep="\t", header=None)[0].astype(str).tolist()
    barcodes = pd.read_csv(expr_dir / "barcodes.tsv", sep="\t", header=None)[0].astype(str).tolist()
    return mat, genes, barcodes


def _normalize_log(mat: sparse.csc_matrix) -> sparse.csc_matrix:
    totals = np.asarray(mat.sum(axis=0)).ravel().astype(np.float32)
    totals[totals <= 0] = 1.0
    norm = mat @ sparse.diags((10000.0 / totals).astype(np.float32), format="csc")
    norm.data = np.log1p(norm.data)
    return norm


def _mean_nearest5(target_xy: np.ndarray, tumor_xy: np.ndarray) -> np.ndarray:
    diff = target_xy[:, None, :] - tumor_xy[None, :, :]
    dist = np.sqrt(np.sum(diff * diff, axis=2))
    k = min(5, dist.shape[1])
    return np.partition(dist, kth=k - 1, axis=1)[:, :k].mean(axis=1)


def _running_es(stats: pd.Series, gene_set: set[str]) -> tuple[float, int, pd.DataFrame]:
    stats = stats.replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False)
    hits = stats.index.to_series().isin(gene_set).to_numpy()
    n_hits = int(hits.sum())
    n_miss = int(len(stats) - n_hits)
    if n_hits == 0 or n_miss == 0:
        raise ValueError("invalid gene set overlap for running ES")
    weights = np.abs(stats.to_numpy(dtype=np.float64))
    hit_norm = float(weights[hits].sum())
    if hit_norm <= EPS:
        increments = np.where(hits, 1.0 / n_hits, -1.0 / n_miss)
    else:
        increments = np.where(hits, weights / hit_norm, -1.0 / n_miss)
    running = np.cumsum(increments)
    peak_idx = int(np.nanargmax(running))
    curve = pd.DataFrame(
        {
            "rank": np.arange(1, len(stats) + 1),
            "gene": stats.index.to_numpy(),
            "stat": stats.to_numpy(dtype=np.float64),
            "hit": hits,
            "running_es": running,
        }
    )
    return float(running[peak_idx]), int(peak_idx + 1), curve


def _null_nes(stats: pd.Series, n_genes: int, observed_es: float, nperm: int, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    genes = stats.replace([np.inf, -np.inf], np.nan).dropna().index.to_numpy()
    if n_genes >= len(genes):
        return float("nan"), float("nan")
    null = []
    for _ in range(nperm):
        selected = set(rng.choice(genes, size=n_genes, replace=False).tolist())
        es, _, _ = _running_es(stats, selected)
        null.append(es)
    null_arr = np.asarray(null, dtype=np.float64)
    denom = float(np.mean(np.abs(null_arr))) + EPS
    nes = observed_es / denom
    pval = (float(np.sum(null_arr >= observed_es)) + 1.0) / (float(nperm) + 1.0)
    return nes, pval


def main() -> int:
    parser = argparse.ArgumentParser(description="Compute Fig.2c/Fig.2d-style enrichment without R/fgsea.")
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--cytospace_dir", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--tumor_labels", default="Melanoma|Melanoma cells")
    parser.add_argument("--slide_label", default="sample")
    parser.add_argument("--gene_set", default="data/raw/cytospace_fig2c_melanoma/prepared/gene_sets/tcell_exhaustion_zheng_cell2017.txt")
    parser.add_argument("--gene_set_name", default="Exhaustion")
    parser.add_argument("--cell_types", default="CD4 T cells|CD8 T cells")
    parser.add_argument("--nperm", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    cytospace_dir = Path(args.cytospace_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    mat, genes, barcodes = _read_assigned_expression(cytospace_dir)
    norm = _normalize_log(mat)
    loc = pd.read_csv(cytospace_dir / "assigned_locations.csv")
    loc["UniqueCID"] = loc["UniqueCID"].astype(str)
    loc = loc.set_index("UniqueCID").reindex(pd.Index(barcodes))
    loc["CellType"] = loc["CellType"].astype(str)
    coord_cols = ["X", "Y"] if {"X", "Y"}.issubset(loc.columns) else ["row", "col"]

    gene_set_path = Path(args.gene_set)
    if not gene_set_path.is_absolute():
        gene_set_path = root / gene_set_path
    signature = set(g for g in gene_set_path.read_text(encoding="utf-8-sig").splitlines() if g.strip())
    gene_index = pd.Index(genes)
    signature = set(gene_index.intersection(pd.Index(list(signature))).tolist())
    if len(signature) < 5:
        raise ValueError(f"Too few signature genes available: {len(signature)}")

    tumor_labels = set(x for x in args.tumor_labels.split("|") if x)
    tumor = loc[loc["CellType"].isin(tumor_labels)]
    if len(tumor) < 5:
        raise ValueError(f"Need at least five tumor mapped cells; found {len(tumor)}")
    tumor_xy = tumor[coord_cols].to_numpy(dtype=np.float64)

    rows = []
    for offset, cell_type in enumerate([x for x in args.cell_types.split("|") if x]):
        target = loc[loc["CellType"].eq(cell_type)].copy()
        if len(target) < 20:
            raise ValueError(f"Too few mapped cells for {cell_type}: {len(target)}")
        target_xy = target[coord_cols].to_numpy(dtype=np.float64)
        dist = _mean_nearest5(target_xy, tumor_xy)
        cutoff = float(np.median(dist))
        close_mask = dist <= cutoff
        far_mask = ~close_mask
        target_positions = np.where(loc.index.to_series().isin(target.index).to_numpy())[0]
        close_cols = target_positions[close_mask]
        far_cols = target_positions[far_mask]
        close_mean = np.asarray(norm[:, close_cols].mean(axis=1)).ravel()
        far_mean = np.asarray(norm[:, far_cols].mean(axis=1)).ravel()
        stats = pd.Series(close_mean - far_mean, index=genes)
        es, peak_rank, curve = _running_es(stats, signature)
        nes, pval = _null_nes(stats, len(signature), es, args.nperm, args.seed + offset)

        prefix = "".join(ch.lower() if ch.isalnum() else "_" for ch in cell_type).strip("_")
        pd.DataFrame(
            {
                "UniqueCID": target.index.astype(str),
                "OriginalCID": target["OriginalCID"].astype(str).to_numpy(),
                "CellType": target["CellType"].astype(str).to_numpy(),
                "SpotID": target["SpotID"].astype(str).to_numpy(),
                "X": target[coord_cols[0]].to_numpy(),
                "Y": target[coord_cols[1]].to_numpy(),
                "mean_distance_to_nearest5_melanoma": dist,
                "distance_group": np.where(close_mask, "close", "far"),
            }
        ).to_csv(out_dir / f"{prefix}_distance_groups.csv", index=False)
        curve["CellType"] = cell_type
        curve.to_csv(out_dir / f"{prefix}_enrichment_curve.csv", index=False)
        rows.append(
            {
                "method": "CytoSPACE",
                "slide": args.slide_label,
                "cell_type": cell_type,
                "pathway": args.gene_set_name,
                "ES": es,
                "NES": nes,
                "pval": pval,
                "padj": pval,
                "peak_rank": peak_rank,
                "n_mapped_cells": int(len(target)),
                "n_close": int(close_mask.sum()),
                "n_far": int(far_mask.sum()),
                "n_exhaustion_genes_used": int(len(signature)),
                "backend": "python_permutation_gsea",
                "nperm": int(args.nperm),
            }
        )

    out = pd.DataFrame(rows)
    out.to_csv(out_dir / "fig2c_official_enrichment_summary.csv", index=False)
    print(f"[OK] wrote: {out_dir / 'fig2c_official_enrichment_summary.csv'}")
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
