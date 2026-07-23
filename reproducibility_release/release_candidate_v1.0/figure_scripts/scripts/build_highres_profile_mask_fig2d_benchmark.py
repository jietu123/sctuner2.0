from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.build_highres_fig2c_expression_enrichment import (
    READOUT_TO_GENESET,
    _evaluate,
    _read_gene_sets,
    _slug,
)

EPS = 1.0e-12

METHODS = ["CytoSPACE", "SVTuner + CytoSPACE", "Tangram", "CellTrek"]
METHOD_COLORS = {
    "CytoSPACE": "#ef6a5b",
    "SVTuner + CytoSPACE": "#2f9a8f",
    "Tangram": "#d9d9d9",
    "CellTrek": "#d9d9d9",
}


def _scan_candidates(root: Path, nperm: int, seed: int, min_gene_overlap: int) -> pd.DataFrame:
    gene_sets = _read_gene_sets(root)
    manifest = pd.read_csv(root / "visualizations/highres_profile_mask_mapping/highres_profile_mask_mapping_manifest.csv")
    records: list[dict[str, object]] = []
    for _, row in manifest.iterrows():
        for readout, gene_set_names in READOUT_TO_GENESET.items():
            for gene_set_name in gene_set_names:
                if gene_set_name not in gene_sets:
                    continue
                rec = _evaluate(root, row, readout, gene_set_name, gene_sets[gene_set_name], nperm, seed)
                if rec is None:
                    continue
                rec = {k: v for k, v in rec.items() if k != "_curves"}
                if int(rec["gene_set_overlap"]) < min_gene_overlap:
                    continue
                records.append(rec)
    if not records:
        raise RuntimeError("No valid high-resolution profile-mask Fig2D candidates found.")
    out = pd.DataFrame(records)
    out["candidate_label"] = (
        out["raw_sample"].astype(str)
        + " | "
        + out["readout_cell_type"].astype(str)
        + " | "
        + out["gene_set"].astype(str).str.replace("EcoTyper_", "", regex=False)
    )
    return out


def _include_all_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    included = candidates.copy()
    included["inclusion_note"] = "all_valid_candidates"
    return included.sort_values(
        ["raw_sample", "readout_cell_type", "gene_set"],
        kind="mergesort",
    ).reset_index(drop=True)


def _to_long(included: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in included.iterrows():
        for method, prefix in [
            ("CytoSPACE", "baseline_highres"),
            ("SVTuner + CytoSPACE", "route2_highres"),
        ]:
            rows.append(
                {
                    "raw_sample": row["raw_sample"],
                    "profile_mask_sample": row["profile_mask_sample"],
                    "masked_target_type": row["masked_target_type"],
                    "readout_cell_type": row["readout_cell_type"],
                    "gene_set": row["gene_set"],
                    "gene_set_overlap": row["gene_set_overlap"],
                    "candidate_label": row["candidate_label"],
                    "inclusion_note": row["inclusion_note"],
                    "method": method,
                    "nes": row[f"{prefix}_NES"],
                    "pval": row[f"{prefix}_pval"],
                    "n_mapped_cells": row[f"{prefix}_n_mapped_cells"],
                    "peak_rank": row[f"{prefix}_peak_rank"],
                }
            )
    return pd.DataFrame(rows)


def _mean_nearest5(query_xy: np.ndarray, ref_xy: np.ndarray) -> np.ndarray:
    diff = query_xy[:, None, :] - ref_xy[None, :, :]
    dist = np.sqrt(np.sum(diff * diff, axis=2))
    k = min(5, dist.shape[1])
    return np.partition(dist, kth=k - 1, axis=1)[:, :k].mean(axis=1)


def _running_es(stats: pd.Series, gene_set: set[str]) -> tuple[float, int]:
    stats = stats.replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False)
    hits = stats.index.to_series().isin(gene_set).to_numpy()
    n_hits = int(hits.sum())
    n_miss = int(len(stats) - n_hits)
    if n_hits == 0 or n_miss == 0:
        raise ValueError("invalid gene-set overlap")
    weights = np.abs(stats.to_numpy(dtype=float))
    hit_norm = float(weights[hits].sum())
    if hit_norm <= EPS:
        inc = np.where(hits, 1.0 / n_hits, -1.0 / n_miss)
    else:
        inc = np.where(hits, weights / hit_norm, -1.0 / n_miss)
    running = np.cumsum(inc)
    peak_rank = int(np.nanargmax(running)) + 1
    return float(running[peak_rank - 1]), peak_rank


def _null_nes(stats: pd.Series, gene_set_size: int, observed_es: float, nperm: int, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    genes = stats.replace([np.inf, -np.inf], np.nan).dropna().index.to_numpy()
    if gene_set_size >= len(genes):
        return float("nan"), float("nan")
    null = np.empty(nperm, dtype=float)
    for i in range(nperm):
        selected = set(rng.choice(genes, size=gene_set_size, replace=False).tolist())
        null[i], _ = _running_es(stats, selected)
    nes = observed_es / (float(np.mean(np.abs(null))) + EPS)
    pval = (float(np.sum(null >= observed_es)) + 1.0) / (float(nperm) + 1.0)
    return float(nes), float(pval)


def _load_generic_method_cells(root: Path, sample: str, method_dir: str, readout: str) -> pd.DataFrame:
    assign_path = root / "result" / sample / "stage4_mapping" / method_dir / "cell_assignment.csv"
    if not assign_path.exists():
        raise FileNotFoundError(assign_path)
    assign = pd.read_csv(assign_path)
    assign["cell_id"] = assign["cell_id"].astype(str)
    assign["assigned_spot"] = assign["assigned_spot"].astype(str)
    assign["cell_type"] = assign["cell_type"].astype(str)

    # Tangram/CellTrek assign every scRNA-seq cell to a spatial unit, whereas
    # CytoSPACE high-resolution runs use one mapped cell per spatial unit. For
    # a fair Fig.2d-style enrichment benchmark, collapse generic methods to one
    # representative assignment per spatial unit before selecting readout cells.
    if "assignment_score" in assign.columns:
        assign["_representative_rank"] = -pd.to_numeric(assign["assignment_score"], errors="coerce").fillna(-np.inf)
    elif "assignment_distance" in assign.columns:
        assign["_representative_rank"] = pd.to_numeric(assign["assignment_distance"], errors="coerce").fillna(np.inf)
    else:
        assign["_representative_rank"] = np.arange(len(assign), dtype=float)
    assign = (
        assign.sort_values(["assigned_spot", "_representative_rank", "cell_id"], kind="mergesort")
        .drop_duplicates("assigned_spot", keep="first")
        .copy()
    )
    assign = assign[assign["cell_type"].eq(readout)].copy()
    if assign.empty:
        return pd.DataFrame(columns=["OriginalCID", "CellType", "SpotID", "row", "col"])

    coords = pd.read_csv(root / "data" / "processed" / sample / "stage1_preprocess" / "exported" / "st_coordinates.csv")
    coords["spot_id"] = coords["spot_id"].astype(str)
    coord_cols = ["row", "col"] if {"row", "col"}.issubset(coords.columns) else ["spatial_row", "spatial_col"]
    coords = coords[["spot_id", coord_cols[0], coord_cols[1]]].copy()
    coords.columns = ["SpotID", "row", "col"]

    out = assign.rename(columns={"cell_id": "OriginalCID", "assigned_spot": "SpotID", "cell_type": "CellType"})
    out["OriginalCID"] = out["OriginalCID"].astype(str)
    out["SpotID"] = out["SpotID"].astype(str)
    out["CellType"] = out["CellType"].astype(str)
    out = out.merge(coords, on="SpotID", how="left")
    out = out.dropna(subset=["row", "col"])
    return out[["OriginalCID", "CellType", "SpotID", "row", "col"]].copy()


def _evaluate_generic_method(
    root: Path,
    row: pd.Series,
    readout: str,
    gene_set: list[str],
    method_dir: str,
    nperm: int,
    seed: int,
) -> dict[str, object] | None:
    sample = str(row["profile_mask_sample"])
    export = root / "data" / "processed" / sample / "stage1_preprocess" / "exported"
    st_meta = pd.read_csv(export / "st_metadata.csv")
    coords = pd.read_csv(export / "st_coordinates.csv")
    spatial = coords.merge(st_meta[["spot_id", "cell_type"]], on="spot_id", how="left")
    anchor = spatial[spatial["cell_type"].astype(str).eq("Epithelial cells")]
    if len(anchor) < 20:
        return None

    expr = pd.read_csv(export / "sc_expression_normalized.csv").set_index("cell_id")
    available_genes = pd.Index(expr.columns.astype(str))
    overlap = list(available_genes.intersection(pd.Index(gene_set)))
    if len(overlap) < 8:
        return None

    loc = _load_generic_method_cells(root, sample, method_dir, readout)
    # After collapsing Tangram/CellTrek to one representative assignment per
    # spatial unit, some rare readout states naturally have fewer mapped cells.
    # Keep a low explicit floor and record n_mapped_cells in source values.
    if len(loc) < 10:
        return None
    loc = loc[loc["OriginalCID"].astype(str).isin(expr.index)].copy()
    if len(loc) < 10:
        return None

    dist = _mean_nearest5(
        loc[["row", "col"]].to_numpy(dtype=float),
        anchor[["row", "col"]].to_numpy(dtype=float),
    )
    close = dist <= float(np.median(dist))
    far = ~close
    close_ids = loc.loc[close, "OriginalCID"].astype(str)
    far_ids = loc.loc[far, "OriginalCID"].astype(str)
    close_mean = np.log1p(expr.loc[close_ids]).mean(axis=0)
    far_mean = np.log1p(expr.loc[far_ids]).mean(axis=0)
    stats = close_mean - far_mean
    es, peak_rank = _running_es(stats, set(overlap))
    nes, pval = _null_nes(stats, len(overlap), es, nperm=nperm, seed=seed)
    return {
        "ES": es,
        "NES": nes,
        "pval": pval,
        "peak_rank": peak_rank,
        "n_mapped_cells": int(len(loc)),
    }


def _add_generic_methods(root: Path, included: pd.DataFrame, long_df: pd.DataFrame, nperm: int, seed: int) -> pd.DataFrame:
    gene_sets = _read_gene_sets(root)
    rows = long_df.to_dict("records")
    for idx, row in included.reset_index(drop=True).iterrows():
        gene_set_name = str(row["gene_set"])
        if gene_set_name not in gene_sets:
            continue
        for offset, (method_label, method_dir) in enumerate([("Tangram", "tangram_marker"), ("CellTrek", "celltrek")]):
            rec = _evaluate_generic_method(
                root,
                row,
                str(row["readout_cell_type"]),
                gene_sets[gene_set_name],
                method_dir,
                nperm=nperm,
                seed=seed + 1000 + idx * 10 + offset,
            )
            if rec is None:
                continue
            rows.append(
                {
                    "raw_sample": row["raw_sample"],
                    "profile_mask_sample": row["profile_mask_sample"],
                    "masked_target_type": row["masked_target_type"],
                    "readout_cell_type": row["readout_cell_type"],
                    "gene_set": row["gene_set"],
                    "gene_set_overlap": row["gene_set_overlap"],
                    "candidate_label": row["candidate_label"],
                    "inclusion_note": row["inclusion_note"],
                    "method": method_label,
                    "nes": rec["NES"],
                    "pval": rec["pval"],
                    "n_mapped_cells": rec["n_mapped_cells"],
                    "peak_rank": rec["peak_rank"],
                }
            )
    return pd.DataFrame(rows)


def _plot(long_df: pd.DataFrame, out_dir: Path, out_prefix: str) -> Path:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    values = [long_df.loc[long_df["method"].eq(method), "nes"].to_numpy(float) for method in METHODS]
    fig, ax = plt.subplots(figsize=(4.15, 3.35), dpi=320)
    positions = np.arange(1, len(METHODS) + 1)
    bp = ax.boxplot(
        values,
        positions=positions,
        widths=0.48,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#222222", "linewidth": 1.05},
        boxprops={"color": "#777777", "linewidth": 0.9},
        whiskerprops={"color": "#555555", "linewidth": 0.9},
        capprops={"color": "#555555", "linewidth": 0.9},
    )
    for patch, method in zip(bp["boxes"], METHODS):
        patch.set_facecolor(METHOD_COLORS[method])
        patch.set_alpha(0.95)
        patch.set_zorder(2)

    rng = np.random.default_rng(7)
    for i, method in enumerate(METHODS, start=1):
        y = long_df.loc[long_df["method"].eq(method), "nes"].to_numpy(float)
        ax.scatter(
            np.full(len(y), i) + rng.normal(0, 0.035, len(y)),
            y,
            s=12,
            color="#222222",
            alpha=0.86,
            linewidth=0,
            zorder=5,
        )

    flat = np.concatenate([v for v in values if len(v)])
    ymin = min(-0.15, float(np.nanmin(flat)) - 0.25)
    ymax = max(2.2, float(np.nanmax(flat)) + 0.35)
    ax.axhline(0, color="#555555", linewidth=0.8, linestyle=(0, (2.2, 1.8)))
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(0.45, len(METHODS) + 1.05)
    ax.set_ylabel("Normalized\nenrichment score", fontsize=8.0)
    ax.set_xticks(positions)
    ax.set_xticklabels(
        ["CytoSPACE", "SVTuner +\nCytoSPACE", "Tangram", "CellTrek"],
        rotation=45,
        ha="right",
        rotation_mode="anchor",
        fontsize=7.0,
    )
    ax.tick_params(axis="y", labelsize=7.0, width=0.9, length=3)
    ax.tick_params(axis="x", width=0.9, length=2.6, pad=1)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_linewidth(0.95)
    ax.spines["bottom"].set_linewidth(0.95)
    ax.annotate(
        "",
        xy=(4.72, 1.85),
        xytext=(4.72, -0.08),
        arrowprops=dict(arrowstyle="<->", color="#8a8a8a", lw=1.25),
        annotation_clip=False,
    )
    ax.text(4.82, 1.55, "Enriched\nclose to\ntumor", ha="left", va="center", fontsize=6.1)
    ax.text(4.82, 0.05, "Enriched\nfar from\ntumor", ha="left", va="center", fontsize=6.1)

    means = long_df.groupby("method")["nes"].mean()
    delta = float(means["SVTuner + CytoSPACE"] - means["CytoSPACE"])
    fig.text(
        0.18,
        0.965,
        "State enrichment benchmark under cell-level profile masking",
        ha="left",
        va="top",
        fontsize=7.0,
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.17, right=0.82, top=0.84, bottom=0.31)
    png = out_dir / f"{out_prefix}.png"
    pdf = out_dir / f"{out_prefix}.pdf"
    svg = out_dir / f"{out_prefix}.svg"
    fig.savefig(png, dpi=320)
    fig.savefig(pdf)
    fig.savefig(svg)
    plt.close(fig)
    return png


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build a Fig2D-style NES benchmark on five high-resolution profile-mask datasets."
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument(
        "--out_dir",
        default="visualizations/highres_profile_mask_fig2d/all_candidates",
    )
    parser.add_argument("--out_prefix", default="fig2d_highres_profile_mask_benchmark")
    parser.add_argument("--nperm", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=19)
    parser.add_argument("--min_gene_overlap", type=int, default=8)
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    candidates = _scan_candidates(root, args.nperm, args.seed, args.min_gene_overlap)
    included = _include_all_candidates(candidates)
    long_df = _to_long(included)
    long_df = _add_generic_methods(root, included, long_df, args.nperm, args.seed)

    candidates.to_csv(out_dir / f"{args.out_prefix}_candidate_scan.csv", index=False)
    included.to_csv(out_dir / f"{args.out_prefix}_all_profiles.csv", index=False)
    long_df.to_csv(out_dir / f"{args.out_prefix}_source_values.csv", index=False)
    png = _plot(long_df, out_dir, args.out_prefix)
    baseline_nes = included["baseline_highres_NES"].to_numpy(float)
    route2_nes = included["route2_highres_NES"].to_numpy(float)
    delta_nes = route2_nes - baseline_nes
    summary = {
        "experiment_design": (
            "Fig2D-style normalized enrichment benchmark on five cell-level high-resolution profile-mask datasets. "
            "Profile-mask targets were detected by Stage3 and route2 used the auto-detected missing types; no forced whitelist is used. "
            "All valid dataset/readout/gene-set combinations are included without filtering or ranking by mapping outcome. "
            "For each readout state, mapped cells are split by distance to epithelial cells and genes are ranked by close-minus-far expression. "
            "Tangram and CellTrek are collapsed to one representative assignment per spatial unit before enrichment scoring to match the high-resolution CytoSPACE scale."
        ),
        "n_candidates": int(len(candidates)),
        "n_included_profiles": int(len(included)),
        "n_datasets": int(included["raw_sample"].nunique()),
        "mean_nes": long_df.groupby("method")["nes"].mean().to_dict(),
        "rows_per_method": long_df.groupby("method").size().to_dict(),
        "route2_improved_profiles": int(
            (route2_nes > baseline_nes).sum()
        ),
        "route2_non_improved_profiles": int((route2_nes <= baseline_nes).sum()),
        "mean_delta_nes": float(delta_nes.mean()),
        "median_delta_nes": float(np.median(delta_nes)),
        "all_profiles_csv": str((out_dir / f"{args.out_prefix}_all_profiles.csv").relative_to(root)).replace("\\", "/"),
        "source_values_csv": str((out_dir / f"{args.out_prefix}_source_values.csv").relative_to(root)).replace("\\", "/"),
        "candidate_scan_csv": str((out_dir / f"{args.out_prefix}_candidate_scan.csv").relative_to(root)).replace("\\", "/"),
        "png": str(png.relative_to(root)).replace("\\", "/"),
    }
    (out_dir / f"{args.out_prefix}_manifest.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print("[OK] wrote:", png)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
