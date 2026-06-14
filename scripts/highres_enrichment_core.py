from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap


EPS = 1.0e-12
CYTOSPACE_METHODS = [
    ("CytoSPACE", "baseline_highres"),
    ("SVTuner + CytoSPACE", "route2_highres"),
]
BENCHMARK_METHODS = ["CytoSPACE", "SVTuner + CytoSPACE", "Tangram", "CellTrek"]
BENCHMARK_COLORS = {
    "CytoSPACE": "#ef6a5b",
    "SVTuner + CytoSPACE": "#2f9a8f",
    "Tangram": "#d9d9d9",
    "CellTrek": "#d9d9d9",
}


def slug(value: str) -> str:
    return "".join(ch.lower() if ch.isalnum() else "_" for ch in str(value)).strip("_")


def read_gene_sets(root: Path) -> dict[str, list[str]]:
    path = root / "data/raw/cytospace_fig2c_melanoma/prepared/gene_sets/table_s6_gene_sets_long.csv"
    df = pd.read_csv(path)
    return {
        str(name): sub.sort_values("rank")["gene"].astype(str).tolist()
        for name, sub in df.groupby("gene_set", sort=False)
    }


def _mean_nearest5(query_xy: np.ndarray, ref_xy: np.ndarray) -> np.ndarray:
    diff = query_xy[:, None, :] - ref_xy[None, :, :]
    dist = np.sqrt(np.sum(diff * diff, axis=2))
    k = min(5, dist.shape[1])
    return np.partition(dist, kth=k - 1, axis=1)[:, :k].mean(axis=1)


def _running_es(stats: pd.Series, gene_set: set[str]) -> tuple[float, int, pd.DataFrame]:
    stats = stats.replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False)
    hits = stats.index.to_series().isin(gene_set).to_numpy()
    n_hits = int(hits.sum())
    n_miss = int(len(stats) - n_hits)
    if n_hits == 0 or n_miss == 0:
        raise ValueError("invalid gene-set overlap")
    weights = np.abs(stats.to_numpy(dtype=np.float64))
    hit_norm = float(weights[hits].sum())
    if hit_norm <= EPS:
        increments = np.where(hits, 1.0 / n_hits, -1.0 / n_miss)
    else:
        increments = np.where(hits, weights / hit_norm, -1.0 / n_miss)
    running = np.cumsum(increments)
    peak_rank = int(np.nanargmax(running)) + 1
    curve = pd.DataFrame(
        {
            "rank": np.arange(1, len(stats) + 1),
            "gene": stats.index.to_numpy(),
            "stat": stats.to_numpy(dtype=float),
            "hit": hits,
            "running_es": running,
        }
    )
    return float(running[peak_rank - 1]), peak_rank, curve


def _null_nes(
    stats: pd.Series,
    gene_set_size: int,
    observed_es: float,
    nperm: int,
    seed: int,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    genes = stats.replace([np.inf, -np.inf], np.nan).dropna().index.to_numpy()
    if gene_set_size >= len(genes):
        return float("nan"), float("nan")
    null = np.empty(nperm, dtype=float)
    for index in range(nperm):
        selected = set(rng.choice(genes, size=gene_set_size, replace=False).tolist())
        null[index], _, _ = _running_es(stats, selected)
    nes = observed_es / (float(np.mean(np.abs(null))) + EPS)
    pval = (float(np.sum(null >= observed_es)) + 1.0) / (float(nperm) + 1.0)
    return float(nes), float(pval)


def _load_cytospace_cells(root: Path, sample: str, suffix: str, readout: str) -> pd.DataFrame:
    path = (
        root
        / "result"
        / sample
        / f"stage4_cytospace_{suffix}"
        / "cytospace_output"
        / "assigned_locations.csv"
    )
    locations = pd.read_csv(path)
    return locations[locations["CellType"].astype(str).eq(readout)].copy()


def evaluate_cytospace_pair(
    root: Path,
    mapping_row: pd.Series,
    readout: str,
    gene_set_name: str,
    gene_set: list[str],
    nperm: int,
    seed: int,
) -> dict[str, object] | None:
    sample = str(mapping_row["profile_mask_sample"])
    export = root / "data" / "processed" / sample / "stage1_preprocess" / "exported"
    st_meta = pd.read_csv(export / "st_metadata.csv")
    coords = pd.read_csv(export / "st_coordinates.csv")
    spatial = coords.merge(st_meta[["spot_id", "cell_type"]], on="spot_id", how="left")
    anchor = spatial[spatial["cell_type"].astype(str).eq("Epithelial cells")]
    if len(anchor) < 20:
        return None

    expression = pd.read_csv(export / "sc_expression_normalized.csv").set_index("cell_id")
    overlap = list(pd.Index(expression.columns.astype(str)).intersection(pd.Index(gene_set)))
    if len(overlap) < 8:
        return None

    result: dict[str, object] = {
        "raw_sample": str(mapping_row["raw_sample"]),
        "profile_mask_sample": sample,
        "masked_target_type": str(mapping_row["masked_target_type"]),
        "anchor_cell_type": "Epithelial cells",
        "readout_cell_type": readout,
        "gene_set": gene_set_name,
        "gene_set_overlap": int(len(overlap)),
    }
    curves: dict[str, pd.DataFrame] = {}
    for offset, (_, suffix) in enumerate(CYTOSPACE_METHODS):
        locations = _load_cytospace_cells(root, sample, suffix, readout)
        if len(locations) < 40:
            return None
        locations = locations[locations["OriginalCID"].astype(str).isin(expression.index)].copy()
        if len(locations) < 40:
            return None
        distance = _mean_nearest5(
            locations[["row", "col"]].to_numpy(dtype=float),
            anchor[["row", "col"]].to_numpy(dtype=float),
        )
        close = distance <= float(np.median(distance))
        close_ids = locations.loc[close, "OriginalCID"].astype(str)
        far_ids = locations.loc[~close, "OriginalCID"].astype(str)
        stats = np.log1p(expression.loc[close_ids]).mean(axis=0) - np.log1p(
            expression.loc[far_ids]
        ).mean(axis=0)
        es, peak_rank, curve = _running_es(stats, set(overlap))
        nes, pval = _null_nes(stats, len(overlap), es, nperm, seed + offset * 100)
        result[f"{suffix}_ES"] = es
        result[f"{suffix}_NES"] = nes
        result[f"{suffix}_pval"] = pval
        result[f"{suffix}_peak_rank"] = peak_rank
        result[f"{suffix}_n_mapped_cells"] = int(len(locations))
        result[f"{suffix}_n_close"] = int(close.sum())
        result[f"{suffix}_n_far"] = int((~close).sum())
        curves[suffix] = curve
    result["delta_NES"] = float(result["route2_highres_NES"]) - float(
        result["baseline_highres_NES"]
    )
    result["route2_peak_frac"] = float(result["route2_highres_peak_rank"]) / float(
        len(curves["route2_highres"])
    )
    result["baseline_peak_frac"] = float(result["baseline_highres_peak_rank"]) / float(
        len(curves["baseline_highres"])
    )
    result["_curves"] = curves
    return result


def _format_p(value: float) -> str:
    return "0.001" if value < 0.0015 else f"{value:.3f}"


def _smooth_curve(values: np.ndarray, window: int = 61) -> np.ndarray:
    if len(values) < 7:
        return values
    width = min(window, len(values) - (1 - len(values) % 2))
    if width % 2 == 0:
        width -= 1
    if width < 5:
        return values
    padding = width // 2
    return np.convolve(
        np.pad(values, (padding, padding), mode="edge"),
        np.ones(width) / width,
        mode="valid",
    )


def _gradient_line(ax, x: np.ndarray, y: np.ndarray) -> None:
    cmap = LinearSegmentedColormap.from_list(
        "rank_gradient",
        ["#f03b20", "#ffb000", "#78b95b", "#2b8cbe"],
    )
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    collection = LineCollection(
        segments,
        cmap=cmap,
        norm=plt.Normalize(x.min(), x.max()),
    )
    collection.set_array(x)
    collection.set_linewidth(1.8)
    collection.set_zorder(4)
    ax.add_collection(collection)


def _draw_enrichment_panel(
    ax,
    curve: pd.DataFrame,
    label: str,
    nes: float,
    pval: float,
) -> None:
    curve = curve.sort_values("rank")
    x = curve["rank"].to_numpy(dtype=float)
    x = x / float(x.max())
    y = _smooth_curve(curve["running_es"].to_numpy(dtype=float))
    hit_indices = np.where(curve["hit"].astype(bool).to_numpy())[0]
    cmap = LinearSegmentedColormap.from_list(
        "rank_gradient_hits",
        ["#f03b20", "#ffb000", "#78b95b", "#2b8cbe"],
    )
    hit_colors = cmap(x[hit_indices])
    _gradient_line(ax, x, y)
    ax.vlines(x[hit_indices], 0, y[hit_indices], color=hit_colors, linewidth=0.72, alpha=0.95, zorder=3)
    ax.scatter(x[hit_indices], y[hit_indices], s=8.2, color=hit_colors, edgecolor="none", alpha=0.95, zorder=5)
    ax.axhline(0, color="#222222", linewidth=0.9)
    y_low = min(0.0, float(y.min()))
    y_high = max(0.0, float(y.max()))
    span = max(y_high - y_low, 0.2)
    ax.set_ylim(y_low - span * 0.10, y_high + span * 0.18)
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.tick_params(axis="y", labelsize=6.4, width=0.9, length=2.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_linewidth(0.9)
    ax.text(0.02, 1.03, label, transform=ax.transAxes, ha="left", va="bottom", fontsize=7.7, weight="bold")
    ax.text(0.70, 0.76, f"NES = {nes:.2f}\nP = {_format_p(pval)}", transform=ax.transAxes, ha="left", va="center", fontsize=6.4)


def plot_fig2c(result: dict[str, object], out_dir: Path, out_prefix: str) -> tuple[Path, Path]:
    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
    fig = plt.figure(figsize=(3.05, 3.85), dpi=320)
    grid = fig.add_gridspec(3, 1, height_ratios=[1, 1, 0.24], hspace=0.28, left=0.20, right=0.98, top=0.90, bottom=0.12)
    baseline_ax = fig.add_subplot(grid[0, 0])
    route2_ax = fig.add_subplot(grid[1, 0], sharex=baseline_ax)
    distance_ax = fig.add_subplot(grid[2, 0], sharex=baseline_ax)
    _draw_enrichment_panel(baseline_ax, result["_curves"]["baseline_highres"], "CytoSPACE", float(result["baseline_highres_NES"]), float(result["baseline_highres_pval"]))
    _draw_enrichment_panel(route2_ax, result["_curves"]["route2_highres"], "SVTuner + CytoSPACE", float(result["route2_highres_NES"]), float(result["route2_highres_pval"]))
    fig.text(0.08, 0.57, "Running enrichment score", rotation=90, va="center", ha="center", fontsize=8.2)
    title_readout = str(result["readout_cell_type"]).replace("Monocytes and Macrophages", "Mono/Mac")
    title_geneset = (
        str(result["gene_set"])
        .replace("EcoTyper_Monocytes_and_Macrophages_", "EcoTyper Mono/Mac ")
        .replace("EcoTyper_CD8_T_cells_", "EcoTyper CD8 ")
        .replace("EcoTyper_CD4_T_cells_", "EcoTyper CD4 ")
        .replace("_", " ")
    )
    fig.text(0.58, 0.975, f"{title_readout} {title_geneset}", ha="center", va="top", fontsize=7.8, linespacing=0.9)
    distance_ax.axis("off")
    distance_ax.annotate("", xy=(0.86, 0.54), xytext=(0.20, 0.54), arrowprops=dict(arrowstyle="->", lw=1.8, color="#8a8a8a"))
    distance_ax.text(0.01, 0.54, "Increasing distance\nfrom epithelial cells", ha="left", va="center", fontsize=6.0)
    png = out_dir / f"{out_prefix}.png"
    pdf = out_dir / f"{out_prefix}.pdf"
    svg = out_dir / f"{out_prefix}.svg"
    fig.savefig(png, dpi=320)
    fig.savefig(pdf)
    fig.savefig(svg)
    plt.close(fig)
    return png, pdf


def to_long(selected: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, row in selected.iterrows():
        for method, prefix in CYTOSPACE_METHODS:
            rows.append(
                {
                    "raw_sample": row["raw_sample"],
                    "profile_mask_sample": row["profile_mask_sample"],
                    "masked_target_type": row["masked_target_type"],
                    "readout_cell_type": row["readout_cell_type"],
                    "gene_set": row["gene_set"],
                    "gene_set_overlap": row["gene_set_overlap"],
                    "candidate_label": row["candidate_label"],
                    "method": method,
                    "nes": row[f"{prefix}_NES"],
                    "pval": row[f"{prefix}_pval"],
                    "n_mapped_cells": row[f"{prefix}_n_mapped_cells"],
                    "peak_rank": row[f"{prefix}_peak_rank"],
                }
            )
    return pd.DataFrame(rows)


def _load_generic_cells(root: Path, sample: str, method_dir: str, readout: str) -> pd.DataFrame:
    assignment = pd.read_csv(root / "result" / sample / "stage4_mapping" / method_dir / "cell_assignment.csv")
    assignment["cell_id"] = assignment["cell_id"].astype(str)
    assignment["assigned_spot"] = assignment["assigned_spot"].astype(str)
    assignment["cell_type"] = assignment["cell_type"].astype(str)
    if "assignment_score" in assignment.columns:
        assignment["_representative_rank"] = -pd.to_numeric(assignment["assignment_score"], errors="coerce").fillna(-np.inf)
    elif "assignment_distance" in assignment.columns:
        assignment["_representative_rank"] = pd.to_numeric(assignment["assignment_distance"], errors="coerce").fillna(np.inf)
    else:
        assignment["_representative_rank"] = np.arange(len(assignment), dtype=float)
    assignment = assignment.sort_values(
        ["assigned_spot", "_representative_rank", "cell_id"],
        kind="mergesort",
    ).drop_duplicates("assigned_spot", keep="first")
    assignment = assignment[assignment["cell_type"].eq(readout)].copy()
    if assignment.empty:
        return pd.DataFrame(columns=["OriginalCID", "CellType", "SpotID", "row", "col"])
    coords = pd.read_csv(root / "data" / "processed" / sample / "stage1_preprocess" / "exported" / "st_coordinates.csv")
    coords["spot_id"] = coords["spot_id"].astype(str)
    coord_columns = ["row", "col"] if {"row", "col"}.issubset(coords.columns) else ["spatial_row", "spatial_col"]
    coords = coords[["spot_id", *coord_columns]].copy()
    coords.columns = ["SpotID", "row", "col"]
    output = assignment.rename(columns={"cell_id": "OriginalCID", "assigned_spot": "SpotID", "cell_type": "CellType"})
    output["SpotID"] = output["SpotID"].astype(str)
    return output.merge(coords, on="SpotID", how="left").dropna(subset=["row", "col"])[
        ["OriginalCID", "CellType", "SpotID", "row", "col"]
    ]


def _evaluate_generic(
    root: Path,
    row: pd.Series,
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
    expression = pd.read_csv(export / "sc_expression_normalized.csv").set_index("cell_id")
    overlap = list(pd.Index(expression.columns.astype(str)).intersection(pd.Index(gene_set)))
    locations = _load_generic_cells(root, sample, method_dir, str(row["readout_cell_type"]))
    locations = locations[locations["OriginalCID"].astype(str).isin(expression.index)].copy()
    if len(anchor) < 20 or len(overlap) < 8 or len(locations) < 10:
        return None
    distance = _mean_nearest5(
        locations[["row", "col"]].to_numpy(dtype=float),
        anchor[["row", "col"]].to_numpy(dtype=float),
    )
    close = distance <= float(np.median(distance))
    stats = np.log1p(expression.loc[locations.loc[close, "OriginalCID"].astype(str)]).mean(axis=0) - np.log1p(
        expression.loc[locations.loc[~close, "OriginalCID"].astype(str)]
    ).mean(axis=0)
    es, peak_rank, _ = _running_es(stats, set(overlap))
    nes, pval = _null_nes(stats, len(overlap), es, nperm, seed)
    return {"NES": nes, "pval": pval, "peak_rank": peak_rank, "n_mapped_cells": int(len(locations))}


def add_generic_methods(
    root: Path,
    selected: pd.DataFrame,
    long_df: pd.DataFrame,
    nperm: int,
    seed: int,
) -> pd.DataFrame:
    gene_sets = read_gene_sets(root)
    rows = long_df.to_dict("records")
    for index, row in selected.reset_index(drop=True).iterrows():
        gene_set = gene_sets[str(row["gene_set"])]
        for offset, (method, method_dir) in enumerate([("Tangram", "tangram_marker"), ("CellTrek", "celltrek")]):
            result = _evaluate_generic(root, row, gene_set, method_dir, nperm, seed + 1000 + index * 10 + offset)
            if result is None:
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
                    "method": method,
                    "nes": result["NES"],
                    "pval": result["pval"],
                    "n_mapped_cells": result["n_mapped_cells"],
                    "peak_rank": result["peak_rank"],
                }
            )
    return pd.DataFrame(rows)


def plot_fig2d(long_df: pd.DataFrame, out_dir: Path, out_prefix: str) -> Path:
    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
    values = [long_df.loc[long_df["method"].eq(method), "nes"].to_numpy(float) for method in BENCHMARK_METHODS]
    fig, ax = plt.subplots(figsize=(4.15, 3.35), dpi=320)
    positions = np.arange(1, len(BENCHMARK_METHODS) + 1)
    boxplot = ax.boxplot(
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
    for patch, method in zip(boxplot["boxes"], BENCHMARK_METHODS):
        patch.set_facecolor(BENCHMARK_COLORS[method])
        patch.set_alpha(0.95)
        patch.set_zorder(2)
    rng = np.random.default_rng(7)
    for position, method in enumerate(BENCHMARK_METHODS, start=1):
        y = long_df.loc[long_df["method"].eq(method), "nes"].to_numpy(float)
        ax.scatter(np.full(len(y), position) + rng.normal(0, 0.035, len(y)), y, s=12, color="#222222", alpha=0.86, linewidth=0, zorder=5)
    flat = np.concatenate([value for value in values if len(value)])
    ax.axhline(0, color="#555555", linewidth=0.8, linestyle=(0, (2.2, 1.8)))
    ax.set_ylim(min(-0.15, float(np.nanmin(flat)) - 0.25), max(2.2, float(np.nanmax(flat)) + 0.35))
    ax.set_xlim(0.45, len(BENCHMARK_METHODS) + 1.05)
    ax.set_ylabel("Normalized\nenrichment score", fontsize=8.0)
    ax.set_xticks(positions)
    ax.set_xticklabels(["CytoSPACE", "SVTuner +\nCytoSPACE", "Tangram", "CellTrek"], rotation=45, ha="right", rotation_mode="anchor", fontsize=7.0)
    ax.tick_params(axis="y", labelsize=7.0, width=0.9, length=3)
    ax.tick_params(axis="x", width=0.9, length=2.6, pad=1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.95)
    ax.spines["bottom"].set_linewidth(0.95)
    ax.annotate("", xy=(4.72, 1.85), xytext=(4.72, -0.08), arrowprops=dict(arrowstyle="<->", color="#8a8a8a", lw=1.25), annotation_clip=False)
    ax.text(4.82, 1.55, "Enriched\nclose to\ntumor", ha="left", va="center", fontsize=6.1)
    ax.text(4.82, 0.05, "Enriched\nfar from\ntumor", ha="left", va="center", fontsize=6.1)
    fig.text(0.18, 0.965, "State enrichment benchmark under cell-level profile masking", ha="left", va="top", fontsize=7.0, fontweight="bold")
    fig.subplots_adjust(left=0.17, right=0.82, top=0.84, bottom=0.31)
    png = out_dir / f"{out_prefix}.png"
    pdf = out_dir / f"{out_prefix}.pdf"
    svg = out_dir / f"{out_prefix}.svg"
    fig.savefig(png, dpi=320)
    fig.savefig(pdf)
    fig.savefig(svg)
    plt.close(fig)
    return png
