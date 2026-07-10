#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.spatial import cKDTree


SAMPLE = "cytospace_fig2d_tme_brca_her2_ffpe_profile_mask_t_cells_sc_missing_plasma_cells"
GROUP = "biological_application"
OUT_DIR = "visualizations/biological_application_brca_her2_joint_stage3ab"
DEFAULT_PAIRS = ["C3_CD81", "APP_CD74", "TIMP1_CD63", "VIM_CD44"]
IMMUNE_TYPES = [
    "T-cells",
    "CD4 T cells",
    "CD8 T cells",
    "B cells",
    "NK cells",
    "Monocytes and Macrophages",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot spatial hotspots for representative LR pairs in the joint BRCA HER2 Stage3A/Stage3B scene."
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--sample", default=SAMPLE)
    parser.add_argument("--storage_group", default=GROUP)
    parser.add_argument("--out_dir", default=OUT_DIR)
    parser.add_argument("--pairs", default=",".join(DEFAULT_PAIRS))
    parser.add_argument("--target_quantile", type=float, default=0.85)
    parser.add_argument("--neighbors", type=int, default=6)
    return parser.parse_args()


def hotspot_cmap() -> LinearSegmentedColormap:
    colors = [
        (0.00, "#D8D8D8"),
        (0.14, "#9E9E9E"),
        (0.28, "#6E6E6E"),
        (0.43, "#7A1F78"),
        (0.62, "#D13F5C"),
        (0.80, "#F98E5D"),
        (1.00, "#FFF2B2"),
    ]
    return LinearSegmentedColormap.from_list("gray_to_hotspot", colors)


def read_coords(path: Path) -> pd.DataFrame:
    coords = pd.read_csv(path)
    id_col = "spot_id" if "spot_id" in coords.columns else coords.columns[0]
    coords[id_col] = coords[id_col].astype(str)
    coords = coords.drop_duplicates(id_col).set_index(id_col)
    lower = {c.lower(): c for c in coords.columns}
    if "col" in lower and "row" in lower:
        x_col, y_col = lower["col"], lower["row"]
    elif "x" in lower and "y" in lower:
        x_col, y_col = lower["x"], lower["y"]
    else:
        numeric = [c for c in coords.columns if pd.api.types.is_numeric_dtype(coords[c])]
        if len(numeric) < 2:
            raise ValueError(f"Cannot infer coordinate columns from {path}")
        y_col, x_col = numeric[:2]
    result = coords[[x_col, y_col]].apply(pd.to_numeric, errors="coerce").dropna()
    result.columns = ["x", "y"]
    result["y_plot"] = -result["y"]
    return result


def read_expression(path: Path) -> pd.DataFrame:
    expr = pd.read_csv(path, index_col=0)
    expr.index = expr.index.astype(str)
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    arr = expr.to_numpy(dtype=float)
    lo = np.nanmin(arr, axis=0)
    hi = np.nanmax(arr, axis=0)
    denom = np.where((hi - lo) > 1e-12, hi - lo, 1.0)
    scaled = (arr - lo) / denom
    return pd.DataFrame(scaled, index=expr.index, columns=expr.columns)


def marker_percentile(expr: pd.DataFrame, marker_genes: list[str]) -> pd.Series:
    present = [gene for gene in marker_genes if gene in expr.columns]
    if not present:
        raise ValueError("No Plasma marker genes are present in ST expression.")
    return expr[present].mean(axis=1).rank(method="average", pct=True)


def fractional_path(root: Path, sample: str, suffix: str) -> Path:
    return root / "result" / sample / f"stage4_cytospace{suffix}" / "cytospace_output" / "fractional_abundances_by_spot.csv"


def read_fractional(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, index_col=0)
    frame.index = frame.index.astype(str)
    frame = frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    frame = frame.drop(columns=["Unknown_sc_only"], errors="ignore")
    totals = frame.sum(axis=1)
    nonzero = totals > 0
    frame.loc[nonzero] = frame.loc[nonzero].div(totals.loc[nonzero], axis=0)
    frame["Abstained"] = np.where(totals <= 1e-12, 1.0, 0.0)
    return frame


def align_fractionals(frames: list[pd.DataFrame]) -> list[pd.DataFrame]:
    columns = sorted(set().union(*(set(frame.columns) for frame in frames)))
    return [frame.reindex(columns=columns, fill_value=0.0) for frame in frames]


def abundance_sum(frame: pd.DataFrame, cell_types: list[str]) -> pd.Series:
    present = [cell_type for cell_type in cell_types if cell_type in frame.columns]
    if not present:
        return pd.Series(0.0, index=frame.index)
    return frame[present].sum(axis=1)


def neighborhood_mean(values: np.ndarray, coords: pd.DataFrame, neighbors: int) -> np.ndarray:
    xy = coords[["x", "y"]].to_numpy()
    tree = cKDTree(xy)
    k = min(neighbors + 1, len(coords))
    _, indices = tree.query(xy, k=k)
    if indices.ndim == 1:
        indices = indices[:, None]
    return values[indices].mean(axis=1)


def draw_panel(
    ax: plt.Axes,
    coords: pd.DataFrame,
    score: pd.Series,
    core_mask: pd.Series,
    title: str,
    norm: Normalize,
) -> None:
    size = 6.6
    ax.scatter(coords["x"], coords["y_plot"], s=size, c="#E4E0D8", marker="h", linewidths=0, alpha=0.75)
    sc = ax.scatter(
        coords["x"],
        coords["y_plot"],
        s=size,
        c=score.loc[coords.index],
        cmap=hotspot_cmap(),
        norm=norm,
        marker="h",
        linewidths=0,
        alpha=0.98,
    )
    core = coords.loc[core_mask.loc[coords.index]]
    ax.scatter(
        core["x"],
        core["y_plot"],
        s=size * 2.35,
        facecolors="none",
        edgecolors="#00CFE8",
        marker="h",
        linewidths=0.5,
        alpha=1.0,
        rasterized=False,
    )
    ax.set_title(title, fontsize=9.3, fontweight="bold", pad=4)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal", adjustable="box")
    for spine in ax.spines.values():
        spine.set_visible(False)
    return sc


def draw_core_panel(
    ax: plt.Axes,
    coords: pd.DataFrame,
    score: pd.Series,
    core_mask: pd.Series,
    title: str,
    norm: Normalize,
    cmap: LinearSegmentedColormap,
) -> None:
    size = 7.4
    core_bool = core_mask.loc[coords.index]
    background = coords.loc[~core_bool]
    core = coords.loc[core_bool]
    ax.scatter(
        background["x"],
        background["y_plot"],
        s=size,
        c="#E8E8E8",
        marker="h",
        linewidths=0,
        alpha=0.42,
    )
    sc = ax.scatter(
        core["x"],
        core["y_plot"],
        s=size * 1.65,
        c=score.loc[core.index],
        cmap=cmap,
        norm=norm,
        marker="h",
        linewidths=0,
        alpha=0.98,
    )
    ax.scatter(
        core["x"],
        core["y_plot"],
        s=size * 2.3,
        facecolors="none",
        edgecolors="#00CFE8",
        marker="h",
        linewidths=0.45,
        alpha=0.95,
    )
    ax.set_title(title, fontsize=8.9, fontweight="bold", pad=4)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal", adjustable="box")
    for spine in ax.spines.values():
        spine.set_visible(False)
    return sc


def main() -> None:
    args = parse_args()
    root = Path(args.project_root).resolve()
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    validation = json.loads((out_dir / "joint_stage3ab_validation_summary.json").read_text(encoding="utf-8"))

    pairs = [pair.strip() for pair in args.pairs.split(",") if pair.strip()]
    requested_single_pair = len(pairs) == 1
    requested_pair_name = pairs[0].lower() if requested_single_pair else ""
    out_stem = (
        f"joint_stage3ab_representative_lr_hotspot_{requested_pair_name}"
        if requested_single_pair
        else "joint_stage3ab_representative_lr_hotspots"
    )
    exported = root / "data" / "processed" / args.storage_group / args.sample / "stage1_preprocess" / "exported"
    coords = read_coords(exported / "st_coordinates.csv")
    expr = read_expression(exported / "st_expression_normalized.csv")
    baseline = read_fractional(fractional_path(root, args.sample, "_bioapp_joint_dropout_baseline"))
    svtuner = read_fractional(fractional_path(root, args.sample, "_bioapp_joint_svtuner"))
    baseline, svtuner = align_fractionals([baseline, svtuner])

    common = coords.index.intersection(expr.index).intersection(baseline.index).intersection(svtuner.index)
    coords = coords.loc[common]
    expr = expr.loc[common]
    baseline = baseline.loc[common]
    svtuner = svtuner.loc[common]
    plasma = marker_percentile(expr, validation["stage3b"]["marker_genes_used"])
    core_mask = plasma >= plasma.quantile(args.target_quantile)

    baseline_immune = abundance_sum(baseline, IMMUNE_TYPES)
    svtuner_immune = abundance_sum(svtuner, IMMUNE_TYPES)
    baseline_neighbor = pd.Series(neighborhood_mean(baseline_immune.to_numpy(), coords, args.neighbors), index=coords.index)
    svtuner_neighbor = pd.Series(neighborhood_mean(svtuner_immune.to_numpy(), coords, args.neighbors), index=coords.index)

    score_rows: list[dict[str, object]] = []
    pair_scores: dict[str, tuple[pd.Series, pd.Series]] = {}
    for pair in pairs:
        ligand, receptor = pair.split("_", 1)
        if ligand not in expr.columns or receptor not in expr.columns:
            continue
        expression_product = expr[ligand] * neighborhood_mean(expr[receptor].to_numpy(), coords, args.neighbors)
        baseline_score = expression_product * baseline_immune * baseline_neighbor
        svtuner_score = expression_product * svtuner_immune * svtuner_neighbor
        pair_scores[pair] = (baseline_score, svtuner_score)
        for method, values in [("CytoSPACE baseline", baseline_score), ("SVTuner Stage3A+B", svtuner_score)]:
            score_rows.append(
                {
                    "lr_pair": pair,
                    "method": method,
                    "mean_all": float(values.mean()),
                    "mean_plasma_core": float(values.loc[core_mask].mean()),
                    "top95_score": float(values.quantile(0.95)),
                }
            )
    if not pair_scores:
        raise ValueError("No requested LR pairs are present in the ST expression matrix.")

    score_summary = pd.DataFrame(score_rows)
    summary_name = (
        f"joint_stage3ab_lr_hotspot_summary_{requested_pair_name}.csv"
        if requested_single_pair
        else "joint_stage3ab_lr_hotspot_summary.csv"
    )
    score_summary.to_csv(out_dir / summary_name, index=False)
    hotspot_frame = pd.DataFrame(index=coords.index)
    hotspot_frame["x"] = coords["x"]
    hotspot_frame["y"] = coords["y"]
    hotspot_frame["plasma_core"] = core_mask
    for pair, (baseline_score, svtuner_score) in pair_scores.items():
        hotspot_frame[f"{pair}__cytospace_baseline"] = baseline_score
        hotspot_frame[f"{pair}__svtuner_stage3ab"] = svtuner_score
    spatial_name = (
        f"joint_stage3ab_lr_hotspot_spatial_scores_{requested_pair_name}.csv"
        if requested_single_pair
        else "joint_stage3ab_lr_hotspot_spatial_scores.csv"
    )
    hotspot_frame.to_csv(out_dir / spatial_name)

    vmax = max(float(max(score.max() for pair in pair_scores.values() for score in pair)), 1e-12)
    norm = Normalize(vmin=0, vmax=vmax)
    single_pair = len(pair_scores) == 1
    fig_size = (8.8, 4.15) if single_pair else (8.9, 2.55 * len(pair_scores))
    fig, axes = plt.subplots(len(pair_scores), 2, figsize=fig_size, dpi=230)
    if len(pair_scores) == 1:
        axes = np.array([axes])
    last_sc = None
    for row, (pair, (baseline_score, svtuner_score)) in enumerate(pair_scores.items()):
        baseline_mean = float(baseline_score.loc[core_mask].mean())
        svtuner_mean = float(svtuner_score.loc[core_mask].mean())
        reduction = baseline_mean - svtuner_mean
        reduction_pct = 100.0 * reduction / baseline_mean if baseline_mean > 1e-12 else 0.0
        last_sc = draw_panel(
            axes[row, 0],
            coords,
            baseline_score,
            core_mask,
            f"{pair}\nCytoSPACE core mean={baseline_mean:.3f}",
            norm,
        )
        last_sc = draw_panel(
            axes[row, 1],
            coords,
            svtuner_score,
            core_mask,
            f"{pair}\nSVTuner core mean={svtuner_mean:.3f}\nreduction={reduction_pct:.0f}%",
            norm,
        )
    fig.suptitle(
        "Representative LR-pair spatial hotspot around the Plasma-rich unsupported core"
        if single_pair
        else "Representative LR-pair spatial hotspots around the Plasma-rich unsupported core",
        fontsize=10.8,
        fontweight="bold",
        y=0.985,
    )
    fig.text(
        0.5,
        0.01,
        "cyan outline: same Plasma marker-defined top15 core in both columns; core mean values are reported in titles.",
        ha="center",
        fontsize=8.5,
    )
    cbar = fig.colorbar(last_sc, ax=axes.ravel().tolist(), fraction=0.024, pad=0.045)
    cbar.set_label("Spatial LR hotspot score", fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    fig.savefig(out_dir / f"{out_stem}.png", bbox_inches="tight", pad_inches=0.05)
    fig.savefig(out_dir / f"{out_stem}.pdf", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


if __name__ == "__main__":
    main()
