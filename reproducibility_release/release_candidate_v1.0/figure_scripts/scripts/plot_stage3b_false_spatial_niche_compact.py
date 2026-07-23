#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from plot_stage3b_false_spatial_niche_ncem_upgrade import (  # noqa: E402
    CONDITIONS,
    STATE_COLORS,
    abbreviate,
    prepare_experiment,
)

COMPACT_LABELS = {
    "Endothelial cells": "Endothelial",
    "Epithelial cells": "Epithelial",
    "Monocytes and Macrophages": "Mono/Macro",
    "Plasma cells": "Plasma",
}


def compact_label(label: str) -> str:
    return COMPACT_LABELS.get(label, abbreviate(label))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create the compact Stage3B false-niche figure containing only "
            "reference-relative neighborhood changes and focal-neighbor coupling."
        )
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument(
        "--metadata_csv",
        default=(
            "visualizations/stage3b_realdata_candidate_scan/spatial_9x2/"
            "stage3b_reference_dropout_spatial_stack_recommended_6x2_metadata.csv"
        ),
    )
    parser.add_argument(
        "--sample",
        default="cytospace_fig2d_tme_brca_her2_ffpe_sc_missing_plasma_cells",
    )
    parser.add_argument("--target_quantile", type=float, default=0.85)
    parser.add_argument("--neighbors", type=int, default=6)
    parser.add_argument("--states", type=int, default=4)
    parser.add_argument("--bootstrap", type=int, default=250)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--out_dir",
        default=(
            "visualizations/stage3b_false_spatial_niche/"
            "brca_her2_plasma_compact"
        ),
    )
    return parser.parse_args()


def active_rows(
    matrices: list[np.ndarray],
    labels: list[str],
    keep_abstained: bool,
) -> list[int]:
    combined = np.maximum.reduce([np.abs(matrix) for matrix in matrices])
    activity = combined.max(axis=1)
    keep = np.flatnonzero(activity > 1e-10).tolist()
    if keep_abstained and "Abstained" in labels:
        abstained = labels.index("Abstained")
        if abstained not in keep:
            keep.append(abstained)
    return sorted(keep)


def active_coupling_types(
    matrices: list[np.ndarray],
    labels: list[str],
) -> list[int]:
    biological = [
        index for index, label in enumerate(labels) if label != "Abstained"
    ]
    combined = np.maximum.reduce([np.abs(matrix) for matrix in matrices])
    activity = combined.sum(axis=0) + combined.sum(axis=1)
    maximum = float(activity[biological].max()) if biological else 0.0
    cutoff = maximum * 0.10
    return [index for index in biological if activity[index] >= cutoff]


def draw_composition_heatmap(
    axis: plt.Axes,
    matrix: np.ndarray,
    row_indices: list[int],
    labels: list[str],
    limit: float,
    show_ylabels: bool,
) -> None:
    shown = matrix[row_indices, :]
    sns.heatmap(
        shown,
        ax=axis,
        cmap="Reds",
        vmin=0,
        vmax=limit,
        cbar=False,
        linewidths=0.35,
        linecolor="white",
        xticklabels=[str(index + 1) for index in range(shown.shape[1])],
        yticklabels=(
            [abbreviate(labels[index]) for index in row_indices]
            if show_ylabels
            else False
        ),
    )
    axis.set_xlabel("Neighborhood state", fontsize=8.2)
    axis.set_ylabel("Neighbor type" if show_ylabels else "", fontsize=8.2)
    axis.tick_params(axis="x", labelrotation=0, labelsize=7.3, length=0)
    axis.tick_params(axis="y", labelrotation=0, labelsize=7.3, length=0)


def draw_coupling_heatmap(
    axis: plt.Axes,
    matrix: np.ndarray,
    indices: list[int],
    labels: list[str],
    limit: float,
    title: str,
    show_ylabels: bool,
) -> None:
    shown = matrix[np.ix_(indices, indices)]
    sns.heatmap(
        shown,
        ax=axis,
        cmap="bwr",
        center=0,
        vmin=-limit,
        vmax=limit,
        cbar=False,
        linewidths=0.35,
        linecolor="white",
        xticklabels=[compact_label(labels[index]) for index in indices],
        yticklabels=(
            [compact_label(labels[index]) for index in indices]
            if show_ylabels
            else False
        ),
        square=True,
    )
    axis.set_title(title, fontsize=10.0, fontweight="bold", pad=7)
    axis.set_xlabel("Neighboring assigned cell type", fontsize=8.2)
    axis.set_ylabel(
        "Focal assigned cell type" if show_ylabels else "",
        fontsize=8.2,
    )
    axis.tick_params(axis="x", labelrotation=52, labelsize=7.0, length=0)
    axis.tick_params(axis="y", labelrotation=0, labelsize=7.0, length=0)


def save_plot_data(
    experiment: dict[str, object],
    out_dir: Path,
    file_stem: str,
) -> None:
    labels = list(experiment["labels"])
    compositions = experiment["state_compositions"]
    rows: list[dict[str, object]] = []
    for condition in CONDITIONS:
        matrix = np.asarray(compositions[condition])
        for type_index, cell_type in enumerate(labels):
            for state in range(matrix.shape[1]):
                rows.append(
                    {
                        "condition": condition,
                        "neighbor_type": cell_type,
                        "neighborhood_state": state + 1,
                        "mean_neighbor_fraction": matrix[type_index, state],
                    }
                )
    pd.DataFrame(rows).to_csv(
        out_dir / f"{file_stem}_panel_a_neighborhood_composition.csv",
        index=False,
    )

    coupling_rows: list[dict[str, object]] = []
    for method in ("CytoSPACE dropout", "SVTuner"):
        matrix = np.asarray(experiment["differences"][method])
        for focal_index, focal_type in enumerate(labels):
            for neighbor_index, neighbor_type in enumerate(labels):
                coupling_rows.append(
                    {
                        "method": method,
                        "focal_assigned_type": focal_type,
                        "neighboring_assigned_type": neighbor_type,
                        "coupling_change_vs_ideal_abstention": matrix[
                            focal_index, neighbor_index
                        ],
                    }
                )
    pd.DataFrame(coupling_rows).to_csv(
        out_dir / f"{file_stem}_panel_b_coupling_changes.csv",
        index=False,
    )


def plot_compact(
    experiment: dict[str, object],
    out_dir: Path,
    file_stem: str,
) -> tuple[Path, Path]:
    sns.set_theme(style="white", context="paper")
    labels = list(experiment["labels"])
    compositions = experiment["state_compositions"]
    composition_matrices = {
        condition: np.asarray(compositions[condition])
        for condition in CONDITIONS
    }
    neighborhood_rows = active_rows(
        list(composition_matrices.values()),
        labels,
        keep_abstained=True,
    )
    neighborhood_limit = max(
        float(matrix[neighborhood_rows, :].max())
        for matrix in composition_matrices.values()
    )
    neighborhood_limit = max(neighborhood_limit, 1e-6)

    coupling_changes = {
        method: np.asarray(experiment["differences"][method])
        for method in ("CytoSPACE dropout", "SVTuner")
    }
    coupling_indices = active_coupling_types(
        list(coupling_changes.values()),
        labels,
    )
    coupling_limit = max(
        float(
            np.abs(
                matrix[np.ix_(coupling_indices, coupling_indices)]
            ).max()
        )
        for matrix in coupling_changes.values()
    )
    coupling_limit = max(coupling_limit, 1e-6)

    fig = plt.figure(figsize=(12.0, 9.2), dpi=350, facecolor="white")
    outer = fig.add_gridspec(
        2,
        1,
        height_ratios=[1.08, 0.92],
        hspace=0.46,
        left=0.13,
        right=0.91,
        top=0.90,
        bottom=0.15,
    )
    panel_a = outer[0].subgridspec(
        2,
        3,
        height_ratios=[0.88, 1.12],
        hspace=0.13,
        wspace=0.12,
    )
    panel_b = outer[1].subgridspec(1, 2, wspace=0.10)
    scatter_axes = [fig.add_subplot(panel_a[0, index]) for index in range(3)]
    heat_axes = [fig.add_subplot(panel_a[1, index]) for index in range(3)]
    axes_b = [fig.add_subplot(panel_b[0, index]) for index in range(2)]

    condition_titles = (
        "Full reference control",
        "CytoSPACE dropout",
        "SVTuner with Stage3B abstention",
    )
    state_labels = np.asarray(experiment["state_labels"])
    for column, (axis, condition, title) in enumerate(
        zip(scatter_axes, CONDITIONS, condition_titles)
    ):
        embedding = np.asarray(experiment["embeddings"][condition])
        for state in np.unique(state_labels):
            selected = state_labels == state
            axis.scatter(
                embedding[selected, 0],
                embedding[selected, 1],
                s=8.5,
                color=STATE_COLORS[int(state) % len(STATE_COLORS)],
                linewidths=0,
                alpha=0.88,
                rasterized=True,
            )
        axis.set_title(title, fontsize=9.2, fontweight="bold", pad=4)
        axis.set_xticks([])
        axis.set_yticks([])
        sns.despine(ax=axis, left=True, bottom=True)

    for index, (axis, condition) in enumerate(zip(heat_axes, CONDITIONS)):
        draw_composition_heatmap(
            axis,
            composition_matrices[condition],
            neighborhood_rows,
            labels,
            neighborhood_limit,
            show_ylabels=index == 0,
        )

    method_titles = ("CytoSPACE dropout", "SVTuner with Stage3B abstention")
    methods = ("CytoSPACE dropout", "SVTuner")
    for index, (axis, method, title) in enumerate(
        zip(axes_b, methods, method_titles)
    ):
        draw_coupling_heatmap(
            axis,
            coupling_changes[method],
            coupling_indices,
            labels,
            coupling_limit,
            title,
            show_ylabels=index == 0,
        )

    scalar_a = plt.cm.ScalarMappable(
        norm=plt.Normalize(0, neighborhood_limit),
        cmap="Reds",
    )
    colorbar_a = fig.colorbar(
        scalar_a,
        ax=[*scatter_axes, *heat_axes],
        fraction=0.025,
        pad=0.025,
        shrink=0.62,
    )
    colorbar_a.set_label(
        "Mean neighborhood fraction",
        fontsize=7.3,
    )
    colorbar_a.ax.tick_params(labelsize=6.6)

    scalar_b = plt.cm.ScalarMappable(
        norm=plt.Normalize(-coupling_limit, coupling_limit),
        cmap="bwr",
    )
    colorbar_b = fig.colorbar(
        scalar_b,
        ax=axes_b,
        fraction=0.025,
        pad=0.025,
        shrink=0.72,
    )
    colorbar_b.set_label(
        "Focal-neighbor coupling change\nrelative to ideal abstention",
        fontsize=7.3,
    )
    colorbar_b.ax.tick_params(labelsize=6.6)

    fig.canvas.draw()
    top_a = scatter_axes[0].get_position().y1
    top_b = axes_b[0].get_position().y1
    fig.text(
        0.035,
        top_a + 0.062,
        "A",
        fontsize=14,
        fontweight="bold",
        va="top",
    )
    fig.text(
        0.063,
        top_a + 0.062,
        "Neighborhood states and spatial composition",
        fontsize=11.0,
        fontweight="bold",
        va="top",
    )
    fig.text(
        0.035,
        top_b + 0.062,
        "B",
        fontsize=14,
        fontweight="bold",
        va="top",
    )
    fig.text(
        0.063,
        top_b + 0.062,
        "Focal-neighbor coupling deviations in the Plasma-cell marker-rich region",
        fontsize=11.0,
        fontweight="bold",
        va="top",
    )

    source_label = str(experiment["row"]["pair_id"]).replace("_", " ")
    target_label = str(experiment["row"]["target_type"])
    fig.text(
        0.13,
        0.035,
        (
            f"{source_label}; {target_label} removed from the SC reference. "
            "Panel A state colors are fixed from the full-reference control; similar "
            "embedding geometry indicates similar global neighborhood structure, "
            "not mapping correctness.\n"
            "Panel B uses a shared scale: positive values indicate excess focal-neighbor "
            "coupling and negative values indicate depletion. Evaluation markers are "
            "not used by Stage3B."
        ),
        fontsize=7.2,
        color="#444444",
    )

    png = out_dir / f"{file_stem}_compact_neighborhood_coupling.png"
    pdf = out_dir / f"{file_stem}_compact_neighborhood_coupling.pdf"
    fig.savefig(png, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    plt.close(fig)
    return png, pdf


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(root / args.metadata_csv)
    experiment = prepare_experiment(
        root=root,
        metadata=metadata,
        sample=args.sample,
        target_quantile=args.target_quantile,
        neighbors=args.neighbors,
        states=args.states,
        repeats=args.bootstrap,
        seed=args.seed,
    )
    file_stem = str(experiment["row"]["sample"])
    save_plot_data(experiment, out_dir, file_stem)
    png, _ = plot_compact(experiment, out_dir, file_stem)
    print(f"[done] {png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
