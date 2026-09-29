#!/usr/bin/env python
"""Reproduce the frozen C1 paired dot-box figure and export editable SVG."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "visualizations" / "method_comparison" / "c1_multiseed"
SEED_SCORES = ROOT / "result" / "c1_seed_scores.csv"
HISTORICAL = (
    ROOT
    / "visualizations"
    / "method_comparison"
    / "composite_no_noise"
    / "composition_recovery_7mapping_methods_composite_no_noise_abstention_aware_scenario.csv"
)

SCENARIOS = [
    "BC-1", "BC-2", "BC-3",
    "Lung-1", "Lung-2", "Lung-3",
    "Brain-1", "Brain-2", "Brain-3",
]

HISTORICAL_SAMPLES = {
    "real_brca7_endothelial_marker_control_sc_missing_endothelial_cells": "BC-1",
    "real_brca7_endothelial_marker_missing_epithelial_cells_sc_missing_endothelial_cells": "BC-2",
    "real_brca7_endothelial_marker_missing_epithelial_cells_pcs_sc_missing_endothelial_cells": "BC-3",
    "human_lung_5loc_fine9_clustered_sim_sc_missing_b_cell": "Lung-1",
    "human_lung_5loc_fine9_clustered_sim_missing_at2_sc_missing_b_cell": "Lung-2",
    "human_lung_5loc_fine9_clustered_sim_missing_at2_fibroblast_sc_missing_b_cell": "Lung-3",
    "mouse_brain_refined7_balanced_clustered_sim_sc_missing_ext_l56": "Brain-1",
    "mouse_brain_refined7_balanced_clustered_sim_missing_micro_fill_inh_pvalb_sc_missing_ext_l56": "Brain-2",
    "mouse_brain_refined7_balanced_clustered_sim_missing_micro_oligo_2_fill_inh_pvalb_sc_missing_ext_l56": "Brain-3",
}

METHODS = ["CytoSPACE", "SVTuner"]
COLORS = {"CytoSPACE": "#6E838D", "SVTuner": "#20B5B1"}
FILLS = {"CytoSPACE": "#AEBCC2", "SVTuner": "#8EDAD7"}


def load_complete_scores() -> pd.DataFrame:
    repeated = pd.read_csv(SEED_SCORES)
    repeated = repeated.loc[repeated["simulation_seed"].between(43, 51)].copy()

    historical = pd.read_csv(HISTORICAL)
    historical = historical.loc[historical["method"].isin(METHODS)].copy()
    historical["scenario_id"] = historical["sample"].map(HISTORICAL_SAMPLES)
    historical = historical.dropna(subset=["scenario_id"])
    seed42 = historical.pivot(index="scenario_id", columns="method", values="composition_recovery").reset_index()
    seed42 = seed42.rename(columns={"CytoSPACE": "cytospace_score", "SVTuner": "svtuner_score"})
    seed42["simulation_seed"] = 42
    seed42["sample"] = "historical_seed42"

    combined = pd.concat(
        [repeated, seed42[["scenario_id", "simulation_seed", "sample", "cytospace_score", "svtuner_score"]]],
        ignore_index=True,
    )
    combined = combined.loc[combined["scenario_id"].isin(SCENARIOS)].copy()
    counts = combined.groupby("scenario_id")["simulation_seed"].nunique()
    if not counts.reindex(SCENARIOS).eq(10).all() or len(combined) != 90:
        raise ValueError("C1 frozen dataset is not 9 scenarios x 10 seeds")
    return combined


def build_long_table(combined: pd.DataFrame) -> pd.DataFrame:
    long = combined.melt(
        id_vars=["scenario_id", "simulation_seed"],
        value_vars=["cytospace_score", "svtuner_score"],
        var_name="method_key",
        value_name="score",
    )
    long["method"] = long["method_key"].map(
        {"cytospace_score": "CytoSPACE", "svtuner_score": "SVTuner"}
    )
    overall = (
        long.groupby(["simulation_seed", "method"], as_index=False)["score"]
        .mean()
        .assign(scenario_id="Overall")
    )
    return pd.concat(
        [long[["scenario_id", "simulation_seed", "method", "score"]], overall],
        ignore_index=True,
    )


def main() -> int:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.edgecolor": "#9EA8AD",
            "axes.linewidth": 0.75,
        }
    )
    scores = build_long_table(load_complete_scores())
    order = [*SCENARIOS, "Overall"]
    x = np.array([0.0, 1.0, 2.0, 3.25, 4.25, 5.25, 6.50, 7.50, 8.50, 9.82])
    offsets = {"CytoSPACE": -0.115, "SVTuner": 0.115}

    fig, ax = plt.subplots(figsize=(7.068, 3.236), dpi=500)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    for method in METHODS:
        grouped = [
            scores.loc[(scores["scenario_id"] == scenario) & (scores["method"] == method), "score"].to_numpy()
            for scenario in order
        ]
        positions = x + offsets[method]
        bp = ax.boxplot(
            grouped,
            positions=positions,
            widths=0.24,
            patch_artist=True,
            showfliers=False,
            manage_ticks=False,
            boxprops={"facecolor": FILLS[method], "edgecolor": COLORS[method], "linewidth": 1.15, "alpha": 0.42},
            medianprops={"color": COLORS[method], "linewidth": 2.2},
            whiskerprops={"color": COLORS[method], "linewidth": 1.15},
            capprops={"color": COLORS[method], "linewidth": 1.15},
        )
        for patch in bp["boxes"]:
            patch.set_alpha(0.42)
        for xpos, values in zip(positions, grouped):
            ax.scatter(
                np.full(len(values), xpos),
                values,
                s=8.0,
                color=COLORS[method],
                edgecolor="white",
                linewidth=0.28,
                alpha=0.95,
                zorder=4,
            )

    for separator in [(x[2] + x[3]) / 2, (x[5] + x[6]) / 2, (x[8] + x[9]) / 2]:
        ax.axvline(separator, color="#E5E9EA", linewidth=0.75, zorder=0)

    group_centers = {
        "BRCA": x[:3].mean(),
        "Lung": x[3:6].mean(),
        "Brain": x[6:9].mean(),
        "Overall": x[9],
    }
    for label, center in group_centers.items():
        ax.text(center, 1.035, label, transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=11.2, fontweight="bold", color="#394247")

    handles = [
        Line2D([0], [0], marker="o", linestyle="none", markersize=6.2,
               markerfacecolor=COLORS[m], markeredgecolor="none", label=m)
        for m in METHODS
    ]
    ax.legend(
        handles=handles,
        loc="lower left",
        bbox_to_anchor=(-0.015, 1.125),
        frameon=False,
        ncol=2,
        columnspacing=1.15,
        handletextpad=0.45,
        fontsize=10.0,
    )

    ax.set_xlim(-0.55, 10.35)
    ax.set_ylim(0.185, 0.97)
    ax.set_yticks(np.arange(0.2, 1.0, 0.1))
    ax.set_xticks(x)
    ax.set_xticklabels(order, fontsize=8.8)
    ax.set_ylabel("Abstention-aware recovery score", fontsize=10.2, labelpad=8)
    ax.tick_params(axis="y", labelsize=8.3, width=0.7, length=3)
    ax.tick_params(axis="x", width=0.7, length=3, pad=5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)

    fig.subplots_adjust(left=0.095, right=0.992, bottom=0.16, top=0.79)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    svg = OUT_DIR / "c1_paired_difference_dotwhisker.svg"
    fig.savefig(svg, facecolor="white")
    plt.close(fig)
    print(f"[done] {svg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
