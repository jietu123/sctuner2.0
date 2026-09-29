from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualizations" / "c5_decomposition"
OUT.mkdir(parents=True, exist_ok=True)

STAGE3A = ROOT / "result" / "c5_stage3a_decomposition" / "c5_stage3a_exclusion_by_dataset.csv"
STAGE3B = ROOT / "result" / "c5_decomposition" / "c5_stage3b_decomposition_by_dataset.csv"
FINAL = ROOT / "result" / "c5_decomposition" / "c5_final_quantitative_source.csv"

COLORS = {
    "teal": "#168C8C",
    "teal_dark": "#116B70",
    "blue": "#557A95",
    "orange": "#D48245",
    "gray": "#7A8793",
    "light_gray": "#D5DCE1",
    "text": "#263238",
}

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.facecolor": "white",
    }
)


def finish(ax, grid_axis="x"):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#8F989F")
    ax.spines["bottom"].set_color("#8F989F")
    ax.tick_params(colors=COLORS["text"], width=0.8)
    ax.grid(axis=grid_axis, color="#E7EBEE", linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def compact_label(dataset, target):
    special = {
        "adult_mouse_kidney_real_profile_mask_endo": "Mouse kidney — Endo",
        "ffpe_mouse_brain_sagittal_real_profile_mask_microglia": "Mouse brain FFPE — Microglia",
        "human_breast_cancer_real_profile_mask_basal_cell": "Breast cancer — Basal cell",
        "human_breast_cancer_visium_ff_wta_real_profile_mask_macrophage": "Breast cancer FF — Macrophage",
        "human_breast_cancer_wta_120_real_profile_mask_endothelial_cell": "Breast cancer FFPE — Endothelial",
        "human_cervical_cancer_real_profile_mask_epithelial_cell": "Cervical cancer — Epithelial",
        "human_heart_ff_real_profile_mask_endothelial_cell": "Human heart — Endothelial",
        "human_intestine_cancer_real_profile_mask_endothelial_cell": "Intestine cancer — Endothelial",
        "human_lymph_node_real_profile_mask_b_cell": "Lymph node — B cell",
        "mouse_embryo_real_profile_mask_erythroid": "Mouse embryo — Erythroid",
        "highres_humanbreastcancerpatient1_profile_mask_monocytes_and_macrophages": "Breast cancer P1 — Mono./macro.",
        "highres_humancoloncancerpatient1_profile_mask_fibroblasts": "Colon cancer P1 — Fibroblasts",
        "highres_humanlungcancerpatient1_profile_mask_plasma_cells": "Lung cancer P1 — Plasma cells",
        "highres_humanmelanomapatient1_profile_mask_fibroblasts": "Melanoma P1 — Fibroblasts",
        "highres_humanmelanomapatient2_profile_mask_b_cells": "Melanoma P2 — B cells",
    }
    return special.get(dataset, f"{dataset} — {target}")


# Figure 1: Stage3A exclusion.
a = pd.read_csv(STAGE3A)
a["display_label"] = [compact_label(d, t) for d, t in zip(a["dataset"], a["target"])]
a["target_exclusion_percent"] = 100 * a["target_exclusion_fraction"].astype(float)
a["detection_summary"] = "15/15 = 100%"
a_source = a[
    [
        "dataset",
        "target",
        "resolution",
        "capacity",
        "target_total_cells",
        "target_dropped_cells",
        "target_exclusion_fraction",
        "target_exclusion_percent",
        "target_detected_missing",
        "off_target_exclusion",
        "off_target_details",
        "display_label",
        "detection_summary",
    ]
].copy()
a_source.to_csv(OUT / "c5_fig1_stage3a_exclusion_source.csv", index=False)

fig, ax = plt.subplots(figsize=(7.2, 5.2))
y = np.arange(len(a))
x = a["target_exclusion_percent"].to_numpy()
for yi, xi in zip(y, x):
    ax.hlines(yi, 0, xi, color="#CDD5DA", linewidth=1.0, zorder=1)
ax.scatter(x, y, s=30, color=COLORS["teal"], edgecolor="white", linewidth=0.6, zorder=3)
ax.set_yticks(y, a["display_label"])
ax.invert_yaxis()
ax.set_xlim(0, 105)
ax.set_xticks([0, 25, 50, 75, 100], ["0", "25", "50", "75", "100"])
ax.set_xlabel("Target reference cells excluded by Stage3A (%)")
ax.text(
    0.985,
    0.035,
    "Target detected in 15/15 experiments (100%)",
    transform=ax.transAxes,
    ha="right",
    va="bottom",
    color=COLORS["teal_dark"],
    fontsize=9,
    fontweight="bold",
)
idx = a.index[a["dataset"] == "highres_humanmelanomapatient2_profile_mask_b_cells"][0]
ax.annotate(
    "Additional off-target exclusion: NK cells",
    xy=(x[idx], idx),
    xytext=(64, idx - 0.55),
    ha="left",
    va="center",
    fontsize=7.5,
    color="#6A5151",
    arrowprops=dict(arrowstyle="-", color="#A08A8A", lw=0.8),
)
finish(ax, "x")
fig.subplots_adjust(left=0.38, right=0.98, bottom=0.12, top=0.98)
save(fig, "c5_fig1_stage3a_target_exclusion")


# Figure 2: standalone vs sequential Stage3B.
b = pd.read_csv(STAGE3B)
b["standalone_percent"] = 100 * b["standalone_fraction"].astype(float)
b["sequential_percent"] = 100 * b["sequential_fraction"].astype(float)
b["display_label"] = [compact_label(d, t) for d, t in zip(b["dataset"], b["target"])]
b_source = b[
    [
        "dataset",
        "target",
        "resolution",
        "capacity",
        "total_units",
        "standalone_withheld",
        "standalone_fraction",
        "standalone_percent",
        "sequential_withheld",
        "sequential_fraction",
        "sequential_percent",
        "delta_withheld",
        "delta_fraction",
        "mask_intersection",
        "mask_union",
        "jaccard",
        "display_label",
    ]
].copy()
b_source.to_csv(OUT / "c5_fig2_stage3b_decomposition_source.csv", index=False)

fig, ax = plt.subplots(figsize=(5.8, 4.7))
route_x = np.array([0.0, 1.0])
resolution_colors = {"spot-resolution": COLORS["blue"], "MERSCOPE": COLORS["orange"]}
for _, row in b.iterrows():
    vals = 100 * np.array([row["standalone_fraction"], row["sequential_fraction"]], dtype=float)
    color = resolution_colors[row["resolution"]]
    ax.plot(route_x, vals, color=color, linewidth=1.0, alpha=0.68, zorder=2)
    ax.scatter(route_x, vals, s=22, color=color, edgecolor="white", linewidth=0.45, zorder=3)
ax.set_xlim(-0.25, 1.25)
ax.set_xticks(route_x, ["Standalone Stage3B", "Sequential Stage3B"])
ax.set_ylabel("Spatial units withheld (%)")
ax.set_ylim(bottom=0)
ax.legend(
    handles=[
        Line2D([0], [0], marker="o", color=COLORS["blue"], lw=1.2, markersize=5, label="Spot-resolution (n=10)"),
        Line2D([0], [0], marker="o", color=COLORS["orange"], lw=1.2, markersize=5, label="MERSCOPE (n=5)"),
    ],
    loc="upper right",
    frameon=False,
    fontsize=8,
)
finish(ax, "y")
fig.tight_layout(pad=0.8)
save(fig, "c5_fig2_stage3b_standalone_vs_sequential")


# Figure 3: MERSCOPE four-route descriptive readout.
c = pd.read_csv(FINAL)
c = c[c["resolution"].eq("MERSCOPE")].copy()
c["display_label"] = [compact_label(d, t).replace(" — ", "\n") for d, t in zip(c["dataset"], c["target"])]
route_cols = {
    "Baseline": "baseline_peak_es",
    "Stage3A-only": "stage3a_only_peak_es",
    "Stage3B-only": "stage3b_only_peak_es",
    "Full": "full_peak_es",
}
c_long = c.melt(
    id_vars=["dataset", "target", "resolution", "capacity", "display_label", "mapping_quality_metric", "mapping_quality_historical_script"],
    value_vars=list(route_cols.values()),
    var_name="route_column",
    value_name="peak_es_low_masked_support_enrichment",
)
c_long["route"] = c_long["route_column"].map({v: k for k, v in route_cols.items()})
c_long["route"] = pd.Categorical(c_long["route"], categories=list(route_cols), ordered=True)
c_long = c_long.sort_values(["dataset", "route"]).drop(columns="route_column")
c_long.to_csv(OUT / "c5_fig3_merscope_four_route_source.csv", index=False)

fig, ax = plt.subplots(figsize=(7.2, 4.1))
base_x = np.arange(len(c))
offsets = np.array([-0.24, -0.08, 0.08, 0.24])
route_colors = ["#818B95", "#557A95", "#D49450", "#168C8C"]
route_markers = ["o", "s", "^", "D"]
for offset, (route, col), color, marker in zip(offsets, route_cols.items(), route_colors, route_markers):
    vals = c[col].astype(float).clip(lower=0).to_numpy()
    ax.scatter(
        base_x + offset,
        vals,
        s=34,
        marker=marker,
        color=color,
        edgecolor="white",
        linewidth=0.55,
        label=route,
        zorder=3,
    )
ax.set_xticks(base_x, c["display_label"])
ax.set_ylabel("Peak ES")
ax.set_ylim(bottom=0)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=4, frameon=False, fontsize=8, handletextpad=0.35, columnspacing=1.1)
ax.text(
    0.5,
    -0.24,
    "Reused high-resolution profile-mask descriptive readout; not ground-truth accuracy.",
    transform=ax.transAxes,
    ha="center",
    va="top",
    fontsize=7.7,
    color="#65717A",
    style="italic",
)
finish(ax, "y")
fig.subplots_adjust(left=0.10, right=0.99, bottom=0.30, top=0.91)
save(fig, "c5_fig3_merscope_four_route_peak_es")

print("SUCCESS")
