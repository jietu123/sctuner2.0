from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Patch, Polygon


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualizations" / "c5_decomposition"
OUT.mkdir(parents=True, exist_ok=True)
STAGE3A = ROOT / "result" / "c5_stage3a_decomposition" / "c5_stage3a_exclusion_by_dataset.csv"
STAGE3B = ROOT / "result" / "c5_decomposition" / "c5_stage3b_decomposition_by_dataset.csv"

TEXT = "#263238"
SPOT = "#557A95"
SPOT_LIGHT = "#E7EEF3"
MER = "#168C8C"
MER_LIGHT = "#E3F1F0"
FLAG = "#B95C62"

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 9, "axes.labelsize": 10,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "savefig.facecolor": "white",
})

LABELS = {
    "adult_mouse_kidney_real_profile_mask_endo": ("Mouse kidney", "Endo"),
    "ffpe_mouse_brain_sagittal_real_profile_mask_microglia": ("Mouse brain FFPE", "Microglia"),
    "human_breast_cancer_real_profile_mask_basal_cell": ("Breast cancer", "Basal cell"),
    "human_breast_cancer_visium_ff_wta_real_profile_mask_macrophage": ("Breast cancer FF", "Macrophage"),
    "human_breast_cancer_wta_120_real_profile_mask_endothelial_cell": ("Breast cancer FFPE", "Endothelial cell"),
    "human_cervical_cancer_real_profile_mask_epithelial_cell": ("Cervical cancer", "Epithelial cell"),
    "human_heart_ff_real_profile_mask_endothelial_cell": ("Human heart", "Endothelial cell"),
    "human_intestine_cancer_real_profile_mask_endothelial_cell": ("Intestine cancer", "Endothelial cell"),
    "human_lymph_node_real_profile_mask_b_cell": ("Lymph node", "B cell"),
    "mouse_embryo_real_profile_mask_erythroid": ("Mouse embryo", "Erythroid"),
    "highres_humanbreastcancerpatient1_profile_mask_monocytes_and_macrophages": ("Breast cancer P1", "Mono./macro."),
    "highres_humancoloncancerpatient1_profile_mask_fibroblasts": ("Colon cancer P1", "Fibroblasts"),
    "highres_humanlungcancerpatient1_profile_mask_plasma_cells": ("Lung cancer P1", "Plasma cells"),
    "highres_humanmelanomapatient1_profile_mask_fibroblasts": ("Melanoma P1", "Fibroblasts"),
    "highres_humanmelanomapatient2_profile_mask_b_cells": ("Melanoma P2", "B cells"),
}


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


# Figure 1: 3 x 5 experiment tile matrix.
a = pd.read_csv(STAGE3A)
a["short_dataset"] = a["dataset"].map(lambda x: LABELS[x][0])
a["short_target"] = a["dataset"].map(lambda x: LABELS[x][1])
a["target_exclusion_percent"] = 100 * a["target_exclusion_fraction"].astype(float)
a["tile_group"] = np.where(a["resolution"].eq("MERSCOPE"), "MERSCOPE (capacity=1)", "Spot-resolution")
a["tile_row"] = list(np.repeat([0, 1], 5)) + [2] * 5
a["tile_col"] = list(np.tile(np.arange(5), 3))
a[["dataset", "target", "short_dataset", "short_target", "resolution", "capacity",
   "target_total_cells", "target_dropped_cells", "target_exclusion_fraction",
   "target_exclusion_percent", "target_detected_missing", "off_target_exclusion",
   "off_target_details", "tile_group", "tile_row", "tile_col"]].to_csv(
       OUT / "c5_fig1_stage3a_experiment_tile_matrix_source.csv", index=False)

fig, ax = plt.subplots(figsize=(10.0, 5.4))
ax.set_xlim(0, 5)
ax.set_ylim(0, 3.55)
ax.axis("off")
tile_w, tile_h = 0.89, 0.76
for _, row in a.iterrows():
    col, grid_row = int(row["tile_col"]), int(row["tile_row"])
    x, y = col + 0.055, 2.48 - grid_row * 0.91
    is_mer = row["resolution"] == "MERSCOPE"
    flagged = str(row["off_target_exclusion"]).lower() == "true"
    face, edge = (MER_LIGHT, MER) if is_mer else (SPOT_LIGHT, SPOT)
    if flagged:
        edge = FLAG
    tile = FancyBboxPatch((x, y), tile_w, tile_h,
                          boxstyle="round,pad=0.012,rounding_size=0.035",
                          facecolor=face, edgecolor=edge, linewidth=2.2 if flagged else 1.05)
    ax.add_patch(tile)
    ax.text(x + 0.07, y + 0.59, row["short_dataset"], ha="left", va="center",
            fontsize=8.2, color=TEXT, fontweight="bold")
    ax.text(x + 0.07, y + 0.40, f"Masked: {row['short_target']}", ha="left", va="center",
            fontsize=7.5, color=TEXT)
    ax.text(x + 0.07, y + 0.18, "100% target exclusion", ha="left", va="center",
            fontsize=8.0, color=edge, fontweight="bold")
    if flagged:
        corner = Polygon([[x + tile_w - 0.19, y + tile_h], [x + tile_w, y + tile_h],
                          [x + tile_w, y + tile_h - 0.19]], closed=True,
                         facecolor=FLAG, edgecolor=FLAG)
        ax.add_patch(corner)
        ax.text(x + 0.50, y + 0.055, "Off-target: NK cells (n=66)", ha="center", va="bottom",
                fontsize=6.8, color=FLAG, fontweight="bold")

ax.text(0.055, 3.35, "Spot-resolution experiments (n=10)", ha="left", va="center",
        fontsize=9.2, color=SPOT, fontweight="bold")
ax.text(0.055, 1.53, "MERSCOPE experiments (n=5; capacity=1)", ha="left", va="center",
        fontsize=9.2, color=MER, fontweight="bold")
ax.text(4.945, 3.35, "15/15 masked targets detected and fully excluded", ha="right", va="center",
        fontsize=10, color="#116B70", fontweight="bold")
ax.legend(handles=[Patch(facecolor=SPOT_LIGHT, edgecolor=SPOT, label="Target-only exclusion"),
                   Patch(facecolor=MER_LIGHT, edgecolor=FLAG, linewidth=2, label="Additional off-target exclusion")],
          loc="lower left", bbox_to_anchor=(0.0, -0.025), frameon=False, ncol=2,
          fontsize=8, handlelength=1.5, columnspacing=1.5)
fig.subplots_adjust(left=0.015, right=0.995, top=0.98, bottom=0.07)
save(fig, "c5_fig1_stage3a_experiment_tile_matrix")


# Figure 2: diverging delta bars.
b = pd.read_csv(STAGE3B)
b["short_dataset"] = b["dataset"].map(lambda x: f"{LABELS[x][0]} — {LABELS[x][1]}")
b["delta_pp"] = 100 * (b["sequential_fraction"].astype(float) - b["standalone_fraction"].astype(float))
b_plot = pd.concat([b[b["resolution"].eq(r)].sort_values("delta_pp")
                    for r in ["spot-resolution", "MERSCOPE"]], ignore_index=True)
b_plot["plot_order"] = np.arange(len(b_plot))
b_plot[["dataset", "target", "short_dataset", "resolution", "capacity", "total_units",
        "standalone_withheld", "standalone_fraction", "sequential_withheld", "sequential_fraction",
        "delta_withheld", "delta_fraction", "delta_pp", "plot_order"]].to_csv(
            OUT / "c5_fig2_stage3b_incremental_diverging_bar_source.csv", index=False)

fig, ax = plt.subplots(figsize=(8.2, 5.9))
y = np.arange(len(b_plot), dtype=float)
y[10:] += 0.9
colors = np.where(b_plot["resolution"].eq("MERSCOPE"), MER, SPOT)
bars = ax.barh(y, b_plot["delta_pp"], height=0.56, color=colors, edgecolor="white", linewidth=0.5, zorder=3)
ax.axvline(0, color="#59636A", linewidth=0.9, zorder=2)
for yi, value in zip(y, b_plot["delta_pp"]):
    label = f"{value:+.2f}".replace("-", "−")
    if abs(value) < 0.005:
        ax.vlines(0, yi - 0.28, yi + 0.28, color=colors[int(np.where(y == yi)[0][0])], linewidth=2.4, zorder=4)
    ax.text(value + (0.17 if value >= 0 else -0.17), yi, label,
            ha="left" if value >= 0 else "right", va="center", fontsize=7.7, color=TEXT)
ax.axhline((y[9] + y[10]) / 2, color="#CDD5DA", linewidth=0.9)
ax.set_yticks(y, b_plot["short_dataset"])
ax.set_xlabel("Change in withheld spatial units (percentage points)")
ax.set_ylim(y[-1] + 0.65, -0.65)
xmin, xmax = min(-11.4, b_plot["delta_pp"].min() - 0.8), max(5.3, b_plot["delta_pp"].max() + 0.8)
ax.set_xlim(xmin, xmax)
ax.grid(axis="x", color="#E7EBEE", linewidth=0.7, zorder=0)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_color("#919AA0")
ax.tick_params(axis="y", length=0, pad=5)
ax.tick_params(axis="x", colors=TEXT)
spot_mean = b.loc[b["resolution"].eq("spot-resolution"), "delta_pp"].mean()
mer_mean = b.loc[b["resolution"].eq("MERSCOPE"), "delta_pp"].mean()
ax.text(0.01, 0.985, f"Spot-resolution mean Δ = {spot_mean:+.2f} pp", transform=ax.transAxes,
        ha="left", va="top", color=SPOT, fontsize=8.8, fontweight="bold")
ax.text(0.01, 0.315, f"MERSCOPE mean Δ = {mer_mean:+.2f} pp", transform=ax.transAxes,
        ha="left", va="top", color=MER, fontsize=8.8, fontweight="bold")
fig.subplots_adjust(left=0.37, right=0.985, bottom=0.11, top=0.98)
save(fig, "c5_fig2_stage3b_incremental_diverging_bar")

print("SUCCESS")
