from pathlib import Path
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualizations" / "c5_decomposition"
OUT.mkdir(parents=True, exist_ok=True)
STAGE3A = ROOT / "result" / "c5_stage3a_decomposition" / "c5_stage3a_exclusion_by_dataset.csv"
STAGE3B = ROOT / "result" / "c5_decomposition" / "c5_stage3b_decomposition_by_dataset.csv"

TEXT = "#263238"
TARGET = "#2A7F83"
OFF_TARGET = "#C15B64"
GROUP_SPOT = "#557A95"
GROUP_MER = "#C77C3D"

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 9, "axes.labelsize": 10,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "savefig.facecolor": "white",
})

LABELS = {
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


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def off_target_count(value):
    if pd.isna(value):
        return 0
    hits = re.findall(r"(?:n\s*=\s*|:)\s*(\d+)", str(value), flags=re.I)
    return sum(int(x) for x in hits)


# Figure 1: excluded-cell composition.
a = pd.read_csv(STAGE3A)
a["display_label"] = a["dataset"].map(LABELS)
a["target_excluded_cells"] = a["target_dropped_cells"].astype(int)
a["off_target_excluded_cells"] = a["off_target_details"].map(off_target_count).astype(int)
a["total_excluded_cells"] = a["target_excluded_cells"] + a["off_target_excluded_cells"]
a["off_target_population"] = np.where(a["off_target_excluded_cells"].gt(0), "NK cells", "")
a[["dataset", "target", "display_label", "resolution", "capacity", "target_total_cells",
   "target_excluded_cells", "off_target_population", "off_target_excluded_cells",
   "total_excluded_cells", "target_detected_missing"]].to_csv(
       OUT / "c5_fig1_stage3a_excluded_cell_composition_source.csv", index=False)

fig, ax = plt.subplots(figsize=(8.3, 5.8))
y = np.arange(len(a), dtype=float)
y[10:] += 0.8
target_n = a["target_excluded_cells"].to_numpy()
off_n = a["off_target_excluded_cells"].to_numpy()
ax.barh(y, target_n, height=0.57, color=TARGET, edgecolor="white", linewidth=0.45,
        label="Masked target cells", zorder=3)
ax.barh(y, off_n, left=target_n, height=0.57, color=OFF_TARGET, edgecolor="white", linewidth=0.45,
        label="Off-target cells", zorder=4)
max_total = int(a["total_excluded_cells"].max())
for yi, tn, on in zip(y, target_n, off_n):
    ax.text(tn + max_total * 0.012, yi, f"{tn:,}", ha="left", va="center", fontsize=7.4, color=TEXT)
    if on > 0:
        ax.text(tn + on / 2, yi, "NK\n66", ha="center", va="center", fontsize=6.7,
                color="white", fontweight="bold", linespacing=0.9)
ax.axhline((y[9] + y[10]) / 2, color="#CDD5DA", linewidth=0.9)
ax.set_yticks(y, a["display_label"])
ax.set_xlabel("Reference cells excluded by Stage3A (n)")
ax.set_xlim(0, max_total * 1.16)
ax.set_ylim(y[-1] + 0.65, -0.65)
ax.grid(axis="x", color="#E7EBEE", linewidth=0.7, zorder=0)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_color("#919AA0")
ax.tick_params(axis="y", length=0, pad=5)
ax.text(0.0, 1.015, "Masked target fully excluded in 15/15 experiments",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=10,
        color=TEXT, fontweight="bold")
ax.text(0.01, 0.965, "Spot-resolution (n=10)", transform=ax.transAxes,
        ha="left", va="top", fontsize=8.5, color=GROUP_SPOT, fontweight="bold")
ax.text(0.01, 0.305, "MERSCOPE (n=5; capacity=1)", transform=ax.transAxes,
        ha="left", va="top", fontsize=8.5, color=GROUP_MER, fontweight="bold")
ax.legend(handles=[Patch(facecolor=TARGET, label="Masked target cells"),
                   Patch(facecolor=OFF_TARGET, label="Off-target cells")],
          loc="lower right", frameon=False, fontsize=8, ncol=2,
          handlelength=1.3, columnspacing=1.1)
fig.subplots_adjust(left=0.36, right=0.985, bottom=0.11, top=0.95)
save(fig, "c5_fig1_stage3a_excluded_cell_composition")


# Figure 2: four-column numeric heatmap.
b = pd.read_csv(STAGE3B)
b["display_label"] = b["dataset"].map(LABELS)
b["standalone_percent"] = 100 * b["standalone_fraction"].astype(float)
b["sequential_percent"] = 100 * b["sequential_fraction"].astype(float)
b["delta_pp"] = 100 * (b["sequential_fraction"].astype(float) - b["standalone_fraction"].astype(float))
b["mask_jaccard"] = b["jaccard"].astype(float)
b_plot = pd.concat([b[b["resolution"].eq("spot-resolution")], b[b["resolution"].eq("MERSCOPE")]],
                   ignore_index=True)
b_plot["heatmap_row"] = np.arange(len(b_plot))
b_plot[["dataset", "target", "display_label", "resolution", "capacity", "total_units",
        "standalone_withheld", "standalone_fraction", "standalone_percent",
        "sequential_withheld", "sequential_fraction", "sequential_percent",
        "delta_withheld", "delta_fraction", "delta_pp", "mask_intersection", "mask_union",
        "mask_jaccard", "heatmap_row"]].to_csv(
            OUT / "c5_fig2_stage3b_decomposition_heatmap_source.csv", index=False)

values = b_plot[["standalone_percent", "sequential_percent", "delta_pp", "mask_jaccard"]].to_numpy()
rate_norm = Normalize(vmin=0, vmax=max(values[:, 0].max(), values[:, 1].max()))
delta_lim = max(abs(values[:, 2].min()), abs(values[:, 2].max()))
delta_norm = TwoSlopeNorm(vmin=-delta_lim, vcenter=0, vmax=delta_lim)
jaccard_norm = Normalize(vmin=0, vmax=1)
cmaps = [mpl.colormaps["YlGnBu"], mpl.colormaps["YlGnBu"], mpl.colormaps["RdBu_r"], mpl.colormaps["Greens"]]
norms = [rate_norm, rate_norm, delta_norm, jaccard_norm]

fig, ax = plt.subplots(figsize=(7.2, 6.2))
ypos = np.arange(len(b_plot), dtype=float)
ypos[10:] += 0.7
for i, yi in enumerate(ypos):
    for j in range(4):
        val = values[i, j]
        color = cmaps[j](norms[j](val))
        ax.add_patch(Rectangle((j, yi - 0.43), 0.94, 0.86,
                               facecolor=color, edgecolor="white", linewidth=1.0))
        if j < 2:
            label = f"{val:.2f}"
        elif j == 2:
            label = f"{val:+.2f}".replace("-", "−")
        else:
            label = f"{val:.2f}"
        rgb = np.array(color[:3])
        luminance = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
        ax.text(j + 0.47, yi, label, ha="center", va="center", fontsize=8,
                color="white" if luminance < 0.54 else TEXT,
                fontweight="bold" if luminance < 0.45 else "normal")

ax.set_xlim(0, 3.94)
ax.set_ylim(ypos[-1] + 0.55, -0.7)
ax.set_yticks(ypos, b_plot["display_label"])
ax.set_xticks(np.arange(4) + 0.47,
              ["Standalone\nwithheld (%)", "Sequential\nwithheld (%)",
               "Δ withheld\n(pp)", "Mask\nJaccard"])
ax.xaxis.tick_top()
ax.tick_params(axis="x", length=0, pad=8, colors=TEXT)
ax.tick_params(axis="y", length=0, pad=5, colors=TEXT)
for spine in ax.spines.values():
    spine.set_visible(False)
ax.axhline((ypos[9] + ypos[10]) / 2, color="#AAB4BA", linewidth=1.1)
ax.text(-0.03, 0.985, "Spot-resolution (n=10)", transform=ax.transAxes,
        ha="right", va="top", rotation=90, fontsize=8, color=GROUP_SPOT, fontweight="bold")
ax.text(-0.03, 0.285, "MERSCOPE (n=5; capacity=1)", transform=ax.transAxes,
        ha="right", va="top", rotation=90, fontsize=8, color=GROUP_MER, fontweight="bold")

# Compact scale keys for the three distinct numerical scales.
cax1 = fig.add_axes([0.42, 0.055, 0.16, 0.015])
cax2 = fig.add_axes([0.63, 0.055, 0.12, 0.015])
cax3 = fig.add_axes([0.80, 0.055, 0.12, 0.015])
cb1 = mpl.colorbar.ColorbarBase(cax1, cmap=cmaps[0], norm=rate_norm, orientation="horizontal")
cb2 = mpl.colorbar.ColorbarBase(cax2, cmap=cmaps[2], norm=delta_norm, orientation="horizontal")
cb3 = mpl.colorbar.ColorbarBase(cax3, cmap=cmaps[3], norm=jaccard_norm, orientation="horizontal")
for cb, label in [(cb1, "Withheld (%)"), (cb2, "Δ (pp)"), (cb3, "Jaccard")]:
    cb.ax.tick_params(labelsize=6.5, length=2, pad=1)
    cb.set_label(label, fontsize=7, labelpad=1)
    cb.outline.set_linewidth(0.5)
fig.subplots_adjust(left=0.39, right=0.97, bottom=0.12, top=0.90)
save(fig, "c5_fig2_stage3b_decomposition_heatmap")

print("SUCCESS")
