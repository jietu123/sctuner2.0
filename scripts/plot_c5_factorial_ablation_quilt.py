from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualizations" / "c5_decomposition"
INPUT = ROOT / "result" / "c5_decomposition" / "c5_merscope_four_route_suppression_wide.csv"
OUT.mkdir(parents=True, exist_ok=True)

BACKGROUND = "#FFF9F4"
CARD = "#FFFDFC"
TEXT = "#2F3542"
AUX = "#667085"
STAGE3A = "#8EA7F8"
STAGE3B = "#C7A6E8"
INTERACTION = "#D95776"
BORDER = "#E7DED6"
CMAP = LinearSegmentedColormap.from_list(
    "editorial_suppression", ["#FAF2E8", "#F7D8BE", "#F3BA98", "#E99A82"]
)

SHORT = {
    "highres_humanbreastcancerpatient1_profile_mask_monocytes_and_macrophages": "Breast cancer P1",
    "highres_humancoloncancerpatient1_profile_mask_fibroblasts": "Colon cancer P1",
    "highres_humanlungcancerpatient1_profile_mask_plasma_cells": "Lung cancer P1",
    "highres_humanmelanomapatient1_profile_mask_fibroblasts": "Melanoma P1",
    "highres_humanmelanomapatient2_profile_mask_b_cells": "Melanoma P2",
}

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 9,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "savefig.facecolor": BACKGROUND,
})

df = pd.read_csv(INPUT)
source_cols = [
    "dataset", "target", "capacity", "baseline", "stage3a_only", "stage3b_only", "full",
    "effect_stage3a", "effect_stage3b", "effect_full",
    "increment_stage3b_after_stage3a", "increment_stage3a_after_stage3b",
]
df[source_cols].to_csv(OUT / "c5_fig2_factorial_ablation_quilt_source.csv", index=False)

score_cols = ["baseline", "stage3a_only", "stage3b_only", "full"]
all_scores = df[score_cols].to_numpy(float)
norm = Normalize(vmin=float(all_scores.min()), vmax=float(all_scores.max()))


def rounded(ax, xy, width, height, facecolor, edgecolor=BORDER, linewidth=0.8,
            radius=0.025, zorder=1):
    patch = FancyBboxPatch(
        xy, width, height,
        boxstyle=f"round,pad=0.008,rounding_size={radius}",
        facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth, zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def draw_card(fig, position, row):
    ax = fig.add_axes(position)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    rounded(ax, (0.01, 0.01), 0.98, 0.98, CARD, BORDER, 1.0, 0.035)

    ax.text(0.055, 0.925, SHORT[row["dataset"]], ha="left", va="center",
            fontsize=10.2, color=TEXT, fontweight="bold")
    ax.text(0.055, 0.845, f"Masked target: {row['target']}", ha="left", va="center",
            fontsize=7.5, color=AUX)

    tile_w, tile_h = 0.205, 0.225
    x0, x1 = 0.165, 0.397
    y_top, y_bottom = 0.515, 0.245
    ax.text((x0 + x1 + tile_w) / 2, 0.785, "Stage3A", ha="center", va="center",
            fontsize=7.8, color=STAGE3A, fontweight="bold")
    ax.text(x0 + tile_w / 2, 0.735, "OFF", ha="center", va="center", fontsize=6.8, color=AUX)
    ax.text(x1 + tile_w / 2, 0.735, "ON", ha="center", va="center", fontsize=6.8, color=AUX)
    ax.text(0.075, 0.515, "Stage3B\nOFF", ha="center", va="center",
            fontsize=6.8, color=STAGE3B, fontweight="bold", linespacing=1.0)
    ax.text(0.075, 0.245, "Stage3B\nON", ha="center", va="center",
            fontsize=6.8, color=STAGE3B, fontweight="bold", linespacing=1.0)

    tiles = [
        (x0, y_top, "Baseline", float(row["baseline"])),
        (x1, y_top, "Stage3A-only", float(row["stage3a_only"])),
        (x0, y_bottom, "Stage3B-only", float(row["stage3b_only"])),
        (x1, y_bottom, "Full", float(row["full"])),
    ]
    best_value = max(value for _, _, _, value in tiles)
    for x, y, label, value in tiles:
        is_best_full = label == "Full" and np.isclose(value, best_value)
        edge = INTERACTION if is_best_full else "#F4ECE6"
        rounded(ax, (x, y), tile_w, tile_h, CMAP(norm(value)), edge,
                1.8 if is_best_full else 0.8, 0.018, 2)
        ax.text(x + tile_w / 2, y + 0.135, f"{value:.3f}", ha="center", va="center",
                fontsize=10.1, color=TEXT, fontweight="bold")
        ax.text(x + tile_w / 2, y + 0.055, label, ha="center", va="center",
                fontsize=6.3, color=AUX)
        if is_best_full:
            ax.text(x + tile_w - 0.018, y + tile_h - 0.018, "best", ha="right", va="top",
                    fontsize=5.8, color=INTERACTION, fontweight="bold")

    effect_x, effect_w = 0.655, 0.295
    ax.text(effect_x, 0.755, "Effects", ha="left", va="center",
            fontsize=7.7, color=TEXT, fontweight="bold")
    effects = [
        ("A effect", float(row["effect_stage3a"]), STAGE3A, "#F0F3FF"),
        ("B effect", float(row["effect_stage3b"]), STAGE3B, "#F5EFFA"),
        ("B after A", float(row["increment_stage3b_after_stage3a"]), INTERACTION, "#FFF0F2"),
    ]
    for idx, (label, value, accent, fill) in enumerate(effects):
        y = 0.585 - idx * 0.165
        rounded(ax, (effect_x, y), effect_w, 0.13, fill, "none", 0, 0.015, 1)
        ax.add_patch(FancyBboxPatch((effect_x, y), 0.014, 0.13,
                                   boxstyle="round,pad=0,rounding_size=0.006",
                                   facecolor=accent, edgecolor="none", zorder=2))
        ax.text(effect_x + 0.03, y + 0.088, label, ha="left", va="center",
                fontsize=6.7, color=AUX)
        ax.text(effect_x + 0.03, y + 0.038, f"{value:+.3f}".replace("-", "−"),
                ha="left", va="center", fontsize=8.4, color=TEXT, fontweight="bold")


fig = plt.figure(figsize=(12.0, 7.6), facecolor=BACKGROUND)
fig.text(0.04, 0.974, "Factorial ablation of low-support suppression",
         ha="left", va="top", fontsize=12, color=TEXT, fontweight="bold")
fig.text(0.04, 0.944, "MERSCOPE profile-masking datasets · capacity = 1",
         ha="left", va="top", fontsize=8.3, color=AUX)

positions = [
    [0.035, 0.545, 0.30, 0.375],
    [0.350, 0.545, 0.30, 0.375],
    [0.665, 0.545, 0.30, 0.375],
    [0.190, 0.145, 0.30, 0.375],
    [0.510, 0.145, 0.30, 0.375],
]
for position, (_, row) in zip(positions, df.iterrows()):
    draw_card(fig, position, row)

cax = fig.add_axes([0.365, 0.073, 0.27, 0.018])
cb = mpl.colorbar.ColorbarBase(cax, cmap=CMAP, norm=norm, orientation="horizontal")
cb.ax.tick_params(labelsize=7, length=2, pad=2, colors=AUX)
cb.set_label("Low-support suppression score", fontsize=8, color=TEXT, labelpad=3)
cb.outline.set_edgecolor(BORDER)
cb.outline.set_linewidth(0.7)

stem = "c5_fig2_factorial_ablation_quilt"
fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor=BACKGROUND)
fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor=BACKGROUND)
fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", facecolor=BACKGROUND)
plt.close(fig)
print("SUCCESS")
