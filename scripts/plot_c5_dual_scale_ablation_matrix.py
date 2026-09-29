from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, SymLogNorm
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualizations" / "c5_decomposition"
INPUT = ROOT / "result" / "c5_decomposition" / "c5_merscope_four_route_suppression_wide.csv"
OUT.mkdir(parents=True, exist_ok=True)

BACKGROUND = "#FFF9F4"
TEXT = "#2F3542"
AUX = "#667085"
GRID = "#E7DED6"
STAGE3A = "#8EA7F8"
STAGE3B = "#C7A6E8"
OFF = "#EDE8E1"

ABS_CMAP = LinearSegmentedColormap.from_list(
    "warm_outcome", ["#FFF1C7", "#F8D19A", "#F3B27F", "#ED8F82"]
)
EFFECT_CMAP = LinearSegmentedColormap.from_list(
    "module_effect", ["#8EA7F8", "#FFF6E9", "#D95776"]
)

ORDER = [
    "highres_humanbreastcancerpatient1_profile_mask_monocytes_and_macrophages",
    "highres_humancoloncancerpatient1_profile_mask_fibroblasts",
    "highres_humanlungcancerpatient1_profile_mask_plasma_cells",
    "highres_humanmelanomapatient1_profile_mask_fibroblasts",
    "highres_humanmelanomapatient2_profile_mask_b_cells",
]
SHORT = {
    ORDER[0]: "Breast cancer P1",
    ORDER[1]: "Colon cancer P1",
    ORDER[2]: "Lung cancer P1",
    ORDER[3]: "Melanoma P1",
    ORDER[4]: "Melanoma P2",
}

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 9,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "savefig.facecolor": BACKGROUND,
})

df = pd.read_csv(INPUT).set_index("dataset").loc[ORDER].reset_index()
df["delta_a"] = df["stage3a_only"] - df["baseline"]
df["delta_b"] = df["stage3b_only"] - df["baseline"]
df["delta_b_after_a"] = df["full"] - df["stage3a_only"]
df["delta_full"] = df["full"] - df["baseline"]

source_cols = [
    "dataset", "target", "baseline", "stage3a_only", "stage3b_only", "full",
    "delta_a", "delta_b", "delta_b_after_a", "delta_full",
]
df[source_cols].to_csv(OUT / "c5_fig2_dual_scale_ablation_matrix_source.csv", index=False)

absolute_cols = ["baseline", "stage3a_only", "stage3b_only", "full"]
effect_cols = ["delta_a", "delta_b", "delta_b_after_a", "delta_full"]
absolute_labels = ["Baseline", "Stage3A-only", "Stage3B-only", "Full"]
effect_labels = ["ΔA", "ΔB", "ΔB|A", "ΔFull"]
x_positions = [0.0, 1.0, 2.0, 3.0, 4.45, 5.45, 6.45, 7.45]
cell_w, cell_h = 0.92, 0.82

absolute_values = df[absolute_cols].to_numpy(float)
effect_values = df[effect_cols].to_numpy(float)
abs_norm = Normalize(vmin=float(absolute_values.min()), vmax=float(absolute_values.max()))
effect_limit = float(np.max(np.abs(effect_values)))
effect_norm = SymLogNorm(linthresh=0.002, linscale=0.65, vmin=-effect_limit,
                         vmax=effect_limit, base=10)

fig, ax = plt.subplots(figsize=(11.4, 6.6), facecolor=BACKGROUND)
ax.set_facecolor(BACKGROUND)
ax.axis("off")
ax.set_xlim(-2.55, 8.55)
ax.set_ylim(-0.72, 6.78)

# Unified block headers.
ax.text(1.96, 6.52, "Absolute outcomes", ha="center", va="center",
        fontsize=10.5, color=TEXT, fontweight="bold")
ax.text(6.41, 6.52, "Module effects", ha="center", va="center",
        fontsize=10.5, color=TEXT, fontweight="bold")

# Factorial column annotations.
ax.text(-0.15, 6.08, "Stage3A", ha="right", va="center",
        fontsize=8.0, color=STAGE3A, fontweight="bold")
ax.text(-0.15, 5.73, "Stage3B", ha="right", va="center",
        fontsize=8.0, color=STAGE3B, fontweight="bold")
a_states = [False, True, False, True]
b_states = [False, False, True, True]
for j, x in enumerate(x_positions[:4]):
    for y, state, on_color in [(5.92, a_states[j], STAGE3A), (5.57, b_states[j], STAGE3B)]:
        color = on_color if state else OFF
        ax.add_patch(Rectangle((x + 0.08, y), 0.76, 0.25,
                               facecolor=color, edgecolor=BACKGROUND, linewidth=0.8))
        ax.text(x + 0.46, y + 0.125, "ON" if state else "OFF",
                ha="center", va="center", fontsize=7.0,
                color="white" if state else AUX, fontweight="bold")

# Column labels.
for x, label in zip(x_positions[:4], absolute_labels):
    ax.text(x + cell_w / 2, 5.31, label, ha="center", va="center",
            fontsize=8.0, color=TEXT)
for x, label in zip(x_positions[4:], effect_labels):
    ax.text(x + cell_w / 2, 5.31, label, ha="center", va="center",
            fontsize=8.4, color=TEXT, fontweight="bold")
ax.plot([3.08, 3.84], [5.16, 5.16], color="#D95776", linewidth=1.5)

# Divider between the two numerical scales.
ax.add_patch(Rectangle((4.08, -0.08), 0.12, 5.35,
                       facecolor="#F0E9E3", edgecolor="none", zorder=0))

for i, row in df.iterrows():
    y = 4.22 - i
    ax.text(-2.42, y + 0.50, SHORT[row["dataset"]], ha="left", va="center",
            fontsize=9.1, color=TEXT, fontweight="bold")
    ax.text(-2.42, y + 0.25, str(row["target"]), ha="left", va="center",
            fontsize=7.5, color=AUX)

    for j, col in enumerate(absolute_cols):
        value = float(row[col])
        x = x_positions[j]
        edge = "#D86F70" if col == "full" else BACKGROUND
        ax.add_patch(Rectangle((x, y), cell_w, cell_h,
                               facecolor=ABS_CMAP(abs_norm(value)), edgecolor=edge,
                               linewidth=1.15 if col == "full" else 0.9))
        ax.text(x + cell_w / 2, y + cell_h / 2, f"{value:.3f}",
                ha="center", va="center", fontsize=9.1, color=TEXT, fontweight="bold")

    for j, col in enumerate(effect_cols):
        value = float(row[col])
        x = x_positions[j + 4]
        ax.add_patch(Rectangle((x, y), cell_w, cell_h,
                               facecolor=EFFECT_CMAP(effect_norm(value)),
                               edgecolor=BACKGROUND, linewidth=0.9))
        rgb = np.array(EFFECT_CMAP(effect_norm(value))[:3])
        luminance = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
        ax.text(x + cell_w / 2, y + cell_h / 2,
                f"{value:+.3f}".replace("-", "−"), ha="center", va="center",
                fontsize=9.0, color="white" if luminance < 0.59 else TEXT,
                fontweight="bold")

# Fine matrix guides without turning the figure into a default heatmap.
for i in range(6):
    y_line = 4.18 - i
    ax.plot([-2.45, 8.37], [y_line, y_line], color=GRID, linewidth=0.55, zorder=0)

# Dual legends.
cax_abs = fig.add_axes([0.285, 0.075, 0.235, 0.018])
cax_eff = fig.add_axes([0.635, 0.075, 0.235, 0.018])
cb_abs = mpl.colorbar.ColorbarBase(cax_abs, cmap=ABS_CMAP, norm=abs_norm,
                                   orientation="horizontal")
cb_eff = mpl.colorbar.ColorbarBase(cax_eff, cmap=EFFECT_CMAP, norm=effect_norm,
                                   orientation="horizontal")
cb_abs.set_label("Low-support suppression score", fontsize=8, color=TEXT, labelpad=3)
cb_eff.set_label("Module effect", fontsize=8, color=TEXT, labelpad=3)
for cb in [cb_abs, cb_eff]:
    cb.ax.tick_params(labelsize=6.8, length=2, pad=2, colors=AUX)
    cb.outline.set_edgecolor(GRID)
    cb.outline.set_linewidth(0.6)
cb_eff.set_ticks([-0.1, -0.01, 0, 0.01, 0.1])
cb_eff.set_ticklabels(["−0.10", "−0.01", "0", "+0.01", "+0.10"])

fig.text(0.5, 0.025,
         "Absolute outcomes are shown on the left; route-specific effects are shown on the right.",
         ha="center", va="bottom", fontsize=7.6, color=AUX)
fig.subplots_adjust(left=0.035, right=0.985, top=0.96, bottom=0.13)

stem = "c5_fig2_dual_scale_ablation_matrix"
fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor=BACKGROUND)
fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor=BACKGROUND)
fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", facecolor=BACKGROUND)
plt.close(fig)
print("SUCCESS")
