from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
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
A_ON = "#8EA7F8"
B_ON = "#C7A6E8"
OFF = "#EDE8E1"
ABS_CMAP = LinearSegmentedColormap.from_list(
    "absolute_outcomes", ["#FFF3CF", "#F9D59D", "#F5B47F", "#EF8E7E"]
)
EFFECT_CMAP = LinearSegmentedColormap.from_list(
    "module_effects", ["#879FF0", "#FFF8EE", "#D95776"]
)

ORDER = [
    "highres_humanbreastcancerpatient1_profile_mask_monocytes_and_macrophages",
    "highres_humancoloncancerpatient1_profile_mask_fibroblasts",
    "highres_humanlungcancerpatient1_profile_mask_plasma_cells",
    "highres_humanmelanomapatient1_profile_mask_fibroblasts",
    "highres_humanmelanomapatient2_profile_mask_b_cells",
]
SHORT = {
    ORDER[0]: "Breast cancer P1", ORDER[1]: "Colon cancer P1",
    ORDER[2]: "Lung cancer P1", ORDER[3]: "Melanoma P1", ORDER[4]: "Melanoma P2",
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
source_cols = ["dataset", "target", "baseline", "stage3a_only", "stage3b_only", "full",
               "delta_a", "delta_b", "delta_b_after_a"]
df[source_cols].to_csv(OUT / "c5_fig2_ablation_heatmap_source.csv", index=False)

abs_cols = ["baseline", "stage3a_only", "stage3b_only", "full"]
eff_cols = ["delta_a", "delta_b", "delta_b_after_a"]
abs_labels = ["Baseline", "Stage3A-only", "Stage3B-only", "Full"]
eff_labels = ["ΔA", "ΔB", "ΔB|A"]
abs_x = [0.0, 1.0, 2.0, 3.0]
eff_x = [4.35, 5.20, 6.05]
abs_w, eff_w, cell_h = 0.93, 0.78, 0.82

abs_norm = Normalize(vmin=0.78, vmax=0.96)
effect_values = df[eff_cols].to_numpy(float)
eff_min = min(-0.012, float(effect_values.min()))
eff_max = max(0.105, float(effect_values.max()))
eff_norm = TwoSlopeNorm(vmin=eff_min, vcenter=0.0, vmax=eff_max)

fig, ax = plt.subplots(figsize=(10.0, 5.75), facecolor=BACKGROUND)
ax.set_facecolor(BACKGROUND)
ax.axis("off")
ax.set_xlim(-2.50, 6.95)
ax.set_ylim(-0.32, 6.28)

ax.text(1.965, 6.10, "Absolute outcomes", ha="center", va="center",
        fontsize=9.6, color=TEXT, fontweight="bold")
ax.text(5.595, 6.10, "Module effects", ha="center", va="center",
        fontsize=9.6, color=TEXT, fontweight="bold")

# Lightweight 2 x 4 factorial annotation strips.
ax.text(-0.13, 5.84, "Stage3A", ha="right", va="center",
        fontsize=7.5, color=A_ON, fontweight="bold")
ax.text(-0.13, 5.56, "Stage3B", ha="right", va="center",
        fontsize=7.5, color=B_ON, fontweight="bold")
a_states = [False, True, False, True]
b_states = [False, False, True, True]
for j, x in enumerate(abs_x):
    for y, state, color in [(5.755, a_states[j], A_ON), (5.475, b_states[j], B_ON)]:
        ax.add_patch(Rectangle((x + 0.10, y), 0.73, 0.17,
                               facecolor=color if state else OFF,
                               edgecolor=BACKGROUND, linewidth=0.7))
        ax.text(x + 0.465, y + 0.085, "ON" if state else "OFF",
                ha="center", va="center", fontsize=6.6,
                color="white" if state else AUX, fontweight="bold")

for x, label in zip(abs_x, abs_labels):
    ax.text(x + abs_w / 2, 5.20, label, ha="center", va="center",
            fontsize=7.8, color=TEXT)
for x, label in zip(eff_x, eff_labels):
    ax.text(x + eff_w / 2, 5.20, label, ha="center", va="center",
            fontsize=8.2, color=TEXT, fontweight="bold")
ax.text(eff_x[2] + eff_w / 2, 5.04, "Full − A-only", ha="center", va="center",
        fontsize=6.2, color=AUX)

# Restrained gap/divider between absolute and effect blocks.
ax.add_patch(Rectangle((4.02, -0.05), 0.10, 5.20,
                       facecolor="#F0E9E3", edgecolor="none", zorder=0))

for i, row in df.iterrows():
    y = 4.18 - i
    ax.text(-2.40, y + 0.51, SHORT[row["dataset"]], ha="left", va="center",
            fontsize=8.8, color=TEXT, fontweight="bold")
    ax.text(-2.40, y + 0.26, str(row["target"]), ha="left", va="center",
            fontsize=7.2, color=AUX)

    for j, col in enumerate(abs_cols):
        value = float(row[col])
        x = abs_x[j]
        ax.add_patch(Rectangle((x, y), abs_w, cell_h,
                               facecolor=ABS_CMAP(abs_norm(value)),
                               edgecolor=BACKGROUND, linewidth=1.0))
        ax.text(x + abs_w / 2, y + cell_h / 2, f"{value:.3f}",
                ha="center", va="center", fontsize=8.9, color=TEXT, fontweight="bold")

    for j, col in enumerate(eff_cols):
        value = float(row[col])
        x = eff_x[j]
        color = EFFECT_CMAP(eff_norm(value))
        ax.add_patch(Rectangle((x, y), eff_w, cell_h, facecolor=color,
                               edgecolor=BACKGROUND, linewidth=1.0))
        rgb = np.array(color[:3])
        luminance = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
        ax.text(x + eff_w / 2, y + cell_h / 2,
                f"{value:+.4f}".replace("-", "−"), ha="center", va="center",
                fontsize=8.7, color="white" if luminance < 0.58 else TEXT,
                fontweight="bold")

# Fine row guides preserve a unified scientific matrix field.
for i in range(6):
    y_line = 4.14 - i
    ax.plot([-2.42, 6.83], [y_line, y_line], color=GRID, linewidth=0.48, zorder=0)

cax_abs = fig.add_axes([0.285, 0.057, 0.225, 0.017])
cax_eff = fig.add_axes([0.650, 0.057, 0.225, 0.017])
cb_abs = mpl.colorbar.ColorbarBase(cax_abs, cmap=ABS_CMAP, norm=abs_norm,
                                   orientation="horizontal")
cb_eff = mpl.colorbar.ColorbarBase(cax_eff, cmap=EFFECT_CMAP, norm=eff_norm,
                                   orientation="horizontal")
cb_abs.set_label("Low-support suppression score", fontsize=7.8, color=TEXT, labelpad=3)
cb_eff.set_label("Module effect", fontsize=7.8, color=TEXT, labelpad=3)
cb_abs.set_ticks([0.78, 0.84, 0.90, 0.96])
cb_eff.set_ticks([eff_min, 0.0, 0.05, 0.10])
cb_eff.set_ticklabels([f"{eff_min:.2f}".replace("-", "−"), "0", "+0.05", "+0.10"])
for cb in [cb_abs, cb_eff]:
    cb.ax.tick_params(labelsize=6.7, length=2, pad=2, colors=AUX)
    cb.outline.set_edgecolor(GRID)
    cb.outline.set_linewidth(0.6)

fig.subplots_adjust(left=0.035, right=0.985, top=0.985, bottom=0.105)

stem = "c5_fig2_ablation_heatmap"
fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor=BACKGROUND)
fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor=BACKGROUND)
fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", facecolor=BACKGROUND)
plt.close(fig)
print("SUCCESS")
