from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualizations" / "c5_decomposition"
SOURCE = OUT / "c5_fig3_spatial_transition_map_source.csv"

BACKGROUND = "#FFF9F4"
TEXT = "#2F3542"
AUX = "#667085"
COLORS = {
    "Neither withheld": "#BFB6AC",
    "Withheld by both": "#8EA7F8",
    "Standalone-only withheld": "#F3B27F",
    "Sequential-only withheld": "#D95776",
}
ORDER = [
    "Neither withheld",
    "Withheld by both",
    "Standalone-only withheld",
    "Sequential-only withheld",
]

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 9,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "savefig.facecolor": BACKGROUND,
})

data = pd.read_csv(SOURCE)
counts = data["mask_transition"].value_counts().reindex(ORDER, fill_value=0).astype(int)
expected = {
    "Neither withheld": 2039,
    "Withheld by both": 253,
    "Standalone-only withheld": 269,
    "Sequential-only withheld": 0,
}
if len(data) != 2561 or counts.to_dict() != expected:
    raise RuntimeError("Frozen spatial transition counts do not match the canonical C5 result")

new_source = OUT / "c5_fig3_spatial_transition_map_source.csv"
data.to_csv(new_source, index=False)

fig = plt.figure(figsize=(6.9, 5.25), facecolor=BACKGROUND)
ax = fig.add_axes([0.025, 0.035, 0.68, 0.93])
ax.set_facecolor(BACKGROUND)

draw_specs = {
    "Neither withheld": {"size": 9.0, "alpha": 0.50, "zorder": 1},
    "Standalone-only withheld": {"size": 15.2, "alpha": 0.94, "zorder": 2},
    "Withheld by both": {"size": 15.7, "alpha": 0.96, "zorder": 3},
    "Sequential-only withheld": {"size": 16.0, "alpha": 0.97, "zorder": 4},
}
draw_order = [
    "Neither withheld",
    "Standalone-only withheld",
    "Withheld by both",
    "Sequential-only withheld",
]
for state in draw_order:
    subset = data[data["mask_transition"].eq(state)]
    if subset.empty:
        continue
    spec = draw_specs[state]
    ax.scatter(subset["col"], subset["row"], s=spec["size"],
               color=COLORS[state], edgecolor="none", alpha=spec["alpha"],
               zorder=spec["zorder"], rasterized=False)

ax.set_aspect("equal", adjustable="datalim")
ax.invert_yaxis()
ax.set_xticks([])
ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)

fig.text(0.745, 0.895,
         "Standalone: 522 / 2561 (20.38%)\nSequential: 253 / 2561 (9.88%)",
         ha="left", va="top", fontsize=8.6, color=TEXT, linespacing=1.55)

handles = [
    Patch(facecolor=COLORS[state], edgecolor="none",
          label=f"{state} (n={expected[state]})")
    for state in ORDER
]
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.735, 0.745),
           frameon=False, fontsize=8.0, handlelength=1.25, handleheight=0.85,
           borderaxespad=0, labelspacing=0.9, handletextpad=0.55)

stem = "c5_fig3_spatial_transition_map"
fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor=BACKGROUND)
fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor=BACKGROUND)
fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", facecolor=BACKGROUND)
plt.close(fig)
print("SUCCESS")
