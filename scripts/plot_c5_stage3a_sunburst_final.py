from pathlib import Path
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Wedge
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualizations" / "c5_decomposition"
SOURCE = OUT / "c5_fig1_stage3a_sunburst_source.csv"
data = pd.read_csv(SOURCE)

MODE = os.environ.get("C5_STAGE3A_PALETTE", "C").lower()
if MODE == "c_plus":
    P = {
        "background": "#FFF9F4", "text": "#2F3542", "aux": "#667085",
        "leader": "#B7C0CC", "spot_group": "#8EA7F8", "merscope_group": "#C7A6E8",
        "off": "#D95776", "off_label": "#C94E6E", "off_line": "#DF7890",
        "center_accent": "#F0A85C", "target_legend": "#F6C667",
    }
    SPOT_EXPERIMENT = ["#DCE7FF", "#C9D9FF", "#E8F0FF"]
    MERSCOPE_EXPERIMENT = ["#EBDCF8", "#DCC5F2", "#F3EAFE"]
    SPOT_TARGET = ["#F6C667", "#F7B87A", "#F4C98B"]
    MERSCOPE_TARGET = ["#F3B26B", "#F6C07D", "#EFAF87"]
else:
    P = {
        "background": "#FFF9F3", "text": "#2F3542", "aux": "#616B77",
        "leader": "#D0D5DC", "spot_group": "#79C6F2", "merscope_group": "#C6A5E8",
        "off": "#E66A7A", "off_label": "#E66A7A", "off_line": "#E66A7A",
        "center_accent": "#F2C675", "target_legend": "#F2C675",
    }
    SPOT_EXPERIMENT = ["#DAF0FC"]
    MERSCOPE_EXPERIMENT = ["#EBDDF7"]
    SPOT_TARGET = ["#F2C675"]
    MERSCOPE_TARGET = ["#F2C675"]

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 9,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "savefig.facecolor": P["background"],
})

group_ranges = {}
for group in ["spot-resolution", "MERSCOPE"]:
    subset = data[data["resolution"].eq(group)]
    group_ranges[group] = (subset["experiment_angle_start"].max(),
                           subset["experiment_angle_end"].min())

fig, ax = plt.subplots(figsize=(9.2, 7.7), facecolor=P["background"])
ax.set_facecolor(P["background"])
ax.set_aspect("equal")
ax.axis("off")


def wedge(start, end, radius, width, color, linewidth=1.0):
    ax.add_patch(Wedge((0, 0), radius, end, start, width=width,
                       facecolor=color, edgecolor=P["background"], linewidth=linewidth))


for group in ["spot-resolution", "MERSCOPE"]:
    start, end = group_ranges[group]
    color = P["spot_group"] if group == "spot-resolution" else P["merscope_group"]
    wedge(start, end, 0.97, 0.38, color, 2.0)
    mid = np.deg2rad((start + end) / 2)
    label = "Spot-resolution\n(n=10)" if group == "spot-resolution" else "MERSCOPE\n(n=5, capacity=1)"
    rotation = np.rad2deg(mid) - 90
    if rotation < -90:
        rotation += 180
    if rotation > 90:
        rotation -= 180
    ax.text(0.78 * np.cos(mid), 0.78 * np.sin(mid), label,
            ha="center", va="center", color="white", fontsize=8.5,
            fontweight="bold", rotation=rotation, rotation_mode="anchor")

label_points = []
group_indices = {"spot-resolution": 0, "MERSCOPE": 0}
for _, row in data.iterrows():
    start, end = row["experiment_angle_start"], row["experiment_angle_end"]
    group = row["resolution"]
    group_index = group_indices[group]
    group_indices[group] += 1
    if group == "spot-resolution":
        experiment_color = SPOT_EXPERIMENT[group_index % len(SPOT_EXPERIMENT)]
        target_color = SPOT_TARGET[group_index % len(SPOT_TARGET)]
    else:
        experiment_color = MERSCOPE_EXPERIMENT[group_index % len(MERSCOPE_EXPERIMENT)]
        target_color = MERSCOPE_TARGET[group_index % len(MERSCOPE_TARGET)]
    wedge(start, end, 1.43, 0.44, experiment_color, 1.05)
    target_end = start - row["target_angular_extent"]
    wedge(start, target_end, 1.79, 0.34, target_color, 0.9)
    if row["off_target_excluded_cells"] > 0:
        wedge(target_end, end, 1.79, 0.34, P["off"], 0.9)
    mid = np.deg2rad((start + end) / 2)
    label_points.append({
        "x0": 1.81 * np.cos(mid), "y0": 1.81 * np.sin(mid),
        "side": 1 if np.cos(mid) >= 0 else -1,
        "label": f"{row['short_label']}\nTarget n={int(row['target_excluded_cells']):,}",
    })

for side in [-1, 1]:
    points = [point for point in label_points if point["side"] == side]
    points.sort(key=lambda point: point["y0"])
    lower, upper, gap = -1.73, 1.73, 0.225
    ys = [max(lower, min(upper, point["y0"])) for point in points]
    for i in range(1, len(ys)):
        ys[i] = max(ys[i], ys[i - 1] + gap)
    if ys and ys[-1] > upper:
        shift = ys[-1] - upper
        ys = [value - shift for value in ys]
        for i in range(len(ys) - 2, -1, -1):
            ys[i] = min(ys[i], ys[i + 1] - gap)
    x_text, x_elbow = 2.44 * side, 2.02 * side
    for point, y_text in zip(points, ys):
        ax.plot([point["x0"], x_elbow, x_text - 0.05 * side],
                [point["y0"], y_text, y_text], color=P["leader"], linewidth=0.55)
        ax.text(x_text, y_text, point["label"],
                ha="left" if side == 1 else "right", va="center", fontsize=7.0,
                color=P["text"], fontweight="normal", linespacing=1.05)

flag = data[data["off_target_excluded_cells"].gt(0)].iloc[0]
angle = np.deg2rad(flag["experiment_angle_end"] + flag["off_target_angular_extent"] / 2)
ax.annotate("NK cells\n(n=66)", xy=(1.66 * np.cos(angle), 1.66 * np.sin(angle)),
            xytext=(0.48, 1.92), ha="center", va="bottom", fontsize=7.5,
            color=P["off_label"], fontweight="bold",
            arrowprops=dict(arrowstyle="-", color=P["off_line"], linewidth=0.85,
                            connectionstyle="arc3,rad=-0.08"))

ax.text(0, 0.13, "15 / 15", ha="center", va="center", fontsize=22,
        color=P["text"], fontweight="bold")
ax.text(0, -0.08, "masked targets", ha="center", va="center", fontsize=9.5,
        color=P["aux"])
ax.text(0, -0.24, "fully excluded", ha="center", va="center", fontsize=9.5,
        color=P["center_accent"], fontweight="bold")
ax.legend(handles=[Patch(facecolor=P["target_legend"], edgecolor="none", label="Masked target exclusion"),
                   Patch(facecolor=P["off"], edgecolor="none", label="Additional off-target exclusion")],
          loc="lower center", bbox_to_anchor=(0.5, -0.045), ncol=2,
          frameon=False, fontsize=8, handlelength=1.4, columnspacing=1.6,
          labelcolor=P["text"])
ax.set_xlim(-3.25, 3.25)
ax.set_ylim(-2.30, 2.22)
fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.06)

stem = "c5_fig1_stage3a_sunburst_paletteC_plus" if MODE == "c_plus" else "c5_fig1_stage3a_sunburst_paletteC"
fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor=P["background"])
fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor=P["background"])
fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", facecolor=P["background"])
plt.close(fig)
print("SUCCESS")
