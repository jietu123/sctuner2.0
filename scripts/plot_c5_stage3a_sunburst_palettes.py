from pathlib import Path
import runpy

import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Wedge
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualizations" / "c5_decomposition"
ctx = runpy.run_path(str(ROOT / "scripts" / "plot_c5_stage3a_sunburst.py"))
source = ctx["source"]
group_ranges = ctx["group_ranges"]

PALETTES = {
    "paletteA": {
        "background": "#FFF9F2", "text": "#2E3440", "aux": "#5C6675", "leader": "#AAB2BD",
        "spot_group": "#86BFF2", "merscope_group": "#9FD6B2",
        "spot_experiment": "#DCECFB", "merscope_experiment": "#DDF1E4",
        "target": "#F3B37A", "off": "#D95C75",
    },
    "paletteB": {
        "background": "#FFFDF7", "text": "#2F3440", "aux": "#606A78", "leader": "#B3BBC5",
        "spot_group": "#8FA7F4", "merscope_group": "#F4D77C",
        "spot_experiment": "#DEE5FF", "merscope_experiment": "#FFF0BF",
        "target": "#F29A8A", "off": "#C94F7C",
    },
    "paletteC": {
        "background": "#FFF9F3", "text": "#2F3542", "aux": "#616B77", "leader": "#B4BCC7",
        "spot_group": "#79C6F2", "merscope_group": "#C6A5E8",
        "spot_experiment": "#DAF0FC", "merscope_experiment": "#EBDDF7",
        "target": "#F6BE5A", "off": "#E66A7A",
    },
}


def render(p, suffix):
    fig, ax = plt.subplots(figsize=(9.2, 7.7), facecolor=p["background"])
    ax.set_facecolor(p["background"])
    ax.set_aspect("equal")
    ax.axis("off")

    def wedge(start, end, radius, width, color, lw=1.0):
        ax.add_patch(Wedge((0, 0), radius, end, start, width=width,
                           facecolor=color, edgecolor=p["background"], linewidth=lw))

    for group in ["spot-resolution", "MERSCOPE"]:
        start, end, _ = group_ranges[group]
        color = p["spot_group"] if group == "spot-resolution" else p["merscope_group"]
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
    for _, row in source.reset_index(drop=True).iterrows():
        start, end = row["experiment_angle_start"], row["experiment_angle_end"]
        experiment_color = p["spot_experiment"] if row["resolution"] == "spot-resolution" else p["merscope_experiment"]
        wedge(start, end, 1.43, 0.44, experiment_color, 1.05)
        target_end = start - row["target_angular_extent"]
        wedge(start, target_end, 1.79, 0.34, p["target"], 0.9)
        if row["off_target_excluded_cells"] > 0:
            wedge(target_end, end, 1.79, 0.34, p["off"], 0.9)
        mid = np.deg2rad((start + end) / 2)
        label_points.append({
            "x0": 1.81 * np.cos(mid), "y0": 1.81 * np.sin(mid),
            "side": 1 if np.cos(mid) >= 0 else -1,
            "label": f"{row['short_label']}\nTarget n={row['target_excluded_cells']:,}",
            "flagged": row["off_target_excluded_cells"] > 0,
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
                    [point["y0"], y_text, y_text], color=p["leader"], linewidth=0.65)
            color = p["off"] if point["flagged"] else p["text"]
            ax.text(x_text, y_text, point["label"],
                    ha="left" if side == 1 else "right", va="center", fontsize=7.0,
                    color=color, fontweight="bold" if point["flagged"] else "normal",
                    linespacing=1.05)

    flag = source[source["off_target_excluded_cells"].gt(0)].iloc[0]
    angle = np.deg2rad(flag["experiment_angle_end"] + flag["off_target_angular_extent"] / 2)
    ax.annotate("NK cells\n(n=66)", xy=(1.66 * np.cos(angle), 1.66 * np.sin(angle)),
                xytext=(1.30, -2.02), ha="center", va="top", fontsize=7.5,
                color=p["off"], fontweight="bold",
                arrowprops=dict(arrowstyle="-", color=p["off"], linewidth=0.85,
                                connectionstyle="arc3,rad=-0.15"))

    ax.text(0, 0.13, "15 / 15", ha="center", va="center", fontsize=22,
            color=p["text"], fontweight="bold")
    ax.text(0, -0.08, "masked targets", ha="center", va="center", fontsize=9.5, color=p["aux"])
    ax.text(0, -0.24, "fully excluded", ha="center", va="center", fontsize=9.5,
            color=p["target"], fontweight="bold")
    ax.legend(handles=[Patch(facecolor=p["target"], edgecolor="none", label="Masked target exclusion"),
                       Patch(facecolor=p["off"], edgecolor="none", label="Additional off-target exclusion")],
              loc="lower center", bbox_to_anchor=(0.5, -0.045), ncol=2,
              frameon=False, fontsize=8, handlelength=1.4, columnspacing=1.6,
              labelcolor=p["text"])
    ax.set_xlim(-3.25, 3.25)
    ax.set_ylim(-2.30, 2.22)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.06)
    stem = f"c5_fig1_stage3a_sunburst_{suffix}"
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor=p["background"])
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor=p["background"])
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", facecolor=p["background"])
    plt.close(fig)


for suffix, palette in PALETTES.items():
    render(palette, suffix)
print("PALETTES_SUCCESS")
