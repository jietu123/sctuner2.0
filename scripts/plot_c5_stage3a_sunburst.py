from pathlib import Path
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Wedge
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualizations" / "c5_decomposition"
INPUT = ROOT / "result" / "c5_stage3a_decomposition" / "c5_stage3a_exclusion_by_dataset.csv"
OUT.mkdir(parents=True, exist_ok=True)

TEXT = "#27343B"
GROUP_COLORS = {"spot-resolution": "#557A95", "MERSCOPE": "#218B8D"}
GROUP_LIGHT = {
    "spot-resolution": ["#DDE7ED", "#CCDCE5"],
    "MERSCOPE": ["#D8ECEA", "#C6E2DF"],
}
TARGET_COLOR = "#4B938F"
OFF_COLOR = "#BF5962"

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 9,
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


def parse_off_target(value):
    if pd.isna(value):
        return 0
    hits = re.findall(r"(?:n\s*=\s*|:)\s*(\d+)", str(value), flags=re.I)
    return sum(int(x) for x in hits)


def add_ring_wedge(ax, start, end, r_outer, width, color, lw=1.0):
    ax.add_patch(Wedge((0, 0), r_outer, end, start, width=width,
                       facecolor=color, edgecolor="white", linewidth=lw))


df = pd.read_csv(INPUT)
df["short_label"] = df["dataset"].map(LABELS)
df["target_excluded_cells"] = df["target_dropped_cells"].astype(int)
df["off_target_excluded_cells"] = df["off_target_details"].map(parse_off_target).astype(int)
df["total_excluded_cells"] = df["target_excluded_cells"] + df["off_target_excluded_cells"]
df["group_label"] = np.where(df["resolution"].eq("MERSCOPE"),
                             "MERSCOPE (n=5, capacity=1)", "Spot-resolution (n=10)")

total = df["total_excluded_cells"].sum()
start_angle = 90.0
records = []
group_ranges = {}
experiment_angles = {}
for group in ["spot-resolution", "MERSCOPE"]:
    group_df = df[df["resolution"].eq(group)]
    group_total = group_df["total_excluded_cells"].sum()
    group_extent = 360.0 * group_total / total
    group_end = start_angle - group_extent
    group_ranges[group] = (start_angle, group_end, group_total)
    cursor = start_angle
    for _, row in group_df.iterrows():
        extent = 360.0 * row["total_excluded_cells"] / total
        end = cursor - extent
        experiment_angles[row["dataset"]] = (cursor, end)
        records.append({
            **row.to_dict(),
            "group_total_excluded_cells": int(group_total),
            "experiment_angle_start": cursor,
            "experiment_angle_end": end,
            "experiment_angular_extent": extent,
            "target_angular_extent": extent * row["target_excluded_cells"] / row["total_excluded_cells"],
            "off_target_angular_extent": extent * row["off_target_excluded_cells"] / row["total_excluded_cells"],
        })
        cursor = end
    start_angle = group_end

source = pd.DataFrame(records)
source[["dataset", "target", "short_label", "resolution", "capacity", "group_label",
        "target_total_cells", "target_excluded_cells", "off_target_excluded_cells",
        "total_excluded_cells", "off_target_exclusion", "off_target_details",
        "group_total_excluded_cells", "experiment_angle_start", "experiment_angle_end",
        "experiment_angular_extent", "target_angular_extent", "off_target_angular_extent"]].to_csv(
            OUT / "c5_fig1_stage3a_sunburst_source.csv", index=False)

fig, ax = plt.subplots(figsize=(9.2, 7.7))
ax.set_aspect("equal")
ax.axis("off")

# First ring: experiment groups.
for group in ["spot-resolution", "MERSCOPE"]:
    start, end, _ = group_ranges[group]
    add_ring_wedge(ax, start, end, r_outer=0.97, width=0.38,
                   color=GROUP_COLORS[group], lw=2.0)
    mid = np.deg2rad((start + end) / 2)
    radius = 0.78
    label = "Spot-resolution\n(n=10)" if group == "spot-resolution" else "MERSCOPE\n(n=5, capacity=1)"
    rotation = np.rad2deg(mid) - 90
    if rotation < -90:
        rotation += 180
    if rotation > 90:
        rotation -= 180
    ax.text(radius * np.cos(mid), radius * np.sin(mid), label,
            ha="center", va="center", color="white", fontsize=8.5,
            fontweight="bold", rotation=rotation, rotation_mode="anchor")

# Second and third rings: experiments and excluded-cell composition.
label_points = []
for idx, row in source.reset_index(drop=True).iterrows():
    start, end = row["experiment_angle_start"], row["experiment_angle_end"]
    group = row["resolution"]
    exp_color = GROUP_LIGHT[group][idx % 2]
    add_ring_wedge(ax, start, end, r_outer=1.43, width=0.44, color=exp_color, lw=1.05)

    target_end = start - row["target_angular_extent"]
    add_ring_wedge(ax, start, target_end, r_outer=1.79, width=0.34,
                   color=TARGET_COLOR, lw=0.9)
    if row["off_target_excluded_cells"] > 0:
        add_ring_wedge(ax, target_end, end, r_outer=1.79, width=0.34,
                       color=OFF_COLOR, lw=0.9)

    mid_deg = (start + end) / 2
    mid = np.deg2rad(mid_deg)
    label_points.append({
        "angle": mid,
        "x0": 1.81 * np.cos(mid),
        "y0": 1.81 * np.sin(mid),
        "side": 1 if np.cos(mid) >= 0 else -1,
        "label": f"{row['short_label']}\nTarget n={row['target_excluded_cells']:,}",
        "flagged": row["off_target_excluded_cells"] > 0,
    })

# Distribute external labels separately on each side to prevent overlap.
for side in [-1, 1]:
    points = [p for p in label_points if p["side"] == side]
    points.sort(key=lambda p: p["y0"])
    lower, upper, gap = -1.73, 1.73, 0.225
    ys = [max(lower, min(upper, p["y0"])) for p in points]
    for i in range(1, len(ys)):
        ys[i] = max(ys[i], ys[i - 1] + gap)
    if ys and ys[-1] > upper:
        shift = ys[-1] - upper
        ys = [v - shift for v in ys]
        for i in range(len(ys) - 2, -1, -1):
            ys[i] = min(ys[i], ys[i + 1] - gap)
    x_text = 2.44 * side
    x_elbow = 2.02 * side
    for p, y_text in zip(points, ys):
        ax.plot([p["x0"], x_elbow, x_text - 0.05 * side],
                [p["y0"], y_text, y_text], color="#A8B0B5", linewidth=0.65)
        color = OFF_COLOR if p["flagged"] else TEXT
        ax.text(x_text, y_text, p["label"],
                ha="left" if side == 1 else "right", va="center",
                fontsize=7.0, color=color,
                fontweight="bold" if p["flagged"] else "normal", linespacing=1.05)

# Explicit annotation for the only off-target component.
flag_row = source[source["off_target_excluded_cells"].gt(0)].iloc[0]
off_mid_deg = flag_row["experiment_angle_end"] + flag_row["off_target_angular_extent"] / 2
off_mid = np.deg2rad(off_mid_deg)
xy = (1.66 * np.cos(off_mid), 1.66 * np.sin(off_mid))
ax.annotate("NK cells\n(n=66)", xy=xy, xytext=(1.30, -2.02),
            ha="center", va="top", fontsize=7.5, color=OFF_COLOR, fontweight="bold",
            arrowprops=dict(arrowstyle="-", color=OFF_COLOR, linewidth=0.85,
                            connectionstyle="arc3,rad=-0.15"))

# Center summary.
ax.text(0, 0.13, "15 / 15", ha="center", va="center", fontsize=22,
        color=TEXT, fontweight="bold")
ax.text(0, -0.08, "masked targets", ha="center", va="center", fontsize=9.5, color=TEXT)
ax.text(0, -0.24, "fully excluded", ha="center", va="center", fontsize=9.5,
        color=TARGET_COLOR, fontweight="bold")

ax.legend(handles=[Patch(facecolor=TARGET_COLOR, edgecolor="none", label="Masked target exclusion"),
                   Patch(facecolor=OFF_COLOR, edgecolor="none", label="Additional off-target exclusion")],
          loc="lower center", bbox_to_anchor=(0.5, -0.045), ncol=2,
          frameon=False, fontsize=8, handlelength=1.4, columnspacing=1.6)
ax.set_xlim(-3.25, 3.25)
ax.set_ylim(-2.30, 2.22)
fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.06)

stem = "c5_fig1_stage3a_sunburst"
fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
plt.close(fig)
print("SUCCESS")
