from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualizations" / "c5_decomposition"
OUT.mkdir(parents=True, exist_ok=True)
STAGE3A = ROOT / "result" / "c5_stage3a_decomposition" / "c5_stage3a_exclusion_by_dataset.csv"
STAGE3B = ROOT / "result" / "c5_decomposition" / "c5_stage3b_decomposition_by_dataset.csv"

TEAL, DARK_TEAL = "#168C8C", "#116B70"
BLUE, ORANGE, RED = "#557A95", "#D48245", "#B95C62"
GRAY, TEXT = "#7A8793", "#263238"
mpl.rcParams.update({"font.family": "Arial", "font.size": 9, "axes.labelsize": 10,
                     "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "axes.linewidth": 0.8,
                     "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
                     "savefig.facecolor": "white"})

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


def clean_axis(ax, grid_axis=None):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#90999F")
    ax.spines["bottom"].set_color("#90999F")
    ax.tick_params(colors=TEXT, width=0.8)
    if grid_axis:
        ax.grid(axis=grid_axis, color="#E8ECEF", linewidth=0.7, zorder=0)
        ax.set_axisbelow(True)


# Figure 1: compact Stage3A diagnostic audit.
a = pd.read_csv(STAGE3A)
a["display_label"] = a["dataset"].map(LABELS)
a["target_exclusion_percent"] = 100 * a["target_exclusion_fraction"].astype(float)
off_bool = a["off_target_exclusion"].astype(str).str.lower().eq("true")
a["off_target_label"] = np.where(off_bool, a["off_target_details"].fillna("Detected"), "None")
a["audit_result"] = np.where(
    a["target_detected_missing"].astype(str).str.lower().eq("true")
    & np.isclose(a["target_exclusion_fraction"].astype(float), 1.0),
    "Detected and fully excluded", "Not fully excluded")
a[["dataset", "target", "resolution", "capacity", "target_total_cells", "target_dropped_cells",
   "target_exclusion_fraction", "target_exclusion_percent", "target_detected_missing", "audit_result",
   "off_target_exclusion", "off_target_details", "off_target_label", "display_label"]].to_csv(
       OUT / "c5_fig1_stage3a_diagnostic_audit_source.csv", index=False)

fig, ax = plt.subplots(figsize=(7.3, 5.35))
y = np.arange(len(a))
ax.set_xlim(-0.55, 1.72)
ax.set_ylim(len(a) - 0.35, -1.15)
for yi in y:
    ax.axhline(yi + 0.5, color="#EDF0F2", linewidth=0.65, zorder=0)
ax.scatter(np.zeros(len(a)), y, marker="s", s=130, color=TEAL, edgecolor="white", linewidth=0.7, zorder=3)
for yi in y:
    ax.text(0, yi, "100%", ha="center", va="center", fontsize=7.2, color="white", fontweight="bold", zorder=4)
off = off_bool.to_numpy()
ax.scatter(np.ones((~off).sum()), y[~off], s=53, facecolor="white", edgecolor="#B8C1C7", linewidth=1.0, zorder=3)
for yi in y[~off]:
    ax.text(1.13, yi, "None", ha="left", va="center", fontsize=8, color=GRAY)
for yi, detail in zip(y[off], a.loc[off, "off_target_label"]):
    ax.scatter([1], [yi], s=62, color=RED, edgecolor="white", linewidth=0.7, zorder=3)
    ax.text(1.13, yi, detail, ha="left", va="center", fontsize=8, color=RED, fontweight="bold")
ax.set_xticks([0, 1], ["Masked target exclusion", "Off-target exclusion"])
ax.xaxis.tick_top()
ax.tick_params(axis="x", length=0, pad=8, labelsize=9, colors=TEXT)
ax.set_yticks(y, a["display_label"])
ax.tick_params(axis="y", length=0, pad=5)
for side in ax.spines.values():
    side.set_visible(False)
ax.text(0.0, -0.11, "15/15 masked targets detected and fully excluded", transform=ax.transAxes,
        ha="left", va="top", fontsize=10, color=DARK_TEAL, fontweight="bold")
fig.subplots_adjust(left=0.40, right=0.97, top=0.91, bottom=0.10)
save(fig, "c5_fig1_stage3a_diagnostic_audit")


# Figure 2: dataset-level horizontal dumbbell plot.
b = pd.read_csv(STAGE3B)
b["display_label"] = b["dataset"].map(LABELS)
b["standalone_percent"] = 100 * b["standalone_fraction"].astype(float)
b["sequential_percent"] = 100 * b["sequential_fraction"].astype(float)
b["delta_pp"] = 100 * b["delta_fraction"].astype(float)
b_plot = pd.concat([b[b["resolution"].eq(r)].sort_values("delta_pp", ascending=False)
                    for r in ["spot-resolution", "MERSCOPE"]], ignore_index=True)
b_plot["plot_order"] = np.arange(len(b_plot))
b_plot[["dataset", "target", "resolution", "capacity", "total_units", "display_label",
        "standalone_withheld", "standalone_fraction", "standalone_percent", "sequential_withheld",
        "sequential_fraction", "sequential_percent", "delta_withheld", "delta_fraction", "delta_pp",
        "mask_intersection", "mask_union", "jaccard", "plot_order"]].to_csv(
            OUT / "c5_fig2_stage3b_horizontal_dumbbell_source.csv", index=False)

fig, ax = plt.subplots(figsize=(8.0, 5.8))
y = np.arange(len(b_plot), dtype=float)
y[10:] += 0.8
max_x = max(b_plot["standalone_percent"].max(), b_plot["sequential_percent"].max())
for yi, (_, row) in zip(y, b_plot.iterrows()):
    x0, x1 = row["standalone_percent"], row["sequential_percent"]
    ax.plot([x0, x1], [yi, yi], color="#AEB7BD", linewidth=1.2, zorder=1)
    ax.scatter(x0, yi, s=29, color="#7B8790", edgecolor="white", linewidth=0.5, zorder=3)
    ax.scatter(x1, yi, s=34, color=TEAL, edgecolor="white", linewidth=0.5, zorder=4)
    ax.text(max(x0, x1) + 0.015 * max_x, yi, f"{row['delta_pp']:+.2f} pp", va="center", fontsize=7.4, color=TEXT)
ax.axhline((y[9] + y[10]) / 2, color="#D7DDE1", linewidth=0.9)
ax.set_yticks(y, b_plot["display_label"])
ax.set_xlabel("Spatial units withheld (%)")
ax.set_xlim(0, max_x * 1.27)
ax.set_ylim(y[-1] + 0.65, -0.65)
spot_mean = 100 * b.loc[b["resolution"].eq("spot-resolution"), "delta_fraction"].mean()
mer_mean = 100 * b.loc[b["resolution"].eq("MERSCOPE"), "delta_fraction"].mean()
ax.text(0.01, 0.985, f"Spot-resolution (n=10)  mean Δ = {spot_mean:+.2f} pp", transform=ax.transAxes,
        ha="left", va="top", color=BLUE, fontsize=8.3, fontweight="bold")
ax.text(0.01, 0.315, f"MERSCOPE (n=5)  mean Δ = {mer_mean:+.2f} pp", transform=ax.transAxes,
        ha="left", va="top", color=ORANGE, fontsize=8.3, fontweight="bold")
ax.legend(handles=[Line2D([0], [0], marker="o", linestyle="none", color="#7B8790", markersize=5.3,
                          label="Standalone Stage3B"),
                   Line2D([0], [0], marker="o", linestyle="none", color=TEAL, markersize=5.5,
                          label="Sequential Stage3B")],
          loc="lower right", frameon=False, fontsize=8, handletextpad=0.35)
clean_axis(ax, "x")
fig.subplots_adjust(left=0.37, right=0.98, bottom=0.11, top=0.98)
save(fig, "c5_fig2_stage3b_horizontal_dumbbell")


# Figure 3: representative MERSCOPE spatial mask transition.
dataset = "highres_humanbreastcancerpatient1_profile_mask_monocytes_and_macrophages"
row = b.loc[b["dataset"].eq(dataset)].iloc[0]
coords_path = ROOT / "data" / "processed" / dataset / "stage1_preprocess" / "exported" / "st_coordinates.csv"
s = pd.read_csv(Path(row["standalone_mask_path"]), usecols=["spot_id", "is_unsupported_region"])
q = pd.read_csv(Path(row["sequential_mask_path"]), usecols=["spot_id", "is_unsupported_region"])
coords = pd.read_csv(coords_path, usecols=["spot_id", "row", "col"])
s["standalone_withheld"] = s["is_unsupported_region"].astype(str).str.lower().eq("true")
q["sequential_withheld"] = q["is_unsupported_region"].astype(str).str.lower().eq("true")
m = coords.merge(s[["spot_id", "standalone_withheld"]], on="spot_id", validate="one_to_one")
m = m.merge(q[["spot_id", "sequential_withheld"]], on="spot_id", validate="one_to_one")
conditions = [~m["standalone_withheld"] & ~m["sequential_withheld"],
              m["standalone_withheld"] & m["sequential_withheld"],
              m["standalone_withheld"] & ~m["sequential_withheld"],
              ~m["standalone_withheld"] & m["sequential_withheld"]]
categories = ["Neither withheld", "Withheld by both", "Standalone-only withheld", "Sequential-only withheld"]
m["mask_transition"] = np.select(conditions, categories, default="Unclassified")
if len(m) != 2561 or int(m["standalone_withheld"].sum()) != 522 or int(m["sequential_withheld"].sum()) != 253:
    raise RuntimeError("Representative mask counts do not match the frozen C5 results")
m[["spot_id", "row", "col", "standalone_withheld", "sequential_withheld", "mask_transition"]].to_csv(
    OUT / "c5_fig3_merscope_breast_mask_transition_source.csv", index=False)

palette = {"Neither withheld": "#D9DEE2", "Withheld by both": "#2A7F83",
           "Standalone-only withheld": "#D98245", "Sequential-only withheld": "#B85C79"}
fig, ax = plt.subplots(figsize=(6.6, 5.6))
for category in categories:
    sub = m[m["mask_transition"].eq(category)]
    ax.scatter(sub["col"], sub["row"], s=8.5 if category == "Neither withheld" else 12,
               color=palette[category], edgecolor="none", alpha=0.78 if category == "Neither withheld" else 0.92,
               label=f"{category} (n={len(sub)})")
ax.set_aspect("equal", adjustable="datalim")
ax.invert_yaxis()
ax.set_xticks([])
ax.set_yticks([])
for side in ax.spines.values():
    side.set_visible(False)
ax.legend(handles=[Patch(facecolor=palette[c], edgecolor="none",
                         label=f"{c} (n={(m['mask_transition'] == c).sum()})") for c in categories],
          loc="upper left", bbox_to_anchor=(1.01, 0.86), frameon=False, fontsize=8, borderaxespad=0)
ax.text(1.01, 0.99, "Standalone withheld: 522 / 2561 (20.38%)\nSequential withheld: 253 / 2561 (9.88%)",
        transform=ax.transAxes, ha="left", va="top", fontsize=8.5, color=TEXT, linespacing=1.5)
ax.text(0.0, 1.015, "Stage3B spatial mask transition after Stage3A reference filtering",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=10, color=TEXT)
fig.subplots_adjust(left=0.02, right=0.72, bottom=0.03, top=0.94)
save(fig, "c5_fig3_merscope_breast_mask_transition")

print("SUCCESS")
