from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualizations" / "c5_decomposition"
INPUT = ROOT / "result" / "c5_decomposition" / "c5_stage3b_decomposition_by_dataset.csv"
OUT.mkdir(parents=True, exist_ok=True)

BACKGROUND = "#FFF9F4"
TEXT = "#2F3542"
AUX = "#667085"
BORDER = "#E3DDD7"
COLORS = {
    "Neither withheld": "#F3EEE8",
    "Withheld by both": "#8EA7F8",
    "Standalone-only": "#F3A77E",
    "Sequential-only": "#C7A6E8",
}

LABELS = {
    "adult_mouse_kidney_real_profile_mask_endo": ("Mouse kidney", "Endo"),
    "ffpe_mouse_brain_sagittal_real_profile_mask_microglia": ("Mouse brain FFPE", "Microglia"),
    "human_breast_cancer_real_profile_mask_basal_cell": ("Breast cancer", "Basal cell"),
    "human_breast_cancer_visium_ff_wta_real_profile_mask_macrophage": ("Breast cancer FF", "Macrophage"),
    "human_breast_cancer_wta_120_real_profile_mask_endothelial_cell": ("Breast cancer FFPE", "Endothelial"),
    "human_cervical_cancer_real_profile_mask_epithelial_cell": ("Cervical cancer", "Epithelial"),
    "human_heart_ff_real_profile_mask_endothelial_cell": ("Human heart", "Endothelial"),
    "human_intestine_cancer_real_profile_mask_endothelial_cell": ("Intestine cancer", "Endothelial"),
    "human_lymph_node_real_profile_mask_b_cell": ("Lymph node", "B cell"),
    "mouse_embryo_real_profile_mask_erythroid": ("Mouse embryo", "Erythroid"),
    "highres_humanbreastcancerpatient1_profile_mask_monocytes_and_macrophages": ("Breast cancer P1", "Mono./macro."),
    "highres_humancoloncancerpatient1_profile_mask_fibroblasts": ("Colon cancer P1", "Fibroblasts"),
    "highres_humanlungcancerpatient1_profile_mask_plasma_cells": ("Lung cancer P1", "Plasma cells"),
    "highres_humanmelanomapatient1_profile_mask_fibroblasts": ("Melanoma P1", "Fibroblasts"),
    "highres_humanmelanomapatient2_profile_mask_b_cells": ("Melanoma P2", "B cells"),
}

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 9,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "savefig.facecolor": BACKGROUND,
})


def as_bool(series):
    return series.astype(str).str.lower().eq("true")


summary = pd.read_csv(INPUT)
rows = []
for _, record in summary.iterrows():
    standalone = pd.read_csv(Path(record["standalone_mask_path"]),
                             usecols=["spot_id", "is_unsupported_region"])
    sequential = pd.read_csv(Path(record["sequential_mask_path"]),
                             usecols=["spot_id", "is_unsupported_region"])
    standalone["standalone"] = as_bool(standalone["is_unsupported_region"])
    sequential["sequential"] = as_bool(sequential["is_unsupported_region"])
    joined = standalone[["spot_id", "standalone"]].merge(
        sequential[["spot_id", "sequential"]], on="spot_id", how="inner", validate="one_to_one")
    total = len(joined)
    expected = int(float(record["total_units"]))
    if total != expected or total != len(standalone) or total != len(sequential):
        raise RuntimeError(f"Mask ID alignment failed for {record['dataset']}")
    neither = int((~joined["standalone"] & ~joined["sequential"]).sum())
    both = int((joined["standalone"] & joined["sequential"]).sum())
    standalone_only = int((joined["standalone"] & ~joined["sequential"]).sum())
    sequential_only = int((~joined["standalone"] & joined["sequential"]).sum())
    if neither + both + standalone_only + sequential_only != total:
        raise RuntimeError(f"Transition counts do not sum for {record['dataset']}")
    rows.append({
        "dataset": record["dataset"], "target": record["target"],
        "resolution": record["resolution"], "capacity": int(record["capacity"]),
        "total_units": total, "neither_n": neither, "both_n": both,
        "standalone_only_n": standalone_only, "sequential_only_n": sequential_only,
        "neither_fraction": neither / total, "both_fraction": both / total,
        "standalone_only_fraction": standalone_only / total,
        "sequential_only_fraction": sequential_only / total,
        "jaccard": float(record["jaccard"]),
        "short_dataset": LABELS[record["dataset"]][0],
        "short_target": LABELS[record["dataset"]][1],
    })

source = pd.DataFrame(rows)
source.to_csv(OUT / "c5_fig2_stage3b_transition_fingerprint_source.csv", index=False)


def draw_fingerprint(ax, x, y, size, row):
    stable = row["neither_fraction"] + row["both_fraction"]
    transition = row["standalone_only_fraction"] + row["sequential_only_fraction"]

    if stable > 0:
        stable_w = size * stable
        neither_h = size * row["neither_fraction"] / stable
        both_h = size * row["both_fraction"] / stable
        if row["neither_fraction"] > 0:
            ax.add_patch(Rectangle((x, y + both_h), stable_w, neither_h,
                                   facecolor=COLORS["Neither withheld"], edgecolor=BACKGROUND,
                                   linewidth=1.0))
        if row["both_fraction"] > 0:
            ax.add_patch(Rectangle((x, y), stable_w, both_h,
                                   facecolor=COLORS["Withheld by both"], edgecolor=BACKGROUND,
                                   linewidth=1.0))

    if transition > 0:
        trans_x = x + size * stable
        trans_w = size * transition
        standalone_h = size * row["standalone_only_fraction"] / transition
        sequential_h = size * row["sequential_only_fraction"] / transition
        if row["standalone_only_fraction"] > 0:
            ax.add_patch(Rectangle((trans_x, y + sequential_h), trans_w, standalone_h,
                                   facecolor=COLORS["Standalone-only"], edgecolor=BACKGROUND,
                                   linewidth=1.0))
        if row["sequential_only_fraction"] > 0:
            ax.add_patch(Rectangle((trans_x, y), trans_w, sequential_h,
                                   facecolor=COLORS["Sequential-only"], edgecolor=BACKGROUND,
                                   linewidth=1.0))

    ax.add_patch(Rectangle((x, y), size, size, facecolor="none", edgecolor=BORDER,
                           linewidth=0.9))
    ax.text(x + size / 2, y + size + 0.20, row["short_dataset"], ha="center", va="bottom",
            fontsize=8.0, color=TEXT, fontweight="bold")
    ax.text(x + size / 2, y + size + 0.07, row["short_target"], ha="center", va="bottom",
            fontsize=7.4, color=AUX)
    ax.text(x + size / 2, y - 0.10, f"J = {row['jaccard']:.2f}", ha="center", va="top",
            fontsize=7.0, color=AUX)


fig, ax = plt.subplots(figsize=(10.8, 7.2), facecolor=BACKGROUND)
ax.set_facecolor(BACKGROUND)
ax.axis("off")

tile = 1.15
x_positions = [0.25 + 1.62 * i for i in range(5)]
y_positions_spot = [4.15, 2.35]
y_merscope = 0.22

spot = source[source["resolution"].eq("spot-resolution")].reset_index(drop=True)
merscope = source[source["resolution"].eq("MERSCOPE")].reset_index(drop=True)
for i, (_, row) in enumerate(spot.iterrows()):
    draw_fingerprint(ax, x_positions[i % 5], y_positions_spot[i // 5], tile, row)
for i, (_, row) in enumerate(merscope.iterrows()):
    draw_fingerprint(ax, x_positions[i], y_merscope, tile, row)

ax.text(0.25, 5.88, "Spot-resolution (n=10)", ha="left", va="center",
        fontsize=10, color="#557A95", fontweight="bold")
ax.text(0.25, 1.92, "MERSCOPE (n=5; capacity=1)", ha="left", va="center",
        fontsize=10, color="#9B7EC4", fontweight="bold")

handles = [Patch(facecolor=COLORS[name], edgecolor=BORDER, linewidth=0.5, label=name)
           for name in ["Neither withheld", "Withheld by both", "Standalone-only", "Sequential-only"]]
ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.035), ncol=4,
          frameon=False, fontsize=8.2, handlelength=1.35, columnspacing=1.5)
ax.text(7.95, -0.12, "Tile area represents the fraction of spatial units.",
        ha="right", va="top", fontsize=7.4, color=AUX)
ax.set_xlim(0, 8.25)
ax.set_ylim(-0.42, 6.08)
fig.subplots_adjust(left=0.025, right=0.985, top=0.99, bottom=0.08)

stem = "c5_fig2_stage3b_transition_fingerprint"
fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor=BACKGROUND)
fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor=BACKGROUND)
fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight", facecolor=BACKGROUND)
plt.close(fig)
print("SUCCESS")
