from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualizations" / "c8_direct_paired_comparison"
SOURCE_PATH = OUT / "C8_direct_paired_comparison_source_values.csv"
TABLE_PATH = OUT / "C8_direct_paired_comparison_stats_table.csv"
STATS_PATH = ROOT / "result" / "c8_direct_paired_comparison" / "c8_paired_statistics.csv"

COLORS = {
    "cytospace": "#C66C1B",
    "svtuner": "#006D77",
    "profile": "#7157A6",
    "profile_light": "#E8E2F1",
    "merscope": "#3B8264",
    "merscope_light": "#E1EEE8",
    "neutral": "#A7B0B7",
    "light": "#DDE2E5",
    "grid": "#EFF2F3",
    "row": "#F7F8F8",
    "text": "#202A33",
    "annotation": "#43505A",
    "secondary": "#65717A",
    "negative": "#B86657",
    "white": "#FFFFFF",
}

METRIC_ORDER = [
    "A1_reciprocal_suppression",
    "A2_cosine_similarity",
    "A3_ecotyper_experiment_mean",
    "B1_merscope_suppression",
    "B2_merscope_peak_es",
]

METRIC_BY_ID = {
    "A1_reciprocal_suppression": "Reciprocal suppression",
    "A2_cosine_similarity": "Reconstructed-expression cosine similarity",
    "A3_ecotyper_experiment_mean": "EcoTyper normalized enrichment (experiment mean)",
    "B1_merscope_suppression": "MERSCOPE low-support suppression score",
    "B2_merscope_peak_es": "MERSCOPE peak ES / low-support enrichment",
}

A_TITLES = {
    "A1_reciprocal_suppression": "Reciprocal suppression",
    "A2_cosine_similarity": "Expression cosine similarity",
    "A3_ecotyper_experiment_mean": "EcoTyper enrichment",
}

B_LABELS = {
    "A1_reciprocal_suppression": "Reciprocal\nsuppression",
    "A2_cosine_similarity": "Expression cosine\nsimilarity",
    "A3_ecotyper_experiment_mean": "EcoTyper\nenrichment",
    "B1_merscope_suppression": "MERSCOPE\nlow-support\nsuppression",
    "B2_merscope_peak_es": "MERSCOPE\nPeak ES",
}

MERSCOPE_LABELS = {
    "highres_humanbreastcancerpatient1_profile_mask_monocytes_and_macrophages": "Breast cancer P1 | Mono./macro.",
    "highres_humancoloncancerpatient1_profile_mask_fibroblasts": "Colon cancer P1 | Fibroblasts",
    "highres_humanlungcancerpatient1_profile_mask_plasma_cells": "Lung cancer P1 | Plasma cells",
    "highres_humanmelanomapatient1_profile_mask_fibroblasts": "Melanoma P1 | Fibroblasts",
    "highres_humanmelanomapatient2_profile_mask_b_cells": "Melanoma P2 | B cells",
}

ECOTYPER_LABELS = {
    "brca_er_her2_fresh_frozen": "BRCA ER/HER2 FF",
    "brca_her2_ffpe": "BRCA HER2 FFPE",
    "brca_tnbc_fresh_frozen": "BRCA TNBC FF",
    "crc_fresh_frozen": "CRC FF",
    "melanoma_slide1": "Melanoma S1",
    "melanoma_slide2": "Melanoma S2",
}


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7.5,
        "axes.titlesize": 8.3,
        "axes.labelsize": 7.4,
        "xtick.labelsize": 6.8,
        "ytick.labelsize": 6.8,
        "legend.fontsize": 6.7,
        "axes.edgecolor": COLORS["light"],
        "axes.linewidth": 0.65,
        "axes.labelcolor": COLORS["text"],
        "xtick.color": COLORS["text"],
        "ytick.color": COLORS["text"],
        "text.color": COLORS["text"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def panel_heading(ax: plt.Axes, letter: str, title: str, y: float = 1.16) -> None:
    ax.text(-0.02, y, letter, transform=ax.transAxes, ha="right", va="top",
            fontsize=11.5, fontweight="bold", color=COLORS["text"])
    ax.text(0.02, y, title, transform=ax.transAxes, ha="left", va="top",
            fontsize=8.7, fontweight="semibold", color=COLORS["text"])


def clean_axis(ax: plt.Axes, xgrid: bool = False, ygrid: bool = True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(COLORS["light"])
    ax.spines["bottom"].set_color(COLORS["light"])
    if xgrid:
        ax.grid(axis="x", color=COLORS["grid"], linewidth=0.55, zorder=0)
    if ygrid:
        ax.grid(axis="y", color=COLORS["grid"], linewidth=0.55, zorder=0)
    ax.tick_params(length=2.3, width=0.55)


def paired_values(source: pd.DataFrame, panel: str, metric: str) -> pd.DataFrame:
    part = source.loc[(source["panel"] == panel) & (source["metric"] == metric)].copy()
    wide = part.pivot(index="independent_unit", columns="method", values="value")
    delta = part.drop_duplicates("independent_unit").set_index("independent_unit")["delta"]
    wide["delta"] = delta
    return wide.reset_index()


def draw_panel_a(fig: plt.Figure, spec, source: pd.DataFrame, stats: pd.DataFrame) -> None:
    gs = spec.subgridspec(2, 3, height_ratios=[0.16, 1.0], hspace=0.10, wspace=0.34)
    header = fig.add_subplot(gs[0, :])
    header.axis("off")
    header.text(-0.015, 0.75, "A", ha="right", va="center", fontsize=11.5,
                fontweight="bold", color=COLORS["text"])
    header.text(0.015, 0.75, "Low-resolution paired performance", ha="left", va="center",
                fontsize=8.7, fontweight="semibold", color=COLORS["text"])
    axes = []
    for i, metric_id in enumerate(METRIC_ORDER[:3]):
        ax = fig.add_subplot(gs[1, i])
        axes.append(ax)
        metric = METRIC_BY_ID[metric_id]
        wide = paired_values(source, "A", metric)
        row = stats.loc[stats["metric_id"] == metric_id].iloc[0]

        values = np.r_[wide["CytoSPACE"].to_numpy(), wide["SVTuner"].to_numpy()]
        data_span = max(np.ptp(values), 0.035 * max(abs(values.mean()), 1.0))
        lo = values.min() - 0.15 * data_span
        hi = values.max() + 0.31 * data_span
        for _, pair in wide.iterrows():
            favorable = pair["delta"] > 0
            ax.plot([0, 1], [pair["CytoSPACE"], pair["SVTuner"]],
                    color=COLORS["profile"] if favorable else COLORS["neutral"],
                    alpha=0.48 if favorable else 0.62,
                    linewidth=0.78 if favorable else 0.92, zorder=1)
        ax.scatter(np.zeros(len(wide)), wide["CytoSPACE"], s=26,
                   color=COLORS["cytospace"], edgecolor="white", linewidth=0.55, zorder=3)
        ax.scatter(np.ones(len(wide)), wide["SVTuner"], s=26,
                   color=COLORS["svtuner"], edgecolor="white", linewidth=0.55, zorder=3)
        ax.set_xlim(-0.30, 1.30)
        ax.set_ylim(lo, hi)
        ax.set_xticks([0, 1], ["CytoSPACE", "SVTuner"])
        ax.set_title(A_TITLES[metric_id], pad=7, fontweight="semibold", color=COLORS["text"])

        line1 = f"Mean Δ  {row.mean_delta_raw:+.4f}   95% CI  [{row.mean_delta_ci_low:.4f}, {row.mean_delta_ci_high:.4f}]"
        line2 = f"P = {row.wilcoxon_p_two_sided:.4g}     Favorable  {int(row.wins)}/{int(row.n_pairs)}"
        ax.text(0.5, 0.985, line1 + "\n" + line2, transform=ax.transAxes,
                ha="center", va="top", fontsize=6.25, linespacing=1.30,
                color=COLORS["annotation"], fontweight="medium")
        clean_axis(ax, ygrid=True)

def draw_panel_b(fig: plt.Figure, spec, stats: pd.DataFrame) -> None:
    gs = spec.subgridspec(1, 2, width_ratios=[1.05, 1.75], wspace=0.045)
    ax = fig.add_subplot(gs[0, 0])
    tx = fig.add_subplot(gs[0, 1], sharey=ax)
    ordered = stats.set_index("metric_id").loc[METRIC_ORDER].reset_index()
    y = np.arange(4, -1, -1)

    for yi in y:
        if yi % 2 == 0:
            ax.axhspan(yi - 0.45, yi + 0.45, color=COLORS["row"], zorder=-2)
            tx.axhspan(yi - 0.45, yi + 0.45, color=COLORS["row"], zorder=-2)

    ax.axvline(0, color=COLORS["secondary"], linewidth=0.85, linestyle=(0, (3, 2)), zorder=0)
    colors = [COLORS["profile"]] * 3 + [COLORS["merscope"]] * 2
    ax.scatter(ordered["cohen_dz"], y, s=36, color=colors,
               edgecolor="white", linewidth=0.65, zorder=4)
    ax.set_yticks(y, [B_LABELS[m] for m in ordered["metric_id"]])
    ax.set_xlim(-0.12, max(1.34, ordered["cohen_dz"].max() + 0.13))
    ax.set_ylim(-0.58, 4.62)
    ax.set_xlabel("Standardized paired effect (Cohen $d_z$)", labelpad=4)
    ax.grid(axis="x", color=COLORS["grid"], linewidth=0.6, zorder=-1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color(COLORS["light"])
    ax.tick_params(axis="y", length=0, pad=7)
    ax.tick_params(axis="x", length=2.4, width=0.55)

    tx.set_xlim(0, 1)
    tx.set_ylim(ax.get_ylim())
    tx.axis("off")
    headers = [(0.02, "Mean Δ [95% CI]", "left"), (0.73, "Favorable", "center"), (0.92, "P value", "center")]
    for x, label, align in headers:
        tx.text(x, 4.58, label, ha=align, va="bottom", fontsize=6.4,
                fontweight="semibold", color=COLORS["annotation"])
    tx.plot([0.01, 0.99], [4.48, 4.48], color=COLORS["light"], linewidth=0.7)
    for yi, row in zip(y, ordered.itertuples(index=False)):
        ci = f"{row.mean_delta_raw:+.4f}  [{row.mean_delta_ci_low:.4f}, {row.mean_delta_ci_high:.4f}]"
        tx.text(0.02, yi, ci, ha="left", va="center", fontsize=6.35, color=COLORS["text"])
        tx.text(0.73, yi, f"{int(row.wins)}/{int(row.n_pairs)}", ha="center", va="center",
                fontsize=6.35, color=COLORS["text"])
        tx.text(0.92, yi, f"{row.wilcoxon_p_two_sided:.4g}", ha="center", va="center",
                fontsize=6.35, color=COLORS["text"])

    panel_heading(ax, "B", "Direct effect summary", y=1.18)
    ax.text(0.02, 1.075, "Profile masking", transform=ax.transAxes, color=COLORS["profile"],
            fontsize=6.3, fontweight="semibold", va="top")
    ax.text(0.36, 1.075, "MERSCOPE", transform=ax.transAxes, color=COLORS["merscope"],
            fontsize=6.3, fontweight="semibold", va="top")


def merscope_label(unit: str) -> str:
    dataset = unit.split(" | ")[0]
    label = MERSCOPE_LABELS.get(dataset, dataset)
    return label.replace(" | ", "\n")


def draw_dumbbell(ax: plt.Axes, wide: pd.DataFrame, title: str, favorable: str) -> None:
    wide = wide.sort_values("delta", ascending=True).reset_index(drop=True)
    y = np.arange(len(wide))
    for yi, pair in wide.iterrows():
        ax.plot([pair["CytoSPACE"], pair["SVTuner"]], [yi, yi],
                color=COLORS["light"], linewidth=1.9, zorder=1)
    ax.scatter(wide["CytoSPACE"], y, s=23, color=COLORS["cytospace"],
               edgecolor="white", linewidth=0.45, zorder=3, label="CytoSPACE")
    ax.scatter(wide["SVTuner"], y, s=23, color=COLORS["svtuner"],
               edgecolor="white", linewidth=0.45, zorder=3, label="SVTuner")
    ax.set_yticks(y, [merscope_label(u) for u in wide["independent_unit"]])
    ax.set_title(title, loc="left", pad=5.5, fontsize=7.4, fontweight="semibold")
    ax.text(1.0, 1.01, favorable, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=6.25, color=COLORS["merscope"], fontweight="bold")
    ax.grid(axis="x", color=COLORS["grid"], linewidth=0.55, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color(COLORS["light"])
    ax.tick_params(axis="y", length=0, pad=5, labelsize=5.9)
    ax.tick_params(axis="x", length=2.2, width=0.5, labelsize=6.2)


def draw_panel_c(fig: plt.Figure, spec, source: pd.DataFrame) -> None:
    gs = spec.subgridspec(3, 1, height_ratios=[0.19, 1.0, 1.0], hspace=0.56)
    header = fig.add_subplot(gs[0, 0])
    header.axis("off")
    header.text(-0.055, 0.68, "C", ha="right", va="center", fontsize=11.5,
                fontweight="bold", color=COLORS["text"])
    header.text(0.00, 0.68, "MERSCOPE paired outcomes", ha="left", va="center",
                fontsize=8.7, fontweight="semibold", color=COLORS["text"])
    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="none", markersize=4.8,
               markerfacecolor=COLORS["cytospace"], markeredgecolor="white", markeredgewidth=0.4),
        Line2D([0], [0], marker="o", linestyle="none", markersize=4.8,
               markerfacecolor=COLORS["svtuner"], markeredgecolor="white", markeredgewidth=0.4),
    ]
    header.legend(legend_handles, ["CytoSPACE", "SVTuner"], frameon=False, ncol=2,
                  loc="center right", bbox_to_anchor=(1.0, 0.68), borderaxespad=0,
                  handlelength=0.7, handletextpad=0.25, columnspacing=0.65, fontsize=5.9)
    header.set_xlim(0, 1)
    header.set_ylim(0, 1)
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[2, 0])
    b1 = paired_values(source, "C", METRIC_BY_ID["B1_merscope_suppression"])
    b2 = paired_values(source, "C", METRIC_BY_ID["B2_merscope_peak_es"])
    draw_dumbbell(ax1, b1, "Low-support suppression", "5/5 favorable")
    draw_dumbbell(ax2, b2, "Peak ES", "2/5 favorable")
    ax2.set_xlabel("Peak ES", fontsize=6.5, labelpad=2)
def draw_panel_d(fig: plt.Figure, spec, source: pd.DataFrame) -> None:
    gs = spec.subgridspec(4, 1, height_ratios=[0.16, 1.0, 0.18, 0.11], hspace=0.34)
    header = fig.add_subplot(gs[0, 0])
    header.axis("off")
    header.text(-0.055, 0.68, "D", ha="right", va="center", fontsize=11.5,
                fontweight="bold", color=COLORS["text"])
    header.text(0.00, 0.68, "EcoTyper readout differences", ha="left", va="center",
                fontsize=8.7, fontweight="semibold", color=COLORS["text"])
    ax = fig.add_subplot(gs[1, 0])
    part = source.loc[source["panel"] == "D"].copy()
    part["readout"] = part["metric"].str.extract(r"EcoTyper (.+) normalized enrichment", expand=False)
    delta = part.drop_duplicates(["independent_unit", "readout"])
    pivot = delta.pivot(index="independent_unit", columns="readout", values="delta")
    row_order = [k for k in ECOTYPER_LABELS if k in pivot.index]
    pivot = pivot.loc[row_order, ["CD4 T cells", "CD8 T cells"]]
    vals = pivot.to_numpy(float)
    vmax = max(abs(vals.min()), abs(vals.max()))
    cmap = LinearSegmentedColormap.from_list(
        "ecotyper_v2", [COLORS["negative"], "#F6F4F1", COLORS["profile"]]
    )
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im = ax.imshow(vals, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
    ax.set_xticks(np.arange(2), ["CD4 T cells", "CD8 T cells"])
    ax.set_yticks(np.arange(len(pivot)), [ECOTYPER_LABELS[x] for x in pivot.index])
    ax.tick_params(length=0, pad=5)
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            v = vals[i, j]
            rgba = cmap(norm(v))
            luminance = 0.2126 * rgba[0] + 0.7152 * rgba[1] + 0.0722 * rgba[2]
            ax.text(j, i, f"{v:+.3f}", ha="center", va="center", fontsize=6.35,
                    fontweight="medium", color="white" if luminance < 0.60 else COLORS["text"])
    for edge in np.arange(-0.5, vals.shape[0], 1):
        ax.axhline(edge, color="white", linewidth=1.35)
    for edge in np.arange(-0.5, vals.shape[1], 1):
        ax.axvline(edge, color="white", linewidth=1.35)
    for spine in ax.spines.values():
        spine.set_visible(False)

    note = fig.add_subplot(gs[2, 0])
    note.axis("off")
    note.text(0.5, 0.42, "12 readouts shown descriptively;\npaired inference uses 6 experiment means.",
              ha="center", va="center", fontsize=5.8, linespacing=1.15,
              color=COLORS["annotation"])
    cax = fig.add_subplot(gs[3, 0])
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
    cbar.set_label("SVTuner - CytoSPACE (Δ normalized enrichment)", fontsize=6.2, labelpad=3)
    cbar.ax.tick_params(labelsize=5.8, length=2, width=0.45)
    cbar.outline.set_linewidth(0.45)


def main() -> None:
    set_style()
    source = pd.read_csv(SOURCE_PATH)
    frozen_table = pd.read_csv(TABLE_PATH)
    stats = pd.read_csv(STATS_PATH)

    if len(source) != 101 or len(frozen_table) != 5 or len(stats) != 5:
        raise ValueError("Frozen C8 inputs do not match the finalized input dimensions")
    if stats["metric_id"].tolist() != METRIC_ORDER:
        raise ValueError("Frozen C8 metric order changed")

    fig = plt.figure(figsize=(7.35, 8.15), facecolor="white")
    outer = fig.add_gridspec(
        3, 1, height_ratios=[2.08, 1.53, 2.34],
        left=0.115, right=0.975, top=0.958, bottom=0.082, hspace=0.52
    )
    draw_panel_a(fig, outer[0], source, stats)
    draw_panel_b(fig, outer[1], stats)
    bottom = outer[2].subgridspec(1, 2, width_ratios=[1.42, 0.98], wspace=0.43)
    draw_panel_c(fig, bottom[0], source)
    draw_panel_d(fig, bottom[1], source)

    png = OUT / "C8_direct_paired_comparison_main_v2.png"
    pdf = OUT / "C8_direct_paired_comparison_main_v2.pdf"
    svg = OUT / "C8_direct_paired_comparison_main_v2.svg"
    fig.savefig(png, dpi=600, facecolor="white")
    fig.savefig(pdf, facecolor="white")
    fig.savefig(svg, facecolor="white")
    plt.close(fig)

    notes = (
        "C8 figure V2 refinement notes\n"
        "\n"
        "- This revision is visual refinement only.\n"
        "- No statistical result, scientific definition, or panel content was changed.\n"
        "- Panel A was strengthened as the primary visual focus with clearer paired observations and aligned annotations.\n"
        "- Panel B was reformatted as an integrated forest-style effect summary with aligned confidence-interval, favorable-pair, and P-value columns.\n"
        "- Panel C received harmonized labels, spacing, method encoding, and favorable-pair annotations.\n"
        "- Panel D was enlarged and refined with clearer tiles, typography, color scale, and inferential-unit note.\n"
        "- Typography, panel hierarchy, margins, line weights, and color contrast were unified across the composite figure.\n"
        "- B2 remains visible and unchanged; A3 and B1 are not labeled statistically significant.\n"
    )
    (OUT / "C8_direct_paired_comparison_v2_notes.txt").write_text(notes, encoding="utf-8")
    print(png)
    print(pdf)
    print(svg)


if __name__ == "__main__":
    main()
