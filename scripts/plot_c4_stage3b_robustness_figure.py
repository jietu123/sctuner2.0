#!/usr/bin/env python
"""Create the frozen C4 Stage3B supplementary robustness figure."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
CALIBRATION_DIR = ROOT / "result" / "c4_stage3b_calibration"
PERTURBATION_CSV = (
    ROOT
    / "result"
    / "c4_stage3b_technical_perturbation"
    / "c4_2_technical_perturbation_pilot_summary.csv"
)
ROBUSTNESS_DIR = ROOT / "result" / "c4_stage3b_statistical_robustness"
OUTPUT_DIR = ROOT / "visualizations" / "c4_stage3b_calibration"

CASES = (
    (
        "cytospace_fig2d_tme_brca_tnbc_fresh_frozen_sc_missing_plasma_cells",
        "BRCA TNBC\nPlasma dropout",
        "BRCA TNBC",
    ),
    (
        "cytospace_fig2d_tme_crc_fresh_frozen_sc_missing_b_cells",
        "CRC\nB-cell dropout",
        "CRC",
    ),
    (
        "cytospace_fig2d_tme_brca_her2_ffpe_sc_missing_plasma_cells",
        "BRCA HER2 FFPE\nPlasma dropout",
        "BRCA HER2\nFFPE",
    ),
)

FEATURES = (
    ("relative_reconstruction_error", "Reconstruction\nerror"),
    ("cosine_deficit", "Cosine\ndeficit"),
    ("positive_residual_fraction", "Positive residual\nfraction"),
    ("residual_concentration", "Residual\nconcentration"),
)

TEAL = "#177E83"
TEAL_LIGHT = "#77B8B7"
BLUE_GRAY = "#667D8A"
LIGHT_GRAY = "#B9C0C3"
GRID = "#E7EAEB"
TEXT = "#273238"


def ecdf(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.sort(np.asarray(values, dtype=float))
    y = np.arange(1, len(x) + 1, dtype=float) / len(x)
    return x, y


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#7A858A")
    ax.spines["bottom"].set_color("#7A858A")
    ax.spines["left"].set_linewidth(0.65)
    ax.spines["bottom"].set_linewidth(0.65)
    ax.tick_params(colors=TEXT, width=0.65, length=2.5)
    ax.grid(axis="y", color=GRID, linewidth=0.55, zorder=0)


def load_calibration(sample: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    observed_path = (
        ROOT
        / "data"
        / "processed"
        / "stage3b_realdata_reference_dropout"
        / sample
        / "stage3b_st_unsupported"
        / "spot_unsupported_scores.csv"
    )
    null_path = CALIBRATION_DIR / sample / "c4_1_raw_null_statistics.csv"
    observed = pd.read_csv(observed_path)
    null = pd.read_csv(null_path)
    columns = [feature for feature, _ in FEATURES]
    observed = observed.loc[:, columns].apply(pd.to_numeric, errors="raise")
    null = null.loc[:, columns].apply(pd.to_numeric, errors="raise")
    if observed.isna().any().any() or null.isna().any().any():
        raise ValueError(f"NaN in C4.1 calibration data for {sample}")
    return observed, null


def load_robustness(sample: str) -> pd.DataFrame:
    if sample == CASES[0][0]:
        path = ROBUSTNESS_DIR / "c4_3_bh_by_summary.csv"
    else:
        path = ROBUSTNESS_DIR / sample / "c4_3_bh_by_summary.csv"
    frame = pd.read_csv(path).set_index("method")
    if set(frame.index) != {"BH", "BY"}:
        raise ValueError(f"Unexpected C4.3 methods in {path}")
    return frame


def add_panel_heading(fig: plt.Figure, x: float, letter: str, title: str) -> None:
    fig.text(x, 0.968, letter, fontsize=10.5, fontweight="bold", color=TEXT, va="top")
    fig.text(x + 0.014, 0.968, title, fontsize=9.2, color=TEXT, va="top")


def main() -> int:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7.5,
            "axes.labelsize": 7.5,
            "axes.titlesize": 8.0,
            "xtick.labelsize": 6.7,
            "ytick.labelsize": 6.7,
            "legend.fontsize": 6.6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    calibration = {sample: load_calibration(sample) for sample, _, _ in CASES}
    perturbation = pd.read_csv(PERTURBATION_CSV).set_index("condition")
    robustness = {sample: load_robustness(sample) for sample, _, _ in CASES}

    fig = plt.figure(figsize=(13.8, 5.15), dpi=300, facecolor="white")
    outer = fig.add_gridspec(
        1,
        3,
        width_ratios=(2.55, 0.93, 1.12),
        left=0.052,
        right=0.985,
        bottom=0.16,
        top=0.89,
        wspace=0.28,
    )

    panel_a = outer[0, 0].subgridspec(3, 4, wspace=0.34, hspace=0.34)
    for row, (sample, row_label, _) in enumerate(CASES):
        observed, null = calibration[sample]
        for column, (feature, feature_label) in enumerate(FEATURES):
            ax = fig.add_subplot(panel_a[row, column])
            null_x, null_y = ecdf(null[feature].to_numpy(dtype=float))
            observed_x, observed_y = ecdf(observed[feature].to_numpy(dtype=float))
            ax.step(
                null_x,
                null_y,
                where="post",
                color=LIGHT_GRAY,
                linewidth=0.9,
                label="Pseudo-null",
                zorder=2,
            )
            ax.step(
                observed_x,
                observed_y,
                where="post",
                color=TEAL,
                linewidth=1.15,
                label="Observed",
                zorder=3,
            )
            ax.set_ylim(0, 1.01)
            ax.set_yticks((0, 0.5, 1.0))
            if column:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel(f"{row_label}\nECDF", labelpad=4.0)
            if row == 0:
                ax.set_title(feature_label, pad=3.5, fontweight="normal")
            if row == len(CASES) - 1:
                ax.set_xlabel("Statistic value", labelpad=2.5)
            style_axis(ax)
            if row == 0 and column == 0:
                ax.legend(
                    loc="lower right",
                    frameon=False,
                    handlelength=1.5,
                    handletextpad=0.45,
                    borderaxespad=0.25,
                )

    ax_b = fig.add_subplot(outer[0, 1])
    b_order = [
        "unperturbed_control",
        "library_size_75",
        "library_size_50",
        "random_dropout_10",
        "random_dropout_20",
    ]
    b_labels = [
        "Baseline",
        "Lib. size\n75%",
        "Lib. size\n50%",
        "Dropout\n10%",
        "Dropout\n20%",
    ]
    b_values = perturbation.loc[b_order, "withheld_rate"].to_numpy(dtype=float) * 100.0
    x_b = np.arange(len(b_order))
    colors = [BLUE_GRAY, TEAL_LIGHT, TEAL_LIGHT, TEAL_LIGHT, TEAL_LIGHT]
    ax_b.bar(x_b, b_values, width=0.62, color=colors, edgecolor="white", linewidth=0.4, zorder=2)
    ax_b.scatter(x_b, b_values, s=13, color=[BLUE_GRAY, TEAL, TEAL, TEAL, TEAL], zorder=3)
    for x, value in zip(x_b, b_values):
        ax_b.text(x, value + 0.22, f"{value:.2f}%", ha="center", va="bottom", fontsize=6.3, color=TEXT)
    ax_b.set_ylabel("Withheld rate (%)")
    ax_b.set_xticks(x_b, b_labels)
    ax_b.set_ylim(0, 8.0)
    ax_b.set_yticks(np.arange(0, 9, 2))
    style_axis(ax_b)

    ax_c = fig.add_subplot(outer[0, 2])
    x_c = np.arange(len(CASES))
    bh = np.asarray(
        [robustness[sample].loc["BH", "significant_fraction"] * 100 for sample, _, _ in CASES]
    )
    by = np.asarray(
        [robustness[sample].loc["BY", "significant_fraction"] * 100 for sample, _, _ in CASES]
    )
    width = 0.31
    ax_c.bar(x_c - width / 2, bh, width=width, color=BLUE_GRAY, label="BH", zorder=2)
    ax_c.bar(x_c + width / 2, by, width=width, color=TEAL, label="BY", zorder=2)
    for positions, values in ((x_c - width / 2, bh), (x_c + width / 2, by)):
        for x, value in zip(positions, values):
            ax_c.text(x, value + 0.34, f"{value:.2f}%", ha="center", va="bottom", fontsize=6.2, color=TEXT)
    ax_c.set_ylabel("Candidate rate (%)")
    ax_c.set_xticks(x_c, [case[2] for case in CASES])
    ax_c.set_ylim(0, 25.5)
    ax_c.set_yticks(np.arange(0, 26, 5))
    ax_c.legend(frameon=False, loc="upper left", ncol=2, handlelength=1.2, columnspacing=1.0)
    ax_c.text(
        0.98,
        0.965,
        "200 $\\rightarrow$ 1000 spatial permutations:\n"
        "final withheld-mask Jaccard = 1.000\nin all three datasets",
        transform=ax_c.transAxes,
        ha="right",
        va="top",
        fontsize=6.2,
        color="#5D686D",
        linespacing=1.25,
    )
    style_axis(ax_c)

    add_panel_heading(fig, 0.010, "A", "Observed vs pseudo-null calibration")
    add_panel_heading(fig, 0.600, "B", "Technical perturbation control")
    add_panel_heading(fig, 0.790, "C", "Statistical robustness")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    pdf_path = OUTPUT_DIR / "c4_stage3b_robustness_figure.pdf"
    png_path = OUTPUT_DIR / "c4_stage3b_robustness_figure.png"
    fig.savefig(pdf_path, facecolor="white")
    fig.savefig(png_path, dpi=300, facecolor="white")
    plt.close(fig)
    print(pdf_path.resolve())
    print(png_path.resolve())
    print("figure_size_inches=13.8x5.15")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
