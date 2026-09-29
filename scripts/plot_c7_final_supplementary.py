"""Generate the final minimal C7 supplementary PR curve and metrics table."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
C7_DIR = ROOT / "visualizations" / "bioapp_experiment" / "C7_cta_extended_performance"
OUT_DIR = C7_DIR / "final_supplementary"

PR_CSV = C7_DIR / "c7_pr_curve_points.csv"
FIXED_CSV = C7_DIR / "c7_fixed_precision_summary.csv"
METRICS_CSV = C7_DIR / "c7_metrics_summary.csv"

TEXT = "#27333A"
SECONDARY = "#667078"
TEAL = "#087F8C"
BASELINE = "#A9AFB3"
ACCENT_30 = "#D99520"
ACCENT_25 = "#B76756"
RULE = "#D6DADD"

PR_WIDTH = 4.65
PR_HEIGHT = 3.45
TABLE_WIDTH = 7.08
TABLE_HEIGHT = 3.80


TABLE_ROWS = [
    ("AUROC", "0.8305", "continuous score"),
    ("Average precision", "0.2837", "continuous score"),
    ("Positive prevalence", "0.0710", "continuous score"),
    ("Frozen binary precision", "0.6122", "frozen binary decision"),
    ("Frozen binary recall", "0.4478", "frozen binary decision"),
    ("Frozen binary F1", "0.5172", "frozen binary decision"),
    ("Specificity", "0.9783", "frozen binary decision"),
    ("Recall at precision >= 0.25", "0.6045", "continuous score"),
    ("Recall at precision >= 0.30", "0.5299", "continuous score"),
    (
        "Precision >= 0.50",
        "Not achieved by any finite withheld_score threshold",
        "continuous score",
    ),
]


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8.0,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 7.4,
            "ytick.labelsize": 7.4,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def read_frozen_inputs() -> tuple[pd.DataFrame, pd.Series]:
    pr = pd.read_csv(PR_CSV)
    pd.read_csv(FIXED_CSV)  # Required frozen C7 input; values below remain fixed.
    metrics = pd.read_csv(METRICS_CSV).iloc[0]
    if list(pr.columns) != ["threshold", "precision", "recall"]:
        raise ValueError("Unexpected C7 PR-curve schema")
    checks = {
        "positive_prevalence": 0.07097457627118645,
        "binary_precision": 0.6122448979591837,
        "binary_recall": 0.44776119402985076,
        "binary_f1": 0.5172413793103449,
        "specificity": 0.9783352337514253,
        "average_precision": 0.2836551170002947,
    }
    for field, expected in checks.items():
        if not np.isclose(float(metrics[field]), expected, rtol=0, atol=1e-12):
            raise ValueError(f"Frozen C7 value changed: {field}")
    return pr, metrics


def save_pr_curve(fig: mpl.figure.Figure) -> None:
    stem = OUT_DIR / "C7_Supp_PR_curve"
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")


def draw_pr_curve(pr: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(PR_WIDTH, PR_HEIGHT))
    fig.subplots_adjust(left=0.145, right=0.970, bottom=0.205, top=0.965)

    ax.plot(
        pr["recall"],
        pr["precision"],
        color=TEAL,
        lw=1.65,
        solid_capstyle="round",
        zorder=3,
    )
    ax.axhline(0.0710, color=BASELINE, lw=0.8, ls=(0, (4, 3)), zorder=1)

    p30 = (0.5299, 0.3008)
    p25 = (0.6045, 0.2500)
    ax.scatter(*p30, s=40, color=ACCENT_30, edgecolor="white", lw=0.7, zorder=5)
    ax.scatter(*p25, s=31, color=ACCENT_25, edgecolor="white", lw=0.7, zorder=5)

    ax.annotate(
        "Precision >= 0.30\nP = 0.3008 · R = 0.5299",
        xy=p30,
        xytext=(0.345, 0.445),
        color=ACCENT_30,
        fontsize=7.2,
        fontweight="bold",
        linespacing=1.23,
        ha="left",
        va="center",
        arrowprops=dict(arrowstyle="-", color=ACCENT_30, lw=0.65, alpha=0.72),
    )
    ax.annotate(
        "Precision >= 0.25\nP = 0.2500 · R = 0.6045",
        xy=p25,
        xytext=(0.665, 0.287),
        color=ACCENT_25,
        fontsize=6.9,
        linespacing=1.23,
        ha="left",
        va="center",
        arrowprops=dict(arrowstyle="-", color=ACCENT_25, lw=0.6, alpha=0.70),
    )

    ax.text(
        0.965,
        0.94,
        "AP = 0.2837",
        transform=ax.transAxes,
        ha="right",
        va="top",
        color=TEXT,
        fontsize=7.7,
        fontweight="bold",
    )
    ax.text(
        0.985,
        0.0710 + 0.012,
        "Prevalence = 0.0710",
        transform=ax.get_yaxis_transform(),
        ha="right",
        va="bottom",
        color=SECONDARY,
        fontsize=6.6,
    )

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 0.54)
    ax.set_xlabel("Recall", color=TEXT)
    ax.set_ylabel("Precision", color=TEXT)
    ax.set_xticks(np.linspace(0, 1, 6))
    ax.set_yticks(np.arange(0, 0.6, 0.1))
    ax.tick_params(colors=SECONDARY, width=0.6, length=3)
    ax.grid(axis="y", color="#EBEEF0", lw=0.5, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#7D858A")
    ax.spines["bottom"].set_color("#7D858A")

    fig.text(
        0.56,
        0.045,
        "No finite score threshold achieved precision >= 0.50.",
        color=SECONDARY,
        fontsize=6.7,
        ha="center",
        va="center",
    )
    save_pr_curve(fig)
    plt.close(fig)


def write_metrics_csv() -> None:
    path = OUT_DIR / "C7_Supp_metrics_table.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Metric", "Value", "Evaluation type"])
        writer.writerows(TABLE_ROWS)


def draw_metrics_table() -> None:
    fig, ax = plt.subplots(figsize=(TABLE_WIDTH, TABLE_HEIGHT))
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    left, right = 0.055, 0.955
    col_metric, col_value, col_type = 0.060, 0.430, 0.775
    top, header_bottom, bottom = 0.935, 0.860, 0.080
    row_height = (header_bottom - bottom) / len(TABLE_ROWS)

    ax.hlines(top, left, right, colors=TEXT, lw=1.0)
    ax.hlines(header_bottom, left, right, colors=TEXT, lw=0.75)
    ax.text(col_metric, 0.898, "Metric", color=TEXT, fontsize=8.1, fontweight="bold", va="center")
    ax.text(col_value, 0.898, "Value", color=TEXT, fontsize=8.1, fontweight="bold", va="center")
    ax.text(col_type, 0.898, "Evaluation type", color=TEXT, fontsize=8.1, fontweight="bold", va="center")

    for index, (metric, value, eval_type) in enumerate(TABLE_ROWS):
        y = header_bottom - (index + 0.5) * row_height
        ax.text(col_metric, y, metric, color=TEXT, fontsize=7.6, va="center", ha="left")
        ax.text(
            col_value,
            y,
            value,
            color=TEXT,
            fontsize=6.9 if index == len(TABLE_ROWS) - 1 else 7.6,
            va="center",
            ha="left",
        )
        ax.text(col_type, y, eval_type, color=SECONDARY, fontsize=7.2, va="center", ha="left")
        if index in (2, 6):
            sep_y = header_bottom - (index + 1) * row_height
            ax.hlines(sep_y, left, right, colors=RULE, lw=0.65)

    ax.hlines(bottom, left, right, colors=TEXT, lw=1.0)
    stem = OUT_DIR / "C7_Supp_metrics_table"
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".png"), dpi=600, facecolor="white")
    plt.close(fig)


def validate_csv() -> None:
    table = pd.read_csv(OUT_DIR / "C7_Supp_metrics_table.csv")
    if list(table.columns) != ["Metric", "Value", "Evaluation type"]:
        raise ValueError("Unexpected supplementary-table schema")
    if len(table) != 10 or table.isna().any().any():
        raise ValueError("Supplementary metrics table is incomplete")
    if table.iloc[-1]["Value"] != "Not achieved by any finite withheld_score threshold":
        raise ValueError("Frozen precision >= 0.50 result changed")


def main() -> None:
    configure_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pr, _ = read_frozen_inputs()
    draw_pr_curve(pr)
    write_metrics_csv()
    draw_metrics_table()
    validate_csv()
    print(f"Wrote final C7 supplementary outputs to {OUT_DIR}")


if __name__ == "__main__":
    main()
