#!/usr/bin/env python3
"""Create the two frozen C7 CTA performance figure drafts."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, PathPatch, Rectangle
from matplotlib.path import Path as MplPath
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "visualizations" / "bioapp_experiment" / "C7_cta_extended_performance"
OUT_DIR = DATA_DIR / "figures"
FIGSIZE = (7.6, 4.35)
BACKGROUND = "#FFF9F4"
TEXT = "#2F3542"
SECONDARY = "#667085"
GRID = "#DED8D1"


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.5,
            "axes.labelcolor": TEXT,
            "axes.edgecolor": "#9AA3AE",
            "xtick.color": SECONDARY,
            "ytick.color": SECONDARY,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": BACKGROUND,
        }
    )


def export(fig: plt.Figure, stem: str) -> list[Path]:
    paths = []
    for suffix, kwargs in (
        (".png", {"dpi": 300}),
        (".pdf", {}),
        (".svg", {}),
    ):
        path = OUT_DIR / f"{stem}{suffix}"
        fig.savefig(path, facecolor=BACKGROUND, **kwargs)
        paths.append(path)
    return paths


def label_point(ax, x, y, text, xytext, color, marker, size=48) -> None:
    ax.scatter(
        [x], [y], s=size, marker=marker, color=color, edgecolor=BACKGROUND,
        linewidth=1.0, zorder=7,
    )
    ax.annotate(
        text,
        xy=(x, y),
        xytext=xytext,
        textcoords="data",
        color=TEXT,
        fontsize=7.2,
        linespacing=1.25,
        ha="left",
        va="center",
        arrowprops={"arrowstyle": "-", "color": color, "lw": 0.9, "alpha": 0.9},
        bbox={"boxstyle": "round,pad=0.28", "fc": "#FFFCF8", "ec": color, "lw": 0.55},
        zorder=8,
    )


def figure_pr(summary: pd.Series, pr: pd.DataFrame, fixed: pd.DataFrame) -> list[Path]:
    fig = plt.figure(figsize=FIGSIZE, facecolor=BACKGROUND)
    grid = fig.add_gridspec(
        1, 2, width_ratios=[3.75, 1.45], left=0.085, right=0.975,
        bottom=0.15, top=0.83, wspace=0.17,
    )
    ax = fig.add_subplot(grid[0, 0])
    info = fig.add_subplot(grid[0, 1])
    ax.set_facecolor(BACKGROUND)
    info.set_facecolor(BACKGROUND)

    recall = pr["recall"].to_numpy(float)
    precision = pr["precision"].to_numpy(float)
    order = np.argsort(recall, kind="stable")
    curve_color = "#187E82"
    ax.fill_between(
        recall[order], precision[order], summary["positive_prevalence"],
        color="#88CAC5", alpha=0.10, zorder=1,
    )
    ax.plot(
        recall, precision, color=curve_color, lw=1.9, solid_capstyle="round",
        label="Continuous withheld score", zorder=3,
    )
    prevalence = float(summary["positive_prevalence"])
    ax.axhline(prevalence, color="#A7A6A4", lw=1.0, ls=(0, (3, 3)), zorder=0)
    ax.text(
        0.985, prevalence + 0.012, f"Prevalence = {prevalence:.4f}",
        color="#777A80", fontsize=7.2, ha="right", va="bottom",
    )

    label_point(
        ax, 0.6044776119, 0.2500,
        "P >= 0.25\nR = 0.6045  |  t = 0.5017",
        (0.655, 0.205), "#F0A35E", "o", 42,
    )
    label_point(
        ax, 0.5298507463, 0.3008474576,
        "P >= 0.30\nR = 0.5299  |  t = 0.5281",
        (0.585, 0.365), "#6D8FD5", "s", 42,
    )
    label_point(
        ax, 0.2089552239, 0.4912280702,
        "Maximum finite-threshold precision\nP = 0.4912  |  R = 0.2090\nt = 0.6233",
        (0.045, 0.565), "#9A73C6", "h", 52,
    )

    frozen_precision = float(summary["binary_precision"])
    frozen_recall = float(summary["binary_recall"])
    ax.scatter(
        [frozen_recall], [frozen_precision], s=92, marker="D", color="#D64F70",
        edgecolor="#FFFDF9", linewidth=1.3, zorder=9,
    )
    ax.annotate(
        "Frozen composite decision\nP = 0.6122  |  R = 0.4478\nF1 = 0.5172",
        xy=(frozen_recall, frozen_precision), xytext=(0.54, 0.625),
        color=TEXT, fontsize=7.5, fontweight="bold", linespacing=1.25,
        ha="left", va="center",
        arrowprops={"arrowstyle": "-", "color": "#D64F70", "lw": 1.15},
        bbox={"boxstyle": "round,pad=0.34", "fc": "#FFF4F6", "ec": "#D64F70", "lw": 0.8},
        zorder=10,
    )

    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 0.68)
    ax.set_xlabel("Recall", fontsize=9)
    ax.set_ylabel("Precision", fontsize=9)
    ax.set_xticks(np.linspace(0, 1, 6))
    ax.set_yticks(np.arange(0, 0.7, 0.1))
    ax.grid(axis="y", color=GRID, lw=0.6, alpha=0.65)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(0.8)
    ax.tick_params(labelsize=7.5, length=3)

    info.axis("off")
    panel = FancyBboxPatch(
        (0.01, 0.03), 0.98, 0.94,
        boxstyle="round,pad=0.022,rounding_size=0.025",
        transform=info.transAxes, facecolor="#FFFDF9", edgecolor="#E4DCD5", lw=0.8,
    )
    info.add_patch(panel)
    info.text(0.10, 0.90, "Continuous score", transform=info.transAxes,
              color=SECONDARY, fontsize=7.5, fontweight="bold")
    info.text(0.10, 0.815, f"AP  {float(summary['average_precision']):.4f}",
              transform=info.transAxes, color=curve_color, fontsize=13, fontweight="bold")
    info.text(0.10, 0.755, f"Prevalence baseline  {prevalence:.4f}",
              transform=info.transAxes, color=TEXT, fontsize=7.4)
    info.plot([0.10, 0.90], [0.69, 0.69], transform=info.transAxes, color="#E6DED7", lw=0.7)

    info.text(0.10, 0.62, "Frozen composite decision", transform=info.transAxes,
              color="#B93F60", fontsize=7.5, fontweight="bold")
    info.text(0.10, 0.535, "P  0.6122     R  0.4478", transform=info.transAxes,
              color=TEXT, fontsize=8.6, fontweight="bold")
    info.text(0.10, 0.475, "F1  0.5172", transform=info.transAxes,
              color=TEXT, fontsize=8.6, fontweight="bold")

    info.plot([0.10, 0.90], [0.41, 0.41], transform=info.transAxes, color="#E6DED7", lw=0.7)
    info.text(0.10, 0.365, "Not a one-dimensional\nthreshold", transform=info.transAxes,
              color="#B93F60", fontsize=7.15, fontweight="bold", linespacing=1.15)
    info.text(
        0.10, 0.235,
        "No rule of the form\nwithheld score >= t reproduces\nthe frozen decision.",
        transform=info.transAxes, color=TEXT, fontsize=7.25, linespacing=1.35,
    )
    info.text(
        0.10, 0.145, "Best match still differs\nin 93 binary decisions.",
        transform=info.transAxes, color="#B93F60", fontsize=7.3, fontweight="bold",
        linespacing=1.3,
    )
    info.text(
        0.10, 0.055, "Finite thresholds: P >= 0.50,\n0.70, and 0.80 not achieved.",
        transform=info.transAxes, color=SECONDARY, fontsize=6.7, linespacing=1.25,
    )

    fig.text(0.085, 0.94, "CTA endpoint discrimination", color=TEXT,
             fontsize=13, fontweight="bold", ha="left")
    fig.text(
        0.085, 0.885,
        "Continuous withheld score and the frozen composite Stage3B decision",
        color=SECONDARY, fontsize=8.5, ha="left",
    )
    return export(fig, "C7_Fig1_PR_performance_panel")


def flow_patch(ax, x0, x1, source, target, color, alpha, zorder) -> None:
    s_top, s_bottom = source
    t_top, t_bottom = target
    dx = (x1 - x0) * 0.43
    vertices = [
        (x0, s_top), (x0 + dx, s_top), (x1 - dx, t_top), (x1, t_top),
        (x1, t_bottom), (x1 - dx, t_bottom), (x0 + dx, s_bottom), (x0, s_bottom),
        (x0, s_top),
    ]
    codes = [
        MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
        MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
        MplPath.CLOSEPOLY,
    ]
    ax.add_patch(PathPatch(MplPath(vertices, codes), facecolor=color, edgecolor="none",
                           alpha=alpha, zorder=zorder))


def metric_badge(ax, x, label, value, color) -> None:
    width, height = 0.165, 0.105
    badge = FancyBboxPatch(
        (x, 0.015), width, height, boxstyle="round,pad=0.014,rounding_size=0.016",
        transform=ax.transAxes, facecolor="#FFFDF9", edgecolor="#E2D9D1", lw=0.75,
        zorder=10,
    )
    ax.add_patch(badge)
    ax.text(x + 0.018, 0.087, label, transform=ax.transAxes, color=SECONDARY,
            fontsize=6.7, fontweight="bold", va="center", zorder=11)
    ax.text(x + 0.018, 0.042, value, transform=ax.transAxes, color=color,
            fontsize=10.1, fontweight="bold", va="center", zorder=11)


def figure_decomposition(summary: pd.Series) -> list[Path]:
    tp, fp, tn, fn = (int(summary[key]) for key in ("tp", "fp", "tn", "fn"))
    total = tp + fp + tn + fn
    positive, negative = tp + fn, tn + fp
    withheld, retained = tp + fp, tn + fn

    fig, ax = plt.subplots(figsize=FIGSIZE, facecolor=BACKGROUND)
    ax.set_facecolor(BACKGROUND)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    x0, x1 = 0.22, 0.73
    node_width = 0.028
    scale = 0.56 / total
    left_top = 0.74
    right_top = 0.74
    gap = 0.055

    pos_interval = (left_top, left_top - positive * scale)
    neg_top = pos_interval[1] - gap
    neg_interval = (neg_top, neg_top - negative * scale)
    withheld_interval = (right_top, right_top - withheld * scale)
    retained_top = withheld_interval[1] - gap
    retained_interval = (retained_top, retained_top - retained * scale)

    pos_withheld_source = (pos_interval[0], pos_interval[0] - tp * scale)
    pos_retained_source = (pos_withheld_source[1], pos_interval[1])
    neg_withheld_source = (neg_interval[0], neg_interval[0] - fp * scale)
    neg_retained_source = (neg_withheld_source[1], neg_interval[1])
    pos_withheld_target = (withheld_interval[0], withheld_interval[0] - tp * scale)
    neg_withheld_target = (pos_withheld_target[1], withheld_interval[1])
    pos_retained_target = (retained_interval[0], retained_interval[0] - fn * scale)
    neg_retained_target = (pos_retained_target[1], retained_interval[1])

    # Large retained flow first; informative minority flows remain on top.
    flow_patch(ax, x0 + node_width, x1, neg_retained_source, neg_retained_target,
               "#AFC2E0", 0.47, 1)
    flow_patch(ax, x0 + node_width, x1, pos_retained_source, pos_retained_target,
               "#F3BA98", 0.50, 2)
    flow_patch(ax, x0 + node_width, x1, neg_withheld_source, neg_withheld_target,
               "#6D8FD5", 0.82, 4)
    flow_patch(ax, x0 + node_width, x1, pos_withheld_source, pos_withheld_target,
               "#D95776", 0.88, 5)

    node_specs = [
        (x0, pos_interval, "#E98477"),
        (x0, neg_interval, "#7897D3"),
        (x1, withheld_interval, "#C74667"),
        (x1, retained_interval, "#BBCBE0"),
    ]
    for x, interval, color in node_specs:
        top, bottom = interval
        ax.add_patch(Rectangle((x, bottom), node_width, top - bottom,
                               facecolor=color, edgecolor="#FFFDF9", lw=0.7, zorder=7))

    ax.text(x0 - 0.025, np.mean(pos_interval), "CTA-positive\n134", ha="right", va="center",
            color="#B95248", fontsize=8.3, fontweight="bold", linespacing=1.25)
    ax.text(x0 - 0.025, np.mean(neg_interval), "CTA-negative\n1,754", ha="right", va="center",
            color="#4D69A3", fontsize=8.3, fontweight="bold", linespacing=1.25)
    ax.text(x1 + node_width + 0.025, np.mean(withheld_interval), "Final withheld\n98",
            ha="left", va="center", color="#B4395A", fontsize=8.3,
            fontweight="bold", linespacing=1.25)
    ax.text(x1 + node_width + 0.025, np.mean(retained_interval),
            "Final forced call / retained\n1,790", ha="left", va="center",
            color="#526C98", fontsize=8.3, fontweight="bold", linespacing=1.25)

    ax.annotate(
        "60 of 134 CTA-positive\nspots withheld",
        xy=(0.52, 0.735), xytext=(0.46, 0.825), ha="center", va="center",
        fontsize=7.4, color="#B43D5B", fontweight="bold",
        arrowprops={"arrowstyle": "-", "color": "#D95776", "lw": 0.9},
        bbox={"boxstyle": "round,pad=0.28", "fc": "#FFF4F5", "ec": "#E8A5B5", "lw": 0.6},
        zorder=12,
    )
    ax.text(0.43, 0.655, "74 CTA-positive spots\nremained forced calls",
            ha="center", va="center", fontsize=7.0, color="#A65D4D", fontweight="bold",
            bbox={"boxstyle": "round,pad=0.25", "fc": "#FFF6EF", "ec": "#EBC7B7", "lw": 0.55},
            zorder=12)
    ax.text(0.60, 0.685, "38 CTA-negative\nspots withheld", ha="center", va="center",
            fontsize=6.9, color="#4C69A5", fontweight="bold",
            bbox={"boxstyle": "round,pad=0.25", "fc": "#F4F7FC", "ec": "#B9C7E2", "lw": 0.55},
            zorder=12)
    ax.text(0.49, 0.40, "1,716 CTA-negative spots retained", ha="center", va="center",
            fontsize=7.0, color="#667A9E", fontweight="bold", alpha=0.95, zorder=8)

    ax.text(0.02, 0.965, "Frozen CTA outcome decomposition", transform=ax.transAxes,
            color=TEXT, fontsize=13, fontweight="bold", ha="left", va="top")
    ax.text(0.02, 0.905, "Composite Stage3B decision across 1,888 formal CTA analysis spots",
            transform=ax.transAxes, color=SECONDARY, fontsize=8.5, ha="left", va="top")
    ax.text(0.205, 0.775, "CTA endpoint", transform=ax.transAxes, color=SECONDARY,
            fontsize=7.2, fontweight="bold", ha="center")
    ax.text(0.755, 0.775, "Final decision", transform=ax.transAxes, color=SECONDARY,
            fontsize=7.2, fontweight="bold", ha="center")

    metric_badge(ax, 0.10, "Precision", f"{float(summary['binary_precision']):.4f}", "#C74667")
    metric_badge(ax, 0.30, "Recall", f"{float(summary['binary_recall']):.4f}", "#E07868")
    metric_badge(ax, 0.50, "F1", f"{float(summary['binary_f1']):.4f}", "#8D6AC1")
    metric_badge(ax, 0.70, "Specificity", f"{float(summary['specificity']):.4f}", "#5876B5")

    return export(fig, "C7_Fig2_outcome_decomposition_panel")


def main() -> None:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(DATA_DIR / "c7_metrics_summary.csv").iloc[0]
    pr = pd.read_csv(DATA_DIR / "c7_pr_curve_points.csv")
    fixed = pd.read_csv(DATA_DIR / "c7_fixed_precision_summary.csv")

    assert int(summary["n_total"]) == 1888
    assert int(summary["n_positive"]) == 134
    assert int(summary["n_negative"]) == 1754
    assert fixed["status"].eq("not achieved").all()
    assert set(np.round(fixed["target_precision"].to_numpy(float), 2)) == {0.5, 0.7, 0.8}

    figure1 = figure_pr(summary, pr, fixed)
    figure2 = figure_decomposition(summary)
    report_path = OUT_DIR / "c7_figure_draft_report.txt"
    report_path.write_text(
        "C7 figure drafts\n"
        f"Figure 1: {figure1[0]} | {figure1[1]} | {figure1[2]}\n"
        f"Figure 2: {figure2[0]} | {figure2[1]} | {figure2[2]}\n"
        f"Figure dimensions: {FIGSIZE[0]:.1f} x {FIGSIZE[1]:.2f} inches\n"
        "Implementation notes: Figures use only the frozen C7 summary, PR-curve, and fixed-precision CSV outputs. "
        "The frozen binary decision is displayed separately from finite withheld-score threshold points.\n",
        encoding="utf-8",
    )
    print(f"Figure 1: {figure1}")
    print(f"Figure 2: {figure2}")
    print(f"Report: {report_path}")


if __name__ == "__main__":
    main()
