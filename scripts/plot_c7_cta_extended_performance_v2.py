"""Create the frozen C7-V2 CTA performance visualizations.

This script is intentionally visualization-only. It reads the frozen C7 tables and
spot-level endpoint/score contracts, checks their frozen counts, and exports two
publication figures without modifying any input or recomputing model outputs.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba
from scipy.stats import gaussian_kde


ROOT = Path(__file__).resolve().parents[1]
C7_DIR = ROOT / "visualizations" / "bioapp_experiment" / "C7_cta_extended_performance"
FIG_DIR = C7_DIR / "figures"
ENDPOINT_CSV = (
    ROOT
    / "visualizations"
    / "bioapp_experiment"
    / "bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping"
    / "spot_level_endpoint_freeze.csv"
)
SCORE_CSV = (
    ROOT
    / "visualizations"
    / "bioapp_experiment"
    / "bioapp_phase7_svtuner_aware_execution_against_frozen_cta_immune_endpoint"
    / "svtuner_immune_all_dropout"
    / "svtuner_immune_all_dropout_spot_level_raw_output_contract.csv"
)

COLORS = {
    "teal": "#087F8C",
    "warm_gray": "#A8A29E",
    "amber": "#D99520",
    "violet": "#8064B0",
    "magenta": "#C23B64",
    "coral": "#D65F4A",
    "blue": "#477DB3",
    "text": "#26323A",
    "secondary": "#667078",
    "grid": "#E9EDF0",
    "neutral_rug": "#C8CDD1",
}

WIDTH_IN = 7.08  # 179.8 mm
FIG1_HEIGHT_IN = 3.72
FIG2_HEIGHT_IN = 3.80


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8.0,
            "axes.labelsize": 8.6,
            "xtick.labelsize": 7.6,
            "ytick.labelsize": 7.6,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def save_figure(fig: mpl.figure.Figure, stem: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "svg"):
        fig.savefig(FIG_DIR / f"{stem}.{ext}", facecolor="white")
    fig.savefig(FIG_DIR / f"{stem}.png", dpi=600, facecolor="white")


def clean_axes(ax: mpl.axes.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#7E878E")
    ax.spines["bottom"].set_color("#7E878E")
    ax.tick_params(colors=COLORS["secondary"])
    ax.xaxis.label.set_color(COLORS["text"])
    ax.yaxis.label.set_color(COLORS["text"])


def plot_pr_landscape() -> None:
    pr = pd.read_csv(C7_DIR / "c7_pr_curve_points.csv")
    metrics = pd.read_csv(C7_DIR / "c7_metrics_summary.csv").iloc[0]
    required = {"threshold", "precision", "recall"}
    if not required.issubset(pr.columns):
        raise ValueError(f"Missing PR columns: {sorted(required - set(pr.columns))}")

    fig, ax = plt.subplots(figsize=(WIDTH_IN, FIG1_HEIGHT_IN))
    fig.subplots_adjust(left=0.105, right=0.965, bottom=0.205, top=0.955)

    ax.plot(
        pr["recall"],
        pr["precision"],
        color=COLORS["teal"],
        lw=1.8,
        solid_capstyle="round",
        zorder=3,
    )
    prevalence = float(metrics["positive_prevalence"])
    ax.axhline(
        prevalence,
        color=COLORS["warm_gray"],
        lw=0.9,
        ls=(0, (4, 3)),
        zorder=1,
    )

    landmarks = {
        "p25": (0.6044776119, 0.2500000000),
        "p30": (0.5298507463, 0.3008474576),
        "maxp": (0.2089552239, 0.4912280702),
        "binary": (0.4477611940, 0.6122448980),
    }
    ax.scatter(*landmarks["p25"], s=28, facecolor="white", edgecolor=COLORS["amber"], lw=1.2, zorder=6)
    ax.scatter(*landmarks["p30"], s=47, facecolor=COLORS["amber"], edgecolor="white", lw=0.8, zorder=7)
    ax.scatter(*landmarks["maxp"], s=42, facecolor=COLORS["violet"], edgecolor="white", lw=0.8, zorder=7)
    ax.scatter(
        *landmarks["binary"],
        s=78,
        marker="D",
        facecolor=COLORS["magenta"],
        edgecolor="white",
        lw=0.9,
        zorder=8,
    )

    ax.annotate(
        "P >= 0.25\nP 0.250 · R 0.604",
        xy=landmarks["p25"],
        xytext=(0.676, 0.226),
        color=COLORS["amber"],
        fontsize=7.2,
        linespacing=1.22,
        ha="left",
        va="center",
        arrowprops=dict(arrowstyle="-", color=to_rgba(COLORS["amber"], 0.62), lw=0.65),
    )
    ax.annotate(
        "P >= 0.30\nP 0.301 · R 0.530",
        xy=landmarks["p30"],
        xytext=(0.625, 0.338),
        color=COLORS["amber"],
        fontsize=7.6,
        fontweight="bold",
        linespacing=1.22,
        ha="left",
        va="center",
        arrowprops=dict(arrowstyle="-", color=to_rgba(COLORS["amber"], 0.72), lw=0.75),
    )
    ax.annotate(
        "Max score-threshold precision\nP = 0.491",
        xy=landmarks["maxp"],
        xytext=(0.070, 0.542),
        color=COLORS["violet"],
        fontsize=7.5,
        linespacing=1.22,
        ha="left",
        va="center",
        arrowprops=dict(arrowstyle="-", color=to_rgba(COLORS["violet"], 0.65), lw=0.7),
    )
    ax.annotate(
        "Composite Stage3B\nP 0.612 · R 0.448\nF1 0.517",
        xy=landmarks["binary"],
        xytext=(0.520, 0.600),
        color=COLORS["magenta"],
        fontsize=8.1,
        fontweight="bold",
        linespacing=1.18,
        ha="left",
        va="center",
        arrowprops=dict(arrowstyle="-", color=to_rgba(COLORS["magenta"], 0.72), lw=0.75),
    )
    ax.text(
        0.520,
        0.537,
        "not score-thresholded",
        color=COLORS["secondary"],
        fontsize=6.9,
        fontstyle="italic",
        ha="left",
        va="top",
    )
    ax.text(
        0.975,
        0.93,
        "AP 0.284\nbaseline 0.071",
        transform=ax.transAxes,
        color=COLORS["text"],
        fontsize=7.6,
        linespacing=1.35,
        ha="right",
        va="top",
    )

    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 0.68)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_xticks(np.linspace(0, 1, 6))
    ax.set_yticks(np.arange(0, 0.7, 0.1))
    ax.grid(axis="y", color=COLORS["grid"], lw=0.55, zorder=0)
    clean_axes(ax)
    fig.text(
        0.53,
        0.050,
        "No finite score threshold reached precision >= 0.50.",
        color=COLORS["secondary"],
        fontsize=7.2,
        ha="center",
        va="center",
    )

    save_figure(fig, "C7_Fig1_PR_landscape_v2")
    plt.close(fig)


def parse_binary(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series
    values = series.astype(str).str.strip().str.lower()
    parsed = values.map({"true": True, "false": False, "1": True, "0": False})
    if parsed.isna().any():
        raise ValueError("Unrecognized withheld_binary value")
    return parsed.astype(bool)


def load_formal_spots() -> pd.DataFrame:
    endpoint = pd.read_csv(ENDPOINT_CSV, usecols=["barcode", "primary_endpoint_status"])
    endpoint = endpoint[endpoint["primary_endpoint_status"].isin(["positive", "negative"])]
    scores = pd.read_csv(SCORE_CSV, usecols=["barcode", "withheld_score", "withheld_binary"])
    merged = endpoint.merge(scores, on="barcode", how="inner", validate="one_to_one")
    merged["withheld_binary"] = parse_binary(merged["withheld_binary"])

    counts = merged["primary_endpoint_status"].value_counts().to_dict()
    if len(merged) != 1888 or counts != {"negative": 1754, "positive": 134}:
        raise ValueError(f"Frozen formal-set mismatch: n={len(merged)}, counts={counts}")
    tab = pd.crosstab(merged["primary_endpoint_status"], merged["withheld_binary"])
    frozen = (
        int(tab.loc["positive", True]),
        int(tab.loc["positive", False]),
        int(tab.loc["negative", True]),
        int(tab.loc["negative", False]),
    )
    if frozen != (60, 74, 38, 1716):
        raise ValueError(f"Frozen binary-count mismatch: {frozen}")
    return merged


def scaled_kde(values: np.ndarray, grid: np.ndarray, height: float) -> np.ndarray:
    density = gaussian_kde(values)(grid)
    return height * density / density.max()


def plot_score_decision_landscape() -> None:
    data = load_formal_spots()
    positive = data[data["primary_endpoint_status"] == "positive"].copy()
    negative = data[data["primary_endpoint_status"] == "negative"].copy()

    score_min = float(data["withheld_score"].min())
    score_max = float(data["withheld_score"].max())
    span = score_max - score_min
    xlo = score_min - 0.045 * span
    xhi = score_max + 0.045 * span
    grid = np.linspace(xlo, xhi, 600)

    fig, ax = plt.subplots(figsize=(WIDTH_IN, FIG2_HEIGHT_IN))
    fig.subplots_adjust(left=0.115, right=0.965, bottom=0.175, top=0.965)
    ax.set_facecolor("white")

    rows = [
        (positive, 1.10, COLORS["coral"], "CTA positive · n = 134", "60 withheld · 74 retained"),
        (negative, 0.25, COLORS["blue"], "CTA negative · n = 1754", "38 withheld · 1716 retained"),
    ]
    for frame, base, color, row_label, count_label in rows:
        values = frame["withheld_score"].to_numpy(float)
        dens = scaled_kde(values, grid, height=0.47)
        ax.fill_between(grid, base, base + dens, color=color, alpha=0.17, linewidth=0, zorder=1)
        ax.plot(grid, base + dens, color=color, lw=1.45, zorder=3)
        ax.hlines(base, xlo, xhi, color=to_rgba(color, 0.38), lw=0.6, zorder=1)

        retained = frame.loc[~frame["withheld_binary"], "withheld_score"].to_numpy(float)
        withheld = frame.loc[frame["withheld_binary"], "withheld_score"].to_numpy(float)
        rug_top = base - 0.035
        ax.vlines(
            retained,
            rug_top - 0.055,
            rug_top,
            color=to_rgba(color, 0.24 if len(frame) < 500 else 0.16),
            lw=0.42,
            zorder=4,
        )
        ax.vlines(
            withheld,
            rug_top - 0.074,
            rug_top + 0.006,
            color=to_rgba(COLORS["magenta"], 0.93),
            lw=0.70,
            zorder=5,
        )
        ax.text(
            xlo,
            base + 0.515,
            row_label,
            color=color,
            fontsize=8.3,
            fontweight="bold",
            ha="left",
            va="bottom",
        )
        ax.text(
            xlo,
            base + 0.472,
            count_label,
            color=COLORS["secondary"],
            fontsize=7.2,
            ha="left",
            va="top",
        )

    threshold = 0.6233151196
    ax.axvline(threshold, color="#AAB0B5", lw=0.75, ls=(0, (3, 3)), zorder=0)
    ax.text(
        threshold + 0.006 * span,
        1.640,
        "t = 0.6233",
        color=COLORS["secondary"],
        fontsize=6.8,
        ha="left",
        va="top",
    )

    ax.text(
        xlo + 0.66 * (xhi - xlo),
        0.825,
        "Composite decision cannot be reduced\nto a single score cutoff",
        color=COLORS["text"],
        fontsize=8.0,
        fontweight="bold",
        linespacing=1.25,
        ha="left",
        va="center",
    )
    ax.text(
        xlo + 0.66 * (xhi - xlo),
        0.710,
        "Best 1D match: 93 discordant spots",
        color=COLORS["secondary"],
        fontsize=7.1,
        ha="left",
        va="center",
    )

    # Small direct key for the rug encoding; no boxed legend.
    key_x = xlo + 0.015 * (xhi - xlo)
    ax.vlines(key_x, -0.005, 0.075, color=COLORS["magenta"], lw=1.1)
    ax.text(key_x + 0.012 * span, 0.035, "withheld", color=COLORS["magenta"], fontsize=6.8, va="center")
    key_x2 = key_x + 0.14 * span
    ax.vlines(key_x2, 0.005, 0.065, color=COLORS["neutral_rug"], lw=0.75)
    ax.text(key_x2 + 0.012 * span, 0.035, "retained", color=COLORS["secondary"], fontsize=6.8, va="center")

    ax.set_xlim(xlo, xhi)
    ax.set_ylim(-0.06, 1.74)
    ax.set_xlabel("Continuous withheld score")
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_color("#7E878E")
    ax.tick_params(axis="x", colors=COLORS["secondary"], width=0.6, length=3)
    ax.xaxis.label.set_color(COLORS["text"])

    save_figure(fig, "C7_Fig2_score_decision_landscape_v2")
    plt.close(fig)


def main() -> None:
    set_style()
    plot_pr_landscape()
    plot_score_decision_landscape()
    print(f"Wrote C7-V2 figures to {FIG_DIR}")


if __name__ == "__main__":
    main()
