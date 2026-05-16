from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
import json


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCORE_DIR = PROJECT_ROOT / "result" / "real_profile_mask_fig2d_score"
VIS_DIR = PROJECT_ROOT / "visualizations" / "simulations" / "real_profile_mask_fig2d_only"

BASELINE = "CytoSPACE"
ROUTE2 = "SVTuner + CytoSPACE"
METHOD_ORDER = [BASELINE, ROUTE2, "Tangram"]
COLORS = {
    BASELINE: "#d81b60",
    ROUTE2: "#2a9d8f",
    "Tangram": "#6a994e",
}
EDGE = "#2d2d2d"
GRID = "#d9d9d9"
BG = "white"


def _sig_label(pvalue: float) -> str:
    if pvalue < 0.001:
        return "***"
    if pvalue < 0.01:
        return "**"
    if pvalue < 0.05:
        return "*"
    return "n.s."


def main() -> int:
    VIS_DIR.mkdir(parents=True, exist_ok=True)

    long_df = pd.read_csv(SCORE_DIR / "fig2d_score_long.csv")
    config = json.loads((SCORE_DIR / "fig2d_score_config.json").read_text(encoding="utf-8"))
    pivot = long_df.pivot_table(index="scenario_label", columns="method", values="value", aggfunc="first")
    pivot = pivot[METHOD_ORDER].dropna()
    low_mask = (pivot[METHOD_ORDER].min(axis=1) < 0.80)
    plot_pivot = pivot.loc[~low_mask].copy()
    low_pivot = pivot.loc[low_mask].copy()

    fig = plt.figure(figsize=(4.35, 4.25), dpi=300)
    fig.patch.set_facecolor(BG)
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[4.0, 0.65], hspace=0.06)
    ax = fig.add_subplot(gs[0])
    ax_low = fig.add_subplot(gs[1], sharex=ax)
    ax.set_facecolor(BG)
    ax_low.set_facecolor(BG)
    ax.grid(axis="y", color=GRID, linestyle="--", linewidth=0.7, alpha=0.9)
    ax.grid(axis="x", visible=False)
    ax_low.grid(axis="y", color=GRID, linestyle="--", linewidth=0.7, alpha=0.9)
    ax_low.grid(axis="x", visible=False)

    positions = np.arange(1, len(METHOD_ORDER) + 1)
    values = [plot_pivot[m].to_numpy(dtype=float) for m in METHOD_ORDER]
    bp = ax.boxplot(
        values,
        positions=positions,
        widths=0.46,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": EDGE, "linewidth": 1.2},
        boxprops={"edgecolor": EDGE, "linewidth": 1.0},
        whiskerprops={"color": EDGE, "linewidth": 0.9},
        capprops={"color": EDGE, "linewidth": 0.9},
    )
    for patch, method in zip(bp["boxes"], METHOD_ORDER):
        patch.set_facecolor(COLORS[method])
        patch.set_alpha(0.85)

    rng = np.random.default_rng(42)
    for i, method in enumerate(METHOD_ORDER, start=1):
        y = plot_pivot[method].to_numpy(dtype=float)
        x = np.full(len(y), i, dtype=float) + rng.normal(0, 0.035, len(y))
        ax.scatter(x, y, s=10, color="black", alpha=0.8, zorder=3, linewidths=0)

    if not low_pivot.empty:
        for _, row in low_pivot.iterrows():
            ax_low.plot(positions, row[METHOD_ORDER].to_numpy(dtype=float), color="#9a9a9a", linewidth=0.65, alpha=0.65, zorder=1)
        for i, method in enumerate(METHOD_ORDER, start=1):
            y = low_pivot[method].to_numpy(dtype=float)
            x = np.full(len(y), i, dtype=float) + rng.normal(0, 0.015, len(y))
            ax_low.scatter(x, y, s=10, color="black", alpha=0.8, zorder=3, linewidths=0)

    all_vals = np.concatenate(values)
    lo = float(np.min(all_vals))
    hi = float(np.max(all_vals))
    span = max(hi - lo, 1e-6)
    ax.set_ylim(lo - span * 0.08, hi + span * 0.20)
    if not low_pivot.empty:
        low_vals = low_pivot[METHOD_ORDER].to_numpy(dtype=float).ravel()
        low_lo = float(low_vals.min())
        low_hi = float(low_vals.max())
        low_span = max(low_hi - low_lo, 1e-6)
        ax_low.set_ylim(low_lo - low_span * 0.45, low_hi + low_span * 0.70)

    yb = hi + span * 0.08
    bracket_specs = [(2, ROUTE2, yb), (3, "Tangram", yb + span * 0.075)]
    for x2, method, y_bracket in bracket_specs:
        paired = pivot[[BASELINE, method]].dropna()
        pvalue = float(
            wilcoxon(
                paired[method].to_numpy(dtype=float),
                paired[BASELINE].to_numpy(dtype=float),
                alternative="greater",
                zero_method="wilcox",
            ).pvalue
        )
        sig = _sig_label(pvalue)
        ax.plot(
            [1, 1, x2, x2],
            [y_bracket - span * 0.02, y_bracket, y_bracket, y_bracket - span * 0.02],
            color=EDGE,
            linewidth=0.85,
        )
        ax.text((1 + x2) / 2, y_bracket + span * 0.010, sig, ha="center", va="bottom", fontsize=8.0, weight="bold", color=EDGE)

    ax.set_xticks(positions)
    ax.tick_params(axis="x", labelbottom=False, bottom=False)
    ax_low.set_xticks(positions)
    ax_low.set_xticklabels(["CytoSPACE", "SVTuner +\nCytoSPACE", "Tangram"], fontsize=8.0, rotation=42, ha="right", rotation_mode="anchor")
    ax.set_ylabel("Suppression enrichment index", fontsize=8.8)
    ax.tick_params(axis="y", labelsize=8.0)
    ax_low.tick_params(axis="y", labelsize=7.0)
    ax.set_xlim(0.55, len(METHOD_ORDER) + 0.45)

    fig.text(0.045, 0.972, "d", fontsize=12.0, weight="bold", va="top")
    fig.text(
        0.185,
        0.972,
        "Benchmarking masked-target suppression\nacross 10 real profile-mask scenarios",
        fontsize=7.0,
        weight="bold",
        va="top",
    )

    for spine in ["top", "right", "bottom"]:
        ax.spines[spine].set_visible(False)
    for spine in ["top", "right"]:
        ax_low.spines[spine].set_visible(False)
    ax.spines["left"].set_color(EDGE)
    ax_low.spines["left"].set_color(EDGE)
    ax_low.spines["bottom"].set_color(EDGE)

    if not low_pivot.empty:
        kwargs = dict(marker=[(-1, -0.55), (1, 0.55)], markersize=5, linestyle="none", color=EDGE, mec=EDGE, mew=0.8, clip_on=False)
        ax.plot([0], [0], transform=ax.transAxes, **kwargs)
        ax_low.plot([0], [1], transform=ax_low.transAxes, **kwargs)

    long_df.to_csv(VIS_DIR / "fig2_panel_d_real_profile_mask_source_values.csv", index=False)
    fig.subplots_adjust(left=0.17, right=0.98, top=0.82, bottom=0.20)

    out_png = VIS_DIR / "fig2_panel_d_real_profile_mask.png"
    out_pdf = VIS_DIR / "fig2_panel_d_real_profile_mask.pdf"
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"[OK] wrote: {out_png}")
    print(f"[OK] wrote: {out_pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
