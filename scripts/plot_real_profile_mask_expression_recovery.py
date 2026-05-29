from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


METHODS = ["CytoSPACE", "SVTuner + CytoSPACE"]
COLORS = {
    "CytoSPACE": "#ee6a5a",
    "SVTuner + CytoSPACE": "#2f9d94",
}


def _sem(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if values.size <= 1:
        return 0.0
    return float(stats.sem(values, nan_policy="omit"))


def _format_p(pvalue: float | None) -> str:
    if pvalue is None or not np.isfinite(pvalue):
        return "P = n/a"
    if pvalue < 1e-3:
        return f"P = {pvalue:.1e}"
    return f"P = {pvalue:.3f}"


def _add_bracket(ax: plt.Axes, x1: float, x2: float, y: float, text: str) -> None:
    h = 0.009
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="black", lw=1.5, clip_on=False)
    ax.text((x1 + x2) / 2, y + h + 0.005, text, ha="center", va="bottom", fontsize=10.5)


def plot_expression_recovery(project_root: Path) -> Path:
    result_dir = project_root / "result" / "real_profile_mask_expression_recovery"
    vis_dir = project_root / "visualizations" / "simulations" / "real_profile_mask_expression_recovery"
    vis_dir.mkdir(parents=True, exist_ok=True)

    by_scenario = pd.read_csv(result_dir / "expression_recovery_by_scenario.csv")
    config_path = result_dir / "expression_recovery_config.json"
    config = json.loads(config_path.read_text(encoding="utf-8")) if config_path.exists() else {}

    x = np.arange(len(METHODS), dtype=float)
    means = np.array([by_scenario[m].mean() for m in METHODS], dtype=float)
    sems = np.array([_sem(by_scenario[m].to_numpy(dtype=float)) for m in METHODS], dtype=float)

    fig, ax = plt.subplots(figsize=(4.15, 4.55), dpi=300)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    bar_width = 0.46
    ax.bar(
        x,
        means,
        width=bar_width,
        color=[COLORS[m] for m in METHODS],
        edgecolor="#555555",
        linewidth=1.2,
        alpha=0.88,
        zorder=2,
    )
    ax.errorbar(
        x,
        means,
        yerr=sems,
        fmt="none",
        ecolor="#333333",
        elinewidth=1.35,
        capsize=5,
        capthick=1.35,
        zorder=4,
    )

    rng = np.random.default_rng(13)
    jitter = rng.normal(0.0, 0.042, size=len(by_scenario))
    for i, row in by_scenario.iterrows():
        j = float(jitter[i])
        for xi, method in zip(x, METHODS, strict=True):
            ax.scatter(xi + j, float(row[method]), s=20, color="#242424", alpha=0.82, zorder=5, linewidth=0)

    pvalue = config.get("wilcoxon_pvalue_route2_gt_baseline")
    y_top = float(np.nanmax(by_scenario[METHODS].to_numpy(dtype=float)))
    _add_bracket(ax, x[0], x[1], y_top + 0.018, _format_p(pvalue))

    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(0.0, min(0.75, y_top + 0.075))
    ax.set_xticks(x)
    ax.set_xticklabels(["CytoSPACE", "SVTuner +\nCytoSPACE"], fontsize=11)
    ax.set_ylabel("Cosine similarity", fontsize=12)
    ax.set_title("Masked-ST expression recovery", fontsize=12.5, fontweight="bold", pad=12)

    # Keep the panel close to the paper-style benchmark barplot; details stay in the CSV/config.

    ax.grid(axis="y", color="#dfdfdf", linestyle="--", linewidth=0.8, alpha=0.85, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.6)
    ax.spines["bottom"].set_linewidth(1.6)
    ax.tick_params(axis="both", width=1.4, length=5)

    fig.subplots_adjust(left=0.25, right=0.97, top=0.83, bottom=0.25)

    out_png = vis_dir / "expression_recovery_cosine_summary_bar.png"
    out_pdf = vis_dir / "expression_recovery_cosine_summary_bar.pdf"
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    return out_png


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project_root", default=".", type=Path)
    args = parser.parse_args()
    out = plot_expression_recovery(args.project_root.resolve())
    print(f"[OK] wrote: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
