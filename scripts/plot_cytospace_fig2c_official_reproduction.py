from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap


PUBLISHED = {
    "CD4 T cells": {"nes": 2.38191397683252, "p": 0.000129433083096039},
    "CD8 T cells": {"nes": 1.81251071297423, "p": 0.000136500136500137},
}


def _format_p(p: float) -> str:
    return "0.0001" if p < 0.00015 else f"{p:.4f}"


def _load_label_values(summary_csv: Path | None) -> dict[str, dict[str, float]]:
    if summary_csv is None:
        return PUBLISHED
    df = pd.read_csv(summary_csv)
    values: dict[str, dict[str, float]] = {}
    for _, row in df.iterrows():
        values[str(row["cell_type"])] = {"nes": float(row["NES"]), "p": float(row["pval"])}
    return values


def _gradient_line(ax, x: np.ndarray, y: np.ndarray, cmap) -> None:
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=cmap, norm=plt.Normalize(x.min(), x.max()))
    lc.set_array(x)
    lc.set_linewidth(1.25)
    lc.set_zorder(4)
    ax.add_collection(lc)


def _plot_panel(ax, curve: pd.DataFrame, title: str, ymax: float, mid_tick: float, label_values: dict[str, dict[str, float]]) -> None:
    curve = curve.sort_values("rank")
    x = curve["rank"].to_numpy(dtype=float)
    y = curve["running_es"].to_numpy(dtype=float)
    x = x / x.max()

    cmap = LinearSegmentedColormap.from_list(
        "cytospace_fig2c_rank",
        ["#f03b20", "#ffd400", "#63b96b", "#2b8cbe"],
    )

    _gradient_line(ax, x, y, cmap)

    hits = curve[curve["hit"].astype(bool)].copy()
    hit_x = hits["rank"].to_numpy(dtype=float) / curve["rank"].max()
    hit_y = hits["running_es"].to_numpy(dtype=float)
    hit_colors = cmap(hit_x)
    ax.vlines(hit_x, 0, hit_y, color=hit_colors, linewidth=0.6, alpha=0.92, zorder=3)
    ax.scatter(hit_x, hit_y, s=7, color=hit_colors, edgecolor="none", zorder=5)

    ax.axhline(0, color="#222222", linewidth=0.9)
    ax.set_xlim(0, 1.0)
    ax.set_ylim(-0.08, ymax)
    ax.set_xticks([])
    ax.set_yticks([0, mid_tick, ymax])
    ax.set_yticklabels([f"{0:g}", f"{mid_tick:g}", f"{ymax:g}"])
    ax.tick_params(axis="y", labelsize=6, width=0.85, length=3)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_position(("data", 0))
    ax.spines["left"].set_linewidth(0.85)
    ax.spines["bottom"].set_linewidth(0.85)

    ax.text(0.50, ymax * 0.88, title, ha="center", va="center", fontsize=8.5)
    pub = label_values[title]
    ax.text(
        0.76,
        ymax * 0.58,
        f"NES = {pub['nes']:.2f}\n\nP = {_format_p(pub['p'])}",
        ha="left",
        va="center",
        fontsize=5.8,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Reproduce CytoSPACE paper Fig. 2c from official source-data values.")
    parser.add_argument("--project_root", default=".")
    parser.add_argument(
        "--curve_dir",
        default="result/cytospace_fig2c_official_slide1_readme_default/fig2c_enrichment_seurat410_fgsea114",
    )
    parser.add_argument(
        "--out_dir",
        default="visualizations/cytospace_fig2c_melanoma",
    )
    parser.add_argument("--summary_csv", default=None, help="Optional computed enrichment summary with NES/pval columns.")
    parser.add_argument("--out_prefix", default="fig2c_official_reproduction")
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    curve_dir = (root / args.curve_dir).resolve()
    out_dir = (root / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_csv = (root / args.summary_csv).resolve() if args.summary_csv else None
    label_values = _load_label_values(summary_csv)

    cd4 = pd.read_csv(curve_dir / "cd4_t_cells_enrichment_curve.csv")
    cd8 = pd.read_csv(curve_dir / "cd8_t_cells_enrichment_curve.csv")

    plt.rcParams.update({
        "font.family": "Arial",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig = plt.figure(figsize=(2.05, 2.25), dpi=300)
    gs = fig.add_gridspec(
        nrows=3,
        ncols=1,
        height_ratios=[1.0, 1.0, 0.22],
        hspace=0.52,
        left=0.16,
        right=0.97,
        top=0.96,
        bottom=0.16,
    )

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax_arrow = fig.add_subplot(gs[2, 0])

    _plot_panel(ax1, cd4, "CD4 T cells", ymax=1.0, mid_tick=0.5, label_values=label_values)
    _plot_panel(ax2, cd8, "CD8 T cells", ymax=0.8, mid_tick=0.4, label_values=label_values)

    ax_arrow.axis("off")
    ax_arrow.annotate(
        "",
        xy=(0.84, 0.52),
        xytext=(0.16, 0.52),
        arrowprops=dict(arrowstyle="->", lw=1.8, color="#8a8a8a"),
    )
    ax_arrow.text(0.02, 0.52, "Increasing distance\nfrom tumor", ha="left", va="center", fontsize=6.2)

    png = out_dir / f"{args.out_prefix}.png"
    pdf = out_dir / f"{args.out_prefix}.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    print(f"[OK] wrote: {png}")
    print(f"[OK] wrote: {pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
