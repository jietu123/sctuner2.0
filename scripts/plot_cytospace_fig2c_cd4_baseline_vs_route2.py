from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap


METHODS = [
    ("CytoSPACE", "baseline", "#d81b60"),
    ("SVTuner + CytoSPACE", "route2", "#2a9d8f"),
]


def _format_p(p: float) -> str:
    return "0.0001" if p < 0.00015 else f"{p:.4f}"


def _load_summary(path: Path) -> dict[str, float]:
    df = pd.read_csv(path)
    row = df[df["cell_type"].astype(str).eq("CD4 T cells")].iloc[0]
    return {"nes": float(row["NES"]), "p": float(row["pval"])}


def _gradient_line(ax, x: np.ndarray, y: np.ndarray, cmap, lw: float = 1.2) -> None:
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=cmap, norm=plt.Normalize(x.min(), x.max()))
    lc.set_array(x)
    lc.set_linewidth(lw)
    lc.set_zorder(4)
    ax.add_collection(lc)


def _draw_curve(ax, curve: pd.DataFrame, y_offset: float, label: str, nes: float, pval: float) -> None:
    curve = curve.sort_values("rank")
    x = curve["rank"].to_numpy(dtype=float)
    y = curve["running_es"].to_numpy(dtype=float)
    x = x / x.max()

    cmap = LinearSegmentedColormap.from_list(
        f"rank_{label}",
        ["#f03b20", "#ffd400", "#63b96b", "#2b8cbe"],
    )
    _gradient_line(ax, x, y + y_offset, cmap)

    hits = curve[curve["hit"].astype(bool)].copy()
    hit_x = hits["rank"].to_numpy(dtype=float) / curve["rank"].max()
    hit_y = hits["running_es"].to_numpy(dtype=float) + y_offset
    hit_colors = cmap(hit_x)
    ax.vlines(hit_x, y_offset, hit_y, color=hit_colors, linewidth=0.55, alpha=0.9, zorder=3)
    ax.scatter(hit_x, hit_y, s=6.5, color=hit_colors, edgecolor="none", zorder=5)
    ax.axhline(y_offset, color="#222222", linewidth=0.85)

    ax.text(0.05, y_offset + 0.86, label, ha="left", va="center", fontsize=7.5, weight="bold")
    ax.text(
        0.72,
        y_offset + 0.64,
        f"NES = {nes:.2f}\nP = {_format_p(pval)}",
        ha="left",
        va="center",
        fontsize=6.0,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot CD4 Fig.2c-style baseline vs route2 comparison.")
    parser.add_argument("--project_root", default=".")
    parser.add_argument(
        "--baseline_dir",
        default=(
            "result/cytospace_fig2c_melanoma_mel1_rep2_force_bcell_unsupported/"
            "fig2c_enrichment_baseline_bcell_unsupported_seurat410_fgsea114"
        ),
    )
    parser.add_argument(
        "--route2_dir",
        default=(
            "result/cytospace_fig2c_melanoma_mel1_rep2_force_bcell_unsupported/"
            "fig2c_enrichment_route2_bcell_unsupported_seurat410_fgsea114"
        ),
    )
    parser.add_argument("--out_dir", default="visualizations/cytospace_fig2c_melanoma")
    parser.add_argument("--out_prefix", default="fig2c_cd4_baseline_vs_route2")
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    baseline_dir = (root / args.baseline_dir).resolve()
    route2_dir = (root / args.route2_dir).resolve()
    out_dir = (root / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    baseline_curve = pd.read_csv(baseline_dir / "cd4_t_cells_enrichment_curve.csv")
    route2_curve = pd.read_csv(route2_dir / "cd4_t_cells_enrichment_curve.csv")
    baseline_summary = _load_summary(baseline_dir / "fig2c_official_enrichment_summary.csv")
    route2_summary = _load_summary(route2_dir / "fig2c_official_enrichment_summary.csv")

    plt.rcParams.update({
        "font.family": "Arial",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig = plt.figure(figsize=(2.55, 2.65), dpi=300)
    gs = fig.add_gridspec(
        nrows=2,
        ncols=1,
        height_ratios=[1.0, 0.22],
        hspace=0.16,
        left=0.16,
        right=0.98,
        top=0.93,
        bottom=0.14,
    )
    ax = fig.add_subplot(gs[0, 0])
    ax_arrow = fig.add_subplot(gs[1, 0])

    _draw_curve(ax, baseline_curve, 1.10, "CytoSPACE", baseline_summary["nes"], baseline_summary["p"])
    _draw_curve(ax, route2_curve, 0.00, "SVTuner + CytoSPACE", route2_summary["nes"], route2_summary["p"])

    ax.set_xlim(0, 1.0)
    ax.set_ylim(-0.08, 2.18)
    ax.set_xticks([])
    ax.set_yticks([0, 0.5, 1.0, 1.1, 1.6, 2.1])
    ax.set_yticklabels(["0", "0.5", "1", "0", "0.5", "1"])
    ax.tick_params(axis="y", labelsize=5.7, width=0.8, length=2.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_linewidth(0.85)
    ax.text(0.50, 2.15, "CD4 T cells", ha="center", va="top", fontsize=9.0)

    delta = route2_summary["nes"] - baseline_summary["nes"]
    ax.text(
        0.02,
        -0.02,
        f"Delta NES (route2 - baseline) = {delta:+.2f}",
        ha="left",
        va="top",
        fontsize=5.6,
        color="#2a9d8f" if delta > 0 else "#555555",
    )

    ax_arrow.axis("off")
    ax_arrow.annotate(
        "",
        xy=(0.86, 0.54),
        xytext=(0.20, 0.54),
        arrowprops=dict(arrowstyle="->", lw=1.8, color="#8a8a8a"),
    )
    ax_arrow.text(0.01, 0.54, "Increasing distance\nfrom tumor", ha="left", va="center", fontsize=6.0)

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
