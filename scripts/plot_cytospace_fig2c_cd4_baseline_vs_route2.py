from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap

EPS = 1.0e-8


METHODS = [
    ("CytoSPACE", "baseline", "#d81b60"),
    ("SVTuner + CytoSPACE", "route2", "#2a9d8f"),
]


def _format_p(p: float) -> str:
    return "0.0001" if p < 0.00015 else f"{p:.4f}"


def _prefix_for_cell_type(cell_type: str) -> str:
    return "".join(ch.lower() if ch.isalnum() else "_" for ch in cell_type).strip("_")


def _load_summary(path: Path, cell_type: str) -> dict[str, float]:
    df = pd.read_csv(path)
    row = df[df["cell_type"].astype(str).eq(cell_type)].iloc[0]
    return {"nes": float(row["NES"]), "p": float(row["pval"])}


def _gradient_line(ax, x: np.ndarray, y: np.ndarray, cmap, lw: float = 1.6) -> None:
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=cmap, norm=plt.Normalize(x.min(), x.max()))
    lc.set_array(x)
    lc.set_linewidth(lw)
    lc.set_zorder(4)
    ax.add_collection(lc)


def _smooth_curve(y: np.ndarray, window: int = 121) -> np.ndarray:
    if len(y) < 7:
        return y.copy()
    win = min(window, len(y) - (1 - len(y) % 2))
    if win < 5:
        return y.copy()
    if win % 2 == 0:
        win -= 1
    pad = win // 2
    padded = np.pad(y, (pad, pad), mode="edge")
    kernel = np.ones(win, dtype=float) / float(win)
    return np.convolve(padded, kernel, mode="valid")


def _nice_limits(y: np.ndarray) -> tuple[float, float]:
    y_min = float(np.min(y))
    y_max = float(np.max(y))
    lo = min(0.0, y_min)
    hi = max(0.0, y_max)
    span = max(hi - lo, 0.18)
    pad = span * 0.14
    return lo - pad, hi + pad


def _draw_curve(
    ax,
    curve: pd.DataFrame,
    label: str,
    nes: float,
    pval: float,
    panel_normalize_display: bool = False,
) -> None:
    curve = curve.sort_values("rank")
    x = curve["rank"].to_numpy(dtype=float)
    y_raw = curve["running_es"].to_numpy(dtype=float)
    if panel_normalize_display:
        denom = max(float(np.max(np.abs(y_raw))), EPS)
        y = y_raw / denom
    else:
        y = y_raw
    y_smooth = _smooth_curve(y)
    x = x / x.max()

    cmap = LinearSegmentedColormap.from_list(
        f"rank_{label}",
        ["#f03b20", "#ffd400", "#63b96b", "#2b8cbe"],
    )
    _gradient_line(ax, x, y_smooth, cmap)

    hits = curve[curve["hit"].astype(bool)].copy()
    hit_x = hits["rank"].to_numpy(dtype=float) / curve["rank"].max()
    hit_y = hits["running_es"].to_numpy(dtype=float)
    if panel_normalize_display:
        hit_y = hit_y / max(float(np.max(np.abs(y_raw))), EPS)
    hit_colors = cmap(hit_x)
    ax.vlines(hit_x, 0.0, hit_y, color=hit_colors, linewidth=0.65, alpha=0.95, zorder=3)
    ax.scatter(hit_x, hit_y, s=7.5, color=hit_colors, edgecolor="none", zorder=5)
    ax.axhline(0.0, color="#222222", linewidth=0.9)

    y_lo, y_hi = _nice_limits(y)
    ax.set_ylim(y_lo, y_hi)

    ax.text(
        0.02,
        1.03,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.4,
        weight="bold",
        clip_on=False,
    )
    ax.text(
        0.72,
        y_hi - (y_hi - y_lo) * 0.23,
        f"NES = {nes:.2f}\nP = {_format_p(pval)}",
        ha="left",
        va="center",
        fontsize=6.2,
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
    parser.add_argument("--cell_type", default="CD4 T cells")
    parser.add_argument("--panel_normalize_display", action="store_true")
    parser.add_argument("--show_delta_label", action="store_true")
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    baseline_dir = (root / args.baseline_dir).resolve()
    route2_dir = (root / args.route2_dir).resolve()
    out_dir = (root / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    cell_type = str(args.cell_type)
    curve_prefix = _prefix_for_cell_type(cell_type)
    baseline_curve = pd.read_csv(baseline_dir / f"{curve_prefix}_enrichment_curve.csv")
    route2_curve = pd.read_csv(route2_dir / f"{curve_prefix}_enrichment_curve.csv")
    baseline_summary = _load_summary(baseline_dir / "fig2c_official_enrichment_summary.csv", cell_type)
    route2_summary = _load_summary(route2_dir / "fig2c_official_enrichment_summary.csv", cell_type)

    plt.rcParams.update({
        "font.family": "Arial",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig = plt.figure(figsize=(2.85, 3.65), dpi=300)
    gs = fig.add_gridspec(
        nrows=3,
        ncols=1,
        height_ratios=[1.0, 1.0, 0.24],
        hspace=0.28,
        left=0.18,
        right=0.98,
        top=0.93,
        bottom=0.12,
    )
    ax_top = fig.add_subplot(gs[0, 0])
    ax_bottom = fig.add_subplot(gs[1, 0], sharex=ax_top)
    ax_arrow = fig.add_subplot(gs[2, 0], sharex=ax_top)

    _draw_curve(
        ax_top,
        baseline_curve,
        "CytoSPACE",
        baseline_summary["nes"],
        baseline_summary["p"],
        panel_normalize_display=bool(args.panel_normalize_display),
    )
    _draw_curve(
        ax_bottom,
        route2_curve,
        "SVTuner + CytoSPACE",
        route2_summary["nes"],
        route2_summary["p"],
        panel_normalize_display=bool(args.panel_normalize_display),
    )

    for panel_ax in (ax_top, ax_bottom):
        panel_ax.set_xlim(0, 1.0)
        panel_ax.tick_params(axis="y", labelsize=6.0, width=0.9, length=2.8)
        panel_ax.spines["top"].set_visible(False)
        panel_ax.spines["right"].set_visible(False)
        panel_ax.spines["bottom"].set_visible(False)
        panel_ax.spines["left"].set_linewidth(0.9)
        panel_ax.set_xticks([])

    y_label = "Panel-scaled enrichment score" if args.panel_normalize_display else "Running enrichment score"
    fig.text(0.07, 0.56, y_label, rotation=90, va="center", ha="center", fontsize=8.0)
    fig.text(0.55, 0.965, cell_type, ha="center", va="top", fontsize=8.8)

    delta = route2_summary["nes"] - baseline_summary["nes"]
    if args.show_delta_label:
        ax_bottom.text(
            0.02,
            ax_bottom.get_ylim()[0] + (ax_bottom.get_ylim()[1] - ax_bottom.get_ylim()[0]) * 0.08,
            f"Delta NES (route2 - baseline) = {delta:+.2f}",
            ha="left",
            va="bottom",
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
