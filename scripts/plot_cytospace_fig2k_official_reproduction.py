from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr


def _load_fig2k_data(xlsx_path: Path) -> pd.DataFrame:
    raw = pd.read_excel(xlsx_path, sheet_name="Table S9", header=None)
    df = raw.iloc[8:, 0:4].copy()
    df.columns = ["known_rank", "state", "cytospace_nes", "merscope_nes"]
    df = df.dropna(subset=["known_rank", "state"]).copy()
    for col in ["known_rank", "cytospace_nes", "merscope_nes"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["known_rank", "cytospace_nes", "merscope_nes"]).copy()
    df["known_rank"] = df["known_rank"].astype(int)
    df["cytospace_pred_rank"] = df["cytospace_nes"].rank(method="average", ascending=True)
    df["merscope_pred_rank"] = df["merscope_nes"].rank(method="average", ascending=True)
    return df.sort_values("known_rank").reset_index(drop=True)


def _format_p(pval: float) -> str:
    if pval >= 0.001:
        return rf"$P$ = {pval:.3f}".rstrip("0").rstrip(".")
    mantissa, exponent = f"{pval:.1e}".split("e")
    return rf"$P$ = {mantissa} x 10$^{{{int(exponent)}}}$"


def _fit_line(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    slope, intercept = np.polyfit(x, y, deg=1)
    x_grid = np.linspace(1, 23, 100)
    return x_grid, slope * x_grid + intercept


def _panel(ax: plt.Axes, df: pd.DataFrame, pred_col: str, title: str, colors: np.ndarray) -> tuple[float, float]:
    x = df["known_rank"].to_numpy(dtype=float)
    y = df[pred_col].to_numpy(dtype=float)
    r, pval = pearsonr(x, y)
    x_grid, y_grid = _fit_line(x, y)

    ax.scatter(x, y, s=7.5, c=colors, edgecolor="none", zorder=3)
    ax.plot(x_grid, y_grid, color="#2a2a2a", linewidth=0.75, zorder=2)
    ax.set_title(title, fontsize=5.9, pad=3)
    ax.set_xlim(0.5, 23.5)
    ax.set_ylim(0.5, 23.5)
    ax.set_aspect("equal", adjustable="box")

    ax.set_xticks(df["known_rank"])
    ax.set_xticklabels(df["state"], rotation=90, fontsize=2.0)
    y_order = df.sort_values(pred_col)["state"].to_list()
    ax.set_yticks(np.arange(1, 24))
    ax.set_yticklabels(y_order, fontsize=2.0)
    ax.tick_params(axis="both", length=1.0, width=0.45, pad=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.65)
    ax.spines["bottom"].set_linewidth(0.65)

    ax.text(0.06, 0.92, rf"$r$ = {r:.2f}" + "\n" + _format_p(pval), transform=ax.transAxes,
            ha="left", va="top", fontsize=5.1)
    return float(r), float(pval)


def main() -> int:
    parser = argparse.ArgumentParser(description="Reproduce CytoSPACE paper Fig. 2k from Supplementary Table S9.")
    parser.add_argument("--project_root", default=".")
    parser.add_argument(
        "--source_xlsx",
        default="data/raw/cytospace_fig2c_melanoma/41587_2023_1697_MOESM3_ESM.xlsx",
    )
    parser.add_argument("--out_dir", default="visualizations/cytospace_fig2k_tcell_states")
    parser.add_argument("--out_prefix", default="fig2k_official_reproduction")
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    out_dir = (root / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    df = _load_fig2k_data((root / args.source_xlsx).resolve())
    df.to_csv(out_dir / f"{args.out_prefix}_source_values.csv", index=False)

    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
    cmap = plt.get_cmap("viridis")
    colors = cmap(np.linspace(0.05, 0.95, len(df)))

    fig, axes = plt.subplots(1, 2, figsize=(3.65, 2.35), dpi=300)
    stats = {}
    stats["cytospace"] = _panel(
        axes[0],
        df,
        "cytospace_pred_rank",
        "scRNA-seq mapped to MERSCOPE\n(whole transcriptome)",
        colors,
    )
    stats["merscope"] = _panel(
        axes[1],
        df,
        "merscope_pred_rank",
        "MERSCOPE\n($n$ = 500 genes)",
        colors,
    )

    axes[0].set_ylabel("Predicted tumor enrichment rank", fontsize=5.7, labelpad=4)
    fig.text(0.56, 0.215, "Known tumor enrichment rank\n(Zheng et al.)", ha="center", va="center", fontsize=5.7)
    fig.suptitle("Tumor/normal enrichment of CD4 T cell states ($n$ = 23)", fontsize=6.4, y=0.965)

    overlay = fig.add_axes([0, 0, 1, 1], frameon=False)
    overlay.set_axis_off()
    overlay.annotate("", xy=(0.105, 0.68), xytext=(0.105, 0.48), xycoords="figure fraction",
                     arrowprops=dict(arrowstyle="<->", color="#8a8a8a", lw=0.8), annotation_clip=False)
    fig.text(0.048, 0.70, "Higher in\ntumor", ha="center", va="bottom", fontsize=5.0)
    fig.text(0.048, 0.46, "Higher in\nadjacent\nnormal", ha="center", va="top", fontsize=5.0)
    overlay.annotate("", xy=(0.82, 0.15), xytext=(0.38, 0.15), xycoords="figure fraction",
                     arrowprops=dict(arrowstyle="<->", color="#8a8a8a", lw=0.8), annotation_clip=False)
    fig.text(0.33, 0.105, "Higher in\nadjacent\nnormal", ha="center", va="top", fontsize=5.0)
    fig.text(0.85, 0.105, "Higher in\ntumor", ha="center", va="top", fontsize=5.0)
    fig.text(0.02, 0.97, "k", fontsize=10.5, fontweight="bold", ha="left", va="top")

    fig.subplots_adjust(left=0.22, right=0.98, top=0.78, bottom=0.36, wspace=0.42)
    png = out_dir / f"{args.out_prefix}.png"
    pdf = out_dir / f"{args.out_prefix}.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)

    print(f"[OK] wrote: {png}")
    print(f"[OK] wrote: {pdf}")
    for name, (r, pval) in stats.items():
        print(f"[INFO] {name}: n={len(df)} r={r:.4f} p={pval:.4g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
