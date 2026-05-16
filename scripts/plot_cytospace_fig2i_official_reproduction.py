from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, t


X_COL = "Ground truth distance to state 32 (Fig. 2i)c"
Y_COL = "Mean Euclidean distance to state 32 by CytoSPACE (base of inner medulla)"
TYPE_COL = "Epithelial type"


def _load_table_s8(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="Table S8", header=7)
    required = {TYPE_COL, X_COL, Y_COL, "Cell type/state numbera", "Epithelial cell statea", "Zone"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Table S8 missing columns: {sorted(missing)}")
    df = df.dropna(subset=[TYPE_COL, X_COL, Y_COL]).copy()
    df["known_distance_rank"] = pd.to_numeric(df[X_COL])
    df["predicted_mean_euclidean_distance"] = pd.to_numeric(df[Y_COL])
    return df


def _linear_fit_with_ci(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    slope, intercept = np.polyfit(x, y, deg=1)
    x_grid = np.linspace(0, 15, 150)
    y_hat = slope * x_grid + intercept

    n = len(x)
    y_fit = slope * x + intercept
    resid = y - y_fit
    dof = max(n - 2, 1)
    s_err = np.sqrt(np.sum(resid**2) / dof)
    x_mean = np.mean(x)
    sxx = np.sum((x - x_mean) ** 2)
    se_mean = s_err * np.sqrt(1 / n + (x_grid - x_mean) ** 2 / sxx)
    tcrit = t.ppf(0.975, dof)
    ci = tcrit * se_mean
    return x_grid, y_hat, ci


def _format_p(pval: float) -> str:
    mantissa, exponent = f"{pval:.1e}".split("e")
    exp = int(exponent)
    return rf"$P$ = {mantissa} x 10$^{{{exp}}}$"


def _plot_panel(ax: plt.Axes, df: pd.DataFrame, epithelial_type: str) -> None:
    sub = df[df[TYPE_COL].eq(epithelial_type)].copy()
    x = sub["known_distance_rank"].to_numpy(dtype=float)
    y = sub["predicted_mean_euclidean_distance"].to_numpy(dtype=float)
    r, pval = pearsonr(x, y)
    x_grid, y_hat, ci = _linear_fit_with_ci(x, y)

    ax.fill_between(x_grid, y_hat - ci, y_hat + ci, color="#c9c9c9", alpha=0.55, linewidth=0)
    ax.plot(x_grid, y_hat, color="#2f65c8", linewidth=1.25)
    ax.scatter(x, y, s=8, color="#111111", linewidth=0, zorder=3)

    ax.set_title(f"{epithelial_type} epithelium", fontsize=8.0, pad=5)
    ax.set_xlim(0, 15.5)
    ax.set_ylim(0, 50)
    ax.set_xticks([0, 5, 10, 15])
    ax.set_yticks([0, 10, 20, 30, 40, 50])
    ax.tick_params(axis="both", labelsize=7.0, width=0.8, length=3.0, pad=1.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.text(
        0.05,
        0.92,
        rf"$r$ = {r:.2f}" + "\n" + _format_p(pval),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.0,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Reproduce CytoSPACE paper Fig. 2i from Supplementary Table S8.")
    parser.add_argument("--project_root", default=".")
    parser.add_argument(
        "--source_xlsx",
        default="data/raw/cytospace_fig2c_melanoma/41587_2023_1697_MOESM3_ESM.xlsx",
    )
    parser.add_argument("--out_dir", default="visualizations/cytospace_fig2i_mouse_kidney")
    parser.add_argument("--out_prefix", default="fig2i_official_reproduction")
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    out_dir = (root / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    df = _load_table_s8((root / args.source_xlsx).resolve())
    df.to_csv(out_dir / f"{args.out_prefix}_source_values.csv", index=False)

    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(4.9, 2.75), dpi=300, sharey=True)
    _plot_panel(axes[0], df, "Nephron")
    _plot_panel(axes[1], df, "Ureteric")
    axes[0].set_ylabel("Predicted distance\n(mean Euclidean)", fontsize=8.0)
    fig.text(0.55, 0.235, "Known distance (rank)", ha="center", va="center", fontsize=8.0)

    fig.suptitle(
        "Predicted vs. known distances of epithelial states\nto the base of the inner medulla",
        fontsize=8.0,
        y=0.98,
    )

    overlay = fig.add_axes([0, 0, 1, 1], frameon=False)
    overlay.set_axis_off()
    overlay.annotate(
        "",
        xy=(0.08, 0.61),
        xytext=(0.08, 0.43),
        xycoords="figure fraction",
        arrowprops=dict(arrowstyle="->", color="#9a9a9a", lw=1.2),
        annotation_clip=False,
    )
    fig.text(0.08, 0.65, "Cortex", ha="center", va="bottom", fontsize=7.5)
    fig.text(0.08, 0.36, "Base\nof inner\nmedulla", ha="center", va="top", fontsize=7.5)
    overlay.annotate(
        "",
        xy=(0.79, 0.19),
        xytext=(0.43, 0.19),
        xycoords="figure fraction",
        arrowprops=dict(arrowstyle="<->", color="#9a9a9a", lw=1.2),
        annotation_clip=False,
    )
    fig.text(0.42, 0.12, "Base\nof inner\nmedulla", ha="center", va="top", fontsize=7.0)
    fig.text(0.83, 0.12, "Cortex", ha="center", va="top", fontsize=7.0)

    fig.subplots_adjust(left=0.23, right=0.96, top=0.74, bottom=0.34, wspace=0.5)
    png = out_dir / f"{args.out_prefix}.png"
    pdf = out_dir / f"{args.out_prefix}.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)

    print(f"[OK] wrote: {png}")
    print(f"[OK] wrote: {pdf}")
    for epithelial_type in ["Nephron", "Ureteric"]:
        sub = df[df[TYPE_COL].eq(epithelial_type)]
        r, pval = pearsonr(sub["known_distance_rank"], sub["predicted_mean_euclidean_distance"])
        print(f"[INFO] {epithelial_type}: n={len(sub)} r={r:.4f} p={pval:.4g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
