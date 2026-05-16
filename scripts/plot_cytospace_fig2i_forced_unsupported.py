from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import pearsonr, t


X_COL = "Ground truth distance to state 32 (Fig. 2i)c"
STATE_COL = "Cell type/state numbera"
TYPE_COL = "Epithelial type"

METHODS = [
    ("stage4_cytospace_baseline", "CytoSPACE", "#ef6a5b", -0.08, "o"),
    ("stage4_cytospace_route2", "SVTuner + CytoSPACE", "#2f9a8f", 0.08, "o"),
]


def _format_p(pval: float) -> str:
    mantissa, exponent = f"{pval:.1e}".split("e")
    return rf"$P$ = {mantissa} x 10$^{{{int(exponent)}}}$"


def _linear_fit_with_ci(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    slope, intercept = np.polyfit(x, y, deg=1)
    x_grid = np.linspace(0, 15, 150)
    y_hat = slope * x_grid + intercept
    n = len(x)
    y_fit = slope * x + intercept
    dof = max(n - 2, 1)
    s_err = np.sqrt(np.sum((y - y_fit) ** 2) / dof)
    x_mean = np.mean(x)
    sxx = np.sum((x - x_mean) ** 2)
    se_mean = s_err * np.sqrt(1 / n + (x_grid - x_mean) ** 2 / sxx)
    ci = t.ppf(0.975, dof) * se_mean
    return x_grid, y_hat, ci


def _state_list(label: object) -> list[str]:
    return [f"State_{x.strip()}" for x in str(label).split(",")]


def _load_coords(root: Path) -> pd.DataFrame:
    path = root / "data" / "raw" / "low_resolution" / "Adult Mouse Kidney" / "spatial" / "tissue_positions_list.csv"
    coords = pd.read_csv(
        path,
        header=None,
        names=["SpotID", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"],
    )
    return coords[["SpotID", "array_row", "array_col"]]


def _load_assigned(root: Path, sample: str, stage4_dir: str) -> pd.DataFrame:
    path = root / "result" / sample / stage4_dir / "cytospace_output" / "assigned_locations.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    required = {"CellType", "SpotID"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {sorted(missing)}")
    coords = _load_coords(root)
    df = df.drop(columns=[c for c in ["row", "col"] if c in df.columns])
    df = df.merge(coords, on="SpotID", how="left").rename(columns={"array_row": "row", "array_col": "col"})
    if df[["row", "col"]].isna().any().any():
        raise ValueError(f"Missing array coordinates after merging {path}")
    return df


def _compute_distances(root: Path, sample: str, table_s8: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for stage4_dir, method, _, _, _ in METHODS:
        assigned = _load_assigned(root, sample, stage4_dir)
        anchor = assigned[assigned["CellType"].astype(str).eq("State_32")][["row", "col"]].to_numpy(dtype=float)
        if len(anchor) < 5:
            raise ValueError(f"{method} has too few State_32 assignments: {len(anchor)}")
        for _, row in table_s8.iterrows():
            states = _state_list(row[STATE_COL])
            xy = assigned[assigned["CellType"].astype(str).isin(states)][["row", "col"]].to_numpy(dtype=float)
            if len(xy) == 0:
                mean_distance = np.nan
                n_assignments = 0
            else:
                dmat = cdist(xy, anchor)
                k = min(5, dmat.shape[1])
                mean_distance = float(np.sort(dmat, axis=1)[:, :k].mean())
                n_assignments = int(len(xy))
            rows.append(
                {
                    "method": method,
                    "stage4_dir": stage4_dir,
                    "epithelial_type": row[TYPE_COL],
                    "state_label": str(row[STATE_COL]),
                    "states": ";".join(states),
                    "known_distance_rank": float(row[X_COL]),
                    "predicted_mean_nearest5_distance_to_state32": mean_distance,
                    "n_assignments": n_assignments,
                }
            )
    return pd.DataFrame(rows)


def _plot_panel(ax: plt.Axes, df: pd.DataFrame, epithelial_type: str) -> None:
    sub_type = df[df["epithelial_type"].eq(epithelial_type)].copy()
    label_lines = []
    for _, method, color, jitter, marker in METHODS:
        sub = sub_type[sub_type["method"].eq(method)].dropna(subset=["predicted_mean_nearest5_distance_to_state32"])
        x = sub["known_distance_rank"].to_numpy(dtype=float)
        y = sub["predicted_mean_nearest5_distance_to_state32"].to_numpy(dtype=float)
        r, pval = pearsonr(x, y)
        label_lines.append((method, r, pval))
        x_grid, y_hat, ci = _linear_fit_with_ci(x, y)
        ax.fill_between(x_grid, y_hat - ci, y_hat + ci, color="#c7c7c7", alpha=0.42, linewidth=0)
        ax.plot(x_grid, y_hat, color=color, linewidth=1.25, alpha=0.95, label=method)
        ax.scatter(
            x + jitter,
            y,
            s=10,
            marker=marker,
            color="#111111",
            edgecolor="#111111",
            linewidth=0,
            alpha=0.90,
            zorder=3,
        )
    ax.set_title(f"{epithelial_type} epithelium", fontsize=7.8, pad=5)
    text = "\n".join([f"{m}: " + rf"$r$ = {r:.2f}" for m, r, _ in label_lines])
    ax.text(0.05, 0.93, text, transform=ax.transAxes, fontsize=5.4, va="top", ha="left")


def _plot_single_method_panel(ax: plt.Axes, df: pd.DataFrame, epithelial_type: str, method: str) -> None:
    method_meta = {m: (color, marker) for _, m, color, _, marker in METHODS}
    color, marker = method_meta[method]
    sub = df[
        df["epithelial_type"].eq(epithelial_type) & df["method"].eq(method)
    ].dropna(subset=["predicted_mean_nearest5_distance_to_state32"])
    x = sub["known_distance_rank"].to_numpy(dtype=float)
    y = sub["predicted_mean_nearest5_distance_to_state32"].to_numpy(dtype=float)
    r, pval = pearsonr(x, y)
    x_grid, y_hat, ci = _linear_fit_with_ci(x, y)
    ax.fill_between(x_grid, y_hat - ci, y_hat + ci, color="#c7c7c7", alpha=0.42, linewidth=0)
    ax.plot(x_grid, y_hat, color=color, linewidth=1.35, alpha=0.95)
    ax.scatter(
        x,
        y,
        s=10,
        marker=marker,
        color="#111111",
        edgecolor="#111111",
        linewidth=0,
        alpha=0.90,
        zorder=3,
    )
    ax.set_title(f"{epithelial_type}\n{method}", fontsize=6.8, pad=4)
    ax.text(
        0.06,
        0.90,
        rf"$r$ = {r:.2f}",
        transform=ax.transAxes,
        fontsize=5.8,
        va="top",
        ha="left",
    )


def _apply_axis_style(ax: plt.Axes) -> None:
    ax.set_xlim(0, 15.5)
    ax.set_ylim(0, 32)
    ax.set_xticks([0, 5, 10, 15])
    ax.set_yticks([0, 10, 20, 30])
    ax.tick_params(axis="both", labelsize=6.4, width=0.8, length=2.8, pad=1.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot Fig.2i-style forced unsupported kidney benchmark.")
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--sample", default="cytospace_fig2i_mouse_kidney_forced_unsupported")
    parser.add_argument(
        "--source_xlsx",
        default="data/raw/cytospace_fig2c_melanoma/41587_2023_1697_MOESM3_ESM.xlsx",
    )
    parser.add_argument("--out_dir", default="visualizations/cytospace_fig2i_mouse_kidney_forced_unsupported")
    parser.add_argument("--out_prefix", default="fig2i_forced_unsupported_baseline_vs_route2")
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    out_dir = (root / args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    table_s8 = pd.read_excel(root / args.source_xlsx, sheet_name="Table S8", header=7)
    distances = _compute_distances(root, args.sample, table_s8)
    distances.to_csv(out_dir / f"{args.out_prefix}_source_values.csv", index=False)

    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(5.2, 2.65), dpi=300, sharey=True)
    _plot_panel(axes[0], distances, "Nephron")
    _plot_panel(axes[1], distances, "Ureteric")
    for ax in axes:
        _apply_axis_style(ax)
    axes[0].set_ylabel("Predicted distance to State 32\n(mean nearest-5 Euclidean)", fontsize=7.0)
    axes[1].legend(frameon=False, fontsize=5.8, loc="lower right", handlelength=1.7)
    fig.text(0.56, 0.18, "Known distance (rank)", ha="center", va="center", fontsize=7.2)
    fig.suptitle(
        "Kidney epithelial spatial ordering under unsupported-state perturbation",
        fontsize=8.0,
        y=0.98,
    )
    overlay = fig.add_axes([0, 0, 1, 1], frameon=False)
    overlay.set_axis_off()
    overlay.annotate(
        "",
        xy=(0.045, 0.62),
        xytext=(0.045, 0.38),
        xycoords="figure fraction",
        arrowprops=dict(arrowstyle="<->", color="#8a8a8a", lw=1.1),
        annotation_clip=False,
    )
    fig.text(0.045, 0.64, "Cortex", ha="center", va="bottom", fontsize=6.2)
    fig.text(0.045, 0.31, "Base\nof inner\nmedulla", ha="center", va="top", fontsize=6.2)
    overlay.annotate(
        "",
        xy=(0.78, 0.115),
        xytext=(0.30, 0.115),
        xycoords="figure fraction",
        arrowprops=dict(arrowstyle="<->", color="#8a8a8a", lw=1.1),
        annotation_clip=False,
    )
    fig.text(0.25, 0.085, "Base\nof inner\nmedulla", ha="center", va="top", fontsize=6.0)
    fig.text(0.82, 0.106, "Cortex", ha="center", va="center", fontsize=6.0)
    fig.subplots_adjust(left=0.20, right=0.98, top=0.78, bottom=0.37, wspace=0.34)

    png = out_dir / f"{args.out_prefix}.png"
    pdf = out_dir / f"{args.out_prefix}.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)

    fig4, axes4 = plt.subplots(2, 2, figsize=(4.9, 3.9), dpi=300, sharex=True, sharey=True)
    panel_specs = [
        ("Nephron", "CytoSPACE"),
        ("Nephron", "SVTuner + CytoSPACE"),
        ("Ureteric", "CytoSPACE"),
        ("Ureteric", "SVTuner + CytoSPACE"),
    ]
    for ax, (etype, method) in zip(axes4.ravel(), panel_specs):
        _plot_single_method_panel(ax, distances, etype, method)
        _apply_axis_style(ax)
    axes4[0, 0].set_ylabel("Predicted distance to State 32\n(mean nearest-5 Euclidean)", fontsize=6.8)
    axes4[1, 0].set_ylabel("Predicted distance to State 32\n(mean nearest-5 Euclidean)", fontsize=6.8)
    fig4.text(0.55, 0.075, "Known distance (rank)", ha="center", va="center", fontsize=7.2)
    fig4.suptitle(
        "Kidney epithelial spatial ordering under unsupported-state perturbation",
        fontsize=8.0,
        y=0.985,
    )
    fig4.subplots_adjust(left=0.17, right=0.98, top=0.86, bottom=0.15, wspace=0.23, hspace=0.46)
    png4 = out_dir / f"{args.out_prefix}_four_panel.png"
    pdf4 = out_dir / f"{args.out_prefix}_four_panel.pdf"
    fig4.savefig(png4, dpi=300)
    fig4.savefig(pdf4)
    plt.close(fig4)

    print(f"[OK] wrote: {png}")
    print(f"[OK] wrote: {pdf}")
    print(f"[OK] wrote: {png4}")
    print(f"[OK] wrote: {pdf4}")
    for (method, etype), sub in distances.groupby(["method", "epithelial_type"]):
        sub = sub.dropna(subset=["predicted_mean_nearest5_distance_to_state32"])
        r, pval = pearsonr(sub["known_distance_rank"], sub["predicted_mean_nearest5_distance_to_state32"])
        print(f"[INFO] {method} {etype}: n={len(sub)} r={r:.4f} p={pval:.4g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
