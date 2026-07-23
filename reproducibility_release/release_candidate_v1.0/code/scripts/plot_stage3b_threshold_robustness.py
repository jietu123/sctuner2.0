#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import Normalize
from scipy.spatial import cKDTree


THRESHOLDS = [0.10, 0.15, 0.20, 0.25]
GROUP_COLORS = {
    "Inside target core": "#4C2A85",
    "Near target core": "#2CB7B8",
    "Outside target-associated region": "#D9D9D9",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot Stage3B target-region threshold robustness summary.")
    p.add_argument("--project_root", default=".")
    p.add_argument(
        "--metadata_csv",
        default="visualizations/stage3b_realdata_candidate_scan/spatial_9x2/stage3b_reference_dropout_spatial_stack_recommended_6x2_metadata.csv",
    )
    p.add_argument(
        "--out_dir",
        default="visualizations/stage3b_realdata_candidate_scan/stage3b_threshold_robustness",
    )
    p.add_argument("--neighbor_rings", type=float, default=2.5)
    return p.parse_args()


def read_coords(path: Path) -> pd.DataFrame:
    coords = pd.read_csv(path)
    id_col = "spot_id" if "spot_id" in coords.columns else coords.columns[0]
    coords[id_col] = coords[id_col].astype(str)
    coords = coords.drop_duplicates(id_col).set_index(id_col)
    lower = {c.lower(): c for c in coords.columns}
    if "col" in lower and "row" in lower:
        x_col, y_col = lower["col"], lower["row"]
    elif "x" in lower and "y" in lower:
        x_col, y_col = lower["x"], lower["y"]
    else:
        numeric = [c for c in coords.columns if pd.api.types.is_numeric_dtype(coords[c])]
        if len(numeric) < 2:
            raise ValueError(f"Cannot infer coordinate columns: {path}")
        y_col, x_col = numeric[0], numeric[1]
    out = coords[[x_col, y_col]].apply(pd.to_numeric, errors="coerce").dropna()
    out.columns = ["x", "y"]
    return out


def read_marker_percentile(expr_path: Path, genes: list[str]) -> pd.Series:
    header = pd.read_csv(expr_path, nrows=0).columns.tolist()
    id_col = header[0]
    present = [g for g in genes if g in header]
    if not present:
        raise ValueError(f"No marker genes present in {expr_path}")
    expr = pd.read_csv(expr_path, index_col=0, usecols=[id_col, *present])
    expr.index = expr.index.astype(str)
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    score = expr.mean(axis=1)
    return score.rank(method="average", pct=True)


def read_blank_mask(root: Path, row: pd.Series) -> pd.Series:
    scores_path = (
        root
        / "data"
        / "processed"
        / str(row["stage3b_group"])
        / str(row["sample"])
        / "stage3b_st_unsupported"
        / "spot_unsupported_scores.csv"
    )
    scores = pd.read_csv(scores_path, index_col=0)
    scores.index = scores.index.astype(str)
    blank = scores["is_unsupported_region"]
    if blank.dtype != bool:
        blank = blank.astype(str).str.lower().isin(["true", "1", "yes"])
    return blank


def display_label(row: pd.Series) -> str:
    mapping = {
        "mouse_embryo_real": "Mouse embryo",
        "brca_tnbc_fresh_frozen": "BRCA TNBC",
        "crc_fresh_frozen": "CRC",
        "brca_her2_ffpe": "BRCA HER2 FFPE",
        "human_intestine_cancer_real": "Intestine cancer",
    }
    pair = str(row["pair_id"])
    return f"{mapping.get(pair, pair.replace('_', ' '))} ({row['target_type']})"


def nearest_neighbor_spacing(xy: np.ndarray) -> float:
    tree = cKDTree(xy)
    distances, _ = tree.query(xy, k=2)
    nn = distances[:, 1]
    nn = nn[np.isfinite(nn) & (nn > 0)]
    if nn.size == 0:
        return 1.0
    return float(np.median(nn))


def scene_table(root: Path, row: pd.Series, neighbor_rings: float) -> pd.DataFrame:
    source_export = (
        root
        / "data"
        / "processed"
        / str(row["storage_group"])
        / str(row["source_sample"])
        / "stage1_preprocess"
        / "exported"
    )
    genes = [g for g in str(row["marker_genes_list"]).split(";") if g]
    coords = read_coords(source_export / "st_coordinates.csv")
    percentile = read_marker_percentile(source_export / "st_expression_normalized.csv", genes)
    blank = read_blank_mask(root, row)
    common = coords.index.intersection(percentile.index).intersection(blank.index)
    coords = coords.loc[common]
    percentile = percentile.loc[common]
    blank = blank.loc[common]

    xy = coords[["x", "y"]].to_numpy(dtype=float)
    radius = nearest_neighbor_spacing(xy) * neighbor_rings
    blank_total = int(blank.sum())
    scene = display_label(row)
    rows = []
    for top_fraction in THRESHOLDS:
        quantile = 1.0 - top_fraction
        target_core = percentile >= float(percentile.quantile(quantile))
        target_xy = coords.loc[target_core, ["x", "y"]].to_numpy(dtype=float)
        if len(target_xy) == 0:
            near_core = pd.Series(False, index=coords.index)
        else:
            tree = cKDTree(target_xy)
            dist_to_core, _ = tree.query(xy, k=1)
            near_core = pd.Series(dist_to_core <= radius, index=coords.index)

        inside = int((blank & target_core).sum())
        near = int((blank & ~target_core & near_core).sum())
        outside = int((blank & ~target_core & ~near_core).sum())
        target_core_n = int(target_core.sum())
        associated = inside + near
        rows.append(
            {
                "scene": scene,
                "pair_id": str(row["pair_id"]),
                "target_type": str(row["target_type"]),
                "top_fraction": top_fraction,
                "threshold_label": f"top{int(round(top_fraction * 100))}",
                "n_spots": int(len(coords)),
                "target_core_spots": target_core_n,
                "stage3b_blank_spots": blank_total,
                "inside_target_core_spots": inside,
                "near_target_core_spots": near,
                "outside_target_associated_spots": outside,
                "target_associated_spots": associated,
                "inside_fraction": inside / blank_total if blank_total else 0.0,
                "near_fraction": near / blank_total if blank_total else 0.0,
                "outside_fraction": outside / blank_total if blank_total else 0.0,
                "target_associated_fraction": associated / blank_total if blank_total else 0.0,
                "core_recall": inside / target_core_n if target_core_n else 0.0,
            }
        )
    return pd.DataFrame(rows)


def plot_figure(detail: pd.DataFrame, out_dir: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.edgecolor": "#333333",
            "axes.linewidth": 1.0,
        }
    )
    sns.set_theme(style="white")

    threshold_order = [f"top{int(round(x * 100))}" for x in THRESHOLDS]
    scene_order = detail.drop_duplicates("scene")["scene"].tolist()

    fig = plt.figure(figsize=(14.8, 8.4), dpi=300, constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=[1.38, 1.0], height_ratios=[0.92, 1.08])
    ax_a = fig.add_subplot(grid[:, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 1])

    x_map = {label: i for i, label in enumerate(threshold_order)}
    y_map = {scene: i for i, scene in enumerate(scene_order)}
    plot_df = detail.copy()
    plot_df["x"] = plot_df["threshold_label"].map(x_map)
    plot_df["y"] = plot_df["scene"].map(y_map)
    size_min = 0.80
    size_max = 1.00
    size_scaled = ((plot_df["target_associated_fraction"] - size_min) / (size_max - size_min)).clip(0, 1)
    sizes = 120 + 720 * size_scaled
    sc = ax_a.scatter(
        plot_df["x"],
        plot_df["y"],
        s=sizes,
        c=plot_df["inside_fraction"],
        cmap="viridis",
        norm=Normalize(0, 1),
        edgecolors="#2b2b2b",
        linewidths=0.55,
    )
    ax_a.set_xticks(range(len(threshold_order)))
    ax_a.set_xticklabels([x.replace("top", "top ") + "%" for x in threshold_order], fontsize=10)
    ax_a.set_yticks(range(len(scene_order)))
    ax_a.set_yticklabels(scene_order, fontsize=9.2)
    ax_a.invert_yaxis()
    ax_a.set_xlabel("Marker-defined target-region threshold", fontsize=10.5)
    ax_a.set_ylabel("")
    ax_a.set_title("A. Stage3B blank-region robustness across target definitions", loc="left", fontsize=12, weight="bold", pad=10)
    ax_a.grid(color="#E4E4E4", linewidth=0.8)
    ax_a.set_axisbelow(True)
    cbar = fig.colorbar(sc, ax=ax_a, shrink=0.72, pad=0.02)
    cbar.set_label("Inside-core fraction", fontsize=9.5)
    cbar.ax.tick_params(labelsize=8.5)
    def bubble_size(value: float) -> float:
        scaled = np.clip((value - size_min) / (size_max - size_min), 0.0, 1.0)
        return 120 + 720 * scaled

    handles = [
        ax_a.scatter([], [], s=bubble_size(frac), c="#BDBDBD", edgecolors="#2b2b2b", linewidths=0.55, label=f"{frac:.1f}")
        for frac in [0.8, 0.9, 1.0]
    ]
    ax_a.legend(
        handles=handles,
        title="Target-associated fraction",
        frameon=False,
        fontsize=8.2,
        title_fontsize=8.4,
        loc="upper center",
        bbox_to_anchor=(0.52, -0.15),
        ncol=4,
        handletextpad=0.9,
        columnspacing=1.35,
        borderaxespad=0.0,
    )

    sns.boxplot(
        data=detail,
        x="threshold_label",
        y="target_associated_fraction",
        order=threshold_order,
        ax=ax_b,
        color="#2CB7B8",
        width=0.55,
        fliersize=0,
        linewidth=1.0,
    )
    sns.stripplot(
        data=detail,
        x="threshold_label",
        y="target_associated_fraction",
        order=threshold_order,
        ax=ax_b,
        color="#222222",
        size=4.2,
        jitter=0.18,
        alpha=0.75,
    )
    ax_b.set_ylim(0.70, 1.015)
    ax_b.set_xlabel("")
    ax_b.set_ylabel("Target-associated fraction", fontsize=10)
    ax_b.set_title("B. Stability summary", loc="left", fontsize=11, weight="bold", pad=8)
    ax_b.grid(axis="y", color="#E0E0E0", linewidth=0.8)
    ax_b.set_axisbelow(True)
    ax_b.tick_params(axis="x", labelsize=9.2)
    ax_b.tick_params(axis="y", labelsize=8.5)

    comp = (
        detail.groupby("threshold_label")[["inside_fraction", "near_fraction", "outside_fraction"]]
        .mean()
        .reindex(threshold_order)
        .reset_index()
    )
    x = np.arange(len(comp))
    bottom = np.zeros(len(comp))
    parts = [
        ("inside_fraction", "Inside target core", GROUP_COLORS["Inside target core"]),
        ("near_fraction", "Near target core", GROUP_COLORS["Near target core"]),
        ("outside_fraction", "Outside target-associated region", GROUP_COLORS["Outside target-associated region"]),
    ]
    for col, label, color in parts:
        vals = comp[col].to_numpy(dtype=float)
        ax_c.bar(x, vals, bottom=bottom, color=color, edgecolor="white", linewidth=0.8, label=label)
        bottom += vals
    ax_c.set_xticks(x)
    ax_c.set_xticklabels([x.replace("top", "top ") + "%" for x in threshold_order], fontsize=9.2)
    ax_c.set_ylim(0, 1)
    ax_c.set_ylabel("Mean blank-region composition", fontsize=10)
    ax_c.set_title("C. Where Stage3B blank spots fall", loc="left", fontsize=11, weight="bold", pad=8)
    ax_c.legend(frameon=False, fontsize=7.7, loc="upper center", bbox_to_anchor=(0.50, -0.18), ncol=1)
    ax_c.tick_params(axis="y", labelsize=8.5)
    ax_c.grid(axis="y", color="#E0E0E0", linewidth=0.8)
    ax_c.set_axisbelow(True)

    for ax in [ax_a, ax_b, ax_c]:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    png = out_dir / "stage3b_threshold_robustness_summary.png"
    pdf = out_dir / "stage3b_threshold_robustness_summary.pdf"
    fig.savefig(png, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[done] {png}")
    print(f"[done] {pdf}")


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(root / args.metadata_csv)
    detail = pd.concat([scene_table(root, row, args.neighbor_rings) for _, row in metadata.iterrows()], ignore_index=True)
    summary = (
        detail.groupby("threshold_label")
        .agg(
            scenes=("scene", "count"),
            mean_target_associated_fraction=("target_associated_fraction", "mean"),
            median_target_associated_fraction=("target_associated_fraction", "median"),
            mean_inside_fraction=("inside_fraction", "mean"),
            mean_near_fraction=("near_fraction", "mean"),
            mean_outside_fraction=("outside_fraction", "mean"),
            mean_core_recall=("core_recall", "mean"),
        )
        .reset_index()
    )
    detail_path = out_dir / "stage3b_threshold_robustness_detail.csv"
    summary_path = out_dir / "stage3b_threshold_robustness_summary.csv"
    detail.to_csv(detail_path, index=False)
    summary.to_csv(summary_path, index=False)
    plot_figure(detail, out_dir)
    print(f"[done] {detail_path}")
    print(f"[done] {summary_path}")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
