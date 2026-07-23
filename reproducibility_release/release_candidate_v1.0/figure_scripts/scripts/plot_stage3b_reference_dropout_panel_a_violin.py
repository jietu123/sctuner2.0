#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


PANEL_COLORS = {
    "Marker-defined target region": "#5B3F92",
    "Stage3B blank region": "#00A6B4",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot Panel A violin comparison for Stage3B reference-dropout experiment.")
    p.add_argument("--project_root", default=".")
    p.add_argument(
        "--metadata_csv",
        default="visualizations/stage3b_realdata_candidate_scan/spatial_9x2/stage3b_reference_dropout_spatial_stack_recommended_6x2_metadata.csv",
    )
    p.add_argument(
        "--out_dir",
        default="visualizations/stage3b_realdata_candidate_scan/panel_a_violin",
    )
    p.add_argument("--target_quantile", type=float, default=0.85)
    return p.parse_args()


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


def scene_label(row: pd.Series) -> str:
    mapping = {
        "mouse_embryo_real": "Mouse embryo",
        "brca_tnbc_fresh_frozen": "BRCA TNBC",
        "crc_fresh_frozen": "CRC",
        "brca_her2_ffpe": "BRCA HER2 FFPE",
        "human_intestine_cancer_real": "Intestine cancer",
    }
    return f"{mapping.get(str(row['pair_id']), str(row['pair_id']).replace('_', ' '))}: {row['target_type']}"


def build_long_table(root: Path, metadata: pd.DataFrame, target_quantile: float) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for scene_idx, (_, row) in enumerate(metadata.iterrows(), start=1):
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
        percentile = read_marker_percentile(source_export / "st_expression_normalized.csv", genes)
        blank = read_blank_mask(root, row)
        common = percentile.index.intersection(blank.index)
        percentile = percentile.loc[common]
        blank = blank.loc[common]

        target_mask = percentile >= float(percentile.quantile(target_quantile))
        scene = scene_label(row)
        target_df = pd.DataFrame(
            {
                "scene": scene,
                "scene_index": scene_idx,
                "group": "Marker-defined target region",
                "marker_percentile": percentile.loc[target_mask].to_numpy(),
            }
        )
        blank_df = pd.DataFrame(
            {
                "scene": scene,
                "scene_index": scene_idx,
                "group": "Stage3B blank region",
                "marker_percentile": percentile.loc[blank].to_numpy(),
            }
        )
        rows.extend([target_df, blank_df])
    return pd.concat(rows, ignore_index=True)


def add_violin(ax, values: np.ndarray, pos: float, color: str) -> None:
    parts = ax.violinplot(
        [values],
        positions=[pos],
        widths=0.72,
        showmeans=False,
        showmedians=False,
        showextrema=False,
        bw_method=0.18,
    )
    body = parts["bodies"][0]
    body.set_facecolor(color)
    body.set_edgecolor("#222222")
    body.set_alpha(0.78)
    body.set_linewidth(0.9)

    q1, med, q3 = np.quantile(values, [0.25, 0.5, 0.75])
    ax.plot([pos - 0.18, pos + 0.18], [med, med], color="#111111", lw=1.6, zorder=5)
    ax.plot([pos, pos], [q1, q3], color="#111111", lw=4.2, solid_capstyle="round", zorder=5)


def build_scene_summary(long_df: pd.DataFrame) -> pd.DataFrame:
    return (
        long_df.groupby(["scene_index", "scene", "group"], as_index=False)["marker_percentile"]
        .agg(
            n_spots="count",
            mean_marker_percentile="mean",
            median_marker_percentile="median",
            q25_marker_percentile=lambda x: x.quantile(0.25),
            q75_marker_percentile=lambda x: x.quantile(0.75),
        )
    )


def plot_panel(long_df: pd.DataFrame, out_dir: Path) -> None:
    order = ["Marker-defined target region", "Stage3B blank region"]
    fig, ax = plt.subplots(figsize=(4.4, 4.15), dpi=320)

    for i, group in enumerate(order, start=1):
        values = long_df.loc[long_df["group"].eq(group), "marker_percentile"].to_numpy(dtype=float)
        add_violin(ax, values, i, PANEL_COLORS[group])

    # Scene-level medians make the multi-scene aggregation visible without overcrowding.
    scene_medians = long_df.groupby(["scene_index", "group"], as_index=False)["marker_percentile"].median()
    rng = np.random.default_rng(11)
    for i, group in enumerate(order, start=1):
        vals = scene_medians.loc[scene_medians["group"].eq(group), "marker_percentile"].to_numpy(dtype=float)
        jitter = rng.normal(0, 0.025, size=len(vals))
        ax.scatter(
            np.full(len(vals), i) + jitter,
            vals,
            s=22,
            c="white",
            edgecolors="#222222",
            linewidths=0.65,
            zorder=7,
        )

    ax.set_ylim(0.0, 1.02)
    ax.set_xlim(0.45, 2.55)
    ax.set_ylabel("Target marker score percentile", fontsize=10.2)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Marker-defined\ntarget region", "Stage3B\nblank region"], fontsize=9.2)
    ax.tick_params(axis="y", labelsize=8.8)
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.8, alpha=0.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("Panel A  Stage3B blank regions preserve target-marker enrichment", loc="left", fontsize=10.8, fontweight="bold", pad=8)

    legend_items = [
        Patch(facecolor=PANEL_COLORS[name], edgecolor="#222222", alpha=0.78, label=name)
        for name in order
    ]
    ax.legend(handles=legend_items, frameon=False, fontsize=7.5, loc="lower left", bbox_to_anchor=(0.01, 0.02))

    fig.tight_layout()
    png = out_dir / "stage3b_reference_dropout_panel_a_marker_percentile_violin.png"
    pdf = out_dir / "stage3b_reference_dropout_panel_a_marker_percentile_violin.pdf"
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
    long_df = build_long_table(root, metadata, args.target_quantile)
    long_path = out_dir / "stage3b_reference_dropout_panel_a_marker_percentile_values.csv"
    summary_path = out_dir / "stage3b_reference_dropout_panel_a_marker_percentile_summary.csv"
    scene_summary_path = out_dir / "stage3b_reference_dropout_panel_a_scene_median_summary.csv"
    long_df.to_csv(long_path, index=False)
    summary = long_df.groupby("group")["marker_percentile"].agg(["count", "mean", "median", "std", "min", "max"]).reset_index()
    summary.to_csv(summary_path, index=False)
    scene_summary = build_scene_summary(long_df)
    scene_summary.to_csv(scene_summary_path, index=False)
    plot_panel(long_df, out_dir)
    print(f"[done] {long_path}")
    print(f"[done] {summary_path}")
    print(f"[done] {scene_summary_path}")
    print(scene_summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
