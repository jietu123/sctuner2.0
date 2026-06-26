#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


GROUP_COLORS = {
    "Inside target core": "#4C2A85",
    "Near target core": "#2CB7B8",
    "Outside target-associated region": "#D9D9D9",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot Stage3B blank-region spatial composition across selected real-data cases.")
    p.add_argument("--project_root", default=".")
    p.add_argument(
        "--metadata_csv",
        default="visualizations/stage3b_realdata_candidate_scan/spatial_9x2/stage3b_reference_dropout_spatial_stack_recommended_6x2_metadata.csv",
    )
    p.add_argument(
        "--out_dir",
        default="visualizations/stage3b_realdata_candidate_scan/panel_a_blank_composition",
    )
    p.add_argument("--target_quantile", type=float, default=0.85)
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
    return f"{mapping.get(pair, pair.replace('_', ' '))}\n{row['target_type']}"


def nearest_neighbor_spacing(xy: np.ndarray) -> float:
    tree = cKDTree(xy)
    distances, _ = tree.query(xy, k=2)
    nn = distances[:, 1]
    nn = nn[np.isfinite(nn) & (nn > 0)]
    if nn.size == 0:
        return 1.0
    return float(np.median(nn))


def load_scene(root: Path, row: pd.Series, target_quantile: float, neighbor_rings: float) -> dict[str, object]:
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

    target_core = percentile >= float(percentile.quantile(target_quantile))
    xy = coords[["x", "y"]].to_numpy(dtype=float)
    spacing = nearest_neighbor_spacing(xy)
    radius = spacing * neighbor_rings

    target_xy = coords.loc[target_core, ["x", "y"]].to_numpy(dtype=float)
    if len(target_xy) == 0:
        near_core = pd.Series(False, index=coords.index)
    else:
        tree = cKDTree(target_xy)
        dist_to_core, _ = tree.query(xy, k=1)
        near_core = pd.Series(dist_to_core <= radius, index=coords.index)

    blank_total = int(blank.sum())
    inside = int((blank & target_core).sum())
    near = int((blank & ~target_core & near_core).sum())
    outside = int((blank & ~target_core & ~near_core).sum())
    if blank_total == 0:
        fractions = {"Inside target core": 0.0, "Near target core": 0.0, "Outside target-associated region": 0.0}
    else:
        fractions = {
            "Inside target core": inside / blank_total,
            "Near target core": near / blank_total,
            "Outside target-associated region": outside / blank_total,
        }

    return {
        "scene": display_label(row),
        "pair_id": str(row["pair_id"]),
        "target_type": str(row["target_type"]),
        "n_spots": len(coords),
        "target_core_spots": int(target_core.sum()),
        "stage3b_blank_spots": blank_total,
        "inside_target_core_spots": inside,
        "near_target_core_spots": near,
        "outside_target_associated_spots": outside,
        "neighbor_spacing": spacing,
        "neighbor_radius": radius,
        **{f"fraction_{k.lower().replace(' ', '_').replace('-', '_')}": v for k, v in fractions.items()},
    }


def plot_panel(summary: pd.DataFrame, out_dir: Path) -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42})
    groups = ["Inside target core", "Near target core", "Outside target-associated region"]
    frac_cols = [f"fraction_{g.lower().replace(' ', '_').replace('-', '_')}" for g in groups]
    y = np.arange(len(summary))

    fig, ax = plt.subplots(figsize=(6.4, 4.25), dpi=320)
    left = np.zeros(len(summary), dtype=float)
    for group, col in zip(groups, frac_cols):
        ax.barh(
            y,
            summary[col].to_numpy(dtype=float),
            left=left,
            height=0.62,
            color=GROUP_COLORS[group],
            edgecolor="white",
            linewidth=0.8,
            label=group,
        )
        left += summary[col].to_numpy(dtype=float)

    ax.set_yticks(y)
    ax.set_yticklabels(summary["scene"], fontsize=8.0)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of Stage3B blank spots", fontsize=9.5)
    ax.set_title("Panel A  Stage3B blank regions concentrate near marker-defined target areas", loc="left", fontsize=10.5, fontweight="bold", pad=8)
    ax.grid(axis="x", color="#DADADA", linewidth=0.8, alpha=0.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="x", labelsize=8.2)
    ax.tick_params(axis="y", length=0)
    ax.legend(frameon=False, fontsize=7.7, loc="lower center", bbox_to_anchor=(0.52, -0.29), ncol=3)

    for row_i, (_, row) in enumerate(summary.iterrows()):
        associated = row["fraction_inside_target_core"] + row["fraction_near_target_core"]
        ax.text(
            min(associated + 0.015, 0.94),
            row_i,
            f"{associated:.2f}",
            va="center",
            ha="left",
            fontsize=7.8,
            color="#1F1F1F",
            fontweight="bold",
        )

    fig.tight_layout()
    png = out_dir / "stage3b_reference_dropout_panel_a_blank_region_composition.png"
    pdf = out_dir / "stage3b_reference_dropout_panel_a_blank_region_composition.pdf"
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
    summary = pd.DataFrame([load_scene(root, row, args.target_quantile, args.neighbor_rings) for _, row in metadata.iterrows()])
    summary_path = out_dir / "stage3b_reference_dropout_panel_a_blank_region_composition_summary.csv"
    summary.to_csv(summary_path, index=False)
    plot_panel(summary, out_dir)
    print(f"[done] {summary_path}")
    print(
        summary[
            [
                "scene",
                "stage3b_blank_spots",
                "inside_target_core_spots",
                "near_target_core_spots",
                "outside_target_associated_spots",
                "fraction_inside_target_core",
                "fraction_near_target_core",
                "fraction_outside_target_associated_region",
            ]
        ].to_string(index=False)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
