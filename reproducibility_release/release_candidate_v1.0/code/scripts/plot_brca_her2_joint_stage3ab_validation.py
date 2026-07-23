#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap, Normalize


SAMPLE = "cytospace_fig2d_tme_brca_her2_ffpe_profile_mask_t_cells_sc_missing_plasma_cells"
GROUP = "biological_application"
OUT_DIR = "visualizations/biological_application_brca_her2_joint_stage3ab"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot a lightweight validation panel for the BRCA HER2 joint Stage3A/Stage3B scenario."
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--sample", default=SAMPLE)
    parser.add_argument("--storage_group", default=GROUP)
    parser.add_argument("--out_dir", default=OUT_DIR)
    return parser.parse_args()


def signature_cmap() -> LinearSegmentedColormap:
    base = colormaps["magma"]
    colors = base(np.linspace(0.18, 1.0, 256))
    return LinearSegmentedColormap.from_list("magma_deep_purple", colors)


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
            raise ValueError(f"Cannot infer coordinates from {path}")
        y_col, x_col = numeric[:2]
    result = coords[[x_col, y_col]].apply(pd.to_numeric, errors="coerce").dropna()
    result.columns = ["x", "y"]
    result["y_plot"] = -result["y"]
    return result


def marker_percentile(expr_path: Path, marker_genes: list[str]) -> pd.Series:
    header = pd.read_csv(expr_path, nrows=0).columns.tolist()
    id_col = header[0]
    present = [gene for gene in marker_genes if gene in header]
    if not present:
        raise ValueError("No requested marker genes are present in ST expression.")
    expr = pd.read_csv(expr_path, index_col=0, usecols=[id_col, *present])
    expr.index = expr.index.astype(str)
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    score = expr.mean(axis=1)
    return score.rank(method="average", pct=True)


def read_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().isin(["true", "1", "yes"])


def load_stage4_summary(root: Path, sample: str, suffix: str) -> dict:
    path = root / "result" / sample / f"stage4_cytospace{suffix}" / "cytospace_output" / "stage4_summary.json"
    return json.loads(path.read_text(encoding="utf-8"))


def draw_spatial_panel(ax, frame: pd.DataFrame) -> None:
    cmap = signature_cmap()
    norm = Normalize(0, 1)
    size = 6.0
    ax.scatter(frame["x"], frame["y_plot"], s=size, c="#E4E0D8", marker="h", linewidths=0, alpha=0.7)
    sc = ax.scatter(
        frame["x"],
        frame["y_plot"],
        s=size,
        c=frame["plasma_marker_percentile"],
        cmap=cmap,
        norm=norm,
        marker="h",
        linewidths=0,
        alpha=0.96,
    )
    target = frame[frame["plasma_top15"]]
    blank = frame[frame["stage3b_blank"]]
    ax.scatter(
        target["x"],
        target["y_plot"],
        s=size * 2.4,
        facecolors="none",
        edgecolors="#00CFE8",
        marker="h",
        linewidths=0.55,
        alpha=1.0,
    )
    ax.scatter(
        blank["x"],
        blank["y_plot"],
        s=size * 1.25,
        facecolors="none",
        edgecolors="#111111",
        marker="h",
        linewidths=0.35,
        alpha=0.88,
    )
    ax.set_title("A. Plasma-cell target region and Stage3B blanking", loc="left", fontweight="bold", fontsize=10)
    ax.text(
        0.01,
        0.02,
        "cyan outline: marker top15 region\nblack outline: Stage3B blank spots",
        transform=ax.transAxes,
        fontsize=8,
        va="bottom",
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal", adjustable="box")
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    cbar.ax.tick_params(labelsize=7)
    cbar.set_label("Plasma marker percentile", fontsize=8)


def draw_stage3a_panel(ax, type_support: pd.DataFrame) -> None:
    plot_types = ["T-cells", "B cells", "CD4 T cells", "CD8 T cells", "Epithelial cells", "Fibroblasts"]
    data = type_support[type_support["orig_type"].isin(plot_types)].copy()
    data["dropped"] = data["Action"].eq("Dropped")
    data = data.set_index("orig_type").loc[[t for t in plot_types if t in set(data["orig_type"])]].reset_index()
    colors = np.where(data["dropped"], "#B2182B", "#4C78A8")
    ax.barh(data["orig_type"], data["support_score"], color=colors, height=0.72)
    ax.axvline(0.7, color="#666666", linestyle="--", linewidth=1.0, alpha=0.75)
    ax.set_xlim(0, 1)
    ax.invert_yaxis()
    ax.set_xlabel("Stage3A support score", fontsize=9)
    ax.set_title("B. Stage3A removes SC-only unsupported T-cells", loc="left", fontweight="bold", fontsize=10)
    for _, row in data.iterrows():
        label = "Dropped" if row["dropped"] else "Kept"
        ax.text(min(float(row["support_score"]) + 0.02, 0.9), row.name, label, va="center", fontsize=8)
    ax.grid(axis="x", color="#DDDDDD", linewidth=0.8)
    ax.set_axisbelow(True)


def draw_mapping_panel(ax, baseline: dict, svtuner: dict, validation: dict) -> None:
    baseline_t = 1.0 if "T-cells" in baseline.get("raw_output_types_by_spot", []) else 0.0
    svtuner_t = 1.0 if "T-cells" in svtuner.get("raw_output_types_by_spot", []) else 0.0
    baseline_blank = 0.0
    svtuner_blank = 1.0 if svtuner.get("stage3b_blank_regions", {}).get("applied_before_mapping") else 0.0
    categories = ["SC-only T-cells\nremoved", "ST-only Plasma region\nblanked before mapping"]
    values = np.array([[1.0 - baseline_t, baseline_blank], [1.0 - svtuner_t, svtuner_blank]])
    x = np.arange(len(categories))
    width = 0.36
    ax.bar(x - width / 2, values[0], width, label="CytoSPACE baseline", color="#A6A6A6")
    ax.bar(x + width / 2, values[1], width, label="SVTuner Stage3A+B", color="#1B9E77")
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=8)
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("Correct handling indicator", fontsize=9)
    ax.set_title("C. Joint mapping outcome", loc="left", fontweight="bold", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(0.0, 1.02))
    for xpos, vals in zip(x, values.T):
        for offset, val in [(-width / 2, vals[0]), (width / 2, vals[1])]:
            if val > 0.5:
                ax.text(xpos + offset, val - 0.08, f"{val:.0f}", ha="center", fontsize=9, fontweight="bold", color="white")
            else:
                ax.text(xpos + offset, val + 0.055, f"{val:.0f}", ha="center", fontsize=9, fontweight="bold")
    ax.text(
        0.02,
        -0.24,
        (
            f"Stage3B Plasma precision={validation['stage3b']['precision']:.2f}, "
            f"recall={validation['stage3b']['recall']:.2f}; "
            f"blank spots={validation['stage3b']['blank_spots']}"
        ),
        transform=ax.transAxes,
        fontsize=8,
        color="#333333",
    )
    ax.grid(axis="y", color="#DDDDDD", linewidth=0.8)
    ax.set_axisbelow(True)


def main() -> None:
    args = parse_args()
    root = Path(args.project_root).resolve()
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    validation_path = out_dir / "joint_stage3ab_validation_summary.json"
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    exported = root / "data" / "processed" / args.storage_group / args.sample / "stage1_preprocess" / "exported"
    stage3b_scores_path = root / "data" / "processed" / args.storage_group / args.sample / "stage3b_st_unsupported" / "spot_unsupported_scores.csv"
    type_support_path = root / "data" / "processed" / args.storage_group / args.sample / "stage3_typematch" / "type_support.csv"

    coords = read_coords(exported / "st_coordinates.csv")
    plasma_score = marker_percentile(exported / "st_expression_normalized.csv", validation["stage3b"]["marker_genes_used"])
    scores = pd.read_csv(stage3b_scores_path, index_col=0)
    scores.index = scores.index.astype(str)
    blank = read_bool(scores["is_unsupported_region"])
    common = coords.index.intersection(plasma_score.index).intersection(blank.index)
    spatial = coords.loc[common].copy()
    spatial["plasma_marker_percentile"] = plasma_score.loc[common]
    spatial["plasma_top15"] = spatial["plasma_marker_percentile"] >= spatial["plasma_marker_percentile"].quantile(0.85)
    spatial["stage3b_blank"] = blank.loc[common].to_numpy()

    type_support = pd.read_csv(type_support_path)
    baseline = load_stage4_summary(root, args.sample, "_bioapp_joint_dropout_baseline")
    svtuner = load_stage4_summary(root, args.sample, "_bioapp_joint_svtuner")

    fig = plt.figure(figsize=(16.2, 4.8), dpi=220)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.0, 1.0], wspace=0.62)
    draw_spatial_panel(fig.add_subplot(gs[0, 0]), spatial)
    draw_stage3a_panel(fig.add_subplot(gs[0, 1]), type_support)
    draw_mapping_panel(fig.add_subplot(gs[0, 2]), baseline, svtuner, validation)
    fig.suptitle("BRCA HER2 FFPE joint Stage3A/Stage3B validation", fontsize=12, fontweight="bold", y=1.04)
    fig.savefig(out_dir / "joint_stage3ab_mechanism_validation.png", bbox_inches="tight")
    fig.savefig(out_dir / "joint_stage3ab_mechanism_validation.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
