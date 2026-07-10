#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial import cKDTree


SAMPLE = "cytospace_fig2d_tme_brca_her2_ffpe_profile_mask_t_cells_sc_missing_plasma_cells"
GROUP = "biological_application"
OUT_DIR = "visualizations/biological_application_brca_her2_joint_stage3ab"
LR_RESOURCE = "data/raw/spatial_communication_reference_resources/liana_1.7.3_omni_resource.csv"

IMMUNE_TYPES = [
    "T-cells",
    "CD4 T cells",
    "CD8 T cells",
    "B cells",
    "NK cells",
    "Monocytes and Macrophages",
]
T_CELL_TYPES = ["T-cells", "CD4 T cells", "CD8 T cells"]
GROUPS = {
    "T/NK": ["T-cells", "CD4 T cells", "CD8 T cells", "NK cells"],
    "B": ["B cells"],
    "Myeloid": ["Monocytes and Macrophages"],
    "Tumor": ["Epithelial cells"],
    "Stromal": ["Fibroblasts", "Endothelial cells", "PVL"],
    "Abstained": ["Abstained"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot putative LR communication changes for the joint BRCA HER2 Stage3A/Stage3B scene."
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--sample", default=SAMPLE)
    parser.add_argument("--storage_group", default=GROUP)
    parser.add_argument("--out_dir", default=OUT_DIR)
    parser.add_argument("--target_quantile", type=float, default=0.85)
    parser.add_argument("--neighbors", type=int, default=6)
    parser.add_argument("--top_lr", type=int, default=80)
    parser.add_argument("--show_lr", type=int, default=18)
    return parser.parse_args()


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
            raise ValueError(f"Cannot infer coordinate columns from {path}")
        y_col, x_col = numeric[:2]
    out = coords[[x_col, y_col]].apply(pd.to_numeric, errors="coerce").dropna()
    out.columns = ["x", "y"]
    return out


def read_expression(path: Path, genes: list[str] | None = None) -> pd.DataFrame:
    if genes is None:
        expr = pd.read_csv(path, index_col=0)
    else:
        header = pd.read_csv(path, nrows=0).columns.tolist()
        id_col = header[0]
        present = [gene for gene in genes if gene in header]
        expr = pd.read_csv(path, index_col=0, usecols=[id_col, *present])
    expr.index = expr.index.astype(str)
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return expr


def marker_percentile(expr: pd.DataFrame, marker_genes: list[str]) -> pd.Series:
    present = [gene for gene in marker_genes if gene in expr.columns]
    if not present:
        raise ValueError("No Plasma marker genes are present in ST expression.")
    score = expr[present].mean(axis=1)
    return score.rank(method="average", pct=True)


def fractional_path(root: Path, sample: str, suffix: str) -> Path:
    return root / "result" / sample / f"stage4_cytospace{suffix}" / "cytospace_output" / "fractional_abundances_by_spot.csv"


def read_fractional(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, index_col=0)
    frame.index = frame.index.astype(str)
    frame = frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    frame = frame.drop(columns=["Unknown_sc_only"], errors="ignore")
    totals = frame.sum(axis=1)
    nonzero = totals > 0
    frame.loc[nonzero] = frame.loc[nonzero].div(totals.loc[nonzero], axis=0)
    frame["Abstained"] = np.where(totals <= 1e-12, 1.0, 0.0)
    return frame


def align_fractionals(frames: list[pd.DataFrame]) -> list[pd.DataFrame]:
    columns = sorted(set().union(*(set(frame.columns) for frame in frames)))
    return [frame.reindex(columns=columns, fill_value=0.0) for frame in frames]


def spatial_edges(coords: pd.DataFrame, neighbors: int) -> tuple[np.ndarray, np.ndarray]:
    xy = coords[["x", "y"]].to_numpy()
    tree = cKDTree(xy)
    k = min(neighbors + 1, len(coords))
    _, indices = tree.query(xy, k=k)
    if indices.ndim == 1:
        indices = indices[:, None]
    source = np.repeat(np.arange(len(coords)), indices.shape[1] - 1)
    destination = indices[:, 1:].reshape(-1)
    valid = source != destination
    return source[valid], destination[valid]


def load_lr_pairs(path: Path, present_genes: set[str]) -> pd.DataFrame:
    lr = pd.read_csv(path)
    lr = lr.loc[lr["resource"].eq("consensus"), ["source_genesymbol", "target_genesymbol"]].drop_duplicates()
    lr.columns = ["ligand", "receptor"]
    simple = ~lr["ligand"].astype(str).str.contains("_", regex=False) & ~lr["receptor"].astype(str).str.contains("_", regex=False)
    lr = lr.loc[simple]
    keep = lr["ligand"].isin(present_genes) & lr["receptor"].isin(present_genes)
    return lr.loc[keep].reset_index(drop=True)


def minmax_expression(expr: pd.DataFrame) -> pd.DataFrame:
    arr = expr.to_numpy(dtype=float)
    lo = np.nanmin(arr, axis=0)
    hi = np.nanmax(arr, axis=0)
    denom = np.where((hi - lo) > 1e-12, hi - lo, 1.0)
    scaled = (arr - lo) / denom
    return pd.DataFrame(scaled, index=expr.index, columns=expr.columns)


def abundance_sum(frame: pd.DataFrame, cell_types: list[str]) -> np.ndarray:
    present = [cell_type for cell_type in cell_types if cell_type in frame.columns]
    if not present:
        return np.zeros(len(frame), dtype=float)
    return frame[present].sum(axis=1).to_numpy(dtype=float)


def lr_edge_weight(expr: pd.DataFrame, lr: pd.DataFrame, source: np.ndarray, destination: np.ndarray) -> np.ndarray:
    weights = []
    for row in lr.itertuples(index=False):
        weights.append(expr[row.ligand].to_numpy()[source] * expr[row.receptor].to_numpy()[destination])
    if not weights:
        return np.zeros(len(source), dtype=float)
    return np.vstack(weights).mean(axis=0)


def lr_pair_scores(
    expr: pd.DataFrame,
    lr: pd.DataFrame,
    source: np.ndarray,
    destination: np.ndarray,
    sender_mass: np.ndarray,
    receiver_mass: np.ndarray,
) -> pd.Series:
    scores: dict[str, float] = {}
    mass = sender_mass[source] * receiver_mass[destination]
    for row in lr.itertuples(index=False):
        product = expr[row.ligand].to_numpy()[source] * expr[row.receptor].to_numpy()[destination]
        scores[f"{row.ligand}_{row.receptor}"] = float(np.mean(product * mass))
    return pd.Series(scores)


def group_matrix(frame: pd.DataFrame, edge_weight: np.ndarray, source: np.ndarray, destination: np.ndarray) -> pd.DataFrame:
    labels = list(GROUPS)
    matrix = np.zeros((len(labels), len(labels)), dtype=float)
    group_values = {label: abundance_sum(frame, types) for label, types in GROUPS.items()}
    for i, sender in enumerate(labels):
        sender_values = group_values[sender][source]
        for j, receiver in enumerate(labels):
            receiver_values = group_values[receiver][destination]
            matrix[i, j] = float(np.mean(edge_weight * sender_values * receiver_values))
    return pd.DataFrame(matrix, index=labels, columns=labels)


def draw_lr_dotplot(ax: plt.Axes, dot_data: pd.DataFrame) -> None:
    methods = ["CytoSPACE baseline", "SVTuner Stage3A+B"]
    pairs = dot_data["lr_pair"].drop_duplicates().tolist()
    pair_to_y = {pair: i for i, pair in enumerate(pairs)}
    method_to_x = {method: i for i, method in enumerate(methods)}
    max_score = max(float(dot_data["score"].max()), 1e-12)
    for row in dot_data.itertuples(index=False):
        ax.scatter(
            method_to_x[row.method],
            pair_to_y[row.lr_pair],
            s=35 + 420 * row.score / max_score,
            c=row.score,
            cmap="Reds",
            vmin=0,
            vmax=max_score,
            edgecolors="#555555",
            linewidths=0.25,
        )
    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels(methods, fontsize=8)
    ax.set_yticks(range(len(pairs)))
    ax.set_yticklabels(pairs, fontsize=7)
    ax.invert_yaxis()
    ax.set_title("A. Core-involved putative LR signal", loc="left", fontweight="bold", fontsize=11)
    ax.set_xlabel("")
    ax.set_ylabel("Ligand_receptor pair", fontsize=9)
    ax.grid(axis="x", color="#E0E0E0", linewidth=0.8)


def draw_group_heatmaps(ax1: plt.Axes, ax2: plt.Axes, baseline_matrix: pd.DataFrame, svtuner_matrix: pd.DataFrame) -> None:
    vmax = max(float(baseline_matrix.max().max()), float(svtuner_matrix.max().max()), 1e-12)
    for ax, matrix, title in [
        (ax1, baseline_matrix, "CytoSPACE baseline"),
        (ax2, svtuner_matrix, "SVTuner Stage3A+B"),
    ]:
        sns.heatmap(
            matrix,
            ax=ax,
            cmap="Reds",
            vmin=0,
            vmax=vmax,
            cbar=ax is ax2,
            cbar_kws={"label": "LR-weighted coupling"},
            linewidths=0.4,
            linecolor="white",
            square=True,
        )
        panel_title = f"B. LR-weighted group coupling\n{title}" if ax is ax1 else title
        ax.set_title(panel_title, fontsize=10, fontweight="bold")
        ax.set_xlabel("Receiver group", fontsize=8)
        ax.set_ylabel("Sender group" if ax is ax1 else "", fontsize=8)
        ax.tick_params(axis="x", labelrotation=45, labelsize=8)
        ax.tick_params(axis="y", labelrotation=0, labelsize=8)


def draw_metric_panel(ax: plt.Axes, metrics: pd.DataFrame) -> None:
    x = np.arange(len(metrics))
    width = 0.36
    ax.bar(x - width / 2, metrics["CytoSPACE baseline"], width, color="#A6A6A6", label="CytoSPACE baseline")
    ax.bar(x + width / 2, metrics["SVTuner Stage3A+B"], width, color="#1B9E77", label="SVTuner Stage3A+B")
    ax.set_xticks(x)
    ax.set_xticklabels(metrics.index, rotation=22, ha="right", fontsize=8)
    ax.set_ylim(0, max(0.01, float(metrics.max().max()) * 1.25))
    ax.set_ylabel("Mean LR-weighted score", fontsize=9)
    ax.set_title("C. Communication interpretation after correction", loc="left", fontweight="bold", fontsize=11)
    ax.legend(frameon=False, fontsize=8)
    for xpos, row in enumerate(metrics.to_numpy()):
        for offset, value in [(-width / 2, row[0]), (width / 2, row[1])]:
            ax.text(xpos + offset, value + float(metrics.max().max()) * 0.035, f"{value:.3f}", ha="center", fontsize=7.5)
    ax.grid(axis="y", color="#DDDDDD", linewidth=0.8)
    ax.set_axisbelow(True)


def main() -> None:
    args = parse_args()
    root = Path(args.project_root).resolve()
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    validation = json.loads((out_dir / "joint_stage3ab_validation_summary.json").read_text(encoding="utf-8"))

    exported = root / "data" / "processed" / args.storage_group / args.sample / "stage1_preprocess" / "exported"
    coords = read_coords(exported / "st_coordinates.csv")
    expr_all = read_expression(exported / "st_expression_normalized.csv")
    marker = marker_percentile(expr_all, validation["stage3b"]["marker_genes_used"])
    target_mask = marker >= marker.quantile(args.target_quantile)

    baseline = read_fractional(fractional_path(root, args.sample, "_bioapp_joint_dropout_baseline"))
    svtuner = read_fractional(fractional_path(root, args.sample, "_bioapp_joint_svtuner"))
    baseline, svtuner = align_fractionals([baseline, svtuner])

    common = coords.index.intersection(expr_all.index).intersection(baseline.index).intersection(svtuner.index)
    coords = coords.loc[common]
    expr_all = expr_all.loc[common]
    baseline = baseline.loc[common]
    svtuner = svtuner.loc[common]
    target_array = target_mask.loc[common].to_numpy(dtype=bool)

    source, destination = spatial_edges(coords, args.neighbors)
    core_edge = target_array[source] | target_array[destination]
    source = source[core_edge]
    destination = destination[core_edge]

    lr = load_lr_pairs(root / LR_RESOURCE, set(expr_all.columns))
    lr_genes = sorted(set(lr["ligand"]).union(set(lr["receptor"])))
    expr = minmax_expression(expr_all[lr_genes])

    edge_expression = []
    for row in lr.itertuples(index=False):
        score = expr[row.ligand].to_numpy()[source] * expr[row.receptor].to_numpy()[destination]
        edge_expression.append(float(score.mean()))
    lr = lr.assign(edge_expression=edge_expression).sort_values("edge_expression", ascending=False).head(args.top_lr).reset_index(drop=True)

    baseline_immune = abundance_sum(baseline, IMMUNE_TYPES)
    svtuner_immune = abundance_sum(svtuner, IMMUNE_TYPES)
    baseline_scores = lr_pair_scores(expr, lr, source, destination, baseline_immune, baseline_immune)
    svtuner_scores = lr_pair_scores(expr, lr, source, destination, svtuner_immune, svtuner_immune)
    score_table = pd.DataFrame(
        {
            "lr_pair": baseline_scores.index,
            "cytospace_baseline": baseline_scores.values,
            "svtuner_stage3ab": svtuner_scores.reindex(baseline_scores.index).to_numpy(),
        }
    )
    score_table["drop_after_svtuner"] = score_table["cytospace_baseline"] - score_table["svtuner_stage3ab"]
    score_table = score_table.sort_values(["drop_after_svtuner", "cytospace_baseline"], ascending=False)
    shown = score_table.head(args.show_lr)
    dot_data = pd.concat(
        [
            shown[["lr_pair", "cytospace_baseline"]].rename(columns={"cytospace_baseline": "score"}).assign(method="CytoSPACE baseline"),
            shown[["lr_pair", "svtuner_stage3ab"]].rename(columns={"svtuner_stage3ab": "score"}).assign(method="SVTuner Stage3A+B"),
        ],
        ignore_index=True,
    )
    dot_data["lr_pair"] = pd.Categorical(dot_data["lr_pair"], categories=shown["lr_pair"].tolist(), ordered=True)
    dot_data = dot_data.sort_values(["lr_pair", "method"])

    edge_weight = lr_edge_weight(expr, lr, source, destination)
    baseline_matrix = group_matrix(baseline, edge_weight, source, destination)
    svtuner_matrix = group_matrix(svtuner, edge_weight, source, destination)

    baseline_t = abundance_sum(baseline, T_CELL_TYPES)
    svtuner_t = abundance_sum(svtuner, T_CELL_TYPES)
    metrics = pd.DataFrame(
        {
            "CytoSPACE baseline": [
                float(np.mean(edge_weight * baseline_immune[source] * baseline_immune[destination])),
                float(np.mean(edge_weight * ((baseline_t[source] * baseline_immune[destination]) + (baseline_immune[source] * baseline_t[destination])))),
                float(np.mean(edge_weight * abundance_sum(baseline, ["Abstained"])[source])),
            ],
            "SVTuner Stage3A+B": [
                float(np.mean(edge_weight * svtuner_immune[source] * svtuner_immune[destination])),
                float(np.mean(edge_weight * ((svtuner_t[source] * svtuner_immune[destination]) + (svtuner_immune[source] * svtuner_t[destination])))),
                float(np.mean(edge_weight * abundance_sum(svtuner, ["Abstained"])[source])),
            ],
        },
        index=[
            "Immune-immune LR\npotential",
            "T-cell-involved LR\npotential",
            "Abstained core-edge\nprotection",
        ],
    )

    score_table.to_csv(out_dir / "joint_stage3ab_communication_lr_pair_scores.csv", index=False)
    dot_data.to_csv(out_dir / "joint_stage3ab_communication_dotplot_data.csv", index=False)
    baseline_matrix.to_csv(out_dir / "joint_stage3ab_communication_group_matrix_cytospace.csv")
    svtuner_matrix.to_csv(out_dir / "joint_stage3ab_communication_group_matrix_svtuner.csv")
    metrics.to_csv(out_dir / "joint_stage3ab_communication_summary_metrics.csv")

    sns.set_theme(style="whitegrid", context="paper")
    fig = plt.figure(figsize=(15.8, 6.1), dpi=220)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.05, 1.05, 0.95], height_ratios=[1.0, 1.0], wspace=0.48, hspace=0.55)
    draw_lr_dotplot(fig.add_subplot(gs[:, 0]), dot_data)
    mid = gs[:, 1].subgridspec(1, 2, wspace=0.18)
    draw_group_heatmaps(fig.add_subplot(mid[0, 0]), fig.add_subplot(mid[0, 1]), baseline_matrix, svtuner_matrix)
    draw_metric_panel(fig.add_subplot(gs[:, 2]), metrics)
    fig.suptitle(
        "BRCA HER2 FFPE communication layer: reference mismatch creates putative immune LR signals",
        fontsize=12,
        fontweight="bold",
        y=1.02,
    )
    fig.savefig(out_dir / "joint_stage3ab_communication_application.png", bbox_inches="tight")
    fig.savefig(out_dir / "joint_stage3ab_communication_application.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
