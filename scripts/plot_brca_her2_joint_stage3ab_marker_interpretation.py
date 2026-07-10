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

SIGNATURES = {
    "Plasma": [],
    "T/NK": ["CD3D", "CD3E", "TRAC", "CD2", "IL7R", "TRBC1", "LTB", "NKG7", "GNLY", "GZMB", "PRF1", "CCL5", "GZMA", "GZMK"],
    "B cell": ["MS4A1", "CD79A", "CD79B", "CD74", "CD37", "HLA-DRA"],
    "Myeloid": ["LYZ", "TYROBP", "C1QA", "C1QB", "C1QC", "FCER1G", "AIF1"],
    "Stromal": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1", "COL6A2"],
    "Epithelial": ["KRT8", "KRT18", "KRT7"],
}

GROUP_TYPES = {
    "T/NK": ["T-cells", "CD4 T cells", "CD8 T cells", "NK cells"],
    "B cell": ["B cells"],
    "Myeloid": ["Monocytes and Macrophages"],
    "Stromal": ["Fibroblasts", "Endothelial cells", "PVL"],
    "Epithelial": ["Epithelial cells"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot marker/pathway interpretation layer for the joint BRCA HER2 Stage3A/Stage3B scene."
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--sample", default=SAMPLE)
    parser.add_argument("--storage_group", default=GROUP)
    parser.add_argument("--out_dir", default=OUT_DIR)
    parser.add_argument("--target_quantile", type=float, default=0.85)
    parser.add_argument("--neighbors", type=int, default=6)
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
    result = coords[[x_col, y_col]].apply(pd.to_numeric, errors="coerce").dropna()
    result.columns = ["x", "y"]
    return result


def read_expression(path: Path) -> pd.DataFrame:
    expr = pd.read_csv(path, index_col=0)
    expr.index = expr.index.astype(str)
    return expr.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def signature_percentiles(expr: pd.DataFrame, signatures: dict[str, list[str]]) -> pd.DataFrame:
    scores: dict[str, pd.Series] = {}
    for name, genes in signatures.items():
        present = [gene for gene in genes if gene in expr.columns]
        if not present:
            scores[name] = pd.Series(np.nan, index=expr.index)
            continue
        raw = expr[present].mean(axis=1)
        scores[name] = raw.rank(method="average", pct=True)
    return pd.DataFrame(scores)


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


def abundance_sum(frame: pd.DataFrame, cell_types: list[str]) -> pd.Series:
    present = [cell_type for cell_type in cell_types if cell_type in frame.columns]
    if not present:
        return pd.Series(0.0, index=frame.index)
    return frame[present].sum(axis=1)


def local_ring(coords: pd.DataFrame, core_mask: pd.Series, neighbors: int) -> pd.Series:
    core_values = core_mask.reindex(coords.index).fillna(False).to_numpy(dtype=bool)
    tree = cKDTree(coords[["x", "y"]].to_numpy())
    k = min(neighbors + 1, len(coords))
    _, indices = tree.query(coords[["x", "y"]].to_numpy(), k=k)
    core_idx = set(np.flatnonzero(core_values).tolist())
    ring = np.zeros(len(coords), dtype=bool)
    for focal, neigh in enumerate(indices[:, 1:]):
        if focal in core_idx:
            continue
        if any(int(n) in core_idx for n in neigh):
            ring[focal] = True
    return pd.Series(ring, index=coords.index)


def draw_signature_heatmap(ax: plt.Axes, region_signature: pd.DataFrame) -> None:
    sns.heatmap(
        region_signature,
        ax=ax,
        cmap="magma",
        vmin=0,
        vmax=1,
        linewidths=0.5,
        linecolor="white",
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "Mean signature percentile"},
    )
    ax.set_title("A. Observed ST marker programs", loc="left", fontweight="bold", fontsize=11)
    ax.set_xlabel("Spatial region", fontsize=9)
    ax.set_ylabel("Marker program", fontsize=9)
    ax.tick_params(axis="x", labelrotation=25, labelsize=8)
    ax.tick_params(axis="y", labelrotation=0, labelsize=8)


def draw_assignment_evidence_panel(ax: plt.Axes, assignment_evidence: pd.DataFrame) -> None:
    melted = assignment_evidence.reset_index(names="program").melt(
        id_vars="program",
        var_name="quantity",
        value_name="value",
    )
    colors = {
        "Observed marker percentile": "#4C78A8",
        "CytoSPACE assigned abundance": "#A6A6A6",
        "SVTuner assigned abundance": "#1B9E77",
    }
    sns.barplot(
        data=melted,
        ax=ax,
        y="program",
        x="value",
        hue="quantity",
        palette=colors,
        orient="h",
    )
    ax.set_xlim(0, 1)
    ax.set_xlabel("Plasma-rich core mean value", fontsize=9)
    ax.set_ylabel("")
    ax.set_title("B. Assignment versus observed marker evidence", loc="left", fontweight="bold", fontsize=11)
    ax.legend(frameon=False, fontsize=7.6, loc="lower center", bbox_to_anchor=(0.5, -0.35), ncol=1)
    ax.grid(axis="x", color="#DDDDDD", linewidth=0.8)
    ax.set_axisbelow(True)


def draw_risk_metrics(ax: plt.Axes, metrics: pd.DataFrame) -> None:
    x = np.arange(len(metrics))
    width = 0.36
    ax.bar(x - width / 2, metrics["CytoSPACE baseline"], width, color="#A6A6A6", label="CytoSPACE baseline")
    ax.bar(x + width / 2, metrics["SVTuner Stage3A+B"], width, color="#1B9E77", label="SVTuner Stage3A+B")
    ax.set_xticks(x)
    ax.set_xticklabels(metrics.index, rotation=23, ha="right", fontsize=8)
    ax.set_ylim(0, max(1.0, float(metrics.max().max()) * 1.25))
    ax.set_ylabel("Mean score in Plasma-rich core", fontsize=9)
    ax.set_title("C. Marker-level interpretation risk", loc="left", fontweight="bold", fontsize=11)
    ax.legend(frameon=False, fontsize=8)
    for xpos, values in enumerate(metrics.to_numpy()):
        for offset, value in [(-width / 2, values[0]), (width / 2, values[1])]:
            ax.text(xpos + offset, value + 0.025, f"{value:.2f}", ha="center", fontsize=8)
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
    expr = read_expression(exported / "st_expression_normalized.csv")
    signatures = dict(SIGNATURES)
    signatures["Plasma"] = validation["stage3b"]["marker_genes_used"]
    sig = signature_percentiles(expr, signatures)

    plasma_core = sig["Plasma"] >= sig["Plasma"].quantile(args.target_quantile)
    common = coords.index.intersection(expr.index).intersection(sig.index)
    coords = coords.loc[common]
    sig = sig.loc[common]
    plasma_core = plasma_core.loc[common]
    ring = local_ring(coords, plasma_core, args.neighbors)
    background = ~(plasma_core | ring)

    baseline = read_fractional(fractional_path(root, args.sample, "_bioapp_joint_dropout_baseline"))
    svtuner = read_fractional(fractional_path(root, args.sample, "_bioapp_joint_svtuner"))
    baseline, svtuner = align_fractionals([baseline, svtuner])
    baseline = baseline.loc[common]
    svtuner = svtuner.loc[common]

    regions = {
        "Background": background,
        "Local ring": ring,
        "Plasma-rich core": plasma_core,
    }
    region_signature = pd.DataFrame(
        {
            region: sig.loc[mask, list(SIGNATURES)].mean()
            for region, mask in regions.items()
        }
    )

    assignment_programs = ["T/NK", "B cell", "Myeloid", "Stromal", "Epithelial"]
    core_mask = plasma_core.to_numpy(dtype=bool)
    assignment_evidence = pd.DataFrame(index=assignment_programs)
    assignment_evidence["Observed marker percentile"] = sig.loc[plasma_core, assignment_programs].mean()
    assignment_evidence["CytoSPACE assigned abundance"] = [
        abundance_sum(baseline, GROUP_TYPES[program]).loc[plasma_core].mean()
        for program in assignment_programs
    ]
    assignment_evidence["SVTuner assigned abundance"] = [
        abundance_sum(svtuner, GROUP_TYPES[program]).loc[plasma_core].mean()
        for program in assignment_programs
    ]

    non_plasma_baseline = 1.0 - baseline.get("Abstained", pd.Series(0.0, index=baseline.index))
    non_plasma_svtuner = 1.0 - svtuner.get("Abstained", pd.Series(0.0, index=svtuner.index))
    tnk_baseline = abundance_sum(baseline, GROUP_TYPES["T/NK"])
    tnk_svtuner = abundance_sum(svtuner, GROUP_TYPES["T/NK"])
    observed_nonplasma = sig[["T/NK", "B cell", "Myeloid", "Stromal", "Epithelial"]].mean(axis=1)
    metrics = pd.DataFrame(
        {
            "CytoSPACE baseline": [
                float((sig["Plasma"] * non_plasma_baseline).loc[plasma_core].mean()),
                float((sig["T/NK"] * tnk_baseline).loc[plasma_core].mean()),
                float((observed_nonplasma * non_plasma_baseline).loc[plasma_core].mean()),
                float(baseline["Abstained"].loc[plasma_core].mean()),
            ],
            "SVTuner Stage3A+B": [
                float((sig["Plasma"] * non_plasma_svtuner).loc[plasma_core].mean()),
                float((sig["T/NK"] * tnk_svtuner).loc[plasma_core].mean()),
                float((observed_nonplasma * non_plasma_svtuner).loc[plasma_core].mean()),
                float(svtuner["Abstained"].loc[plasma_core].mean()),
            ],
        },
        index=[
            "Plasma evidence\nforced into other types",
            "T/NK marker-weighted\nassignment",
            "Supported-marker\nassignment load",
            "Unsupported-core\nabstention",
        ],
    )

    region_signature.to_csv(out_dir / "joint_stage3ab_marker_region_signature_percentiles.csv")
    assignment_evidence.to_csv(out_dir / "joint_stage3ab_marker_assignment_evidence.csv")
    metrics.to_csv(out_dir / "joint_stage3ab_marker_interpretation_metrics.csv")

    sns.set_theme(style="whitegrid", context="paper")
    fig = plt.figure(figsize=(15.4, 4.9), dpi=220)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.14, 1.05], wspace=0.47)
    draw_signature_heatmap(fig.add_subplot(gs[0, 0]), region_signature)
    draw_assignment_evidence_panel(fig.add_subplot(gs[0, 1]), assignment_evidence)
    draw_risk_metrics(fig.add_subplot(gs[0, 2]), metrics)
    fig.suptitle(
        "BRCA HER2 FFPE marker interpretation layer: unsupported Plasma evidence and downstream readout",
        fontsize=12,
        fontweight="bold",
        y=1.04,
    )
    fig.savefig(out_dir / "joint_stage3ab_marker_interpretation_application.png", bbox_inches="tight")
    fig.savefig(out_dir / "joint_stage3ab_marker_interpretation_application.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
