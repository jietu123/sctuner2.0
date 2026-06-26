#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from scipy import sparse

_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from scripts.prepare_cell2location_mouse_brain_stage3b_case import (
    read_10x_h5,
    target_markers_from_counts,
    unique_gene_positions,
)


TARGETS = {
    "thalamic_excitatory": {
        "label": "Thalamic excitatory",
        "drop_types": ["Ext_Thal_1", "Ext_Thal_2"],
        "section": "ST8059051",
    },
    "oligodendrocyte_opc": {
        "label": "Oligodendrocyte/OPC",
        "drop_types": ["Oligo_1", "Oligo_2", "OPC_1", "OPC_2"],
        "section": "ST8059050",
    },
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot marker-defined ST-only target region.")
    p.add_argument("--project_root", default=".")
    p.add_argument("--target", choices=sorted(TARGETS), default="thalamic_excitatory")
    p.add_argument("--section", default=None)
    p.add_argument("--quantile", type=float, default=0.80)
    p.add_argument("--n_markers", type=int, default=30)
    p.add_argument("--point_size", type=float, default=16.0)
    p.add_argument("--out_dir", default="visualizations/cell2location_stage3b_case")
    return p.parse_args()


def log_normalize_rows(values: sparse.spmatrix) -> sparse.csr_matrix:
    values = values.astype(np.float32).tocsr()
    totals = np.asarray(values.sum(axis=1)).ravel()
    totals[totals <= 0] = 1.0
    return values.multiply(10000.0 / totals[:, None]).log1p().tocsr()


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    spec = TARGETS[args.target]
    section = args.section or spec["section"]
    label = spec["label"]
    drop_types = spec["drop_types"]

    source = root / "data" / "raw" / "cell2location_mouse_brain"
    sc_data = ad.read_h5ad(source / "sc.h5ad")
    if sc_data.raw is None:
        raise ValueError("sc.h5ad lacks raw counts")
    labels = sc_data.obs["annotation_1"].astype(str)
    raw_counts = sc_data.raw.X
    if not sparse.issparse(raw_counts):
        raw_counts = sparse.csr_matrix(raw_counts)

    section_dir = source / "mouse_brain_visium_wo_cloupe_data" / "rawdata" / section
    counts, barcodes, st_genes = read_10x_h5(section_dir / "filtered_feature_bc_matrix.h5")
    positions = pd.read_csv(
        section_dir / "spatial" / "tissue_positions_list.csv",
        header=None,
        names=["spot_id", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"],
    )
    positions["spot_id"] = positions["spot_id"].astype(str)
    barcode_pos = {barcode: i for i, barcode in enumerate(barcodes.astype(str))}
    tissue_spots = [
        spot
        for spot in positions.loc[positions["in_tissue"] == 1, "spot_id"].tolist()
        if spot in barcode_pos
    ]
    row_idx = np.array([barcode_pos[spot] for spot in tissue_spots], dtype=int)
    counts = counts[row_idx, :].tocsr()

    st_pos = unique_gene_positions(st_genes)
    markers = target_markers_from_counts(raw_counts, sc_data.raw.var, labels, drop_types, st_pos)
    markers = [gene for gene in markers[: args.n_markers] if gene in st_pos]
    if len(markers) < 8:
        raise ValueError(f"Too few marker genes for {label}: {markers}")

    norm = log_normalize_rows(counts)
    marker_idx = np.array([st_pos[gene] for gene in markers], dtype=int)
    marker_score = np.asarray(norm[:, marker_idx].mean(axis=1)).ravel()
    threshold = float(np.quantile(marker_score, args.quantile))
    target_mask = marker_score >= threshold

    plot_df = positions.set_index("spot_id").loc[tissue_spots].copy()
    plot_df["marker_score"] = marker_score
    plot_df["target_region"] = target_mask

    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"cell2location_{section}_{args.target}_marker_region"

    fig, ax = plt.subplots(figsize=(6.8, 5.6), dpi=300)
    marker_cmap = LinearSegmentedColormap.from_list(
        "target_marker_purple_yellow",
        ["#2b1b6d", "#6f3f97", "#b94c9a", "#f17570", "#f7b267", "#fff2b2"],
    )
    sc_all = ax.scatter(
        plot_df["pxl_col"],
        -plot_df["pxl_row"],
        s=args.point_size,
        c=plot_df["marker_score"],
        cmap=marker_cmap,
        marker="h",
        linewidths=0,
        alpha=0.95,
    )
    ax.set_title(
        f"{label} marker score",
        fontsize=10.5,
        weight="bold",
        pad=8,
    )
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(sc_all, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("target marker score", fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    fig.savefig(out_dir / f"{prefix}.png", bbox_inches="tight")
    fig.savefig(out_dir / f"{prefix}.pdf", bbox_inches="tight")
    plt.close(fig)

    plot_df[["array_row", "array_col", "pxl_row", "pxl_col", "marker_score", "target_region"]].to_csv(
        out_dir / f"{prefix}.csv",
        index_label="spot_id",
    )
    summary = {
        "section": section,
        "target_label": label,
        "drop_types_for_later_reference_dropout": drop_types,
        "target_sc_cells": int(labels.isin(drop_types).sum()),
        "markers_used": markers,
        "region_quantile": args.quantile,
        "st_spots": int(len(plot_df)),
        "target_region_spots": int(target_mask.sum()),
        "target_region_fraction": float(target_mask.mean()),
        "marker_score_threshold": threshold,
        "png": str(out_dir / f"{prefix}.png"),
        "pdf": str(out_dir / f"{prefix}.pdf"),
    }
    (out_dir / f"{prefix}_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
