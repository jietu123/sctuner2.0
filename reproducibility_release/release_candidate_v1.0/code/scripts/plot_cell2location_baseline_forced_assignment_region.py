#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


CATEGORY_ORDER = [
    "Astrocyte",
    "Excitatory neuron",
    "Inhibitory neuron",
    "Neuroblast",
    "Microglia",
    "Endothelial",
    "Oligodendrocyte/OPC",
    "Other",
]

CATEGORY_PALETTE = {
    "Astrocyte": "#1f77b4",
    "Excitatory neuron": "#ff7f0e",
    "Inhibitory neuron": "#2ca02c",
    "Neuroblast": "#d62728",
    "Microglia": "#9467bd",
    "Endothelial": "#8c564b",
    "Oligodendrocyte/OPC": "#e377c2",
    "Other": "#17becf",
}


def broad_cell_type(cell_type: str) -> str:
    cell_type = str(cell_type)
    if cell_type.startswith("Astro"):
        return "Astrocyte"
    if cell_type.startswith("Ext_"):
        return "Excitatory neuron"
    if cell_type.startswith("Inh_"):
        return "Inhibitory neuron"
    if cell_type.startswith("Nb_"):
        return "Neuroblast"
    if cell_type.startswith("Micro"):
        return "Microglia"
    if cell_type.startswith("Endo"):
        return "Endothelial"
    if cell_type.startswith("Oligo") or cell_type.startswith("OPC"):
        return "Oligodendrocyte/OPC"
    return "Other"


def aggregate_to_broad_categories(frac: pd.DataFrame) -> pd.DataFrame:
    broad_columns: dict[str, list[str]] = {}
    for column in frac.columns.astype(str):
        broad_columns.setdefault(broad_cell_type(column), []).append(column)
    broad = pd.DataFrame(index=frac.index)
    for category, columns in broad_columns.items():
        broad[category] = frac[columns].sum(axis=1)
    ordered = [category for category in CATEGORY_ORDER if category in broad.columns]
    return broad.loc[:, ordered]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Plot baseline forced assignments inside a marker-defined ST-only target region."
    )
    p.add_argument("--project_root", default=".")
    p.add_argument(
        "--sample",
        default="cell2loc_scan_st8059051_sc_missing_thalamic_excitatory",
    )
    p.add_argument("--section", default="ST8059051")
    p.add_argument("--target", default="thalamic_excitatory")
    p.add_argument("--target_label", default="Thalamic excitatory")
    p.add_argument("--point_size", type=float, default=16.0)
    p.add_argument("--top_types", type=int, default=8)
    p.add_argument("--out_dir", default="visualizations/cell2location_stage3b_case")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    out_dir = root / args.out_dir
    target_prefix = f"cell2location_{args.section}_{args.target}_marker_region"
    target_path = out_dir / f"{target_prefix}.csv"
    if not target_path.exists():
        raise FileNotFoundError(f"Target marker-region CSV not found: {target_path}")

    mapping_path = (
        root
        / "result"
        / args.sample
        / "stage4_cytospace_baseline"
        / "cytospace_output"
        / "fractional_abundances_by_spot.csv"
    )
    if not mapping_path.exists():
        raise FileNotFoundError(f"Baseline fractional abundance file not found: {mapping_path}")

    target = pd.read_csv(target_path, index_col=0)
    target.index = target.index.astype(str)
    target["target_region"] = target["target_region"].astype(bool)
    frac = pd.read_csv(mapping_path, index_col=0)
    frac.index = frac.index.astype(str)
    frac = frac.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    common = target.index.intersection(frac.index)
    if common.empty:
        raise ValueError("No overlapping spots between target region and baseline mapping")
    target = target.loc[common].copy()
    frac = frac.loc[common].copy()
    region_ids = target.index[target["target_region"]]
    if region_ids.empty:
        raise ValueError("Target marker-defined region is empty")

    broad_frac = aggregate_to_broad_categories(frac)
    dominant = broad_frac.idxmax(axis=1)
    dominant_fraction = broad_frac.max(axis=1)
    target["dominant_type"] = dominant
    target["dominant_fraction"] = dominant_fraction

    region_types = dominant.loc[region_ids].value_counts()
    top_types = [category for category in CATEGORY_ORDER if category in region_types.index]
    top_types = top_types[: args.top_types]
    target["plot_type"] = np.where(
        target["target_region"],
        target["dominant_type"].where(target["dominant_type"].isin(top_types), "Other"),
        "Outside target region",
    )

    palette = {category: CATEGORY_PALETTE[category] for category in top_types}
    palette["Other"] = CATEGORY_PALETTE["Other"]

    fig, ax = plt.subplots(figsize=(6.8, 5.6), dpi=300)
    background = target[~target["target_region"]]
    ax.scatter(
        background["pxl_col"],
        -background["pxl_row"],
        s=args.point_size,
        c="#d9d9d9",
        marker="h",
        linewidths=0,
        alpha=0.72,
    )
    region = target[target["target_region"]]
    for cell_type in top_types + (["Other"] if (region["plot_type"] == "Other").any() else []):
        sub = region[region["plot_type"] == cell_type]
        if sub.empty:
            continue
        ax.scatter(
            sub["pxl_col"],
            -sub["pxl_row"],
            s=args.point_size,
            c=[palette[cell_type]],
            marker="h",
            linewidths=0,
            alpha=0.96,
            label=cell_type,
        )
    ax.set_title(
        f"Baseline assignments in {args.target_label}-enriched region",
        fontsize=10.5,
        weight="bold",
        pad=8,
    )
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.legend(
        bbox_to_anchor=(1.02, 0.5),
        loc="center left",
        frameon=False,
        fontsize=7,
        markerscale=1.0,
        handletextpad=0.2,
        borderaxespad=0,
    )
    prefix = f"cell2location_{args.section}_{args.target}_baseline_forced_assignment_region"
    png = out_dir / f"{prefix}.png"
    pdf = out_dir / f"{prefix}.pdf"
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    composition = (
        dominant.loc[region_ids]
        .value_counts()
        .rename_axis("assigned_type")
        .reset_index(name="spots")
    )
    composition["sort_order"] = composition["assigned_type"].map(
        {category: i for i, category in enumerate(CATEGORY_ORDER)}
    )
    composition = composition.sort_values(["sort_order", "assigned_type"]).drop(columns="sort_order")
    composition["fraction"] = composition["spots"] / composition["spots"].sum()
    composition.to_csv(out_dir / f"{prefix}_composition.csv", index=False)
    target[["pxl_row", "pxl_col", "target_region", "dominant_type", "dominant_fraction", "plot_type"]].to_csv(
        out_dir / f"{prefix}.csv",
        index_label="spot_id",
    )
    summary = {
        "sample": args.sample,
        "section": args.section,
        "target_label": args.target_label,
        "assignment_level": "broad_category",
        "fine_reference_types": int(frac.shape[1]),
        "broad_reference_categories": broad_frac.columns.tolist(),
        "target_region_spots": int(len(region_ids)),
        "assigned_categories": composition.to_dict(orient="records"),
        "png": str(png),
        "pdf": str(pdf),
    }
    (out_dir / f"{prefix}_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
