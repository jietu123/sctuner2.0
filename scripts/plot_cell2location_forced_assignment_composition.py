#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
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


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Plot baseline forced-assignment composition in the ST-only target region."
    )
    p.add_argument("--project_root", default=".")
    p.add_argument("--section", default="ST8059051")
    p.add_argument("--target", default="thalamic_excitatory")
    p.add_argument("--target_label", default="Thalamic excitatory")
    p.add_argument("--top_types", type=int, default=8)
    p.add_argument("--out_dir", default="visualizations/cell2location_stage3b_case")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    out_dir = root / args.out_dir
    spatial_prefix = f"cell2location_{args.section}_{args.target}_baseline_forced_assignment_region"
    comp_path = out_dir / f"{spatial_prefix}_composition.csv"
    if not comp_path.exists():
        raise FileNotFoundError(f"Forced-assignment composition CSV not found: {comp_path}")

    composition = pd.read_csv(comp_path)
    if not {"assigned_type", "spots", "fraction"}.issubset(composition.columns):
        raise ValueError(f"Unexpected composition columns: {composition.columns.tolist()}")
    composition = composition.sort_values("fraction", ascending=False).reset_index(drop=True)
    top = composition.head(args.top_types).copy()
    other_spots = int(composition.iloc[args.top_types:]["spots"].sum())
    other_fraction = float(composition.iloc[args.top_types:]["fraction"].sum())
    if other_spots > 0:
        top = pd.concat(
            [
                top,
                pd.DataFrame(
                    [{"assigned_type": "Other", "spots": other_spots, "fraction": other_fraction}]
                ),
            ],
            ignore_index=True,
        )

    colors = [CATEGORY_PALETTE.get(x, CATEGORY_PALETTE["Other"]) for x in top["assigned_type"]]

    fig, ax = plt.subplots(figsize=(5.6, 4.2), dpi=300)
    bars = ax.barh(
        top["assigned_type"],
        top["fraction"],
        color=colors,
        edgecolor="white",
        linewidth=0.8,
    )
    ax.invert_yaxis()
    ax.set_xlim(0, max(0.5, float(top["fraction"].max()) * 1.18))
    ax.set_xlabel("Fraction of target-region spots", fontsize=9)
    ax.set_ylabel("")
    ax.set_title(
        f"Baseline forced-assignment composition\nin {args.target_label}-enriched region",
        fontsize=10.5,
        weight="bold",
        pad=8,
    )
    ax.tick_params(axis="both", labelsize=8)
    ax.grid(axis="x", color="#d9d9d9", linewidth=0.7)
    ax.set_axisbelow(True)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    total_spots = int(composition["spots"].sum())
    for bar, (_, row) in zip(bars, top.iterrows()):
        ax.text(
            bar.get_width() + 0.006,
            bar.get_y() + bar.get_height() / 2,
            f"{row['fraction'] * 100:.1f}%",
            va="center",
            ha="left",
            fontsize=7.5,
        )

    prefix = f"cell2location_{args.section}_{args.target}_forced_assignment_composition"
    png = out_dir / f"{prefix}.png"
    pdf = out_dir / f"{prefix}.pdf"
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    top.to_csv(out_dir / f"{prefix}.csv", index=False)
    summary = {
        "section": args.section,
        "target_label": args.target_label,
        "assignment_level": "broad_category",
        "target_region_spots": total_spots,
        "top_types": top.to_dict(orient="records"),
        "png": str(png),
        "pdf": str(pdf),
    }
    (out_dir / f"{prefix}_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
