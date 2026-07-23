#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd


GROUP_STYLE = {
    "background": ("Background", "#E8E8E8"),
    "blank_outside_target": ("Blank outside marker region", "#8DA6D4"),
    "hit_target": ("Marker region blanked", "#25A486"),
    "missed_target": ("Marker region missed", "#E36A05"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export the Thalamic Stage3B miss-diagnostics map as a compact editable SVG."
    )
    parser.add_argument("--project_root", default=".")
    parser.add_argument(
        "--spots_csv",
        default=(
            "visualizations/cell2location_stage3b_case/thalamic_top15/"
            "cell2location_ST8059051_thalamic_excitatory_"
            "stage3b_miss_diagnostics_spots.csv"
        ),
    )
    parser.add_argument(
        "--output_svg",
        default=(
            "visualizations/cell2location_stage3b_case/thalamic_top15/"
            "cell2location_ST8059051_thalamic_excitatory_"
            "stage3b_miss_diagnostics_spatial.svg"
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    spots_path = (root / args.spots_csv).resolve()
    output_path = (root / args.output_svg).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    spots = pd.read_csv(spots_path)
    required = {"pxl_col", "pxl_row", "group"}
    missing = sorted(required.difference(spots.columns))
    if missing:
        raise ValueError(f"Missing required columns in {spots_path}: {missing}")
    spots["pxl_col"] = pd.to_numeric(spots["pxl_col"], errors="raise")
    spots["pxl_row"] = pd.to_numeric(spots["pxl_row"], errors="raise")

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "svg.fonttype": "none",
            "svg.image_inline": True,
        }
    )
    fig, ax = plt.subplots(figsize=(7.2, 7.15), dpi=260, facecolor="white")

    for group, (_label, color) in GROUP_STYLE.items():
        selected = spots.loc[spots["group"].astype(str).eq(group)]
        if selected.empty:
            continue
        layer = ax.scatter(
            selected["pxl_col"],
            selected["pxl_row"],
            s=29,
            marker="h",
            c=color,
            edgecolors="white",
            linewidths=0.38,
            rasterized=True,
            zorder=2 if group == "background" else 3,
        )
        layer.set_gid(f"spots_{group}")

    ax.set_aspect("equal", adjustable="box")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(
        "Stage3B blanking vs marker-defined target region",
        fontsize=17.5,
        fontweight="bold",
        pad=16,
    )

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="h",
            linestyle="none",
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=0.45,
            markersize=8.5,
            label=label,
        )
        for label, color in GROUP_STYLE.values()
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.045),
        ncol=2,
        frameon=False,
        fontsize=10.5,
        handletextpad=0.65,
        columnspacing=1.35,
    )
    fig.subplots_adjust(left=0.035, right=0.985, top=0.91, bottom=0.14)
    fig.savefig(output_path, format="svg", facecolor="white")
    plt.close(fig)
    print(f"[done] {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
