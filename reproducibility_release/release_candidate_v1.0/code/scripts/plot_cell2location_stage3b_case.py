#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot the cell2location mouse brain Stage3B case study.")
    p.add_argument("--project_root", default=".")
    p.add_argument(
        "--sample",
        default="cell2loc_mouse_brain_st8059048_sc_missing_oligodendrocyte_opc",
    )
    p.add_argument("--target_label", default="Oligodendrocyte/OPC")
    p.add_argument("--top_types", type=int, default=8)
    p.add_argument("--out_dir", default="visualizations/cell2location_stage3b_case")
    return p.parse_args()


def read_fractional(root: Path, sample: str, suffix: str) -> pd.DataFrame:
    rel = Path(f"stage4_cytospace{suffix}") / "cytospace_output" / "fractional_abundances_by_spot.csv"
    candidates = [
        root / "result" / sample / rel,
        root / "result" / "cell2location_mouse_brain" / sample / rel,
    ]
    path = next((candidate for candidate in candidates if candidate.exists()), candidates[0])
    frame = pd.read_csv(path, index_col=0)
    frame.index = frame.index.astype(str)
    return frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    sample = args.sample
    processed = root / "data" / "processed" / "cell2location_mouse_brain" / sample
    stage1 = processed / "stage1_preprocess"
    export = stage1 / "exported"
    stage3b = processed / "stage3b_st_unsupported"
    result = root / "result" / "cell2location_mouse_brain" / sample

    coords = pd.read_csv(export / "st_coordinates.csv", index_col=0)
    coords.index = coords.index.astype(str)
    st_expr = pd.read_csv(export / "st_expression_normalized.csv", index_col=0)
    st_expr.index = st_expr.index.astype(str)
    markers = [
        x.strip()
        for x in (stage1 / "target_marker_genes.txt").read_text(encoding="utf-8").splitlines()
        if x.strip()
    ]
    markers = [gene for gene in markers[:30] if gene in st_expr.columns]
    if not markers:
        raise ValueError("No target marker genes found in ST expression")
    marker_score = st_expr.loc[:, markers].mean(axis=1).rename("target_marker_score")

    scores = pd.read_csv(stage3b / "spot_unsupported_scores.csv", index_col=0)
    scores.index = scores.index.astype(str)
    scores["stage3b_score"] = -np.log10(scores["spot_pvalue"].clip(lower=np.finfo(float).tiny))
    scores["is_unsupported_region"] = scores["is_unsupported_region"].astype(bool)
    blank_ids = scores.index[scores["is_unsupported_region"]]

    baseline = read_fractional(root, sample, "_baseline")
    svtuner = read_fractional(root, sample, "_stage3b_blank")
    common_blank = blank_ids.intersection(baseline.index).intersection(svtuner.index)
    if common_blank.empty:
        raise ValueError("Stage3B blank spot set is empty or absent from mapping outputs")

    baseline_blank = baseline.loc[common_blank]
    svtuner_blank = svtuner.loc[common_blank]
    baseline_burden = baseline_blank.sum(axis=1)
    svtuner_burden = svtuner_blank.sum(axis=1)

    type_mass = baseline_blank.sum(axis=0).sort_values(ascending=False)
    top = type_mass.head(args.top_types)
    other = type_mass.iloc[args.top_types :].sum()
    if other > 0:
        top.loc["Other"] = other
    top_fraction = top / top.sum() if top.sum() > 0 else top

    plot_df = coords.join(marker_score, how="left").join(scores[["stage3b_score", "is_unsupported_region"]], how="left")
    plot_df["is_unsupported_region"] = plot_df["is_unsupported_region"].fillna(False).astype(bool)

    sns.set_theme(style="whitegrid", context="paper")
    fig = plt.figure(figsize=(14.2, 4.6), dpi=300)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.28, 1.02, 0.9], wspace=0.55)

    ax0 = fig.add_subplot(gs[0, 0])
    base_scatter = ax0.scatter(
        plot_df["pxl_col"],
        -plot_df["pxl_row"],
        c=plot_df["target_marker_score"],
        s=7,
        cmap="magma",
        linewidths=0,
        alpha=0.9,
    )
    blank_df = plot_df[plot_df["is_unsupported_region"]]
    ax0.scatter(
        blank_df["pxl_col"],
        -blank_df["pxl_row"],
        s=20,
        facecolors="none",
        edgecolors="#00A6A6",
        linewidths=0.7,
        label="Stage3B blank",
    )
    ax0.set_title(f"A. {args.target_label} marker field with Stage3B blank mask", fontsize=9, weight="bold")
    ax0.set_xticks([])
    ax0.set_yticks([])
    ax0.set_aspect("equal")
    ax0.legend(frameon=False, loc="lower right", fontsize=7)
    cb = fig.colorbar(base_scatter, ax=ax0, fraction=0.046, pad=0.02)
    cb.set_label("marker score", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    ax1 = fig.add_subplot(gs[0, 1])
    colors = sns.color_palette("tab20", n_colors=len(top_fraction))
    left = 0.0
    for (cell_type, value), color in zip(top_fraction.items(), colors):
        ax1.barh([0], [value], left=left, color=color, edgecolor="white", height=0.46, label=cell_type)
        left += value
    ax1.set_xlim(0, 1)
    ax1.set_yticks([0])
    ax1.set_yticklabels(["CytoSPACE\nbaseline"], fontsize=8)
    ax1.set_xlabel("fraction of assignments inside Stage3B blank spots", fontsize=8)
    ax1.set_title("B. Forced assignment destinations", fontsize=9, weight="bold")
    ax1.legend(
        bbox_to_anchor=(0.5, -0.24),
        loc="upper center",
        frameon=False,
        fontsize=6.2,
        borderaxespad=0,
        ncol=3,
        handlelength=1.2,
        columnspacing=0.8,
    )

    ax2 = fig.add_subplot(gs[0, 2])
    burden_df = pd.DataFrame(
        {
            "method": ["CytoSPACE", "SVTuner + CytoSPACE"],
            "assignment_burden": [baseline_burden.mean(), svtuner_burden.mean()],
            "sem": [baseline_burden.sem(), svtuner_burden.sem()],
        }
    )
    ax2.bar(
        burden_df["method"],
        burden_df["assignment_burden"],
        yerr=burden_df["sem"],
        color=["#CC6677", "#2A9D8F"],
        edgecolor="black",
        linewidth=0.5,
        capsize=3,
    )
    for i, values in enumerate([baseline_burden, svtuner_burden]):
        jitter = np.random.default_rng(42 + i).normal(0, 0.035, size=len(values))
        ax2.scatter(np.full(len(values), i) + jitter, values, s=5, color="black", alpha=0.22, linewidths=0)
    ax2.set_ylim(-0.03, 1.08)
    ax2.set_ylabel("assignment burden\nin blank spots", fontsize=8)
    ax2.set_title("C. Error burden removed by pre-mapping blanking", fontsize=9, weight="bold")
    ax2.tick_params(axis="x", rotation=25, labelsize=8)

    fig.suptitle(
        "Reference-missing case study on cell2location mouse brain Visium",
        fontsize=11,
        weight="bold",
        y=1.02,
    )
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / "cell2location_mouse_brain_stage3b_case_first_version.png"
    pdf = out_dir / "cell2location_mouse_brain_stage3b_case_first_version.pdf"
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    summary = {
        "sample": sample,
        "target_label": args.target_label,
        "markers_used": markers[:30],
        "stage3b_blank_spots": int(len(common_blank)),
        "baseline_mean_assignment_burden": float(baseline_burden.mean()),
        "svtuner_mean_assignment_burden": float(svtuner_burden.mean()),
        "baseline_top_forced_assignment_fraction": top_fraction.to_dict(),
    }
    (out_dir / "cell2location_mouse_brain_stage3b_case_first_version_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    top_fraction.rename("fraction").to_csv(out_dir / "cell2location_mouse_brain_forced_assignment_destinations.csv")
    burden_df.to_csv(out_dir / "cell2location_mouse_brain_assignment_burden.csv", index=False)
    print(f"[saved] {png}")
    print(f"[saved] {pdf}")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
