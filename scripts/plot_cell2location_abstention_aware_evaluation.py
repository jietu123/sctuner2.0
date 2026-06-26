#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Plot abstention-aware evaluation for the cell2location Stage3B case."
    )
    p.add_argument("--project_root", default=".")
    p.add_argument(
        "--sample",
        default="cell2loc_scan_st8059051_sc_missing_thalamic_excitatory",
    )
    p.add_argument("--section", default="ST8059051")
    p.add_argument("--target", default="thalamic_excitatory")
    p.add_argument("--target_label", default="Thalamic excitatory")
    p.add_argument("--n_bins", type=int, default=10)
    p.add_argument("--out_dir", default="visualizations/cell2location_stage3b_case")
    return p.parse_args()


def read_fractional(root: Path, sample: str, suffix: str) -> pd.DataFrame:
    path = (
        root
        / "result"
        / sample
        / f"stage4_cytospace{suffix}"
        / "cytospace_output"
        / "fractional_abundances_by_spot.csv"
    )
    if not path.exists():
        raise FileNotFoundError(f"Missing fractional abundance file: {path}")
    frame = pd.read_csv(path, index_col=0)
    frame.index = frame.index.astype(str)
    return frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    marker_prefix = f"cell2location_{args.section}_{args.target}_marker_region"
    marker_path = out_dir / f"{marker_prefix}.csv"
    if not marker_path.exists():
        raise FileNotFoundError(f"Missing marker-region CSV: {marker_path}")

    blank_path = (
        root
        / "result"
        / args.sample
        / "stage4_cytospace_stage3b_blank"
        / "cytospace_output"
        / "stage3b_blank_spots.csv"
    )
    if not blank_path.exists():
        raise FileNotFoundError(f"Missing Stage3B blank manifest: {blank_path}")
    stage3b_scores_path = (
        root
        / "data"
        / "processed"
        / "cell2location_mouse_brain"
        / args.sample
        / "stage3b_st_unsupported"
        / "spot_unsupported_scores.csv"
    )
    if not stage3b_scores_path.exists():
        raise FileNotFoundError(f"Missing Stage3B spot scores: {stage3b_scores_path}")

    marker = pd.read_csv(marker_path, index_col=0)
    marker.index = marker.index.astype(str)
    marker["marker_score"] = pd.to_numeric(marker["marker_score"], errors="coerce")
    marker["target_region"] = marker["target_region"].astype(bool)

    blank = pd.read_csv(blank_path, index_col=0)
    blank.index = blank.index.astype(str)
    blank_ids = set(blank.index)

    stage3b = pd.read_csv(stage3b_scores_path, index_col=0)
    stage3b.index = stage3b.index.astype(str)
    stage3b["spot_pvalue"] = pd.to_numeric(stage3b["spot_pvalue"], errors="coerce").fillna(1.0)
    stage3b["unsupported_evidence"] = -np.log10(
        stage3b["spot_pvalue"].clip(lower=np.finfo(float).tiny)
    )
    stage3b["is_spot_candidate"] = stage3b["is_spot_candidate"].astype(bool)

    baseline = read_fractional(root, args.sample, "_baseline")
    svtuner = read_fractional(root, args.sample, "_stage3b_blank")

    common = marker.index.intersection(stage3b.index).intersection(baseline.index).intersection(svtuner.index)
    if common.empty:
        raise ValueError("No common spots across marker, Stage3B, baseline, and SVTuner outputs")
    marker = marker.loc[common].copy()
    stage3b = stage3b.loc[common].copy()
    baseline = baseline.loc[common].copy()
    svtuner = svtuner.loc[common].copy()

    marker["marker_percentile"] = marker["marker_score"].rank(method="average", pct=True)
    marker["unsupported_evidence"] = stage3b["unsupported_evidence"]
    marker["unsupported_evidence_percentile"] = marker["unsupported_evidence"].rank(method="average", pct=True)
    marker["bin_id"] = np.minimum(
        np.floor(marker["unsupported_evidence_percentile"] * args.n_bins).astype(int),
        args.n_bins - 1,
    )
    max_bin = int(marker["bin_id"].max())
    marker["bin_label"] = marker["bin_id"].map(lambda x: f"{x + 1}")
    marker["stage3b_blank"] = marker.index.to_series().isin(blank_ids).to_numpy()

    baseline_mapped = baseline.sum(axis=1) > 0
    svtuner_mapped = svtuner.sum(axis=1) > 0

    # Soft abstention-aware score:
    # high Stage3B unsupported-evidence percentile favors blanking;
    # low unsupported-evidence percentile favors keeping mapped.
    # This avoids treating blank rows as ordinary mapping errors and avoids evaluating only blank spots.
    p_blank_desired = marker["unsupported_evidence_percentile"].clip(0, 1)
    baseline_score = np.where(baseline_mapped, 1.0 - p_blank_desired, p_blank_desired)
    svtuner_score = np.where(svtuner_mapped, 1.0 - p_blank_desired, p_blank_desired)
    marker["baseline_score"] = baseline_score
    marker["svtuner_score"] = svtuner_score
    marker["baseline_forced_assignment"] = baseline_mapped.to_numpy()
    marker["svtuner_mapped"] = svtuner_mapped.to_numpy()

    rows = []
    for bin_id, sub in marker.groupby("bin_id", sort=True):
        rows.append(
            {
                "bin_id": int(bin_id),
                "bin_label": f"{int(bin_id) + 1}",
                "n_spots": int(len(sub)),
                "mean_marker_percentile": float(sub["marker_percentile"].mean()),
                "mean_marker_score": float(sub["marker_score"].mean()),
                "mean_unsupported_evidence_percentile": float(
                    sub["unsupported_evidence_percentile"].mean()
                ),
                "mean_unsupported_evidence": float(sub["unsupported_evidence"].mean()),
                "target_region_fraction": float(sub["target_region"].mean()),
                "stage3b_candidate_fraction": float(stage3b.loc[sub.index, "is_spot_candidate"].mean()),
                "baseline_forced_assignment_rate": float(sub["baseline_forced_assignment"].mean()),
                "svtuner_blank_rate": float(sub["stage3b_blank"].mean()),
                "svtuner_preserved_mapping_rate": float(sub["svtuner_mapped"].mean()),
                "baseline_abstention_aware_score": float(sub["baseline_score"].mean()),
                "svtuner_abstention_aware_score": float(sub["svtuner_score"].mean()),
            }
        )
    bins = pd.DataFrame(rows)

    fig, (ax0, ax1) = plt.subplots(
        1,
        2,
        figsize=(10.4, 4.2),
        dpi=300,
        gridspec_kw={"width_ratios": [1.08, 1.0], "wspace": 0.32},
    )

    x = bins["bin_id"].to_numpy()
    labels = bins["bin_label"].tolist()

    ax0.axvspan(max_bin - 1.5, max_bin + 0.5, color="#f2f2f2", zorder=0)
    ax0.plot(
        x,
        bins["baseline_forced_assignment_rate"],
        color="#d95f02",
        marker="o",
        linewidth=2.0,
        markersize=4.2,
        label="CytoSPACE forced assignment",
    )
    ax0.plot(
        x,
        bins["svtuner_blank_rate"],
        color="#1b9e77",
        marker="o",
        linewidth=2.0,
        markersize=4.2,
        label="SVTuner Stage3B blank rate",
    )
    ax0.plot(
        x,
        bins["stage3b_candidate_fraction"],
        color="#666666",
        linestyle=":",
        marker="^",
        linewidth=1.4,
        markersize=3.4,
        label="Stage3B spot-candidate fraction",
    )
    ax0.set_ylim(-0.04, 1.04)
    ax0.set_xticks(x)
    ax0.set_xticklabels(labels)
    ax0.set_xlabel("Stage3B unsupported-evidence decile", fontsize=9)
    ax0.set_ylabel("Fraction of spots", fontsize=9)
    ax0.set_title("A. Selective blanking along unsupported-evidence gradient", fontsize=10, weight="bold")
    ax0.grid(axis="y", color="#d9d9d9", linewidth=0.7)
    ax0.legend(frameon=False, fontsize=7.2, loc="upper left")
    ax0.text(
        max_bin - 1.45,
        0.07,
        "highest\nunsupported\nevidence",
        fontsize=7.2,
        color="#555555",
        va="bottom",
    )

    ax1.axvspan(max_bin - 1.5, max_bin + 0.5, color="#f2f2f2", zorder=0)
    ax1.plot(
        x,
        bins["baseline_abstention_aware_score"],
        color="#d95f02",
        marker="o",
        linewidth=2.0,
        markersize=4.2,
        label="CytoSPACE",
    )
    ax1.plot(
        x,
        bins["svtuner_abstention_aware_score"],
        color="#1b9e77",
        marker="o",
        linewidth=2.0,
        markersize=4.2,
        label="SVTuner + Stage3B blanking",
    )
    ax1.set_ylim(-0.04, 1.04)
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels)
    ax1.set_xlabel("Stage3B unsupported-evidence decile", fontsize=9)
    ax1.set_ylabel("Abstention-aware compatibility score", fontsize=9)
    ax1.set_title("B. Mapping/blanking evaluated by regional support", fontsize=10, weight="bold")
    ax1.grid(axis="y", color="#d9d9d9", linewidth=0.7)
    ax1.legend(frameon=False, fontsize=7.2, loc="upper right")

    fig.suptitle(
        "Abstention-aware evaluation of ST-only unsupported regions",
        fontsize=12,
        weight="bold",
        y=1.03,
    )

    prefix = f"cell2location_{args.section}_{args.target}_abstention_aware_evaluation"
    png = out_dir / f"{prefix}.png"
    pdf = out_dir / f"{prefix}.pdf"
    csv = out_dir / f"{prefix}.csv"
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    bins.to_csv(csv, index=False)
    summary = {
        "sample": args.sample,
        "section": args.section,
        "target_label": args.target_label,
        "n_spots": int(len(marker)),
        "stage3b_blank_spots": int(marker["stage3b_blank"].sum()),
        "marker_defined_target_region_spots": int(marker["target_region"].sum()),
        "stage3b_candidate_spots": int(stage3b["is_spot_candidate"].sum()),
        "baseline_mean_abstention_aware_score": float(marker["baseline_score"].mean()),
        "svtuner_mean_abstention_aware_score": float(marker["svtuner_score"].mean()),
        "baseline_high_unsupported_evidence_score": float(
            bins.loc[bins["bin_id"].isin([max_bin - 1, max_bin]), "baseline_abstention_aware_score"].mean()
        ),
        "svtuner_high_unsupported_evidence_score": float(
            bins.loc[bins["bin_id"].isin([max_bin - 1, max_bin]), "svtuner_abstention_aware_score"].mean()
        ),
        "png": str(png),
        "pdf": str(pdf),
        "csv": str(csv),
    }
    (out_dir / f"{prefix}_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
