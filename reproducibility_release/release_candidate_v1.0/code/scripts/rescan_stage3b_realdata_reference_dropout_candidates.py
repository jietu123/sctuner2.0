#!/usr/bin/env python
from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def slug(value: str) -> str:
    import re

    return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", str(value).lower())).strip("_")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Build, run, and evaluate Stage3B real-data reference-dropout candidates "
            "for marker-region blanking accuracy."
        )
    )
    p.add_argument("--project_root", default=".")
    p.add_argument(
        "--candidate_csv",
        default="visualizations/stage3b_realdata_candidate_scan/fig2d_source_st_target_candidate_scan.csv",
    )
    p.add_argument("--target_group", default="stage3b_realdata_reference_dropout")
    p.add_argument("--top_k_per_pair", type=int, default=0, help="0 means use all candidates.")
    p.add_argument("--overwrite_prepare", action="store_true")
    p.add_argument("--marker_gene_union_only", action="store_true")
    p.add_argument("--rerun_stage3b", action="store_true")
    p.add_argument("--n_spatial_permutations", type=int, default=200)
    p.add_argument("--stage3b_max_genes", type=int, default=0)
    p.add_argument(
        "--run_order",
        choices=["pair_id", "st_spots"],
        default="pair_id",
        help="Order candidates before running Stage3B.",
    )
    p.add_argument("--out_dir", default="visualizations/stage3b_realdata_candidate_scan")
    return p.parse_args()


def read_yaml(path: Path) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def result_dir(root: Path, sample: str, cfg: dict) -> Path:
    group = ((cfg.get("storage") or {}).get("group") or "").strip().strip("/\\")
    return root / "result" / group / sample if group else root / "result" / sample


def processed_dir(root: Path, sample: str, cfg: dict) -> Path:
    group = ((cfg.get("storage") or {}).get("group") or "").strip().strip("/\\")
    return root / "data" / "processed" / group / sample if group else root / "data" / "processed" / sample


def stage1_export_dir(root: Path, sample: str, group: str) -> Path:
    return root / "data" / "processed" / group / sample / "stage1_preprocess" / "exported"


def marker_scores(st_expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in st_expr.columns]
    if not present:
        raise ValueError("No marker genes found in ST expression matrix")
    score = st_expr[present].mean(axis=1)
    return score.astype(float)


def bool_series(frame: pd.DataFrame, column: str) -> pd.Series:
    values = frame[column]
    if values.dtype == bool:
        return values
    return values.astype(str).str.lower().isin(["true", "1", "yes"])


def evaluate_candidate(root: Path, row: pd.Series, target_group: str) -> tuple[list[dict], pd.DataFrame]:
    source_sample = str(row["source_sample"])
    source_group = str(row["storage_group"])
    target_type = str(row["cell_type"])
    sample = f"{source_sample}_sc_missing_{slug(target_type)}"
    cfg = read_yaml(root / "configs" / "datasets" / f"{sample}.yaml")
    scores_path = processed_dir(root, sample, cfg) / "stage3b_st_unsupported" / "spot_unsupported_scores.csv"
    if not scores_path.exists():
        raise FileNotFoundError(f"Missing Stage3B scores: {scores_path}")

    export = stage1_export_dir(root, source_sample, source_group)
    st_expr = pd.read_csv(export / "st_expression_normalized.csv", index_col=0)
    st_expr.index = st_expr.index.astype(str)

    genes = [g for g in str(row["marker_genes_list"]).split(";") if g]
    score = marker_scores(st_expr, genes)
    percentile = score.rank(method="average", pct=True)

    stage3b = pd.read_csv(scores_path, index_col=0)
    stage3b.index = stage3b.index.astype(str)
    common = score.index.intersection(stage3b.index)
    if common.empty:
        raise ValueError(f"No common spots between marker score and Stage3B scores for {sample}")
    score = score.loc[common]
    percentile = percentile.loc[common]
    stage3b = stage3b.loc[common]
    blank = bool_series(stage3b, "is_unsupported_region")

    spot_rows = pd.DataFrame(
        {
            "sample": sample,
            "pair_id": row["pair_id"],
            "source_sample": source_sample,
            "target_type": target_type,
            "marker_score": score,
            "marker_percentile": percentile,
            "stage3b_blank": blank,
        },
        index=common,
    )
    spot_rows.index.name = "spot_id"

    out: list[dict] = []
    for top_label, quantile in [("top15", 0.85), ("top20", 0.80)]:
        threshold = float(score.quantile(quantile))
        target = score >= threshold
        overlap = int((target & blank).sum())
        target_n = int(target.sum())
        blank_n = int(blank.sum())
        recall = overlap / target_n if target_n else 0.0
        precision = overlap / blank_n if blank_n else 0.0
        denom = recall + precision
        f1 = 2 * recall * precision / denom if denom else 0.0
        union = int((target | blank).sum())
        jaccard = overlap / union if union else 0.0
        blank_percentiles = percentile.loc[blank]
        target_percentiles = percentile.loc[target]
        out.append(
            {
                "pair_id": row["pair_id"],
                "source_sample": source_sample,
                "sample": sample,
                "target_type": target_type,
                "top_region": top_label,
                "st_spots": int(len(common)),
                "target_region_spots": target_n,
                "stage3b_blank_spots": blank_n,
                "overlap_spots": overlap,
                "recall": recall,
                "precision": precision,
                "f1": f1,
                "jaccard": jaccard,
                "target_marker_percentile_median": float(target_percentiles.median()) if target_n else math.nan,
                "blank_marker_percentile_median": (
                    float(blank_percentiles.median()) if blank_n else math.nan
                ),
                "blank_marker_percentile_mean": float(blank_percentiles.mean()) if blank_n else math.nan,
                "marker_genes": int(row["marker_genes"]),
                "sc_cells": int(row["sc_cells"]),
                "moran_i": float(row["moran_i"]),
                "candidate_score": float(row["candidate_score"]),
            }
        )
    return out, spot_rows.reset_index()


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    candidates = pd.read_csv(root / args.candidate_csv)
    candidates = candidates[
        ~candidates["cell_type"].astype(str).str.lower().str.contains("unresolved", na=False)
    ].copy()
    if args.run_order == "st_spots" and "st_spots" in candidates.columns:
        candidates = candidates.sort_values(
            ["st_spots", "pair_id", "candidate_score"],
            ascending=[True, True, False],
        )
    else:
        candidates = candidates.sort_values(["pair_id", "candidate_score"], ascending=[True, False])
    if args.top_k_per_pair and args.top_k_per_pair > 0:
        selected = candidates.groupby("pair_id", group_keys=False).head(args.top_k_per_pair).copy()
    else:
        selected = candidates.copy()

    selected_path = out_dir / "stage3b_realdata_reference_dropout_rescan_candidates.csv"
    selected.to_csv(selected_path, index=False)
    print(f"[select] candidates={len(selected)} -> {selected_path}", flush=True)

    prepare_cmd = [
        sys.executable,
        "scripts/prepare_stage3b_realdata_reference_dropout_scenarios.py",
        "--project_root",
        str(root),
        "--selected_csv",
        str(selected_path.relative_to(root)),
        "--target_group",
        args.target_group,
    ]
    if args.overwrite_prepare:
        prepare_cmd.append("--overwrite")
    if args.marker_gene_union_only:
        prepare_cmd.append("--marker_gene_union_only")
    subprocess.run(prepare_cmd, cwd=root, check=True)

    for _, row in selected.iterrows():
        sample = f"{row['source_sample']}_sc_missing_{slug(row['cell_type'])}"
        cfg_path = root / "configs" / "datasets" / f"{sample}.yaml"
        cfg = read_yaml(cfg_path)
        scores_path = processed_dir(root, sample, cfg) / "stage3b_st_unsupported" / "spot_unsupported_scores.csv"
        if scores_path.exists() and not args.rerun_stage3b:
            print(f"[stage3b] skip existing {sample}", flush=True)
            continue
        cmd = [
            sys.executable,
            "-m",
            "src.stages.stage3b_st_unsupported",
            "--project_root",
            str(root),
            "--sample",
            sample,
            "--n_spatial_permutations",
            str(args.n_spatial_permutations),
            "--max_genes",
            str(args.stage3b_max_genes),
        ]
        print(f"[stage3b] run {sample}", flush=True)
        subprocess.run(cmd, cwd=root, check=True)

    metric_rows: list[dict] = []
    spot_frames: list[pd.DataFrame] = []
    for _, row in selected.iterrows():
        rows, spots = evaluate_candidate(root, row, args.target_group)
        metric_rows.extend(rows)
        spot_frames.append(spots)

    metrics = pd.DataFrame(metric_rows)
    metrics_path = out_dir / "stage3b_realdata_reference_dropout_rescan_overlap_metrics.csv"
    metrics.to_csv(metrics_path, index=False)

    spots = pd.concat(spot_frames, ignore_index=True)
    spots_path = out_dir / "stage3b_realdata_reference_dropout_rescan_spot_marker_scores.csv"
    spots.to_csv(spots_path, index=False)

    top15 = metrics[metrics["top_region"].eq("top15")].copy()
    top15["selection_score"] = (
        top15["precision"].fillna(0) * 0.35
        + top15["recall"].fillna(0) * 0.25
        + top15["blank_marker_percentile_median"].fillna(0) * 0.40
    )
    ranked = top15.sort_values(
        ["selection_score", "precision", "blank_marker_percentile_median", "recall"],
        ascending=False,
    )
    ranked_path = out_dir / "stage3b_realdata_reference_dropout_rescan_ranked_top15.csv"
    ranked.to_csv(ranked_path, index=False)

    print(f"[done] metrics: {metrics_path}", flush=True)
    print(f"[done] spot scores: {spots_path}", flush=True)
    print(f"[done] ranked top15: {ranked_path}", flush=True)
    print(ranked.head(20)[[
        "pair_id",
        "target_type",
        "stage3b_blank_spots",
        "overlap_spots",
        "recall",
        "precision",
        "blank_marker_percentile_median",
        "selection_score",
    ]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
