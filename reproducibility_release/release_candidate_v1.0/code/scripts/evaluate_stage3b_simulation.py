#!/usr/bin/env python
"""Evaluate Stage3B detections against simulation truth after detection."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_ROOT = _HERE.parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.stages.storage import processed_dir, read_dataset_config, result_dir
from src.utils.sample_paths import resolve_sample_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--target_type", required=True)
    parser.add_argument("--sim_group", default="real_brca")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = Path(args.project_root).resolve()
    cfg = read_dataset_config(root, args.sample)
    score_path = (
        processed_dir(root, args.sample, cfg)
        / "stage3b_st_unsupported"
        / "spot_unsupported_scores.csv"
    )
    truth_path = (
        resolve_sample_dir(root, args.sample, sim_group=args.sim_group)
        / "sim_truth_spot_type_fraction.csv"
    )
    scores = pd.read_csv(score_path, index_col=0)
    truth = pd.read_csv(truth_path, index_col=0)
    scores.index = scores.index.astype(str)
    truth.index = truth.index.astype(str)
    truth = truth.reindex(scores.index)
    if args.target_type not in truth.columns:
        raise ValueError(f"Target type absent from truth: {args.target_type}")
    fraction = pd.to_numeric(truth[args.target_type], errors="coerce").fillna(0.0)
    detected = scores["is_unsupported_region"].astype(str).str.lower().eq("true")

    dominant = truth.apply(pd.to_numeric, errors="coerce").fillna(0.0).idxmax(axis=1)
    dominant_positive = dominant.eq(args.target_type)
    zero_target = fraction.eq(0.0)
    target_mass_total = float(fraction.sum())
    metrics = {
        "sample": args.sample,
        "target_type": args.target_type,
        "n_spots": int(len(scores)),
        "detected_region_spots": int(detected.sum()),
        "target_positive_spots": int((fraction > 0).sum()),
        "target_dominant_spots": int(dominant_positive.sum()),
        "dominant_region_recall": float(
            (detected & dominant_positive).sum() / max(int(dominant_positive.sum()), 1)
        ),
        "dominant_region_precision": float(
            (detected & dominant_positive).sum() / max(int(detected.sum()), 1)
        ),
        "target_mass_recall": float(
            fraction.loc[detected].sum() / max(target_mass_total, np.finfo(float).eps)
        ),
        "zero_target_false_positive_spots": int((detected & zero_target).sum()),
        "zero_target_false_positive_rate": float(
            (detected & zero_target).sum() / max(int(zero_target.sum()), 1)
        ),
        "mean_target_fraction_detected": float(fraction.loc[detected].mean())
        if detected.any()
        else 0.0,
        "mean_target_fraction_not_detected": float(fraction.loc[~detected].mean())
        if (~detected).any()
        else 0.0,
    }
    output_dir = result_dir(root, args.sample, cfg) / "stage3b_st_unsupported"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "stage3b_truth_evaluation.json"
    output_path.write_text(
        json.dumps(metrics, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(metrics, ensure_ascii=False, indent=2))
    print(f"[DONE] {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
