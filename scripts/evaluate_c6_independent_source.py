#!/usr/bin/env python3
"""Frozen C6 independent-source evaluation (reads completed outputs only)."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


GROUP = "c6_independent_source"
CONTROL = "c6_breast_independent_control"
DROPOUT = "c6_breast_independent_fibroblast_dropout"
TARGET = "Fibroblasts"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project_root", default=".")
    return parser.parse_args()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def object_sha256(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def load_fraction(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, index_col=0)
    frame.index = frame.index.astype(str)
    return frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def normalize_rows(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.div(frame.sum(axis=1).replace(0.0, np.nan), axis=0).fillna(0.0)


def aligned_overlap(prediction: pd.DataFrame, truth: pd.DataFrame) -> pd.Series:
    spots = truth.index.astype(str)
    columns = sorted(set(prediction.columns).union(truth.columns))
    pred = normalize_rows(prediction.reindex(index=spots, columns=columns, fill_value=0.0))
    obs = normalize_rows(truth.reindex(index=spots, columns=columns, fill_value=0.0))
    return pd.Series(np.minimum(pred.to_numpy(), obs.to_numpy()).sum(axis=1), index=spots)


def stage_paths(root: Path, sample: str) -> dict[str, Path]:
    processed = root / "data" / "processed" / "simulation_experiments" / GROUP / sample
    result_grouped = root / "result" / "simulation_experiments" / GROUP / sample
    result_plain = root / "result" / sample
    return {
        "sim": root / "data" / "sim" / GROUP / sample,
        "stage1": processed / "stage1_preprocess" / "exported",
        "stage3": processed / "stage3_typematch",
        "stage3_summary": result_grouped / "stage3_typematch" / "stage3_summary.json",
        "stage3b": processed / "stage3b_st_unsupported",
        "stage3b_summary": result_grouped / "stage3b_st_unsupported" / "stage3b_summary.json",
        "baseline": result_plain / "stage4_cytospace_baseline" / "cytospace_output",
        "svtuner": result_plain / "stage4_cytospace_c6_svtuner" / "cytospace_output",
    }


def stage3_audit(root: Path, sample: str) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    paths = stage_paths(root, sample)
    summary = load_json(paths["stage3_summary"])
    relabel = pd.read_csv(paths["stage3"] / "cell_type_relabel.csv")
    dropped = relabel["action"].astype(str).eq("Dropped")
    detected = [str(value) for value in summary["action_overview"].get("missing_types", [])]
    row = {
        "scenario": sample,
        "input_types": ";".join(sorted(relabel["orig_type"].astype(str).unique())),
        "detected_unsupported_type_count": len(detected),
        "detected_unsupported_types": ";".join(detected),
        "n_input_cells": int(len(relabel)),
        "n_cells_marked_unsupported_or_unknown": int(dropped.sum()),
        "false_exclusion_rate": float(dropped.mean()),
        "false_exclusion": bool(dropped.any()),
    }
    return row, summary, {"relabel": relabel, "dropped": dropped}


def safe_ratio(numerator: int, denominator: int) -> float:
    return float(numerator / denominator) if denominator else 0.0


def main() -> None:
    root = Path(parse_args().project_root).resolve()
    out_dir = root / "result" / "c6_independent_source_evaluation"
    out_dir.mkdir(parents=True, exist_ok=True)

    samples = [CONTROL, DROPOUT]
    paths = {sample: stage_paths(root, sample) for sample in samples}
    truth = {
        sample: load_fraction(paths[sample]["sim"] / "sim_truth_spot_type_fraction.csv")
        for sample in samples
    }
    predictions: dict[str, dict[str, pd.DataFrame]] = {}
    scores: dict[str, pd.DataFrame] = {}
    blank_ids: dict[str, set[str]] = {}
    stage3_rows: list[dict[str, Any]] = []
    stage3_summaries: dict[str, dict[str, Any]] = {}
    stage3_aux: dict[str, dict[str, Any]] = {}
    stage3b_summaries: dict[str, dict[str, Any]] = {}

    for sample in samples:
        audit, summary, aux = stage3_audit(root, sample)
        stage3_rows.append(audit)
        stage3_summaries[sample] = summary
        stage3_aux[sample] = aux
        stage3b_summaries[sample] = load_json(paths[sample]["stage3b_summary"])
        scores[sample] = pd.read_csv(
            paths[sample]["stage3b"] / "spot_unsupported_scores.csv", index_col=0
        )
        scores[sample].index = scores[sample].index.astype(str)
        mask = scores[sample]["is_unsupported_region"].astype(str).str.lower().eq("true")
        blank_ids[sample] = set(scores[sample].index[mask])
        predictions[sample] = {
            "CytoSPACE": load_fraction(paths[sample]["baseline"] / "fractional_abundances_by_spot.csv"),
            "SVTuner": load_fraction(paths[sample]["svtuner"] / "fractional_abundances_by_spot.csv"),
        }

    stage3_audit_frame = pd.DataFrame(stage3_rows)
    stage3_audit_frame.to_csv(out_dir / "c6_stage3a_audit.csv", index=False)

    dominant = normalize_rows(truth[DROPOUT]).idxmax(axis=1)
    truth_positive = dominant.eq(TARGET)
    predicted_positive = pd.Series(
        scores[DROPOUT].index.isin(blank_ids[DROPOUT]), index=scores[DROPOUT].index
    ).reindex(truth[DROPOUT].index, fill_value=False)
    tp = int((truth_positive & predicted_positive).sum())
    fp = int((~truth_positive & predicted_positive).sum())
    fn = int((truth_positive & ~predicted_positive).sum())
    tn = int((~truth_positive & ~predicted_positive).sum())
    pvalues = pd.to_numeric(scores[DROPOUT]["spot_pvalue"], errors="coerce").reindex(
        truth[DROPOUT].index
    )
    evidence = -np.log10(np.clip(pvalues.to_numpy(float), np.finfo(float).tiny, 1.0))
    labels = truth_positive.to_numpy(dtype=int)
    dropout_metrics = {
        "scenario": DROPOUT,
        "target": TARGET,
        "truth_unsupported_spots": int(truth_positive.sum()),
        "predicted_withheld": int(predicted_positive.sum()),
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "TN": tn,
        "precision": safe_ratio(tp, tp + fp),
        "recall": safe_ratio(tp, tp + fn),
        "supported_spot_FPR": safe_ratio(fp, fp + tn),
        "target_dominant_withheld_rate": safe_ratio(tp, tp + fn),
        "non_target_withheld_rate": safe_ratio(fp, fp + tn),
        "AUROC": float(roc_auc_score(labels, evidence)),
        "AP": float(average_precision_score(labels, evidence)),
    }
    pd.DataFrame([dropout_metrics]).to_csv(
        out_dir / "c6_stage3b_dropout_metrics.csv", index=False
    )

    summary_rows: list[dict[str, Any]] = []
    for sample in samples:
        truth_frame = truth[sample]
        target_dominant = normalize_rows(truth_frame).idxmax(axis=1).eq(TARGET)
        for method in ("CytoSPACE", "SVTuner"):
            overlap = aligned_overlap(predictions[sample][method], truth_frame)
            withheld = pd.Series(truth_frame.index.isin(blank_ids[sample]), index=truth_frame.index)
            if method == "CytoSPACE":
                spot_score = overlap
                n_withheld = 0
            elif sample == CONTROL:
                spot_score = overlap.copy()
                spot_score.loc[withheld] = 0.0
                n_withheld = int(withheld.sum())
            else:
                spot_score = overlap.copy()
                spot_score.loc[withheld & target_dominant] = 1.0
                spot_score.loc[withheld & ~target_dominant] = 0.0
                n_withheld = int(withheld.sum())
            summary_rows.append(
                {
                    "scenario": sample,
                    "method": method,
                    "whole_space_recovery": float(spot_score.mean()),
                    "n_spots": int(len(truth_frame)),
                    "n_withheld": n_withheld,
                    "withheld_rate": safe_ratio(n_withheld, len(truth_frame)),
                }
            )
    summary_frame = pd.DataFrame(summary_rows)
    summary_frame.to_csv(out_dir / "c6_summary_metrics.csv", index=False)

    sim_info = {sample: load_json(paths[sample]["sim"] / "sim_info.json") for sample in samples}
    stage3_ids = {
        sample: object_sha256(stage3_summaries[sample]["params"]) for sample in samples
    }
    stage3b_configs = {sample: stage3b_summaries[sample]["config"] for sample in samples}
    stage3b_ids = {sample: object_sha256(stage3b_configs[sample]) for sample in samples}

    score_blank_match: dict[str, bool] = {}
    prediction_id_match: dict[str, dict[str, bool]] = {}
    fraction_rows_valid: dict[str, dict[str, bool]] = {}
    for sample in samples:
        manifest = pd.read_csv(paths[sample]["svtuner"] / "stage3b_blank_spots.csv", index_col=0)
        manifest_ids = set(manifest.index.astype(str))
        score_blank_match[sample] = manifest_ids == blank_ids[sample]
        prediction_id_match[sample] = {}
        fraction_rows_valid[sample] = {}
        for method, frame in predictions[sample].items():
            prediction_id_match[sample][method] = set(frame.index) == set(truth[sample].index)
            sums = frame.sum(axis=1)
            if method == "SVTuner":
                nonblank = ~sums.index.isin(blank_ids[sample])
                fraction_rows_valid[sample][method] = bool((sums.loc[nonblank] > 0).all())
            else:
                fraction_rows_valid[sample][method] = bool((sums > 0).all())

    numeric_metrics = np.asarray(
        summary_frame.select_dtypes(include=[np.number]).to_numpy().ravel().tolist()
        + [value for value in dropout_metrics.values() if isinstance(value, (int, float))],
        dtype=float,
    )
    control_truth_path = paths[CONTROL]["sim"] / "sim_truth_spot_type_fraction.csv"
    dropout_truth_path = paths[DROPOUT]["sim"] / "sim_truth_spot_type_fraction.csv"
    control_st_path = paths[CONTROL]["sim"] / "brca_STdata_GEP.txt"
    dropout_st_path = paths[DROPOUT]["sim"] / "brca_STdata_GEP.txt"
    dropout_meta = pd.read_csv(paths[DROPOUT]["stage1"] / "sc_metadata.csv")
    qc = {
        "truth_spot_universe_identical": set(truth[CONTROL].index) == set(truth[DROPOUT].index),
        "truth_files_sha256_identical": sha256(control_truth_path) == sha256(dropout_truth_path),
        "st_expression_sha256_identical": sha256(control_st_path) == sha256(dropout_st_path),
        "dropout_reference_has_no_fibroblasts": not dropout_meta.astype(str).eq(TARGET).any().any(),
        "prediction_spot_ids_match_truth": prediction_id_match,
        "blank_manifest_matches_scores": score_blank_match,
        "fraction_rows_valid": fraction_rows_valid,
        "metrics_finite": bool(np.isfinite(numeric_metrics).all()),
        "whole_space_scores_in_0_1": bool(summary_frame["whole_space_recovery"].between(0, 1).all()),
        "auroc_ap_in_0_1": bool(
            0 <= dropout_metrics["AUROC"] <= 1 and 0 <= dropout_metrics["AP"] <= 1
        ),
    }
    qc["all_pass"] = bool(
        qc["truth_spot_universe_identical"]
        and qc["truth_files_sha256_identical"]
        and qc["st_expression_sha256_identical"]
        and qc["dropout_reference_has_no_fibroblasts"]
        and all(all(values.values()) for values in prediction_id_match.values())
        and all(score_blank_match.values())
        and all(all(values.values()) for values in fraction_rows_valid.values())
        and qc["metrics_finite"]
        and qc["whole_space_scores_in_0_1"]
        and qc["auroc_ap_in_0_1"]
    )

    provenance = {
        "dataset_A_reference_source": sim_info[CONTROL]["reference_source_sample"],
        "dataset_B_profile_source": sim_info[CONTROL]["profile_source_sample"],
        "independent_profile_source": sim_info[CONTROL]["independent_profile_source"],
        "profile_source_path": sim_info[CONTROL]["profile_source_path"],
        "dropout_target": TARGET,
        "seed": sim_info[CONTROL]["seed"],
        "cell_types": sim_info[CONTROL]["harmonized_types"],
        "three_way_shared_gene_count": sim_info[CONTROL]["three_way_shared_gene_count"],
        "stage3A_config_identity": stage3_ids,
        "stage3A_configs_identical": stage3_ids[CONTROL] == stage3_ids[DROPOUT],
        "stage3B_config_identity": stage3b_ids,
        "stage3B_configs_identical": stage3b_ids[CONTROL] == stage3b_ids[DROPOUT],
        "stage3B_resolved_config": stage3b_configs[CONTROL],
        "stage4_execution_settings": {
            "seed": 42,
            "n_processors": 1,
            "n_subspots": 800,
            "mapping_cells_per_spot": 5,
            "sc_expr_source": "normalized",
            "baseline_filter_mode": "none",
            "svtuner_filter_mode": "plugin_unknown",
            "svtuner_filter_scope": "unsupported_all",
            "svtuner_cell_type_column": "plugin_type",
            "svtuner_stage3b_blank_regions": True,
        },
        "truth_rule": "Fibroblast truth-dominant spot",
        "continuous_evidence_score": "-log10(clip(spot_pvalue, finfo(float).tiny, 1.0))",
        "file_sha256": {
            "control_truth": sha256(control_truth_path),
            "dropout_truth": sha256(dropout_truth_path),
            "control_st_expression": sha256(control_st_path),
            "dropout_st_expression": sha256(dropout_st_path),
        },
        "final_qc": qc,
    }
    with (out_dir / "c6_provenance.json").open("w", encoding="utf-8") as handle:
        json.dump(provenance, handle, indent=2, ensure_ascii=False)

    if not qc["all_pass"]:
        raise RuntimeError("C6 final QC failed; inspect c6_provenance.json")
    print(json.dumps({"status": "PASS", "qc": qc}, ensure_ascii=False))


if __name__ == "__main__":
    main()
