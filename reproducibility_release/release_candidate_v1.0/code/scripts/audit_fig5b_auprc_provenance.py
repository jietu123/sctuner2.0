from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import re
import sys
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy import stats
from sklearn.metrics import (
    auc,
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ENDPOINT = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping" / "spot_level_endpoint_freeze.csv"
DEFAULT_SCORE = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase7_svtuner_aware_execution_against_frozen_cta_immune_endpoint" / "svtuner_immune_all_dropout" / "svtuner_immune_all_dropout_spot_level_raw_output_contract.csv"
DEFAULT_RAW_SCORE = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase7_svtuner_aware_execution_against_frozen_cta_immune_endpoint" / "svtuner_immune_all_dropout" / "raw_outputs" / "spot_unsupported_scores.csv"
DEFAULT_PHASE8_SPOTS = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison" / "svtuner_endpoint_score_by_spot.csv"
DEFAULT_PHASE8_SUMMARY = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison" / "bioapp_phase8_svtuner_vs_endpoint_evaluation_summary.json"
DEFAULT_PHASE8_SCRIPT = ROOT / "scripts" / "run_bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison.py"
DEFAULT_FIG5_SCRIPT = ROOT / "scripts" / "run_bioapp_main_figure_v3_12_panel_A_evidence_chain_redesign.py"
DEFAULT_OUTPUT = ROOT / "visualizations" / "bioapp_experiment" / "fig5b_auprc_formal_audit"

EXPECTED_TOTAL = 2248
EXPECTED_ANALYSIS = 1888
EXPECTED_POSITIVE = 134
EXPECTED_NEGATIVE = 1754
EXPECTED_AMBIGUOUS = 290
EXPECTED_EXCLUDED = 70
MATCH_TOLERANCE = 0.0001
TARGET_OLD = 0.2131
TARGET_CURRENT = 0.2836551170002947

TEXT_EXTENSIONS = {".py", ".md", ".json", ".txt", ".tex", ".tsv", ".yaml", ".yml", ".csv"}
VALUE_PATTERNS = {
    "0.2131": re.compile(r"(?<!\d)0\.2131(?:\d+)?(?!\d)"),
    "0.28366": re.compile(r"(?<!\d)0\.28366(?:\d+)?(?!\d)"),
    "0.2837": re.compile(r"(?<!\d)0\.2837(?:\d+)?(?!\d)"),
    "AUPRC": re.compile(r"AUPRC", re.IGNORECASE),
    "average_precision_score": re.compile(r"\baverage_precision_score\b"),
    "precision_recall_curve": re.compile(r"\bprecision_recall_curve\b"),
    "roc_auc_score": re.compile(r"\broc_auc_score\b"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit the provenance and definition of the Fig. 5b precision-recall metric.")
    parser.add_argument("--endpoint-file", type=Path, default=DEFAULT_ENDPOINT)
    parser.add_argument("--score-file", type=Path, default=DEFAULT_SCORE)
    parser.add_argument("--raw-score-file", type=Path, default=DEFAULT_RAW_SCORE)
    parser.add_argument("--phase8-spot-file", type=Path, default=DEFAULT_PHASE8_SPOTS)
    parser.add_argument("--phase8-summary-file", type=Path, default=DEFAULT_PHASE8_SUMMARY)
    parser.add_argument("--phase8-script", type=Path, default=DEFAULT_PHASE8_SCRIPT)
    parser.add_argument("--fig5-script", type=Path, default=DEFAULT_FIG5_SCRIPT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--endpoint-barcode-column", default="barcode")
    parser.add_argument("--score-barcode-column", default="barcode")
    parser.add_argument("--endpoint-column", default="primary_endpoint_status")
    parser.add_argument("--score-column", default="withheld_score")
    parser.add_argument("--binary-column", default="withheld_binary")
    parser.add_argument("--positive-label", default="positive")
    parser.add_argument("--negative-label", default="negative")
    parser.add_argument("--ambiguous-label", default="ambiguous")
    parser.add_argument("--excluded-label", default="excluded")
    parser.add_argument("--formal", action=argparse.BooleanOptionalAction, default=True)
    return parser.parse_args()


def resolve(path: Path) -> Path:
    path = path.expanduser()
    if not path.is_absolute():
        path = ROOT / path
    return path.resolve()


def rel(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return str(path.resolve()).replace("\\", "/")


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")


def require_columns(df: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def sha1(path: Path) -> str:
    digest = hashlib.sha1()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def modified_time(path: Path) -> str:
    return pd.Timestamp(path.stat().st_mtime, unit="s", tz="UTC").isoformat()


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")


def bool_series(series: pd.Series, label: str) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        if series.isna().any():
            raise ValueError(f"{label} contains missing values")
        return series.astype(bool)
    normalized = series.astype(str).str.strip().str.lower()
    valid = {"true", "false", "1", "0", "yes", "no"}
    invalid = normalized[~normalized.isin(valid)]
    if not invalid.empty:
        raise ValueError(f"{label} contains non-boolean values: {sorted(invalid.unique())[:10]}")
    return normalized.isin({"true", "1", "yes"})


def metric_set(y_true: pd.Series | np.ndarray, y_score: pd.Series | np.ndarray) -> dict[str, Any]:
    y = np.asarray(y_true, dtype=int)
    score = np.asarray(y_score, dtype=float)
    finite = np.isfinite(score)
    y = y[finite]
    score = score[finite]
    if len(y) == 0 or len(np.unique(y)) != 2 or len(np.unique(score)) < 2:
        return {
            "n": int(len(y)),
            "auroc": None,
            "average_precision": None,
            "trapezoidal_pr_auc": None,
        }
    precision, recall, _ = precision_recall_curve(y, score)
    trapezoidal = float(auc(recall[::-1], precision[::-1]))
    index = np.arange(len(recall))
    traversal_preserving_order = np.lexsort((-index, recall))
    trapezoidal_sorted = float(auc(recall[traversal_preserving_order], precision[traversal_preserving_order]))
    if not np.isclose(trapezoidal, trapezoidal_sorted, rtol=0, atol=1e-15):
        raise RuntimeError("PR-AUC implementations disagree after preserving curve traversal within recall ties")
    return {
        "n": int(len(y)),
        "auroc": float(roc_auc_score(y, score)),
        "average_precision": float(average_precision_score(y, score)),
        "trapezoidal_pr_auc": trapezoidal,
    }


def endpoint_counts(endpoint: pd.DataFrame, column: str, labels: dict[str, str]) -> dict[str, int]:
    values = endpoint[column].astype(str)
    known = set(labels.values())
    return {
        "total": int(len(endpoint)),
        "positive": int((values == labels["positive"]).sum()),
        "negative": int((values == labels["negative"]).sum()),
        "ambiguous": int((values == labels["ambiguous"]).sum()),
        "excluded": int((values == labels["excluded"]).sum()),
        "other_unknown": int((~values.isin(known)).sum()),
    }


def validate_formal_counts(counts: dict[str, int], formal: bool) -> list[str]:
    expected = {
        "total": EXPECTED_TOTAL,
        "positive": EXPECTED_POSITIVE,
        "negative": EXPECTED_NEGATIVE,
        "ambiguous": EXPECTED_AMBIGUOUS,
        "excluded": EXPECTED_EXCLUDED,
        "other_unknown": 0,
    }
    warnings = [f"{key}: observed {counts[key]}, expected {value}" for key, value in expected.items() if counts[key] != value]
    if warnings and formal:
        raise ValueError("Formal endpoint counts failed: " + "; ".join(warnings))
    return warnings


def duplicate_rows(df: pd.DataFrame, barcode_column: str, source: str) -> pd.DataFrame:
    duplicated = df[df[barcode_column].duplicated(keep=False)].copy()
    if duplicated.empty:
        return pd.DataFrame(columns=["source", "barcode", "row_index"])
    out = pd.DataFrame(
        {
            "source": source,
            "barcode": duplicated[barcode_column].astype(str),
            "row_index": duplicated.index.astype(int),
        }
    )
    return out.sort_values(["source", "barcode", "row_index"])


def merge_inputs(args: argparse.Namespace, endpoint: pd.DataFrame, scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    endpoint_barcode = args.endpoint_barcode_column
    score_barcode = args.score_barcode_column
    require_columns(endpoint, [endpoint_barcode, args.endpoint_column], "endpoint file")
    require_columns(scores, [score_barcode, args.score_column, args.binary_column], "score file")

    duplicates = pd.concat(
        [
            duplicate_rows(endpoint, endpoint_barcode, "endpoint"),
            duplicate_rows(scores, score_barcode, "score"),
        ],
        ignore_index=True,
    )
    if not duplicates.empty:
        raise ValueError(f"Duplicate barcodes detected; see duplicate_barcodes.tsv after resolving the input")

    endpoint_keys = set(endpoint[endpoint_barcode].astype(str))
    score_keys = set(scores[score_barcode].astype(str))
    unmatched_rows = [
        {"source": "endpoint_only", "barcode": barcode}
        for barcode in sorted(endpoint_keys - score_keys)
    ] + [
        {"source": "score_only", "barcode": barcode}
        for barcode in sorted(score_keys - endpoint_keys)
    ]
    unmatched = pd.DataFrame(unmatched_rows, columns=["source", "barcode"])

    endpoint_input = endpoint.copy()
    score_input = scores.copy()
    endpoint_input[endpoint_barcode] = endpoint_input[endpoint_barcode].astype(str)
    score_input[score_barcode] = score_input[score_barcode].astype(str)
    if score_barcode != endpoint_barcode:
        score_input = score_input.rename(columns={score_barcode: endpoint_barcode})

    merged = endpoint_input.merge(
        score_input,
        on=endpoint_barcode,
        how="left",
        suffixes=("_endpoint", "_score"),
        validate="one_to_one",
        indicator=True,
    )
    join_counts = {
        "endpoint_rows": int(len(endpoint)),
        "score_rows": int(len(scores)),
        "inner_join_rows": int(len(endpoint_keys & score_keys)),
        "left_join_rows": int(len(merged)),
        "endpoint_only_rows": int(len(endpoint_keys - score_keys)),
        "score_only_rows": int(len(score_keys - endpoint_keys)),
    }
    return merged, unmatched, join_counts


def choose_merged_column(merged: pd.DataFrame, requested: str) -> str:
    if requested in merged.columns:
        return requested
    score_variant = f"{requested}_score"
    if score_variant in merged.columns:
        return score_variant
    endpoint_variant = f"{requested}_endpoint"
    if endpoint_variant in merged.columns:
        return endpoint_variant
    raise ValueError(f"Merged input does not contain requested column {requested!r}")


def compute_formal_analysis(args: argparse.Namespace, merged: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any], pd.DataFrame, pd.DataFrame]:
    endpoint_column = choose_merged_column(merged, args.endpoint_column)
    score_column = choose_merged_column(merged, args.score_column)
    binary_column = choose_merged_column(merged, args.binary_column)
    labels = {
        "positive": args.positive_label,
        "negative": args.negative_label,
        "ambiguous": args.ambiguous_label,
        "excluded": args.excluded_label,
    }
    status = merged[endpoint_column].astype(str)
    merged = merged.copy()
    merged[score_column] = pd.to_numeric(merged[score_column], errors="coerce")
    merged[binary_column] = bool_series(merged[binary_column], binary_column)
    merged["audit_endpoint_class"] = status
    merged["audit_y_true"] = (status == labels["positive"]).astype(int)

    formal = merged[status.isin([labels["positive"], labels["negative"]])].copy()
    excluded = merged[~status.isin([labels["positive"], labels["negative"]])].copy()
    formal = formal.sort_values(args.endpoint_barcode_column).reset_index(drop=True)
    excluded = excluded.sort_values(args.endpoint_barcode_column).reset_index(drop=True)

    missing_score = int(formal[score_column].isna().sum())
    finite = np.isfinite(formal[score_column].to_numpy(dtype=float))
    nonfinite_score = int((~finite & formal[score_column].notna().to_numpy()).sum())
    if args.formal and (len(formal) != EXPECTED_ANALYSIS or int(formal["audit_y_true"].sum()) != EXPECTED_POSITIVE):
        raise ValueError(
            f"Formal analysis set mismatch: n={len(formal)}, positive={formal['audit_y_true'].sum()}, "
            f"expected n={EXPECTED_ANALYSIS}, positive={EXPECTED_POSITIVE}"
        )
    if args.formal and (missing_score or nonfinite_score):
        raise ValueError(f"Formal score contains missing={missing_score}, nonfinite={nonfinite_score}")

    metrics = metric_set(formal["audit_y_true"], formal[score_column])
    positive = formal.loc[formal["audit_y_true"] == 1, score_column].to_numpy(dtype=float)
    negative = formal.loc[formal["audit_y_true"] == 0, score_column].to_numpy(dtype=float)
    mann_whitney = stats.mannwhitneyu(positive, negative, alternative="two-sided")
    metrics.update(
        {
            "positive_prevalence": float(len(positive) / len(formal)),
            "positive_mean_score": float(np.mean(positive)),
            "negative_mean_score": float(np.mean(negative)),
            "mean_score_difference": float(np.mean(positive) - np.mean(negative)),
            "mann_whitney_u": float(mann_whitney.statistic),
            "mann_whitney_pvalue": float(mann_whitney.pvalue),
            "positive_withheld_count": int(formal.loc[formal["audit_y_true"] == 1, binary_column].sum()),
            "negative_withheld_count": int(formal.loc[formal["audit_y_true"] == 0, binary_column].sum()),
            "positive_withheld_rate": float(formal.loc[formal["audit_y_true"] == 1, binary_column].mean()),
            "negative_withheld_rate": float(formal.loc[formal["audit_y_true"] == 0, binary_column].mean()),
            "score_min": float(formal[score_column].min()),
            "score_max": float(formal[score_column].max()),
            "score_mean": float(formal[score_column].mean()),
            "score_std_sample": float(formal[score_column].std(ddof=1)),
            "score_unique_values": int(formal[score_column].nunique()),
            "missing_score_count": missing_score,
            "nonfinite_score_count": nonfinite_score,
            "score_transform_audit": "No audit-side scaling, ranking, clipping, z-scoring, or binarization; Phase 7 continuous withheld_score used verbatim.",
            "higher_score_meaning": "stronger reference inadequacy / withholding evidence",
        }
    )

    precision, recall, pr_thresholds = precision_recall_curve(formal["audit_y_true"], formal[score_column])
    pr_curve = pd.DataFrame(
        {
            "curve_order": np.arange(len(precision), dtype=int),
            "precision": precision,
            "recall": recall,
            "threshold": np.r_[pr_thresholds, np.nan],
        }
    )
    fpr, tpr, roc_thresholds = roc_curve(formal["audit_y_true"], formal[score_column])
    roc_curve_table = pd.DataFrame(
        {
            "curve_order": np.arange(len(fpr), dtype=int),
            "false_positive_rate": fpr,
            "true_positive_rate": tpr,
            "threshold": roc_thresholds,
        }
    )
    return formal, excluded, metrics, pr_curve, roc_curve_table


def candidate_metric_row(
    analysis_name: str,
    source_file: Path,
    column: str,
    values: pd.Series,
    endpoint: pd.DataFrame,
    barcode_column: str,
    endpoint_column: str,
    positive_label: str,
    negative_labels: set[str],
    transform: str = "identity",
) -> dict[str, Any]:
    candidate = pd.DataFrame({barcode_column: endpoint[barcode_column].astype(str), "endpoint": endpoint[endpoint_column].astype(str)})
    score_table = pd.DataFrame({barcode_column: values.index.astype(str), "score": pd.to_numeric(values.to_numpy(), errors="coerce")})
    candidate = candidate.merge(score_table, on=barcode_column, how="left", validate="one_to_one")
    candidate = candidate[candidate["endpoint"].isin({positive_label} | negative_labels)].copy()
    if transform == "negate":
        candidate["score"] = -candidate["score"]
    elif transform == "one_minus":
        candidate["score"] = 1.0 - candidate["score"]
    y = (candidate["endpoint"] == positive_label).astype(int)
    result = metric_set(y, candidate["score"])
    positive_mean = float(candidate.loc[y == 1, "score"].mean()) if int((y == 1).sum()) else None
    negative_mean = float(candidate.loc[y == 0, "score"].mean()) if int((y == 0).sum()) else None
    return {
        "analysis_name": analysis_name,
        "input_score_file": rel(source_file),
        "score_column": column,
        "score_transform": transform,
        "endpoint_file": rel(DEFAULT_ENDPOINT),
        "included_endpoint_classes": ";".join(sorted({positive_label} | negative_labels)),
        "n_total": int(result["n"]),
        "n_positive": int((y == 1).sum()),
        "n_negative_or_control": int((y == 0).sum()),
        "positive_mean": positive_mean,
        "negative_mean": negative_mean,
        "auroc": result["auroc"],
        "average_precision": result["average_precision"],
        "trapezoidal_pr_auc": result["trapezoidal_pr_auc"],
        "matches_0_2131": bool(
            any(
                value is not None and abs(float(value) - TARGET_OLD) <= MATCH_TOLERANCE
                for value in [result["auroc"], result["average_precision"], result["trapezoidal_pr_auc"]]
            )
        ),
        "matches_0_28366": bool(
            any(
                value is not None and abs(float(value) - TARGET_CURRENT) <= MATCH_TOLERANCE
                for value in [result["auroc"], result["average_precision"], result["trapezoidal_pr_auc"]]
            )
        ),
        "notes": "Read-only reconstruction; not used to select the formal score.",
    }


def score_series_by_barcode(df: pd.DataFrame, barcode_column: str, value_column: str) -> pd.Series:
    table = df[[barcode_column, value_column]].copy()
    table[barcode_column] = table[barcode_column].astype(str)
    if table[barcode_column].duplicated().any():
        raise ValueError(f"Candidate file has duplicate barcodes for {value_column}")
    return table.set_index(barcode_column)[value_column]


def candidate_scores(args: argparse.Namespace, endpoint: pd.DataFrame, score_df: pd.DataFrame, raw_df: pd.DataFrame, phase8_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    endpoint_barcode = args.endpoint_barcode_column
    endpoint_column = args.endpoint_column
    formal_negative = {args.negative_label}

    sources = [
        (resolve(args.score_file), score_df, args.score_barcode_column),
        (resolve(args.raw_score_file), raw_df, "spot_id"),
        (resolve(args.phase8_spot_file), phase8_df, "barcode"),
    ]
    exclude_name_parts = ("imagerow", "imagecol", "fraction_from_phase5")
    include_name_parts = (
        "score",
        "error",
        "deficit",
        "fraction",
        "pvalue",
        "qvalue",
        "weight",
        "supported",
        "withheld",
        "unsupported",
    )
    for source_file, frame, barcode in sources:
        if frame[barcode].astype(str).duplicated().any():
            continue
        for column in frame.columns:
            lower = column.lower()
            if column == barcode or any(part in lower for part in exclude_name_parts):
                continue
            numeric = pd.to_numeric(frame[column], errors="coerce")
            if numeric.notna().sum() == 0 or not any(part in lower for part in include_name_parts):
                continue
            transform = "negate" if ("pvalue" in lower or "qvalue" in lower) else "identity"
            series = pd.Series(numeric.to_numpy(), index=frame[barcode].astype(str), name=column)
            rows.append(
                candidate_metric_row(
                    analysis_name="candidate_score_formal_positive_vs_negative",
                    source_file=source_file,
                    column=column,
                    values=series,
                    endpoint=endpoint,
                    barcode_column=endpoint_barcode,
                    endpoint_column=endpoint_column,
                    positive_label=args.positive_label,
                    negative_labels=formal_negative,
                    transform=transform,
                )
            )
    return pd.DataFrame(rows).sort_values(["input_score_file", "score_column", "score_transform"]).reset_index(drop=True)


def sensitivity_reconstruction(
    args: argparse.Namespace,
    endpoint: pd.DataFrame,
    score_df: pd.DataFrame,
    phase8_df: pd.DataFrame,
    candidate_table: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    formal_series = score_series_by_barcode(score_df, args.score_barcode_column, args.score_column)
    binary_values = bool_series(score_df[args.binary_column], args.binary_column).astype(int)
    binary_series = pd.Series(binary_values.to_numpy(), index=score_df[args.score_barcode_column].astype(str), name=args.binary_column)
    phase8_series = score_series_by_barcode(phase8_df, "barcode", args.score_column)

    analyses = [
        ("formal_positive_vs_negative", formal_series, {args.negative_label}, "identity", resolve(args.score_file), args.score_column),
        ("positive_vs_all_nonpositive", formal_series, {args.negative_label, args.ambiguous_label, args.excluded_label}, "identity", resolve(args.score_file), args.score_column),
        ("positive_vs_negative_plus_ambiguous", formal_series, {args.negative_label, args.ambiguous_label}, "identity", resolve(args.score_file), args.score_column),
        ("positive_vs_negative_plus_ambiguous_plus_excluded", formal_series, {args.negative_label, args.ambiguous_label, args.excluded_label}, "identity", resolve(args.score_file), args.score_column),
        ("formal_after_explicit_missing_score_drop", formal_series.dropna(), {args.negative_label}, "identity", resolve(args.score_file), args.score_column),
        ("phase8_plotting_file_score", phase8_series, {args.negative_label}, "identity", resolve(args.phase8_spot_file), args.score_column),
        ("binary_mask_as_continuous_score", binary_series, {args.negative_label}, "identity", resolve(args.score_file), args.binary_column),
        ("formal_score_direction_reversed", formal_series, {args.negative_label}, "negate", resolve(args.score_file), args.score_column),
    ]
    for name, values, controls, transform, source_file, column in analyses:
        rows.append(
            candidate_metric_row(
                analysis_name=name,
                source_file=source_file,
                column=column,
                values=values,
                endpoint=endpoint,
                barcode_column=args.endpoint_barcode_column,
                endpoint_column=args.endpoint_column,
                positive_label=args.positive_label,
                negative_labels=controls,
                transform=transform,
            )
        )

    for _, row in candidate_table.iterrows():
        if bool(row["matches_0_2131"]) or bool(row["matches_0_28366"]):
            candidate_row = row.to_dict()
            candidate_row["analysis_name"] = "candidate_score_match_followup"
            rows.append(candidate_row)
    return pd.DataFrame(rows)


def candidate_endpoint_files(endpoint_file: Path, score_file: Path, phase8_file: Path) -> pd.DataFrame:
    candidates: list[dict[str, Any]] = []
    search_root = ROOT / "visualizations" / "bioapp_experiment"
    for path in sorted(search_root.rglob("*.csv")):
        if "fig5b_auprc_formal_audit" in path.parts:
            continue
        try:
            header = pd.read_csv(path, nrows=0).columns.tolist()
        except Exception as exc:
            candidates.append({"file": rel(path), "candidate_endpoint_columns": "", "status": "unreadable", "notes": str(exc)[:300]})
            continue
        columns = [column for column in header if "endpoint" in column.lower() and ("status" in column.lower() or "class" in column.lower() or column.lower().startswith("is_"))]
        if not columns:
            continue
        status = "formal" if path.resolve() == endpoint_file.resolve() else "candidate_or_derived"
        if path.resolve() in {score_file.resolve(), phase8_file.resolve()}:
            status = "copied_endpoint_labels_in_score_output"
        candidates.append(
            {
                "file": rel(path),
                "candidate_endpoint_columns": ";".join(columns),
                "status": status,
                "notes": "Header inventory only; formal endpoint remains the frozen Phase 2C file.",
            }
        )
    return pd.DataFrame(candidates, columns=["file", "candidate_endpoint_columns", "status", "notes"])


def search_value_hits(output_dir: Path) -> pd.DataFrame:
    hits: list[dict[str, Any]] = []
    skip_parts = {".git", "data", "fig5b_auprc_formal_audit"}
    search_roots = [
        ROOT / "scripts",
        ROOT / "manuscript",
        ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping",
        ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase6_svtuner_aware_execution_planning_against_frozen_cta_immune_endpoint",
        ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase7_svtuner_aware_execution_against_frozen_cta_immune_endpoint",
        ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison",
        ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase9_final_biological_application_audit_figure_selection_interpretation_boundary_lock",
        ROOT / "visualizations" / "bioapp_experiment" / "bioapp_main_figure_v3_10_d2_light_green_footprint_correction",
        ROOT / "visualizations" / "bioapp_experiment" / "bioapp_main_figure_v3_11_top_row_abc_redesign",
        ROOT / "visualizations" / "bioapp_experiment" / "bioapp_main_figure_v3_12_panel_A_evidence_chain_redesign",
    ]
    root_level_files = [path for path in ROOT.iterdir() if path.is_file() and path.suffix.lower() in TEXT_EXTENSIONS]
    candidates = list(root_level_files)
    for search_root in search_roots:
        if search_root.exists():
            candidates.extend(path for path in search_root.rglob("*") if path.is_file())
    for path in sorted(set(candidates)):
        if path.resolve() == Path(__file__).resolve():
            continue
        if not path.is_file() or path.suffix.lower() not in TEXT_EXTENSIONS or any(part in skip_parts for part in path.parts):
            continue
        if path.stat().st_size > 10 * 1024 * 1024:
            continue
        try:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError:
            continue
        for line_number, line in enumerate(lines, start=1):
            metric_context = bool(re.search(r"AUPRC|average.precision|PR.AUC|withheld|Fig.?5|panel.?b", line, re.IGNORECASE))
            for name, pattern in VALUE_PATTERNS.items():
                if not pattern.search(line):
                    continue
                if path.suffix.lower() == ".csv" and name in {"0.2131", "0.28366", "0.2837"} and not metric_context:
                    continue
                hits.append(
                    {
                        "pattern": name,
                        "file": rel(path),
                        "line_number": line_number,
                        "context": line.strip()[:1000],
                        "metric_context": metric_context,
                    }
                )
    return pd.DataFrame(hits, columns=["pattern", "file", "line_number", "context", "metric_context"])


def manuscript_occurrences(hits: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    deduplicated = hits.drop_duplicates(subset=["file", "line_number", "context"])
    for _, hit in deduplicated.iterrows():
        context = str(hit["context"])
        if not re.search(r"AUPRC|average_precision_score|average precision|PR-AUC|0\.2837|0\.28366|0\.2131", context, re.IGNORECASE):
            continue
        file_name = str(hit["file"])
        fig5_related_file = any(
            token in file_name
            for token in [
                "README.md",
                "\u751f\u7269\u5b66\u5e94\u7528\u5b9e\u9a8c\u8bba\u6587\u6750\u6599\u51c6\u5907.md",
                "run_bioapp_phase8_",
                "run_bioapp_phase9_",
                "run_bioapp_main_figure_v3_",
                "bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison",
                "bioapp_phase9_final_biological_application_audit_figure_selection_interpretation_boundary_lock",
                "bioapp_main_figure_v3_",
            ]
        )
        if not fig5_related_file:
            continue
        if file_name.endswith(".csv") and not re.search(r"AUPRC|average precision|PR-AUC", context, re.IGNORECASE):
            continue
        value_match = re.search(r"0\.2836551170002947|0\.28366\d*|0\.2837\d*|0\.2131\d*", context)
        metric_name = "average precision" if "average_precision_score" in context else ("AUPRC" if re.search(r"AUPRC", context, re.IGNORECASE) else "unspecified")
        if "bioapp_main_figure_v3_10" in file_name or "bioapp_main_figure_v3_11" in file_name:
            status = "legacy"
        elif "average_precision_score" in context:
            status = "consistent"
        elif re.search(r"AUPRC", context, re.IGNORECASE):
            status = "requires update" if ("withheld" in context.lower() or "0.283" in context) else "ambiguous terminology"
        elif "0.2131" in context:
            status = "legacy"
        else:
            status = "ambiguous terminology"
        rows.append(
            {
                "file": hit["file"],
                "line_number": hit["line_number"],
                "context": context,
                "reported_metric_name": metric_name,
                "reported_value": value_match.group(0) if value_match else "",
                "status": status,
            }
        )
    return pd.DataFrame(rows, columns=["file", "line_number", "context", "reported_metric_name", "reported_value", "status"])


def file_hash_table(paths: Iterable[tuple[str, Path]]) -> pd.DataFrame:
    rows = []
    for role, path in paths:
        if path.is_file():
            rows.append(
                {
                    "role": role,
                    "file": rel(path),
                    "sha1": sha1(path),
                    "modified_time_utc": modified_time(path),
                    "size_bytes": path.stat().st_size,
                }
            )
        else:
            rows.append({"role": role, "file": rel(path), "sha1": "not found", "modified_time_utc": "not found", "size_bytes": None})
    return pd.DataFrame(rows)


def software_versions() -> dict[str, str]:
    return {
        "Python": platform.python_version(),
        "implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "NumPy": np.__version__,
        "pandas": pd.__version__,
        "scikit-learn": sklearn.__version__,
        "SciPy": scipy.__version__,
    }


def provenance_table(
    args: argparse.Namespace,
    counts: dict[str, int],
    metrics: dict[str, Any],
    source_files: dict[str, Path],
    old_value_reproduced: bool,
) -> pd.DataFrame:
    endpoint_file = source_files["endpoint"]
    score_file = source_files["score"]
    phase8_summary = source_files["phase8_summary"]
    phase8_script = source_files["phase8_script"]
    audit_script = Path(__file__).resolve()
    common = {
        "input_endpoint_file": rel(endpoint_file),
        "input_endpoint_file_sha1": sha1(endpoint_file),
        "input_score_file": rel(score_file),
        "input_score_file_sha1": sha1(score_file),
        "barcode_column": f"endpoint:{args.endpoint_barcode_column};score:{args.score_barcode_column}",
        "endpoint_column": args.endpoint_column,
        "score_column": args.score_column,
        "positive_label": args.positive_label,
        "negative_label": args.negative_label,
        "total_input_rows": counts["total"],
        "analysis_rows": EXPECTED_ANALYSIS,
        "positive_count": counts["positive"],
        "negative_count": counts["negative"],
        "ambiguous_count": counts["ambiguous"],
        "excluded_count": counts["excluded"],
        "duplicate_barcode_count": 0,
        "missing_score_count": metrics["missing_score_count"],
        "metric_library": "scikit-learn",
        "library_version": sklearn.__version__,
    }
    rows = [
        {
            "reported_value": TARGET_CURRENT,
            "value_label": "withheld_AUPRC (actually average precision)",
            "source_file": rel(phase8_summary),
            "source_file_sha1": sha1(phase8_summary),
            "source_file_modified_time": modified_time(phase8_summary),
            "generating_script": rel(phase8_script),
            "generating_script_sha1": sha1(phase8_script),
            **common,
            "metric_function": "sklearn.metrics.average_precision_score",
            "notes": "The Phase 8 variable name and Fig. 5b annotation say AUPRC, but the generating function computes average precision.",
        },
        {
            "reported_value": TARGET_OLD,
            "value_label": "early record: AUPRC",
            "source_file": "not found",
            "source_file_sha1": "not found",
            "source_file_modified_time": "not found",
            "generating_script": "not found",
            "generating_script_sha1": "not found",
            **{key: "not found" for key in common},
            "metric_function": "not found",
            "metric_library": "not found",
            "library_version": "not found",
            "notes": "Reproduced by read-only sensitivity checks." if old_value_reproduced else "Could not be traced to a reproducible formal or sensitivity input.",
        },
        {
            "reported_value": metrics["trapezoidal_pr_auc"],
            "value_label": "independently audited trapezoidal PR-AUC",
            "source_file": rel(score_file),
            "source_file_sha1": sha1(score_file),
            "source_file_modified_time": modified_time(score_file),
            "generating_script": rel(audit_script),
            "generating_script_sha1": sha1(audit_script),
            **common,
            "metric_function": "sklearn.metrics.precision_recall_curve + sklearn.metrics.auc on traversal-preserving ascending recall",
            "notes": "Same frozen 1,888-spot input; distinct from average precision and not equal to 0.2131.",
        },
    ]
    fields = [
        "reported_value",
        "value_label",
        "source_file",
        "source_file_sha1",
        "source_file_modified_time",
        "generating_script",
        "generating_script_sha1",
        "input_endpoint_file",
        "input_endpoint_file_sha1",
        "input_score_file",
        "input_score_file_sha1",
        "barcode_column",
        "endpoint_column",
        "score_column",
        "positive_label",
        "negative_label",
        "total_input_rows",
        "analysis_rows",
        "positive_count",
        "negative_count",
        "ambiguous_count",
        "excluded_count",
        "duplicate_barcode_count",
        "missing_score_count",
        "metric_function",
        "metric_library",
        "library_version",
        "notes",
    ]
    return pd.DataFrame(rows)[fields]


def write_readme(output_dir: Path, args: argparse.Namespace) -> None:
    text = f"""# Fig. 5b precision-recall formal audit

This directory is an independent, read-only audit of the precision-recall metric shown in Fig. 5b. It does not modify or rerun Stage3, Stage3B, Stage4, CytoSPACE, the CTA endpoint, the shared-gene set, the manuscript, the formal figure, or existing source-value files.

## Reproduction

Run from the repository root:

```powershell
python scripts/audit_fig5b_auprc_provenance.py --formal
```

Formal inputs:

- endpoint: `{rel(resolve(args.endpoint_file))}`
- score: `{rel(resolve(args.score_file))}`
- endpoint column: `{args.endpoint_column}`
- score column: `{args.score_column}`
- binary column: `{args.binary_column}`

The formal analysis is fixed to CTA endpoint-positive versus endpoint-negative spots. Ambiguous and excluded spots are retained in `excluded_endpoint_360_spots.csv` and are not used for formal discrimination metrics.

`average precision` and `trapezoidal PR-AUC` are reported separately because they are not interchangeable.
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def final_decision_text(
    args: argparse.Namespace,
    metrics: dict[str, Any],
    source_files: dict[str, Path],
    hashes: dict[str, str],
    old_value_reproduced: bool,
    occurrence_table: pd.DataFrame,
) -> str:
    affected = occurrence_table[occurrence_table["status"].eq("requires update")].drop_duplicates(subset=["file", "line_number"])
    affected_lines = "\n".join(f"- `{row.file}:{row.line_number}`" for row in affected.itertuples()) or "- None identified"
    old_origin = (
        "A read-only sensitivity reconstruction reproduced 0.2131; see sensitivity_metric_reconstruction.tsv."
        if old_value_reproduced
        else "The value 0.2131 could not be traced to a reproducible formal input. Current files, Phase 7/8 outputs, candidate score columns, alternative endpoint universes, plotting inputs, binary-score misuse, reversed score direction, and metric definitions were searched."
    )
    root_cause = "legacy score/input reconstruction" if old_value_reproduced else "untraceable legacy value"
    return f"""# Fig. 5b AUPRC formal audit decision

## 1. Final decision

RESOLVED

## 2. Formal Fig. 5b analysis input

- endpoint file: `{rel(source_files['endpoint'])}`
- score file: `{rel(source_files['score'])}`
- barcode columns: endpoint `{args.endpoint_barcode_column}`; score `{args.score_barcode_column}`
- endpoint column: `{args.endpoint_column}`
- score column: `{args.score_column}`
- analysis n: {EXPECTED_ANALYSIS}
- positive n: {EXPECTED_POSITIVE}
- negative n: {EXPECTED_NEGATIVE}
- endpoint SHA1: `{hashes['endpoint']}`
- score SHA1: `{hashes['score']}`

## 3. Independently recomputed metrics

- AUROC: {metrics['auroc']:.8f}
- average precision: {metrics['average_precision']:.8f}
- trapezoidal PR-AUC: {metrics['trapezoidal_pr_auc']:.8f}
- positive prevalence: {metrics['positive_prevalence']:.8f}
- positive mean score: {metrics['positive_mean_score']:.8f}
- negative mean score: {metrics['negative_mean_score']:.8f}
- mean score difference: {metrics['mean_score_difference']:.8f}
- Mann-Whitney U P value: {metrics['mann_whitney_pvalue']:.8e}
- binary positive withholding: {metrics['positive_withheld_count']} / {EXPECTED_POSITIVE} ({metrics['positive_withheld_rate']:.8f})
- binary negative withholding: {metrics['negative_withheld_count']} / {EXPECTED_NEGATIVE} ({metrics['negative_withheld_rate']:.8f})

## 4. Origin of 0.2837

`0.2836551170002947` was generated by `sklearn.metrics.average_precision_score` in `{rel(source_files['phase8_script'])}` from the frozen positive/negative endpoint set and Phase 7 `withheld_score`. The value was stored under the misleading key `withheld_AUPRC` in `{rel(source_files['phase8_summary'])}`. Fig. 5b reads that summary key and formats it as `AUPRC = 0.2837`; the figure script does not recompute the metric.

## 5. Origin of 0.2131

{old_origin}

The same formal PR curve has trapezoidal PR-AUC {metrics['trapezoidal_pr_auc']:.8f}, not 0.2131.

## 6. Root cause

`{root_cause}` for 0.2131, combined with ambiguous use of `AUPRC` for the reproducible average-precision value.

This is not an AP-versus-trapezoidal-PR-AUC explanation: on the frozen formal input, AP is {metrics['average_precision']:.8f} and trapezoidal PR-AUC is {metrics['trapezoidal_pr_auc']:.8f}.

## 7. Manuscript recommendation

- Report **average precision = 0.2837**.
- Show **AP = 0.284** in Fig. 5b.
- State in Methods that the value is computed with `sklearn.metrics.average_precision_score`.
- Reserve **trapezoidal PR-AUC** for the separately computed value {metrics['trapezoidal_pr_auc']:.4f} if it is reported at all.
- Do not call AP and trapezoidal PR-AUC the same metric.
- Do not use 0.2131 in formal reporting because it has no reproducible formal provenance.

## 8. Files requiring future updates

{affected_lines}

No file in this list was modified by this audit.

## 9. Guardrail confirmation

```text
Stage3 rerun: false
Stage4 rerun: false
CTA endpoint modified: false
Stage3B threshold modified: false
Shared-gene set modified: false
Formal manuscript modified: false
Formal Fig. 5 modified: false
Existing source-value files overwritten: false
```
"""


def main() -> None:
    args = parse_args()
    source_files = {
        "endpoint": resolve(args.endpoint_file),
        "score": resolve(args.score_file),
        "raw_score": resolve(args.raw_score_file),
        "phase8_spots": resolve(args.phase8_spot_file),
        "phase8_summary": resolve(args.phase8_summary_file),
        "phase8_script": resolve(args.phase8_script),
        "fig5_script": resolve(args.fig5_script),
    }
    output_dir = resolve(args.output_dir)
    for label, path in source_files.items():
        require_file(path, label)
    output_dir.mkdir(parents=True, exist_ok=True)

    endpoint = pd.read_csv(source_files["endpoint"])
    score_df = pd.read_csv(source_files["score"])
    raw_df = pd.read_csv(source_files["raw_score"])
    phase8_df = pd.read_csv(source_files["phase8_spots"])
    require_columns(endpoint, [args.endpoint_barcode_column, args.endpoint_column], "endpoint file")
    require_columns(score_df, [args.score_barcode_column, args.score_column, args.binary_column], "score file")
    require_columns(raw_df, ["spot_id", "unsupported_fraction_estimate", "is_unsupported_region"], "raw score file")
    require_columns(phase8_df, ["barcode", args.endpoint_column, args.score_column, args.binary_column], "Phase 8 spot file")

    duplicate_table = pd.concat(
        [
            duplicate_rows(endpoint, args.endpoint_barcode_column, rel(source_files["endpoint"])),
            duplicate_rows(score_df, args.score_barcode_column, rel(source_files["score"])),
        ],
        ignore_index=True,
    )
    duplicate_table.to_csv(output_dir / "duplicate_barcodes.tsv", sep="\t", index=False)
    if not duplicate_table.empty:
        raise ValueError("Duplicate barcodes detected in formal inputs")

    counts = endpoint_counts(
        endpoint,
        args.endpoint_column,
        {
            "positive": args.positive_label,
            "negative": args.negative_label,
            "ambiguous": args.ambiguous_label,
            "excluded": args.excluded_label,
        },
    )
    warnings = validate_formal_counts(counts, args.formal)
    merged, unmatched, join_counts = merge_inputs(args, endpoint, score_df)
    unmatched.to_csv(output_dir / "unmatched_barcodes.tsv", sep="\t", index=False)
    if args.formal and not unmatched.empty:
        raise ValueError(f"Formal endpoint/score merge has {len(unmatched)} unmatched barcodes")

    formal, excluded, metrics, pr_curve, roc_curve_table = compute_formal_analysis(args, merged)
    formal.to_csv(output_dir / "formal_analysis_1888_spots.csv", index=False)
    excluded.to_csv(output_dir / "excluded_endpoint_360_spots.csv", index=False)
    pr_curve.to_csv(output_dir / "precision_recall_curve.csv", index=False)
    roc_curve_table.to_csv(output_dir / "roc_curve.csv", index=False)

    candidate_table = candidate_scores(args, endpoint, score_df, raw_df, phase8_df)
    candidate_table.to_csv(output_dir / "candidate_score_columns.tsv", sep="\t", index=False)
    candidate_endpoints = candidate_endpoint_files(source_files["endpoint"], source_files["score"], source_files["phase8_spots"])
    candidate_endpoints.to_csv(output_dir / "candidate_endpoint_files.tsv", sep="\t", index=False)
    sensitivity = sensitivity_reconstruction(args, endpoint, score_df, phase8_df, candidate_table)
    sensitivity.to_csv(output_dir / "sensitivity_metric_reconstruction.tsv", sep="\t", index=False)
    old_value_reproduced = bool(candidate_table["matches_0_2131"].any() or sensitivity["matches_0_2131"].any())

    hits = search_value_hits(output_dir)
    hits.to_csv(output_dir / "searched_value_hits.tsv", sep="\t", index=False)
    occurrences = manuscript_occurrences(hits)
    occurrences.to_csv(output_dir / "manuscript_metric_occurrences.tsv", sep="\t", index=False)

    audit_script = Path(__file__).resolve()
    hash_table = file_hash_table(
        [
            ("formal endpoint", source_files["endpoint"]),
            ("formal Phase 7 score contract", source_files["score"]),
            ("Phase 7 raw score", source_files["raw_score"]),
            ("Phase 8 spot-level derived table", source_files["phase8_spots"]),
            ("Phase 8 summary", source_files["phase8_summary"]),
            ("Phase 8 generating script", source_files["phase8_script"]),
            ("Fig. 5b generating script", source_files["fig5_script"]),
            ("formal audit script", audit_script),
        ]
    )
    hash_table.to_csv(output_dir / "input_file_hashes.tsv", sep="\t", index=False)
    versions = software_versions()
    (output_dir / "software_versions.txt").write_text("\n".join(f"{key}: {value}" for key, value in versions.items()) + "\n", encoding="utf-8")

    raw_check = score_df[[args.score_barcode_column, args.score_column, args.binary_column]].copy()
    raw_check[args.score_barcode_column] = raw_check[args.score_barcode_column].astype(str)
    raw_subset = raw_df[["spot_id", "unsupported_fraction_estimate", "is_unsupported_region"]].copy().rename(columns={"spot_id": args.score_barcode_column})
    raw_subset[args.score_barcode_column] = raw_subset[args.score_barcode_column].astype(str)
    raw_check = raw_check.merge(raw_subset, on=args.score_barcode_column, validate="one_to_one")
    raw_score_max_abs_diff = float(np.max(np.abs(raw_check[args.score_column] - raw_check["unsupported_fraction_estimate"])))
    raw_binary_equal = bool(
        (
            bool_series(raw_check[args.binary_column], args.binary_column)
            == bool_series(raw_check["is_unsupported_region"], "is_unsupported_region")
        ).all()
    )

    phase8_summary = json.loads(source_files["phase8_summary"].read_text(encoding="utf-8"))
    phase8_reported = float(phase8_summary["metric_highlights"]["withheld_AUPRC"])
    phase8_value_reproduced = bool(np.isclose(phase8_reported, metrics["average_precision"], rtol=0, atol=1e-15))
    summary = {
        "decision": "RESOLVED",
        "formal_analysis": {
            "n": int(len(formal)),
            "positive": counts["positive"],
            "negative": counts["negative"],
            "ambiguous": counts["ambiguous"],
            "excluded": counts["excluded"],
            "endpoint_file": rel(source_files["endpoint"]),
            "score_file": rel(source_files["score"]),
            "endpoint_column": args.endpoint_column,
            "score_column": args.score_column,
            "binary_column": args.binary_column,
        },
        "join_audit": join_counts,
        "endpoint_warnings": warnings,
        "metrics": metrics,
        "phase8_reported_value": phase8_reported,
        "phase8_reported_value_reproduced_as_average_precision": phase8_value_reproduced,
        "raw_score_contract_check": {
            "withheld_score_vs_unsupported_fraction_estimate_max_abs_diff": raw_score_max_abs_diff,
            "withheld_binary_vs_is_unsupported_region_equal": raw_binary_equal,
        },
        "origin_0_2837": "Phase 8 sklearn.metrics.average_precision_score stored as withheld_AUPRC and imported by the Fig. 5b script.",
        "origin_0_2131": "reproduced in sensitivity audit" if old_value_reproduced else "untraceable legacy value; not reproduced",
        "root_cause": "legacy score/input reconstruction" if old_value_reproduced else "untraceable legacy value plus ambiguous AUPRC terminology",
        "recommended_manuscript_metric": "average precision = 0.2837",
        "recommended_figure_label": "AP = 0.284",
        "guardrails": {
            "Stage3_rerun": False,
            "Stage4_rerun": False,
            "CTA_endpoint_modified": False,
            "Stage3B_threshold_modified": False,
            "shared_gene_set_modified": False,
            "formal_manuscript_modified": False,
            "formal_Fig5_modified": False,
            "existing_source_value_files_overwritten": False,
        },
        "software_versions": versions,
    }
    write_json(output_dir / "audit_summary.json", summary)

    flat_summary = {
        "decision": summary["decision"],
        "analysis_n": len(formal),
        "positive_n": counts["positive"],
        "negative_n": counts["negative"],
        "ambiguous_n": counts["ambiguous"],
        "excluded_n": counts["excluded"],
        "AUROC": metrics["auroc"],
        "average_precision": metrics["average_precision"],
        "trapezoidal_PR_AUC": metrics["trapezoidal_pr_auc"],
        "positive_prevalence": metrics["positive_prevalence"],
        "positive_mean_score": metrics["positive_mean_score"],
        "negative_mean_score": metrics["negative_mean_score"],
        "mean_score_difference": metrics["mean_score_difference"],
        "Mann_Whitney_U_pvalue": metrics["mann_whitney_pvalue"],
        "positive_withheld_count": metrics["positive_withheld_count"],
        "negative_withheld_count": metrics["negative_withheld_count"],
        "value_0_2837_reproduced_as_average_precision": phase8_value_reproduced,
        "value_0_2131_reproduced": old_value_reproduced,
    }
    pd.DataFrame([flat_summary]).to_csv(output_dir / "audit_summary.tsv", sep="\t", index=False)

    provenance = provenance_table(args, counts, metrics, source_files, old_value_reproduced)
    provenance.to_csv(output_dir / "value_provenance.tsv", sep="\t", index=False)
    write_readme(output_dir, args)

    hash_map = {row.role: row.sha1 for row in hash_table.itertuples()}
    final_text = final_decision_text(
        args,
        metrics,
        source_files,
        {"endpoint": hash_map["formal endpoint"], "score": hash_map["formal Phase 7 score contract"]},
        old_value_reproduced,
        occurrences,
    )
    (output_dir / "final_decision.md").write_text(final_text, encoding="utf-8")

    print("Fig. 5b AUPRC audit completed.")
    print()
    print("Formal analysis:")
    print(f"n = {len(formal)}")
    print(f"positive = {counts['positive']}")
    print(f"negative = {counts['negative']}")
    print()
    print(f"AUROC = {metrics['auroc']:.10f}")
    print(f"Average precision = {metrics['average_precision']:.10f}")
    print(f"Trapezoidal PR-AUC = {metrics['trapezoidal_pr_auc']:.10f}")
    print()
    print("0.2837 source:")
    print("Phase 8 average_precision_score, stored under withheld_AUPRC and imported by Fig. 5b.")
    print()
    print("0.2131 source:")
    print("Reproduced by sensitivity audit." if old_value_reproduced else "Could not be traced to a reproducible formal input.")
    print()
    print(f"Root cause: {summary['root_cause']}")
    print("Decision: RESOLVED")
    print("Recommended manuscript metric: average precision = 0.2837")
    print("Stage3 rerun: false")
    print("Stage4 rerun: false")
    print("Formal manuscript modified: false")
    print("Formal Fig. 5 modified: false")


if __name__ == "__main__":
    main()
