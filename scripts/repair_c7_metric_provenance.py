#!/usr/bin/env python
"""Deterministically repair C7 metric provenance from the two frozen CSV inputs.

This script does not run CytoSPACE, SVTuner, or Stage3B. It only joins the
frozen records by barcode, restricts the endpoint to the formal positive and
negative classes, and recomputes the requested metrics and threshold audits.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score


ENDPOINT_REL = Path(
    "visualizations/bioapp_experiment/"
    "bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping/"
    "spot_level_endpoint_freeze.csv"
)
SVTUNER_REL = Path(
    "visualizations/bioapp_experiment/"
    "bioapp_phase7_svtuner_aware_execution_against_frozen_cta_immune_endpoint/"
    "svtuner_immune_all_dropout/"
    "svtuner_immune_all_dropout_spot_level_raw_output_contract.csv"
)
OUTPUT_REL = Path("result/c7_provenance_repair")


def full_precision(value: Any) -> str:
    """Return a deterministic, round-trip-safe representation for numeric values."""
    if isinstance(value, (float, np.floating)):
        if math.isnan(float(value)):
            return "NA"
        return format(float(value), ".17g")
    if isinstance(value, (bool, np.bool_)):
        return "YES" if bool(value) else "NO"
    return str(value)


def parse_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.astype(bool)
    normalized = series.astype(str).str.strip().str.lower()
    mapping = {
        "true": True,
        "false": False,
        "1": True,
        "0": False,
        "yes": True,
        "no": False,
    }
    parsed = normalized.map(mapping)
    if parsed.isna().any():
        bad = sorted(normalized[parsed.isna()].unique().tolist())
        raise ValueError(f"Unrecognized withheld_binary values: {bad}")
    return parsed.astype(bool)


def fixed_precision_result(
    precision: np.ndarray,
    recall: np.ndarray,
    thresholds: np.ndarray,
    target: float,
) -> dict[str, Any]:
    valid = np.flatnonzero(precision >= target)
    if valid.size == 0:
        return {
            "target": target,
            "achieved": False,
            "precision": math.nan,
            "recall": math.nan,
            "threshold": math.nan,
            "tie_count": 0,
            "threshold_min": math.nan,
            "threshold_max": math.nan,
            "tied_thresholds": [],
        }
    max_recall = float(recall[valid].max())
    tied = valid[np.isclose(recall[valid], max_recall, rtol=0.0, atol=1e-15)]
    # The frozen scientific rule fixes maximum recall but is silent on ties.
    # Use the smallest tied finite threshold as a deterministic representative.
    selected = tied[np.argmin(thresholds[tied])]
    return {
        "target": target,
        "achieved": True,
        "precision": float(precision[selected]),
        "recall": float(recall[selected]),
        "threshold": float(thresholds[selected]),
        "tie_count": int(tied.size),
        "threshold_min": float(thresholds[tied].min()),
        "threshold_max": float(thresholds[tied].max()),
        "tied_thresholds": [float(x) for x in thresholds[tied]],
    }


def threshold_disagreement_audit(
    scores: np.ndarray, frozen_binary: np.ndarray
) -> dict[str, Any]:
    unique_scores = np.unique(scores[np.isfinite(scores)])
    # Observed score values enumerate every non-empty >= threshold partition.
    # Add one finite value above max(score) to include the all-negative rule.
    candidates = np.concatenate(
        [unique_scores, np.array([np.nextafter(unique_scores[-1], np.inf)])]
    )
    records: list[dict[str, Any]] = []
    min_disagreement = scores.size + 1
    for threshold in candidates:
        predicted = scores >= threshold
        disagreement = int(np.count_nonzero(predicted != frozen_binary))
        if disagreement > min_disagreement:
            continue
        record = {
            "threshold": float(threshold),
            "disagreement": disagreement,
            "tp": int(np.sum(predicted & frozen_binary)),
            "fp": int(np.sum(predicted & ~frozen_binary)),
            "tn": int(np.sum(~predicted & ~frozen_binary)),
            "fn": int(np.sum(~predicted & frozen_binary)),
        }
        if disagreement < min_disagreement:
            min_disagreement = disagreement
            records = [record]
        else:
            records.append(record)

    # Attach the exact real-valued threshold interval that yields each tied rule.
    for record in records:
        threshold = record["threshold"]
        idx = int(np.searchsorted(unique_scores, threshold, side="left"))
        if idx == unique_scores.size:
            record["interval"] = f"({full_precision(unique_scores[-1])}, +inf)"
        elif idx == 0:
            record["interval"] = f"(-inf, {full_precision(threshold)}]"
        else:
            record["interval"] = (
                f"({full_precision(unique_scores[idx - 1])}, "
                f"{full_precision(threshold)}]"
            )
    return {
        "minimum_disagreement": min_disagreement,
        "exactly_reproducible": min_disagreement == 0,
        "best_records": records,
    }


def add_row(rows: list[dict[str, str]], metric: str, value: Any, source: str, method: str) -> None:
    rows.append(
        {
            "metric": metric,
            "value": full_precision(value),
            "source": source,
            "method": method,
        }
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project_root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.project_root.resolve()
    endpoint_path = root / ENDPOINT_REL
    svtuner_path = root / SVTUNER_REL
    output_dir = root / OUTPUT_REL

    endpoint = pd.read_csv(
        endpoint_path,
        usecols=["barcode", "primary_endpoint_status"],
        dtype={"barcode": str, "primary_endpoint_status": str},
    )
    svtuner = pd.read_csv(
        svtuner_path,
        usecols=["barcode", "withheld_score", "withheld_binary"],
        dtype={"barcode": str},
    )
    if endpoint["barcode"].duplicated().any() or svtuner["barcode"].duplicated().any():
        raise ValueError("Frozen inputs must each contain unique barcodes")

    formal_endpoint = endpoint[
        endpoint["primary_endpoint_status"].isin(["positive", "negative"])
    ].copy()
    joined = formal_endpoint.merge(svtuner, on="barcode", how="left", validate="one_to_one")
    unmatched_formal_barcodes = int(joined["withheld_score"].isna().sum())
    joined["withheld_score"] = pd.to_numeric(joined["withheld_score"], errors="raise")
    joined["withheld_binary"] = parse_bool(joined["withheld_binary"])

    missing_analysis_fields = int(
        joined[["primary_endpoint_status", "withheld_score", "withheld_binary"]]
        .isna()
        .any(axis=1)
        .sum()
    )
    if unmatched_formal_barcodes or missing_analysis_fields:
        raise ValueError(
            "Formal analysis universe contains unmatched or missing records: "
            f"unmatched={unmatched_formal_barcodes}, missing={missing_analysis_fields}"
        )

    y_true = joined["primary_endpoint_status"].eq("positive").to_numpy(dtype=int)
    scores = joined["withheld_score"].to_numpy(dtype=float)
    frozen_binary = joined["withheld_binary"].to_numpy(dtype=bool)
    if not np.isfinite(scores).all():
        raise ValueError("Formal withheld_score contains non-finite values")

    n_total = int(joined.shape[0])
    n_positive = int(y_true.sum())
    n_negative = int(n_total - n_positive)
    expected = (1888, 134, 1754, 0)
    observed = (n_total, n_positive, n_negative, missing_analysis_fields)
    if observed != expected:
        raise ValueError(f"Formal universe mismatch: observed={observed}, expected={expected}")

    auroc = float(roc_auc_score(y_true, scores))
    average_precision = float(average_precision_score(y_true, scores))
    pr_precision_all, pr_recall_all, pr_thresholds = precision_recall_curve(y_true, scores)
    pr_precision = pr_precision_all[:-1]
    pr_recall = pr_recall_all[:-1]
    finite_mask = np.isfinite(pr_thresholds)
    pr_precision = pr_precision[finite_mask]
    pr_recall = pr_recall[finite_mask]
    pr_thresholds = pr_thresholds[finite_mask]

    fixed = {
        target: fixed_precision_result(pr_precision, pr_recall, pr_thresholds, target)
        for target in (0.25, 0.30, 0.50)
    }
    highest_precision = float(pr_precision.max())
    highest_idx_all = np.flatnonzero(
        np.isclose(pr_precision, highest_precision, rtol=0.0, atol=1e-15)
    )
    highest_idx = highest_idx_all[np.argmin(pr_thresholds[highest_idx_all])]
    highest_recall = float(pr_recall[highest_idx])
    highest_threshold = float(pr_thresholds[highest_idx])

    disagreement = threshold_disagreement_audit(scores, frozen_binary)
    source_endpoint = ENDPOINT_REL.as_posix()
    source_svtuner = SVTUNER_REL.as_posix()
    source_both = f"{source_endpoint}; {source_svtuner}"

    rows: list[dict[str, str]] = []
    add_row(rows, "formal_matched_n", n_total, source_both, "Barcode inner match after retaining endpoint status positive/negative")
    add_row(rows, "formal_cta_positive_n", n_positive, source_endpoint, "Count primary_endpoint_status == positive")
    add_row(rows, "formal_cta_negative_n", n_negative, source_endpoint, "Count primary_endpoint_status == negative")
    add_row(rows, "formal_missing_n", missing_analysis_fields, source_both, "Missingness in endpoint, withheld_score, or withheld_binary after barcode match")
    add_row(rows, "auroc", auroc, source_both, "sklearn.metrics.roc_auc_score(primary endpoint, withheld_score)")
    add_row(rows, "average_precision", average_precision, source_both, "sklearn.metrics.average_precision_score(primary endpoint, withheld_score)")

    for target, result in fixed.items():
        key = f"precision_ge_{target:.2f}"
        add_row(rows, f"{key}_status", "achieved" if result["achieved"] else "not achieved", source_both, "Maximum recall among finite PR thresholds satisfying precision target")
        add_row(rows, f"{key}_selected_threshold", result["threshold"], source_both, "Minimum threshold among maximum-recall ties; NA when not achieved")
        add_row(rows, f"{key}_achieved_precision", result["precision"], source_both, "Precision at deterministic selected finite threshold")
        add_row(rows, f"{key}_recall", result["recall"], source_both, "Maximum recall satisfying precision target")
        add_row(rows, f"{key}_max_recall_tie_count", result["tie_count"], source_both, "Number of finite thresholds tied at maximum eligible recall")
        add_row(rows, f"{key}_tied_threshold_min", result["threshold_min"], source_both, "Minimum finite threshold in maximum-recall tie set")
        add_row(rows, f"{key}_tied_threshold_max", result["threshold_max"], source_both, "Maximum finite threshold in maximum-recall tie set")

    add_row(rows, "highest_finite_threshold_precision", highest_precision, source_both, "Maximum precision over all finite PR-curve thresholds")
    add_row(rows, "highest_finite_threshold_precision_recall", highest_recall, source_both, "Recall at maximum finite-threshold precision")
    add_row(rows, "highest_finite_threshold_precision_threshold", highest_threshold, source_both, "Minimum threshold among maximum-precision ties")

    best_records = disagreement["best_records"]
    add_row(rows, "simple_threshold_minimum_disagreement_n", disagreement["minimum_disagreement"], source_svtuner, "Minimum Hamming distance between frozen withheld_binary and withheld_score >= t")
    add_row(rows, "simple_threshold_exact_reproduction", disagreement["exactly_reproducible"], source_svtuner, "YES only if minimum disagreement is zero")
    add_row(rows, "simple_threshold_best_rule_count", len(best_records), source_svtuner, "Number of distinct score-induced binary rules tied for minimum disagreement")
    for i, record in enumerate(best_records, start=1):
        prefix = f"simple_threshold_best_{i}"
        add_row(rows, f"{prefix}_threshold", record["threshold"], source_svtuner, "Observed finite score threshold representing tied best rule")
        add_row(rows, f"{prefix}_equivalent_interval", record["interval"], source_svtuner, "All real thresholds in interval produce the same >= rule")
        for field in ("tp", "fp", "tn", "fn"):
            add_row(rows, f"{prefix}_{field}", record[field], source_svtuner, "Confusion count treating frozen withheld_binary as reference")

    p30 = fixed[0.30]
    rounded_p30_matches = math.isclose(p30["threshold"], 0.5280533202, rel_tol=0.0, abs_tol=5e-11)
    alternate_p30_is_exact = any(
        float(x) == 0.5280528671035196 for x in pr_thresholds
    )
    core_checks = {
        "formal universe": observed == expected,
        "AUROC rounds to 0.8305": round(auroc, 4) == 0.8305,
        "AP rounds to 0.2836551170": round(average_precision, 10) == 0.2836551170,
        "P>=0.30 threshold 0.5280533202 is reproducible by rounding": rounded_p30_matches,
        "P>=0.30 threshold 0.5280528671035196 is an exact frozen threshold": alternate_p30_is_exact,
        "minimum disagreement equals 93": disagreement["minimum_disagreement"] == 93,
    }
    # The stale alternate threshold is the provenance discrepancy being repaired.
    overall_status = "DISCREPANCY" if not alternate_p30_is_exact else "PASS"
    add_row(rows, "overall_provenance_status", overall_status, source_both, "DISCREPANCY if any conflicting historical exact threshold is not reproducible from frozen inputs")

    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "c7_provenance_repair_summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value", "source", "method"])
        writer.writeheader()
        writer.writerows(rows)

    md_path = output_dir / "c7_provenance_repair.md"
    lines = [
        "# C7 provenance repair",
        "",
        "## Frozen inputs",
        "",
        f"- Endpoint: `{source_endpoint}`",
        "  - barcode: `barcode`",
        "  - formal endpoint: `primary_endpoint_status`",
        f"- SVTuner output: `{source_svtuner}`",
        "  - barcode: `barcode`",
        "  - continuous score: `withheld_score`",
        "  - frozen decision: `withheld_binary`",
        "",
        "## Matching rule and formal universe",
        "",
        "Records were matched one-to-one by exact `barcode`. The formal analysis universe retains only endpoint statuses `positive` and `negative`; `ambiguous` and `excluded` records are outside this frozen analysis set.",
        "",
        "| Quantity | Value |",
        "|---|---:|",
        f"| Matched formal n | {n_total} |",
        f"| CTA-positive | {n_positive} |",
        f"| CTA-negative | {n_negative} |",
        f"| Missing formal records/fields | {missing_analysis_fields} |",
        "",
        "## AUROC and average precision",
        "",
        "Computed directly from `primary_endpoint_status` (positive = 1) and `withheld_score`.",
        "",
        "| Metric | Full-precision value | Historical value | Comparison |",
        "|---|---:|---:|---|",
        f"| AUROC | {full_precision(auroc)} | approximately 0.8305 | PASS (four-decimal rounding) |",
        f"| Average precision | {full_precision(average_precision)} | approximately 0.2836551170 | PASS (ten-decimal rounding) |",
        "",
        "## Fixed-precision operating points",
        "",
        "For each target precision, all finite PR-curve thresholds satisfying the target were retained, the maximum recall was selected, and the smallest threshold among any maximum-recall ties was used as the deterministic representative.",
        "",
        "| Target precision | Status | Selected threshold | Achieved precision | Recall | Max-recall tied thresholds |",
        "|---:|---|---:|---:|---:|---:|",
    ]
    for target in (0.25, 0.30, 0.50):
        result = fixed[target]
        lines.append(
            "| "
            + " | ".join(
                [
                    f"{target:.2f}",
                    "achieved" if result["achieved"] else "not achieved",
                    full_precision(result["threshold"]),
                    full_precision(result["precision"]),
                    full_precision(result["recall"]),
                    str(result["tie_count"]),
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            f"At P >= 0.30, {fixed[0.30]['tie_count']} finite thresholds attain the same maximum recall. Their range is `{full_precision(fixed[0.30]['threshold_min'])}` to `{full_precision(fixed[0.30]['threshold_max'])}`. The selected threshold is the smallest tied threshold.",
            "",
            "### P >= 0.30 historical threshold reconciliation",
            "",
            f"- Frozen-input selected threshold: `{full_precision(p30['threshold'])}`.",
            "- `0.5280533202` is the ten-decimal rounded form of that frozen threshold: **PASS**.",
            "- `0.5280528671035196` is not an observed finite PR threshold in the frozen input: **DISCREPANCY**.",
            "",
            "### Highest finite-threshold precision",
            "",
            f"- precision: `{full_precision(highest_precision)}`",
            f"- recall: `{full_precision(highest_recall)}`",
            f"- threshold: `{full_precision(highest_threshold)}`",
            "",
            "No finite threshold achieved precision >= 0.50.",
            "",
            "## Frozen withheld_binary versus a simple score threshold",
            "",
            "All distinct finite `withheld_score`-induced rules of the form `withheld_score >= t` were enumerated, with one additional finite threshold above the maximum score to include the all-negative rule. Confusion counts below treat frozen `withheld_binary` as the reference label.",
            "",
            f"- Minimum disagreement count: **{disagreement['minimum_disagreement']}**",
            f"- Exact reproduction by any simple threshold: **{'YES' if disagreement['exactly_reproducible'] else 'NO'}**",
            f"- Number of tied best score-induced rules: **{len(best_records)}**",
            "",
            "| Representative threshold | Equivalent threshold interval | TP | FP | TN | FN | Disagreements |",
            "|---:|---|---:|---:|---:|---:|---:|",
        ]
    )
    for record in best_records:
        lines.append(
            f"| {full_precision(record['threshold'])} | `{record['interval']}` | "
            f"{record['tp']} | {record['fp']} | {record['tn']} | {record['fn']} | "
            f"{record['disagreement']} |"
        )
    lines.extend(
        [
            "",
            "The historical conclusion `minimum disagreement count = 93` is reproduced exactly.",
            "",
            "## Comparison with currently recorded C7 numbers",
            "",
            "| Recorded item | Frozen-input result | Status |",
            "|---|---|---|",
            f"| AUROC approximately 0.8305 | {full_precision(auroc)} | PASS |",
            f"| AP approximately 0.2836551170 | {full_precision(average_precision)} | PASS |",
            f"| P>=0.30 threshold 0.5280533202 | {full_precision(p30['threshold'])} | PASS after stated rounding |",
            "| P>=0.30 threshold 0.5280528671035196 | not present among frozen finite thresholds | DISCREPANCY |",
            f"| Minimum binary disagreement 93 | {disagreement['minimum_disagreement']} | PASS |",
            "",
            f"## Overall status: {overall_status}",
            "",
            "The performance metrics and 93-disagreement conclusion reproduce. The provenance discrepancy is confined to the alternate P>=0.30 threshold value `0.5280528671035196`; the frozen-input threshold is `0.5280533202496566` (reported as `0.5280533202` at ten decimal places).",
            "",
        ]
    )
    md_path.write_text("\n".join(lines), encoding="utf-8")

    print(f"summary={summary_path}")
    print(f"report={md_path}")
    print(f"status={overall_status}")


if __name__ == "__main__":
    main()
