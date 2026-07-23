from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import auc, average_precision_score, precision_recall_curve, roc_auc_score, roc_curve


ROOT = Path(__file__).resolve().parents[1]

P2C = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping"
P3R = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase3r_review_resolution_for_cta_endpoint_baseline_feasibility"
P5 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase5_endpoint_specific_evaluation_of_cytospace_baselines_against_frozen_cta_immune_endpoint"
P6 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase6_svtuner_aware_execution_planning_against_frozen_cta_immune_endpoint"
P7 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase7_svtuner_aware_execution_against_frozen_cta_immune_endpoint"

OUT = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison"

BASELINE_ANCHOR = "baseline_immune_all_dropout"
SVTUNER_CONDITION = "svtuner_immune_all_dropout"
SEED = 20260705
BOOTSTRAP_N = 1000


def to_builtin(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): to_builtin(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [to_builtin(v) for v in obj]
    if isinstance(obj, tuple):
        return [to_builtin(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        if np.isnan(obj):
            return None
        return float(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, float) and math.isnan(obj):
        return None
    return obj


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(to_builtin(data), indent=2), encoding="utf-8")


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return str(path).replace("\\", "/")


def ensure_dirs() -> None:
    OUT.mkdir(parents=True, exist_ok=True)


def bool_series(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s.fillna(False)
    return s.astype(str).str.lower().isin(["true", "1", "yes"])


def cohen_d(a: np.ndarray, b: np.ndarray) -> float | None:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return None
    pooled = np.sqrt(((len(a) - 1) * np.var(a, ddof=1) + (len(b) - 1) * np.var(b, ddof=1)) / (len(a) + len(b) - 2))
    if pooled == 0:
        return None
    return float((np.mean(a) - np.mean(b)) / pooled)


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float | None:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) == 0 or len(b) == 0:
        return None
    gt = 0
    lt = 0
    for x in a:
        gt += int(np.sum(x > b))
        lt += int(np.sum(x < b))
    return float((gt - lt) / (len(a) * len(b)))


def bootstrap_delta(a: np.ndarray, b: np.ndarray, n: int = BOOTSTRAP_N, seed: int = SEED) -> tuple[float | None, float | None]:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) == 0 or len(b) == 0:
        return None, None
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n):
        vals.append(float(np.mean(rng.choice(a, size=len(a), replace=True)) - np.mean(rng.choice(b, size=len(b), replace=True))))
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def bootstrap_rate(values: np.ndarray, n: int = BOOTSTRAP_N, seed: int = SEED) -> tuple[float | None, float | None]:
    values = np.asarray(values, dtype=float)
    values = values[~np.isnan(values)]
    if len(values) == 0:
        return None, None
    rng = np.random.default_rng(seed)
    vals = [float(np.mean(rng.choice(values, size=len(values), replace=True))) for _ in range(n)]
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def safe_ratio(num: float, den: float) -> float | None:
    if den == 0 or pd.isna(den):
        return None
    return float(num / den)


def validate_inputs() -> tuple[list[str], dict[str, Any]]:
    errors: list[str] = []
    p2c = read_json(P2C / "bioapp_phase2c_endpoint_freeze_summary.json")
    p5 = read_json(P5 / "bioapp_phase5_endpoint_specific_baseline_evaluation_summary.json")
    p6 = read_json(P6 / "bioapp_phase6_svtuner_execution_planning_summary.json")
    p7 = read_json(P7 / "bioapp_phase7_svtuner_aware_execution_summary.json")

    checks = [
        (p2c.get("decision") == "PASS", "Phase 2C decision is not PASS"),
        (p2c.get("endpoint_frozen") is True, "Phase 2C endpoint is not frozen"),
        (p2c.get("primary_endpoint") == "Immune cells", "Phase 2C primary endpoint is not Immune cells"),
        (p5.get("decision") == "PASS", "Phase 5 decision is not PASS"),
        (p5.get("baseline_vs_endpoint_metric_computed") is True, "Phase 5 baseline metric missing"),
        (p5.get("ready_for_phase6_svtuner_execution") is True, "Phase 5 not ready for Phase 6"),
        (p6.get("decision") == "PASS", "Phase 6 decision is not PASS"),
        (p6.get("baseline_anchor_condition") == BASELINE_ANCHOR, "Phase 6 baseline anchor mismatch"),
        (p6.get("primary_svtuner_condition") == SVTUNER_CONDITION, "Phase 6 SVTuner condition mismatch"),
        (p6.get("threshold_selected_using_endpoint_labels") is False, "Phase 6 threshold policy invalid"),
        (p7.get("decision") == "PASS", "Phase 7 decision is not PASS"),
        (p7.get("SVTuner_run") is True, "Phase 7 SVTuner run missing"),
        (p7.get("SVTuner_Stage3_run") is True, "Phase 7 Stage3 run missing"),
        (p7.get("SVTuner_Stage4_run") is False, "Phase 7 Stage4 unexpectedly run"),
        (p7.get("Stage4_run") is False, "Phase 7 Stage4 flag unexpectedly true"),
        (p7.get("CytoSPACE_rerun") is False, "Phase 7 CytoSPACE rerun unexpectedly true"),
        (p7.get("ready_for_phase8_evaluation") is True, "Phase 7 not ready for Phase 8"),
        (p7.get("spot_level_contract_rows") == 2248, "Phase 7 contract row count mismatch"),
        (p7.get("contract_spot_universe_aligned") is True, "Phase 7 contract universe not aligned"),
        (p7.get("threshold_selected_using_endpoint_labels") is False, "Phase 7 threshold selected using endpoint"),
        (p7.get("prevention_analysis_run") is False, "Phase 7 prevention analysis already run"),
        (p7.get("withheld_enrichment_computed") is False, "Phase 7 withheld enrichment already computed"),
        (p7.get("withheld_AUROC_AUPRC_computed") is False, "Phase 7 AUROC/AUPRC already computed"),
        (
            p7.get("endpoint_positive_nonimmune_assignment_reduction_computed") is False,
            "Phase 7 burden reduction already computed",
        ),
        (p7.get("contradiction_prevention_rate_computed") is False, "Phase 7 prevention rate already computed"),
        (
            p7.get("endpoint_concordant_interpretation_score_computed") is False,
            "Phase 7 endpoint-concordant interpretation already computed",
        ),
        (p7.get("baseline_svtuner_comparison_computed") is False, "Phase 7 baseline/SVTuner comparison already computed"),
    ]
    for ok, msg in checks:
        if not ok:
            errors.append(msg)

    required = [
        P2C / "spot_level_endpoint_freeze.csv",
        P5 / "baseline_endpoint_score_all_runs_by_spot.csv",
        P5 / "baseline_endpoint_metric_summary.csv",
        P5 / "baseline_condition_metric_comparison.csv",
        P6 / "svtuner_metric_alignment_plan.json",
        P6 / "svtuner_threshold_policy.json",
        P6 / "svtuner_output_contract.json",
        P7 / SVTUNER_CONDITION / "svtuner_immune_all_dropout_spot_level_raw_output_contract.csv",
        P7 / SVTUNER_CONDITION / "raw_outputs" / "stage3b_summary.json",
        P7 / SVTUNER_CONDITION / "raw_outputs" / "spot_unsupported_scores.csv",
    ]
    for path in required:
        if not path.exists():
            errors.append(f"Missing required input: {rel(path)}")

    return errors, {"phase2c": p2c, "phase5": p5, "phase6": p6, "phase7": p7}


def gene_count_audit(context: dict[str, Any]) -> dict[str, Any]:
    gene_file = P3R / "formal_gene_intersection.txt"
    raw_text = gene_file.read_text(encoding="utf-8")
    raw_lines = raw_text.splitlines()
    nonempty = [x.strip() for x in raw_lines if x.strip()]
    unique = sorted(set(nonempty))
    first = nonempty[0] if nonempty else ""
    header_detected = first.lower() in {"gene", "genes", "gene_id", "feature", "features"} or first.startswith("#")
    header_aware = len(nonempty) - 1 if header_detected else len(nonempty)
    stage3b = read_json(P7 / SVTUNER_CONDITION / "raw_outputs" / "stage3b_summary.json")
    stage3b_gene_count = stage3b.get("dimensions", {}).get("genes")
    phase_counts = [
        context["phase5"].get("n_gene_overlap_used_in_phase4"),
        context["phase6"].get("n_gene_overlap"),
        context["phase7"].get("n_gene_overlap"),
    ]
    phase_summary_gene_count = next((x for x in phase_counts if x is not None), None)
    discrepancy = len(nonempty) != stage3b_gene_count or len(unique) != stage3b_gene_count
    material = not (stage3b_gene_count == 2000 and len(nonempty) == 2000 and len(unique) == 2000)
    explanation = (
        "formal_gene_intersection.txt contains 2000 non-empty unique gene lines; Stage3B also reports 2000 genes. "
        "The previous 1999 count is attributable to a line-count/read-mode artifact, not a material mismatch."
        if not material
        else "Unresolved difference between gene list counts and Stage3B reported dimensions."
    )
    audit = {
        "source": rel(gene_file),
        "raw_line_count": len(raw_lines),
        "nonempty_line_count": len(nonempty),
        "unique_gene_count": len(unique),
        "header_detected": header_detected,
        "header_aware_gene_count": header_aware,
        "stage3b_reported_gene_count": stage3b_gene_count,
        "phase_summary_gene_count": phase_summary_gene_count,
        "final_effective_gene_count_for_phase8": 2000 if not material else None,
        "count_discrepancy_detected": discrepancy,
        "discrepancy_explanation": explanation,
        "material_mismatch": material,
        "decision_impact": "none" if not material else "REVIEW_REQUIRED",
    }
    pd.DataFrame([audit]).to_csv(OUT / "gene_count_audit.csv", index=False)
    write_json(OUT / "gene_count_audit.json", audit)
    return audit


def load_tables() -> dict[str, pd.DataFrame]:
    endpoint = pd.read_csv(P2C / "spot_level_endpoint_freeze.csv")
    baseline_all = pd.read_csv(P5 / "baseline_endpoint_score_all_runs_by_spot.csv")
    baseline_anchor = baseline_all[baseline_all["baseline_run"] == BASELINE_ANCHOR].copy()
    contract = pd.read_csv(P7 / SVTUNER_CONDITION / "svtuner_immune_all_dropout_spot_level_raw_output_contract.csv")
    raw_scores = pd.read_csv(P7 / SVTUNER_CONDITION / "raw_outputs" / "spot_unsupported_scores.csv")
    return {
        "endpoint": endpoint,
        "baseline_all": baseline_all,
        "baseline_anchor": baseline_anchor,
        "contract": contract,
        "raw_scores": raw_scores,
    }


def endpoint_qc(endpoint: pd.DataFrame) -> dict[str, Any]:
    counts = endpoint["primary_endpoint_status"].value_counts(dropna=False).to_dict()
    qc = {
        "source": rel(P2C / "spot_level_endpoint_freeze.csv"),
        "positive_count": int(counts.get("positive", 0)),
        "negative_count": int(counts.get("negative", 0)),
        "ambiguous_count": int(counts.get("ambiguous", 0)),
        "excluded_count": int(counts.get("excluded", 0)),
        "main_analysis_count": int(counts.get("positive", 0) + counts.get("negative", 0)),
        "used_as_authoritative": True,
        "notes": "spot_level_endpoint_freeze.csv is authoritative; txt list line-count differences are ignored for endpoint status.",
    }
    pd.DataFrame([qc]).to_csv(OUT / "endpoint_analysis_set_qc.csv", index=False)
    return qc


def baseline_import_qc(endpoint: pd.DataFrame, baseline_anchor: pd.DataFrame) -> dict[str, Any]:
    overlap = set(endpoint["barcode"]).intersection(set(baseline_anchor["barcode"]))
    qc = {
        "n_rows": int(len(baseline_anchor)),
        "n_unique_spots": int(baseline_anchor["barcode"].nunique()),
        "n_overlap_with_frozen_endpoint": int(len(overlap)),
        "baseline_immune_score_constant": bool(baseline_anchor["immune_score"].nunique(dropna=True) <= 1),
        "baseline_nonimmune_score_available": bool(baseline_anchor["nonimmune_score"].notna().all()),
        "baseline_anchor_import_status": "ok" if len(overlap) == len(endpoint) and baseline_anchor["barcode"].nunique() == len(endpoint) else "review_required",
        "notes": "baseline anchor imported from Phase 5 without recomputation.",
    }
    pd.DataFrame([qc]).to_csv(OUT / "phase5_baseline_anchor_import_qc.csv", index=False)
    return qc


def phase7_import_qc(endpoint: pd.DataFrame, contract: pd.DataFrame) -> dict[str, Any]:
    overlap = set(endpoint["barcode"]).intersection(set(contract["barcode"]))
    qc = {
        "n_rows": int(len(contract)),
        "n_unique_spots": int(contract["barcode"].nunique()),
        "n_overlap_with_frozen_endpoint": int(len(overlap)),
        "withheld_score_available": bool(contract["withheld_score"].notna().all()),
        "withheld_binary_available": "withheld_binary" in contract.columns,
        "phase7_contract_import_status": "ok" if len(overlap) == len(endpoint) and contract["barcode"].nunique() == len(endpoint) else "review_required",
        "notes": "Phase 7 raw SVTuner contract imported without rerunning SVTuner or Stage3B.",
    }
    pd.DataFrame([qc]).to_csv(OUT / "phase7_svtuner_raw_output_import_qc.csv", index=False)
    return qc


def build_score_table(endpoint: pd.DataFrame, baseline_anchor: pd.DataFrame, contract: pd.DataFrame) -> pd.DataFrame:
    base_cols = [
        "barcode",
        "baseline_run",
        "immune_score",
        "nonimmune_score",
        "dominant_baseline_label",
        "dominant_baseline_fraction",
    ]
    base = baseline_anchor[base_cols].rename(
        columns={
            "baseline_run": "baseline_anchor_run",
            "immune_score": "baseline_immune_score",
            "nonimmune_score": "baseline_nonimmune_score",
            "dominant_baseline_label": "baseline_dominant_label",
            "dominant_baseline_fraction": "baseline_dominant_fraction",
        }
    )
    c = contract[
        [
            "barcode",
            "reference_unrepresented_score",
            "withheld_score",
            "withheld_binary",
            "supported_assignment_score",
            "reference_supported_status",
            "svtuner_immune_score",
            "svtuner_nonimmune_score",
            "svtuner_dominant_label",
            "svtuner_dominant_fraction",
            "score_parse_status",
        ]
    ].copy()
    c["withheld_binary"] = bool_series(c["withheld_binary"])
    merged = endpoint[["barcode", "primary_endpoint_status"]].merge(c, on="barcode", how="left").merge(base, on="barcode", how="left")
    merged["is_endpoint_positive"] = merged["primary_endpoint_status"].eq("positive")
    merged["is_endpoint_negative"] = merged["primary_endpoint_status"].eq("negative")
    merged["is_endpoint_ambiguous"] = merged["primary_endpoint_status"].eq("ambiguous")
    merged["is_endpoint_excluded"] = merged["primary_endpoint_status"].eq("excluded")
    merged["score_source"] = np.where(merged["withheld_score"].notna(), "withheld_score", "reference_unrepresented_score")
    merged["score_parse_status"] = merged["score_parse_status"].fillna("missing_phase7_contract")
    ordered = [
        "barcode",
        "primary_endpoint_status",
        "is_endpoint_positive",
        "is_endpoint_negative",
        "is_endpoint_ambiguous",
        "is_endpoint_excluded",
        "withheld_score",
        "reference_unrepresented_score",
        "withheld_binary",
        "supported_assignment_score",
        "reference_supported_status",
        "baseline_anchor_run",
        "baseline_immune_score",
        "baseline_nonimmune_score",
        "baseline_dominant_label",
        "baseline_dominant_fraction",
        "svtuner_immune_score",
        "svtuner_nonimmune_score",
        "svtuner_dominant_label",
        "svtuner_dominant_fraction",
        "score_source",
        "score_parse_status",
    ]
    merged[ordered].to_csv(OUT / "svtuner_endpoint_score_by_spot.csv", index=False)
    return merged


def compute_withheld_metrics(scores: pd.DataFrame) -> dict[str, Any]:
    main = scores[scores["is_endpoint_positive"] | scores["is_endpoint_negative"]].copy()
    pos = main.loc[main["is_endpoint_positive"], "withheld_score"].astype(float).to_numpy()
    neg = main.loc[main["is_endpoint_negative"], "withheld_score"].astype(float).to_numpy()
    delta = float(np.nanmean(pos) - np.nanmean(neg))
    ci_low, ci_high = bootstrap_delta(pos, neg)
    try:
        mw_p = float(stats.mannwhitneyu(pos, neg, alternative="two-sided").pvalue)
    except Exception:
        mw_p = None
    summary = {
        "n_positive_spots": int(len(pos)),
        "n_negative_spots": int(len(neg)),
        "mean_withheld_score_positive": float(np.nanmean(pos)),
        "mean_withheld_score_negative": float(np.nanmean(neg)),
        "median_withheld_score_positive": float(np.nanmedian(pos)),
        "median_withheld_score_negative": float(np.nanmedian(neg)),
        "delta_mean_positive_minus_negative": delta,
        "fold_change_positive_over_negative": safe_ratio(float(np.nanmean(pos)), float(np.nanmean(neg))),
        "bootstrap_n": BOOTSTRAP_N,
        "bootstrap_95CI_delta_low": ci_low,
        "bootstrap_95CI_delta_high": ci_high,
        "Mann_Whitney_U_pvalue": mw_p,
        "Cohen_d": cohen_d(pos, neg),
        "Cliffs_delta": cliffs_delta(pos, neg),
        "metric_status": "ok",
    }
    pd.DataFrame([summary]).to_csv(OUT / "withheld_enrichment_summary.csv", index=False)

    y = main["is_endpoint_positive"].astype(int).to_numpy()
    x = main["withheld_score"].astype(float).to_numpy()
    if len(np.unique(y)) == 2 and len(np.unique(x[~np.isnan(x)])) > 1:
        auroc = float(roc_auc_score(y, x))
        auprc = float(average_precision_score(y, x))
        fpr, tpr, roc_thresh = roc_curve(y, x)
        precision, recall, pr_thresh = precision_recall_curve(y, x)
        pd.DataFrame({"fpr": fpr, "tpr": tpr, "threshold": np.r_[roc_thresh]}).to_csv(OUT / "withheld_roc_curve.csv", index=False)
        pr_df = pd.DataFrame({"precision": precision, "recall": recall})
        pr_df["threshold"] = np.r_[pr_thresh, np.nan]
        pr_df.to_csv(OUT / "withheld_pr_curve.csv", index=False)
    else:
        auroc = None
        auprc = None
        pd.DataFrame(columns=["fpr", "tpr", "threshold"]).to_csv(OUT / "withheld_roc_curve.csv", index=False)
        pd.DataFrame(columns=["precision", "recall", "threshold"]).to_csv(OUT / "withheld_pr_curve.csv", index=False)
    summary["withheld_AUROC"] = auroc
    summary["withheld_AUPRC"] = auprc
    summary["positive_class_prevalence"] = float(np.mean(y))
    return summary


def compute_binary_metrics(scores: pd.DataFrame) -> dict[str, Any]:
    main = scores[scores["is_endpoint_positive"] | scores["is_endpoint_negative"]].copy()
    pos = main[main["is_endpoint_positive"]]
    neg = main[main["is_endpoint_negative"]]
    pos_with = int(pos["withheld_binary"].sum())
    neg_with = int(neg["withheld_binary"].sum())
    pos_no = int(len(pos) - pos_with)
    neg_no = int(len(neg) - neg_with)
    contingency = pd.DataFrame(
        [
            {"endpoint_status": "positive", "withheld_true": pos_with, "withheld_false": pos_no, "total": len(pos)},
            {"endpoint_status": "negative", "withheld_true": neg_with, "withheld_false": neg_no, "total": len(neg)},
        ]
    )
    contingency.to_csv(OUT / "binary_withheld_endpoint_contingency.csv", index=False)
    try:
        fisher_p = float(stats.fisher_exact([[pos_with, pos_no], [neg_with, neg_no]])[1])
    except Exception:
        fisher_p = None
    pos_rate = pos_with / len(pos) if len(pos) else np.nan
    neg_rate = neg_with / len(neg) if len(neg) else np.nan
    odds_ratio = ((pos_with + 0.5) * (neg_no + 0.5)) / ((pos_no + 0.5) * (neg_with + 0.5))
    summary = {
        "withheld_binary_positive_count": pos_with,
        "withheld_binary_negative_count": neg_with,
        "withheld_binary_positive_rate": float(pos_rate),
        "withheld_binary_negative_rate": float(neg_rate),
        "risk_difference_positive_minus_negative": float(pos_rate - neg_rate),
        "risk_ratio": safe_ratio(float(pos_rate), float(neg_rate)),
        "odds_ratio": float(odds_ratio),
        "Fisher_exact_pvalue": fisher_p,
        "threshold_source": "Phase 7 predefined Stage3B binary output; no Phase 8 threshold tuning",
    }
    pd.DataFrame([summary]).to_csv(OUT / "binary_withheld_enrichment_summary.csv", index=False)
    return summary


def compute_burden_and_prevention(scores: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any], dict[str, Any], dict[str, Any]]:
    comp = scores[
        [
            "barcode",
            "primary_endpoint_status",
            "baseline_nonimmune_score",
            "withheld_score",
            "withheld_binary",
        ]
    ].copy()
    comp["baseline_forced_nonimmune_burden"] = pd.to_numeric(comp["baseline_nonimmune_score"], errors="coerce")
    comp["withheld_score"] = pd.to_numeric(comp["withheld_score"], errors="coerce").clip(0, 1)
    comp["withheld_binary"] = bool_series(comp["withheld_binary"])
    comp["svtuner_retained_forced_nonimmune_burden_binary"] = comp["baseline_forced_nonimmune_burden"] * (~comp["withheld_binary"]).astype(float)
    comp["svtuner_prevented_forced_nonimmune_burden_binary"] = comp["baseline_forced_nonimmune_burden"] * comp["withheld_binary"].astype(float)
    comp["svtuner_retained_forced_nonimmune_burden_continuous"] = comp["baseline_forced_nonimmune_burden"] * (1 - comp["withheld_score"])
    comp["svtuner_prevented_forced_nonimmune_burden_continuous"] = comp["baseline_forced_nonimmune_burden"] * comp["withheld_score"]
    comp["comparison_status"] = "withheld-aware forced-call burden; not post-Stage4 remapping"
    comp.to_csv(OUT / "baseline_svtuner_endpoint_comparison_by_spot.csv", index=False)

    pos = comp[comp["primary_endpoint_status"] == "positive"]
    burden_summary = {
        "n_endpoint_positive_spots": int(len(pos)),
        "mean_baseline_forced_nonimmune_burden": float(pos["baseline_forced_nonimmune_burden"].mean()),
        "mean_retained_burden_binary": float(pos["svtuner_retained_forced_nonimmune_burden_binary"].mean()),
        "mean_prevented_burden_binary": float(pos["svtuner_prevented_forced_nonimmune_burden_binary"].mean()),
        "mean_retained_burden_continuous": float(pos["svtuner_retained_forced_nonimmune_burden_continuous"].mean()),
        "mean_prevented_burden_continuous": float(pos["svtuner_prevented_forced_nonimmune_burden_continuous"].mean()),
        "fraction_burden_prevented_binary": safe_ratio(float(pos["svtuner_prevented_forced_nonimmune_burden_binary"].sum()), float(pos["baseline_forced_nonimmune_burden"].sum())),
        "fraction_burden_prevented_continuous": safe_ratio(float(pos["svtuner_prevented_forced_nonimmune_burden_continuous"].sum()), float(pos["baseline_forced_nonimmune_burden"].sum())),
        "metric_note": "withheld-aware forced-call burden estimate; not a post-Stage4 remapping result",
    }
    pd.DataFrame([burden_summary]).to_csv(OUT / "endpoint_positive_forced_nonimmune_burden_reduction_summary.csv", index=False)

    denom = pos[pos["baseline_forced_nonimmune_burden"] >= 0.5].copy()
    prevented = denom[denom["withheld_binary"]]
    ci_low, ci_high = bootstrap_rate(denom["withheld_binary"].astype(float).to_numpy())
    prevention_summary = {
        "denominator_rule": "CTA Immune-positive spots with baseline_nonimmune_score >= 0.5",
        "denominator_spots": int(len(denom)),
        "prevented_spots": int(len(prevented)),
        "not_prevented_spots": int(len(denom) - len(prevented)),
        "contradiction_prevention_rate": safe_ratio(float(len(prevented)), float(len(denom))),
        "bootstrap_n": BOOTSTRAP_N,
        "bootstrap_95CI_low": ci_low,
        "bootstrap_95CI_high": ci_high,
        "metric_note": "Phase 8 operational contradiction prevention metric; not final biological application claim",
    }
    pd.DataFrame([prevention_summary]).to_csv(OUT / "contradiction_prevention_summary.csv", index=False)
    denom.to_csv(OUT / "contradiction_prevention_spots.csv", index=False)

    interp = comp[comp["primary_endpoint_status"].isin(["positive", "negative"])].copy()
    interp["endpoint_concordant_interpretation_binary"] = np.where(
        interp["primary_endpoint_status"].eq("positive"),
        interp["withheld_binary"].astype(int),
        np.nan,
    )
    interp.to_csv(OUT / "endpoint_concordant_interpretation_by_spot.csv", index=False)
    pos_rate = float(comp.loc[comp["primary_endpoint_status"].eq("positive"), "withheld_binary"].mean())
    neg_rate = float(comp.loc[comp["primary_endpoint_status"].eq("negative"), "withheld_binary"].mean())
    interp_summary = {
        "positive_withheld_rate": pos_rate,
        "negative_withheld_rate": neg_rate,
        "positive_minus_negative_withheld_rate": pos_rate - neg_rate,
        "endpoint_concordant_interpretation_score_positive": pos_rate,
        "specificity_reference_negative_withheld_rate": neg_rate,
    }
    pd.DataFrame([interp_summary]).to_csv(OUT / "endpoint_concordant_interpretation_summary.csv", index=False)
    return comp, burden_summary, prevention_summary, interp_summary


def compute_spatial_overlap(scores: pd.DataFrame, comp: pd.DataFrame) -> pd.DataFrame:
    positive = set(scores.loc[scores["is_endpoint_positive"], "barcode"])
    withheld = set(scores.loc[bool_series(scores["withheld_binary"]), "barcode"])
    baseline_forced = set(comp.loc[comp["baseline_forced_nonimmune_burden"] >= 0.5, "barcode"])
    k = len(withheld)
    top_k = set(scores.sort_values("withheld_score", ascending=False).head(k)["barcode"]) if k > 0 else set()

    rows = []

    def add(name_a: str, a: set[str], name_b: str, b: set[str], notes: str) -> None:
        overlap = a & b
        rows.append(
            {
                "set_A": name_a,
                "set_B": name_b,
                "n_A": len(a),
                "n_B": len(b),
                "n_overlap": len(overlap),
                "overlap_fraction_A": safe_ratio(float(len(overlap)), float(len(a))),
                "overlap_fraction_B": safe_ratio(float(len(overlap)), float(len(b))),
                "jaccard_index": safe_ratio(float(len(overlap)), float(len(a | b))),
                "notes": notes,
            }
        )

    add("CTA Immune-positive", positive, "SVTuner withheld_binary", withheld, "Phase 7 binary withheld set")
    add("CTA Immune-positive", positive, f"high withheld_score top-k (k={k})", top_k, "k fixed to number of withheld_binary spots")
    add("CTA Immune-positive", positive, "baseline forced nonimmune", baseline_forced, "baseline_nonimmune_score >= 0.5")
    add("baseline forced nonimmune in CTA-positive", positive & baseline_forced, "SVTuner withheld_binary", withheld, "denominator subset vs withheld")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "spatial_overlap_summary.csv", index=False)
    return df


def metric_comparison(withheld_summary: dict[str, Any], binary_summary: dict[str, Any], burden: dict[str, Any], prevention: dict[str, Any], interp: dict[str, Any]) -> pd.DataFrame:
    rows = [
        {"metric": "withheld_AUROC", "value": withheld_summary.get("withheld_AUROC"), "source": "SVTuner withheld_score vs endpoint"},
        {"metric": "withheld_AUPRC", "value": withheld_summary.get("withheld_AUPRC"), "source": "SVTuner withheld_score vs endpoint"},
        {"metric": "withheld_enrichment_delta", "value": withheld_summary.get("delta_mean_positive_minus_negative"), "source": "positive minus negative withheld_score"},
        {"metric": "withheld_binary_positive_rate", "value": binary_summary.get("withheld_binary_positive_rate"), "source": "Phase 7 binary withheld output"},
        {"metric": "withheld_binary_negative_rate", "value": binary_summary.get("withheld_binary_negative_rate"), "source": "Phase 7 binary withheld output"},
        {"metric": "contradiction_prevention_rate", "value": prevention.get("contradiction_prevention_rate"), "source": "withheld-aware forced-call burden"},
        {"metric": "fraction_burden_prevented_binary", "value": burden.get("fraction_burden_prevented_binary"), "source": "withheld-aware forced-call burden"},
        {"metric": "fraction_burden_prevented_continuous", "value": burden.get("fraction_burden_prevented_continuous"), "source": "withheld-aware forced-call burden"},
        {"metric": "endpoint_concordant_interpretation_score_positive", "value": interp.get("endpoint_concordant_interpretation_score_positive"), "source": "positive spots withheld rate"},
    ]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "baseline_svtuner_metric_comparison.csv", index=False)
    return df


def plot_spatial(endpoint: pd.DataFrame, scores: pd.DataFrame, comp: pd.DataFrame) -> list[dict[str, str]]:
    df = endpoint[["barcode", "imagecol", "imagerow", "primary_endpoint_status"]].merge(scores[["barcode", "withheld_score", "withheld_binary"]], on="barcode").merge(
        comp[["barcode", "baseline_forced_nonimmune_burden"]], on="barcode"
    )
    pos = df["primary_endpoint_status"].eq("positive")
    fig, axes = plt.subplots(1, 5, figsize=(18, 4), constrained_layout=True)
    for ax in axes:
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])

    status_colors = {"positive": "#d62728", "negative": "#d9d9d9", "ambiguous": "#bdbdbd", "excluded": "#ffffff"}
    axes[0].scatter(df["imagecol"], df["imagerow"], c=df["primary_endpoint_status"].map(status_colors), s=4, linewidths=0)
    axes[0].set_title("A. Frozen CTA Immune endpoint")
    sc = axes[1].scatter(df["imagecol"], df["imagerow"], c=df["baseline_forced_nonimmune_burden"], cmap="magma", s=4, linewidths=0, vmin=0, vmax=1)
    axes[1].set_title("B. Baseline forced nonimmune burden")
    fig.colorbar(sc, ax=axes[1], fraction=0.046)
    sc = axes[2].scatter(df["imagecol"], df["imagerow"], c=df["withheld_score"], cmap="viridis", s=4, linewidths=0, vmin=0, vmax=1)
    axes[2].set_title("C. SVTuner withheld score")
    fig.colorbar(sc, ax=axes[2], fraction=0.046)
    axes[3].scatter(df["imagecol"], df["imagerow"], c=np.where(df["withheld_binary"], "#1f78b4", "#d9d9d9"), s=4, linewidths=0)
    axes[3].set_title("D. SVTuner withheld binary")
    axes[4].scatter(df["imagecol"], df["imagerow"], c="#d9d9d9", s=3, linewidths=0)
    axes[4].scatter(df.loc[pos, "imagecol"], df.loc[pos, "imagerow"], facecolors="none", edgecolors="#d62728", s=18, linewidths=0.6)
    axes[4].scatter(df.loc[df["withheld_binary"], "imagecol"], df.loc[df["withheld_binary"], "imagerow"], c="#1f78b4", s=5, linewidths=0)
    axes[4].set_title("E. Endpoint outline + withheld")
    fig.suptitle("Baseline/SVTuner comparison against frozen CTA Immune endpoint", fontsize=12)
    paths = []
    for ext in ["pdf", "svg"]:
        path = OUT / f"endpoint_baseline_svtuner_three_way_spatial_comparison.{ext}"
        fig.savefig(path)
        paths.append({"file": rel(path), "type": "figure", "description": "Three-way spatial comparison"})
    plt.close(fig)
    return paths


def plot_box(scores: pd.DataFrame) -> list[dict[str, str]]:
    main = scores[scores["primary_endpoint_status"].isin(["positive", "negative"])]
    data = [
        main.loc[main["primary_endpoint_status"].eq("positive"), "withheld_score"].astype(float).to_numpy(),
        main.loc[main["primary_endpoint_status"].eq("negative"), "withheld_score"].astype(float).to_numpy(),
    ]
    fig, ax = plt.subplots(figsize=(4, 4), constrained_layout=True)
    ax.boxplot(data, tick_labels=["Positive", "Negative"], patch_artist=True, boxprops={"facecolor": "#88ccee"})
    ax.set_ylabel("withheld_score")
    ax.set_title("SVTuner withheld score by CTA Immune endpoint status")
    paths = []
    for ext in ["pdf", "svg"]:
        path = OUT / f"withheld_score_positive_vs_negative_boxplot.{ext}"
        fig.savefig(path)
        paths.append({"file": rel(path), "type": "figure", "description": "Withheld score positive vs negative boxplot"})
    plt.close(fig)
    return paths


def plot_roc_pr() -> list[dict[str, str]]:
    paths = []
    roc = pd.read_csv(OUT / "withheld_roc_curve.csv")
    pr = pd.read_csv(OUT / "withheld_pr_curve.csv")
    fig, ax = plt.subplots(figsize=(4, 4), constrained_layout=True)
    if not roc.empty:
        ax.plot(roc["fpr"], roc["tpr"], color="#0072b2")
    ax.plot([0, 1], [0, 1], color="#999999", linestyle="--", linewidth=1)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("Withheld score ROC")
    for ext in ["pdf", "svg"]:
        path = OUT / f"withheld_score_roc_curve.{ext}"
        fig.savefig(path)
        paths.append({"file": rel(path), "type": "figure", "description": "Withheld score ROC curve"})
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4, 4), constrained_layout=True)
    if not pr.empty:
        ax.plot(pr["recall"], pr["precision"], color="#d55e00")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Withheld score PR")
    for ext in ["pdf", "svg"]:
        path = OUT / f"withheld_score_pr_curve.{ext}"
        fig.savefig(path)
        paths.append({"file": rel(path), "type": "figure", "description": "Withheld score PR curve"})
    plt.close(fig)
    return paths


def plot_binary_contingency() -> list[dict[str, str]]:
    cont = pd.read_csv(OUT / "binary_withheld_endpoint_contingency.csv")
    fig, ax = plt.subplots(figsize=(4, 4), constrained_layout=True)
    x = np.arange(len(cont))
    ax.bar(x, cont["withheld_true"], label="withheld", color="#0072b2")
    ax.bar(x, cont["withheld_false"], bottom=cont["withheld_true"], label="not withheld", color="#cccccc")
    ax.set_xticks(x)
    ax.set_xticklabels(cont["endpoint_status"])
    ax.set_ylabel("Spot count")
    ax.set_title("Binary withheld vs CTA endpoint")
    ax.legend(frameon=False)
    paths = []
    for ext in ["pdf", "svg"]:
        path = OUT / f"binary_withheld_endpoint_contingency_plot.{ext}"
        fig.savefig(path)
        paths.append({"file": rel(path), "type": "figure", "description": "Binary withheld endpoint contingency"})
    plt.close(fig)
    return paths


def plot_burden() -> list[dict[str, str]]:
    s = pd.read_csv(OUT / "endpoint_positive_forced_nonimmune_burden_reduction_summary.csv").iloc[0]
    labels = ["Baseline forced", "Retained (binary)", "Prevented (binary)", "Retained (continuous)", "Prevented (continuous)"]
    vals = [
        s["mean_baseline_forced_nonimmune_burden"],
        s["mean_retained_burden_binary"],
        s["mean_prevented_burden_binary"],
        s["mean_retained_burden_continuous"],
        s["mean_prevented_burden_continuous"],
    ]
    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    ax.bar(labels, vals, color=["#999999", "#cc79a7", "#0072b2", "#cc79a7", "#0072b2"])
    ax.set_ylim(0, max(1.0, max(vals) * 1.15))
    ax.set_ylabel("Mean burden in CTA Immune-positive spots")
    ax.set_title("Withheld-aware forced-call burden")
    ax.tick_params(axis="x", rotation=30)
    paths = []
    for ext in ["pdf", "svg"]:
        path = OUT / f"endpoint_positive_forced_nonimmune_burden_reduction.{ext}"
        fig.savefig(path)
        paths.append({"file": rel(path), "type": "figure", "description": "Endpoint-positive forced nonimmune burden reduction"})
    plt.close(fig)
    return paths


def plot_metric_summary(metric_df: pd.DataFrame) -> list[dict[str, str]]:
    plot_df = metric_df[metric_df["metric"].isin(["withheld_AUROC", "withheld_AUPRC", "withheld_enrichment_delta", "contradiction_prevention_rate", "withheld_binary_positive_rate"])].copy()
    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    ax.bar(plot_df["metric"], plot_df["value"], color="#56b4e9")
    ax.set_ylabel("Metric value")
    ax.set_title("Phase 8 SVTuner endpoint metric summary")
    ax.tick_params(axis="x", rotation=35)
    paths = []
    for ext in ["pdf", "svg"]:
        path = OUT / f"phase8_svtuner_endpoint_metric_summary.{ext}"
        fig.savefig(path)
        paths.append({"file": rel(path), "type": "figure", "description": "Phase 8 metric summary"})
    plt.close(fig)
    return paths


def write_readme(decision: str, gene_audit: dict[str, Any]) -> None:
    text = f"""# BioApp Phase 8

Stage: SVTuner-vs-endpoint evaluation and baseline/SVTuner comparison.

This phase reads frozen outputs from Phase 2C, Phase 5, Phase 6, and Phase 7. It does not rerun SVTuner, Stage3, Stage3B, Stage4, CytoSPACE, or any Phase 4 baseline.

Frozen endpoint:
- CTA-defined Immune cells.
- The endpoint is used only for evaluation.
- The endpoint is not redefined in Phase 8.

Baseline anchor:
- {BASELINE_ANCHOR}

SVTuner condition:
- {SVTUNER_CONDITION}

Gene count audit:
- final_effective_gene_count_for_phase8 = {gene_audit.get('final_effective_gene_count_for_phase8')}
- material_mismatch = {gene_audit.get('material_mismatch')}
- explanation = {gene_audit.get('discrepancy_explanation')}

Scores:
- withheld_score and withheld_binary are imported from Phase 7 raw SVTuner contract.
- binary withheld status uses the predefined Phase 7 Stage3B rule; Phase 8 does not tune thresholds.

Endpoint analysis set:
- Main metrics use CTA Immune-positive and CTA Immune-negative spots.
- Ambiguous and excluded spots remain in spot-level tables but are excluded from main metrics.

Operational definitions:
- Burden reduction is a withheld-aware forced-call burden estimate.
- It is not a post-Stage4 remapping result.

Allowed claims:
- SVTuner raw unsupported / withheld outputs were evaluated against the frozen CTA Immune endpoint.
- Baseline and SVTuner-aware outputs were compared under the immune-all-dropout condition.
- Endpoint-specific comparison metrics were generated for final biological-application audit.

Disallowed claims:
- Biological application completed.
- Biological discovery made.
- SVTuner definitively improves biological interpretation.

Decision:
- {decision}

Next:
- BioApp Phase 9 final biological-application audit, figure selection, and interpretation boundary lock.
"""
    (OUT / "README.md").write_text(text, encoding="utf-8")


def write_manifest(files: list[dict[str, str]]) -> None:
    existing = []
    for path in sorted(OUT.rglob("*")):
        if path.is_file():
            existing.append(
                {
                    "file": rel(path),
                    "type": "output",
                    "description": "Phase 8 output",
                    "created_by_phase": "BioApp Phase 8",
                    "status": "generated",
                    "notes": "",
                }
            )
    extra = []
    for item in files:
        extra.append(
            {
                "file": item["file"],
                "type": item.get("type", "figure"),
                "description": item.get("description", ""),
                "created_by_phase": "BioApp Phase 8",
                "status": "generated",
                "notes": "",
            }
        )
    pd.DataFrame(existing + extra).drop_duplicates(subset=["file"]).to_csv(OUT / "manifest.csv", index=False)


def write_decision(decision: str, summary: dict[str, Any]) -> None:
    m = summary["metric_highlights"]
    next_text = (
        "BioApp Phase 9 - final biological-application audit, figure selection, and interpretation boundary lock"
        if decision == "PASS"
        else (
            "Review score parsing, gene-count audit, endpoint analysis set, or baseline/SVTuner comparison before final audit. Do not write final biological application claim."
            if decision == "REVIEW_REQUIRED"
            else "Stop. Fix boundary violation or invalid evaluation before retrying Phase 8."
        )
    )
    text = f"""BioApp Phase 8 - SVTuner-vs-endpoint evaluation and baseline/SVTuner comparison

Decision: {decision}

Input endpoint:
CTA-defined Immune cells

Endpoint status:
Frozen at Phase 2C

Primary endpoint note:
Sparse immune-positive spatial compartments, not a large continuous ROI.

Baseline anchor:
{BASELINE_ANCHOR}

SVTuner condition:
{SVTUNER_CONDITION}

Analysis universe:
2248 frozen endpoint spots
main analysis set = positive + negative = {summary.get('main_analysis_spots')}

Gene count audit:
final_effective_gene_count_for_phase8 = {summary.get('final_effective_gene_count_for_phase8')}
gene_count_material_mismatch = {summary.get('gene_count_material_mismatch')}

SVTuner-vs-endpoint metrics:
mean_withheld_score_positive = {m.get('mean_withheld_score_positive')}
mean_withheld_score_negative = {m.get('mean_withheld_score_negative')}
withheld_enrichment_delta = {m.get('withheld_enrichment_delta')}
withheld_AUROC = {m.get('withheld_AUROC')}
withheld_AUPRC = {m.get('withheld_AUPRC')}
withheld_binary_positive_rate = {m.get('withheld_binary_positive_rate')}
withheld_binary_negative_rate = {m.get('withheld_binary_negative_rate')}

Baseline/SVTuner comparison:
contradiction_prevention_rate = {m.get('contradiction_prevention_rate')}
mean_prevented_forced_nonimmune_burden_positive_binary = {m.get('mean_prevented_forced_nonimmune_burden_positive_binary')}
mean_prevented_forced_nonimmune_burden_positive_continuous = {m.get('mean_prevented_forced_nonimmune_burden_positive_continuous')}

Spatial comparison:
endpoint_baseline_SVTuner_spatial_comparison_available = {summary.get('endpoint_baseline_SVTuner_spatial_comparison_available')}

Boundary checks:
SVTuner rerun: false
Stage3 rerun: false
Stage4 run: false
CytoSPACE rerun: false
Endpoint redefined: false
Threshold selected using endpoint labels: false
Biological application final claim made: false

Allowed claims:
- SVTuner raw unsupported / withheld outputs were evaluated against the frozen CTA Immune endpoint.
- Baseline and SVTuner-aware outputs were compared under the immune-all-dropout condition.
- Endpoint-specific comparison metrics were generated for final biological-application audit.

Disallowed claims:
- Biological application completed.
- Biological discovery made.
- SVTuner definitively improves biological interpretation.

Next:
{next_text}
"""
    (OUT / "decision.txt").write_text(text, encoding="utf-8")


def main() -> int:
    ensure_dirs()
    errors, context = validate_inputs()
    if errors:
        decision = "FAIL"
        summary = {
            "phase": "BioApp Phase 8 SVTuner-vs-endpoint evaluation and baseline/SVTuner comparison",
            "decision": decision,
            "reason": "input validation failed",
            "errors": errors,
            "SVTuner_rerun": False,
            "Stage3_rerun": False,
            "Stage4_run": False,
            "CytoSPACE_rerun": False,
            "endpoint_redefined": False,
            "biological_application_allowed": False,
        }
        write_json(OUT / "bioapp_phase8_svtuner_vs_endpoint_evaluation_summary.json", summary)
        write_json(OUT / "bioapp_phase8_golden_rules_v2_1_check.json", {"decision": decision, "biological_application_allowed": False, "errors": errors})
        (OUT / "decision.txt").write_text("Decision: FAIL\nreason = input validation failed\n", encoding="utf-8")
        return 1

    gene_audit = gene_count_audit(context)
    tables = load_tables()
    endpoint = tables["endpoint"]
    baseline_anchor = tables["baseline_anchor"]
    contract = tables["contract"]
    endpoint_set = endpoint_qc(endpoint)
    baseline_qc = baseline_import_qc(endpoint, baseline_anchor)
    phase7_qc = phase7_import_qc(endpoint, contract)
    pd.DataFrame(
        [
            {"input": "Phase 2C frozen endpoint", "status": "ok", "notes": context["phase2c"].get("decision")},
            {"input": "Phase 5 baseline anchor", "status": baseline_qc["baseline_anchor_import_status"], "notes": BASELINE_ANCHOR},
            {"input": "Phase 6 planning", "status": "ok", "notes": "threshold selected using endpoint labels = false"},
            {"input": "Phase 7 raw SVTuner output", "status": phase7_qc["phase7_contract_import_status"], "notes": SVTUNER_CONDITION},
            {"input": "gene count audit", "status": "ok" if not gene_audit["material_mismatch"] else "review_required", "notes": gene_audit["discrepancy_explanation"]},
        ]
    ).to_csv(OUT / "phase8_input_qc.csv", index=False)

    scores = build_score_table(endpoint, baseline_anchor, contract)
    withheld_summary = compute_withheld_metrics(scores)
    binary_summary = compute_binary_metrics(scores)
    comp, burden_summary, prevention_summary, interp_summary = compute_burden_and_prevention(scores)
    overlap = compute_spatial_overlap(scores, comp)
    metric_df = metric_comparison(withheld_summary, binary_summary, burden_summary, prevention_summary, interp_summary)

    figure_rows: list[dict[str, str]] = []
    figure_rows += plot_spatial(endpoint, scores, comp)
    figure_rows += plot_box(scores)
    figure_rows += plot_roc_pr()
    figure_rows += plot_binary_contingency()
    figure_rows += plot_burden()
    figure_rows += plot_metric_summary(metric_df)
    pd.DataFrame(figure_rows).to_csv(OUT / "figure_manifest.csv", index=False)

    pass_conditions = [
        not gene_audit["material_mismatch"],
        len(scores) == 2248,
        scores["barcode"].nunique() == 2248,
        scores["withheld_score"].notna().all(),
        withheld_summary.get("withheld_AUROC") is not None,
        withheld_summary.get("withheld_AUPRC") is not None,
        len(figure_rows) >= 12,
        baseline_qc["baseline_anchor_import_status"] == "ok",
        phase7_qc["phase7_contract_import_status"] == "ok",
    ]
    decision = "PASS" if all(pass_conditions) else "REVIEW_REQUIRED"
    ready_for_phase9 = decision == "PASS"
    highlights = {
        "mean_withheld_score_positive": withheld_summary.get("mean_withheld_score_positive"),
        "mean_withheld_score_negative": withheld_summary.get("mean_withheld_score_negative"),
        "withheld_enrichment_delta": withheld_summary.get("delta_mean_positive_minus_negative"),
        "withheld_AUROC": withheld_summary.get("withheld_AUROC"),
        "withheld_AUPRC": withheld_summary.get("withheld_AUPRC"),
        "withheld_binary_positive_rate": binary_summary.get("withheld_binary_positive_rate"),
        "withheld_binary_negative_rate": binary_summary.get("withheld_binary_negative_rate"),
        "contradiction_prevention_rate": prevention_summary.get("contradiction_prevention_rate"),
        "mean_prevented_forced_nonimmune_burden_positive_binary": burden_summary.get("mean_prevented_burden_binary"),
        "mean_prevented_forced_nonimmune_burden_positive_continuous": burden_summary.get("mean_prevented_burden_continuous"),
    }
    summary = {
        "phase": "BioApp Phase 8 SVTuner-vs-endpoint evaluation and baseline/SVTuner comparison",
        "decision": decision,
        "stage_type": "SVTuner-vs-endpoint evaluation and baseline/SVTuner comparison / biological application evidence generation",
        "input_phase2c_decision": context["phase2c"].get("decision"),
        "input_phase5_decision": context["phase5"].get("decision"),
        "input_phase6_decision": context["phase6"].get("decision"),
        "input_phase7_decision": context["phase7"].get("decision"),
        "endpoint_frozen": True,
        "primary_endpoint": "Immune cells",
        "primary_endpoint_spatial_note": "sparse immune-positive spatial compartments, not a large continuous ROI",
        "endpoint_positive_spots": endpoint_set["positive_count"],
        "endpoint_negative_spots": endpoint_set["negative_count"],
        "endpoint_ambiguous_spots": endpoint_set["ambiguous_count"],
        "endpoint_excluded_spots": endpoint_set["excluded_count"],
        "main_analysis_spots": endpoint_set["main_analysis_count"],
        "baseline_anchor_condition": BASELINE_ANCHOR,
        "svtuner_condition": SVTUNER_CONDITION,
        "analysis_universe_type": "full_frozen_endpoint_spot_universe",
        "analysis_universe_spots": 2248,
        "gene_count_audit_completed": True,
        "final_effective_gene_count_for_phase8": gene_audit.get("final_effective_gene_count_for_phase8"),
        "gene_count_material_mismatch": gene_audit.get("material_mismatch"),
        "SVTuner_rerun": False,
        "Stage3_rerun": False,
        "Stage4_run": False,
        "CytoSPACE_rerun": False,
        "endpoint_redefined": False,
        "threshold_selected_using_endpoint_labels": False,
        "phase7_raw_output_used": True,
        "spot_level_contract_rows": int(len(contract)),
        "contract_spot_universe_aligned": phase7_qc["phase7_contract_import_status"] == "ok",
        "withheld_enrichment_computed": True,
        "withheld_AUROC_AUPRC_computed": withheld_summary.get("withheld_AUROC") is not None,
        "binary_withheld_enrichment_computed": True,
        "endpoint_positive_nonimmune_assignment_reduction_computed": True,
        "contradiction_prevention_rate_computed": True,
        "endpoint_concordant_interpretation_score_computed": True,
        "baseline_svtuner_comparison_computed": True,
        "metric_highlights": highlights,
        "endpoint_baseline_SVTuner_spatial_comparison_available": bool(len(figure_rows) > 0),
        "biological_application_evidence_chain_candidate_complete": ready_for_phase9,
        "biological_application_final_claim_made": False,
        "biological_application_allowed": False,
        "ready_for_phase9_final_biological_application_audit": ready_for_phase9,
        "allowed_claims": [
            "SVTuner raw unsupported / withheld outputs were evaluated against the frozen CTA Immune endpoint.",
            "Baseline and SVTuner-aware outputs were compared under the immune-all-dropout condition.",
            "Endpoint-specific comparison metrics were generated for final biological-application audit.",
        ],
        "disallowed_claims": [
            "Biological application completed.",
            "Biological discovery made.",
            "SVTuner definitively improves biological interpretation.",
        ],
    }
    golden = {
        "stage_type": summary["stage_type"],
        "biological_question_defined": True,
        "external_endpoint_predefined": True,
        "endpoint_independent_from_SVTuner": True,
        "endpoint_spatially_registered": True,
        "baseline_comparison_available": True,
        "endpoint_specific_improvement_defined": True,
        "endpoint_specific_quantitative_metric_available": True,
        "endpoint_baseline_SVTuner_spatial_comparison_available": summary["endpoint_baseline_SVTuner_spatial_comparison_available"],
        "interpretation_boundary_defined": True,
        "biological_application_evidence_chain_candidate_complete": summary["biological_application_evidence_chain_candidate_complete"],
        "final_biological_application_audit_required": True,
        "biological_application_allowed": False,
        "allowed_claim_level": "SVTuner-vs-endpoint evaluation and baseline/SVTuner comparison only",
        "decision": decision,
    }
    write_json(OUT / "bioapp_phase8_svtuner_vs_endpoint_evaluation_summary.json", summary)
    write_json(OUT / "bioapp_phase8_golden_rules_v2_1_check.json", golden)
    write_decision(decision, summary)
    write_readme(decision, gene_audit)
    write_manifest(figure_rows)

    print("BioApp Phase 8 completed.")
    print()
    print("Decision:")
    print(decision)
    print()
    print("Endpoint:")
    print("CTA-defined Immune cells")
    print()
    print("Endpoint frozen:")
    print("true")
    print()
    print("Baseline anchor:")
    print(BASELINE_ANCHOR)
    print()
    print("SVTuner condition:")
    print(SVTUNER_CONDITION)
    print()
    print("Analysis universe:")
    print("2248 frozen endpoint spots")
    print(f"main analysis spots = {endpoint_set['main_analysis_count']}")
    print()
    print("Gene count audit:")
    print(f"final_effective_gene_count_for_phase8 = {gene_audit.get('final_effective_gene_count_for_phase8')}")
    print(f"material_mismatch = {gene_audit.get('material_mismatch')}")
    print()
    print("Key SVTuner-vs-endpoint metrics:")
    print(f"withheld_AUROC = {highlights['withheld_AUROC']}")
    print(f"withheld_AUPRC = {highlights['withheld_AUPRC']}")
    print(f"withheld_enrichment_delta = {highlights['withheld_enrichment_delta']}")
    print(f"withheld_binary_positive_rate = {highlights['withheld_binary_positive_rate']}")
    print(f"withheld_binary_negative_rate = {highlights['withheld_binary_negative_rate']}")
    print()
    print("Baseline/SVTuner comparison:")
    print(f"contradiction_prevention_rate = {highlights['contradiction_prevention_rate']}")
    print(f"mean_prevented_forced_nonimmune_burden_positive_binary = {highlights['mean_prevented_forced_nonimmune_burden_positive_binary']}")
    print(f"mean_prevented_forced_nonimmune_burden_positive_continuous = {highlights['mean_prevented_forced_nonimmune_burden_positive_continuous']}")
    print()
    print("Boundary checks:")
    print("SVTuner rerun = false")
    print("Stage3 rerun = false")
    print("Stage4 run = false")
    print("CytoSPACE rerun = false")
    print("Endpoint redefined = false")
    print("Threshold selected using endpoint labels = false")
    print("Biological application final claim made = false")
    print()
    print("Ready for Phase 9 final biological-application audit:")
    print(str(ready_for_phase9).lower())
    print()
    print("Next:")
    print(
        "BioApp Phase 9 - final biological-application audit, figure selection, and interpretation boundary lock"
        if decision == "PASS"
        else "Review score parsing, gene-count audit, endpoint analysis set, or baseline/SVTuner comparison before final audit."
    )
    return 0 if decision in {"PASS", "REVIEW_REQUIRED"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
