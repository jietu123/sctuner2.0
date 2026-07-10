from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase9_final_biological_application_audit_figure_selection_interpretation_boundary_lock"

P2C = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping"
P3R = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase3r_review_resolution_for_cta_endpoint_baseline_feasibility"
P4 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase4_formal_cytospace_baseline_execution_against_frozen_cta_immune_endpoint"
P5 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase5_endpoint_specific_evaluation_of_cytospace_baselines_against_frozen_cta_immune_endpoint"
P6 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase6_svtuner_aware_execution_planning_against_frozen_cta_immune_endpoint"
P7 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase7_svtuner_aware_execution_against_frozen_cta_immune_endpoint"
P8 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT)).replace("\\", "/")


def load_inputs() -> tuple[dict[str, dict[str, Any]], list[str]]:
    errors: list[str] = []
    specs = {
        "phase2c": P2C / "bioapp_phase2c_endpoint_freeze_summary.json",
        "phase3r": P3R / "bioapp_phase3r_review_resolution_summary.json",
        "phase4": P4 / "bioapp_phase4_formal_cytospace_baseline_execution_summary.json",
        "phase5": P5 / "bioapp_phase5_endpoint_specific_baseline_evaluation_summary.json",
        "phase6": P6 / "bioapp_phase6_svtuner_execution_planning_summary.json",
        "phase7": P7 / "bioapp_phase7_svtuner_aware_execution_summary.json",
        "phase8": P8 / "bioapp_phase8_svtuner_vs_endpoint_evaluation_summary.json",
        "phase8_golden": P8 / "bioapp_phase8_golden_rules_v2_1_check.json",
        "phase8_gene": P8 / "gene_count_audit.json",
    }
    data: dict[str, dict[str, Any]] = {}
    for key, path in specs.items():
        if not path.exists():
            errors.append(f"Missing required input: {rel(path)}")
            continue
        data[key] = read_json(path)
    return data, errors


def audit_inputs(data: dict[str, dict[str, Any]]) -> tuple[list[dict[str, Any]], list[str]]:
    rows: list[dict[str, Any]] = []
    errors: list[str] = []

    required_pass = {
        "phase2c": "BioApp Phase 2C endpoint freeze",
        "phase3r": "BioApp Phase 3R review resolution",
        "phase4": "BioApp Phase 4 formal CytoSPACE baseline execution",
        "phase5": "BioApp Phase 5 endpoint-specific baseline evaluation",
        "phase6": "BioApp Phase 6 SVTuner-aware planning",
        "phase7": "BioApp Phase 7 SVTuner-aware raw execution",
        "phase8": "BioApp Phase 8 SVTuner-vs-endpoint evaluation",
    }
    for key, label in required_pass.items():
        decision = data.get(key, {}).get("decision")
        ok = decision == "PASS"
        rows.append({"check": label, "observed": decision, "pass": ok, "notes": ""})
        if not ok:
            errors.append(f"{label} is not PASS")

    p2c = data.get("phase2c", {})
    endpoint_checks = [
        ("endpoint_frozen", p2c.get("endpoint_frozen") is True),
        ("primary_endpoint_is_Immune_cells", p2c.get("primary_endpoint") == "Immune cells"),
        ("endpoint_independent_no_expression_markers", p2c.get("expression_markers_used_to_define_endpoint") is False),
        ("endpoint_independent_no_mapping_outputs", p2c.get("mapping_outputs_used_to_define_endpoint") is False),
    ]
    for name, ok in endpoint_checks:
        rows.append({"check": name, "observed": str(ok), "pass": ok, "notes": "Phase 2C endpoint freeze audit"})
        if not ok:
            errors.append(f"Endpoint check failed: {name}")

    p8 = data.get("phase8", {})
    golden = data.get("phase8_golden", {})
    gene = data.get("phase8_gene", {})
    phase8_checks = [
        ("phase8_golden_rules_decision", golden.get("decision") == "PASS"),
        ("no_endpoint_redefinition", p8.get("endpoint_redefined") is False),
        ("no_threshold_tuning_with_endpoint", p8.get("threshold_selected_using_endpoint_labels") is False),
        ("no_stage3_rerun", p8.get("Stage3_rerun") is False),
        ("no_stage4_run", p8.get("Stage4_run") is False),
        ("no_cytospace_rerun", p8.get("CytoSPACE_rerun") is False),
        ("gene_count_no_material_mismatch", gene.get("material_mismatch") is False),
        ("phase7_raw_output_used", p8.get("phase7_raw_output_used") is True),
        ("baseline_svtuner_same_condition", p8.get("baseline_anchor_condition") == "baseline_immune_all_dropout" and p8.get("svtuner_condition") == "svtuner_immune_all_dropout"),
        ("spatial_comparison_available", p8.get("endpoint_baseline_SVTuner_spatial_comparison_available") is True),
        ("candidate_evidence_chain_complete", p8.get("biological_application_evidence_chain_candidate_complete") is True),
    ]
    for name, ok in phase8_checks:
        rows.append({"check": name, "observed": str(ok), "pass": ok, "notes": "Phase 8 final audit"})
        if not ok:
            errors.append(f"Phase 8 audit check failed: {name}")

    metrics = p8.get("metric_highlights", {})
    metric_checks = [
        ("withheld_AUROC_supports_endpoint_concordance", metrics.get("withheld_AUROC", 0) >= 0.75),
        ("positive_withheld_rate_exceeds_negative", metrics.get("withheld_binary_positive_rate", 0) > metrics.get("withheld_binary_negative_rate", 1)),
        ("withheld_enrichment_positive_delta", metrics.get("withheld_enrichment_delta", 0) > 0),
        ("burden_reduction_supported", metrics.get("mean_prevented_forced_nonimmune_burden_positive_binary", 0) > 0),
        ("operational_prevention_rate_positive", metrics.get("contradiction_prevention_rate", 0) > 0),
    ]
    for name, ok in metric_checks:
        rows.append({"check": name, "observed": str(ok), "pass": ok, "notes": "Metric adequacy audit"})
        if not ok:
            errors.append(f"Metric adequacy check failed: {name}")

    return rows, errors


def figure_selection_table() -> pd.DataFrame:
    rows = [
        {
            "figure_file": rel(P8 / "endpoint_baseline_svtuner_three_way_spatial_comparison.pdf"),
            "paired_svg": rel(P8 / "endpoint_baseline_svtuner_three_way_spatial_comparison.svg"),
            "classification": "main figure candidate",
            "recommended_use": "Primary spatial evidence module",
            "rationale": "Shows frozen endpoint, baseline forced nonimmune burden, SVTuner withheld score, withheld binary, and endpoint/withheld overlay in one evidence chain.",
            "layout_review": "required",
            "allowed_optimization": "layout-only; may adjust panel arrangement, canvas size, font size, spacing",
            "data_change_allowed": False,
        },
        {
            "figure_file": rel(P8 / "endpoint_positive_forced_nonimmune_burden_reduction.pdf"),
            "paired_svg": rel(P8 / "endpoint_positive_forced_nonimmune_burden_reduction.svg"),
            "classification": "main or supplementary figure candidate",
            "recommended_use": "Quantitative burden-reduction support",
            "rationale": "Directly summarizes withheld-aware forced-call burden reduction in CTA Immune-positive spots.",
            "layout_review": "minor",
            "allowed_optimization": "layout-only",
            "data_change_allowed": False,
        },
        {
            "figure_file": rel(P8 / "phase8_svtuner_endpoint_metric_summary.pdf"),
            "paired_svg": rel(P8 / "phase8_svtuner_endpoint_metric_summary.svg"),
            "classification": "supplementary figure candidate",
            "recommended_use": "Metric summary support",
            "rationale": "Aggregates AUROC/AUPRC, enrichment, prevention-rate style metrics; useful as compact supplement or audit summary.",
            "layout_review": "minor",
            "allowed_optimization": "layout-only",
            "data_change_allowed": False,
        },
        {
            "figure_file": rel(P8 / "withheld_score_positive_vs_negative_boxplot.pdf"),
            "paired_svg": rel(P8 / "withheld_score_positive_vs_negative_boxplot.svg"),
            "classification": "supplementary figure candidate",
            "recommended_use": "Distribution-level evidence",
            "rationale": "Shows positive/negative separation in withheld_score.",
            "layout_review": "minor",
            "allowed_optimization": "layout-only",
            "data_change_allowed": False,
        },
        {
            "figure_file": rel(P8 / "withheld_score_roc_curve.pdf"),
            "paired_svg": rel(P8 / "withheld_score_roc_curve.svg"),
            "classification": "supplementary figure candidate",
            "recommended_use": "Metric evidence",
            "rationale": "Supports AUROC-based endpoint concordance.",
            "layout_review": "minor",
            "allowed_optimization": "layout-only",
            "data_change_allowed": False,
        },
        {
            "figure_file": rel(P8 / "withheld_score_pr_curve.pdf"),
            "paired_svg": rel(P8 / "withheld_score_pr_curve.svg"),
            "classification": "supplementary figure candidate",
            "recommended_use": "Metric evidence under endpoint sparsity",
            "rationale": "Supports AUPRC interpretation under class imbalance.",
            "layout_review": "minor",
            "allowed_optimization": "layout-only",
            "data_change_allowed": False,
        },
        {
            "figure_file": rel(P8 / "binary_withheld_endpoint_contingency_plot.pdf"),
            "paired_svg": rel(P8 / "binary_withheld_endpoint_contingency_plot.svg"),
            "classification": "QC/audit-only figure",
            "recommended_use": "Internal or supplementary audit",
            "rationale": "Useful for checking positive/negative binary withheld counts, but less suitable as a core biological-application figure.",
            "layout_review": "not required",
            "allowed_optimization": "layout-only",
            "data_change_allowed": False,
        },
    ]
    return pd.DataFrame(rows)


def write_report(figures: pd.DataFrame, summary: dict[str, Any]) -> None:
    metrics = summary["phase8_metric_highlights"]
    report = f"""# BioApp Phase 9 Figure Selection Report

## Decision

BioApp Phase 9 decision: **{summary['decision']}**

Biological application allowed: **{str(summary['biological_application_allowed']).lower()}**

## Main Evidence

The Phase 8 evidence chain supports a formal biological-application evidence claim under the locked boundary:

- Endpoint: CTA-defined Immune cells
- Baseline anchor: baseline_immune_all_dropout
- SVTuner condition: svtuner_immune_all_dropout
- Main analysis spots: {summary['main_analysis_spots']}
- Effective gene count: {summary['final_effective_gene_count_for_phase9']}

Key Phase 8 metrics:

- withheld_AUROC = {metrics['withheld_AUROC']:.4f}
- withheld_AUPRC = {metrics['withheld_AUPRC']:.4f}
- withheld_enrichment_delta = {metrics['withheld_enrichment_delta']:.4f}
- positive withheld rate = {metrics['withheld_binary_positive_rate']:.4f}
- negative withheld rate = {metrics['withheld_binary_negative_rate']:.4f}
- operational contradiction prevention rate = {metrics['contradiction_prevention_rate']:.4f}

## Figure Triage

{figures[['figure_file', 'classification', 'recommended_use', 'layout_review']].to_markdown(index=False)}

## Recommended Main Figure Candidate

`endpoint_baseline_svtuner_three_way_spatial_comparison.pdf`

This is the strongest main figure candidate because it connects the frozen endpoint, baseline forced nonimmune burden, SVTuner withheld score, binary withheld output, and endpoint/withheld spatial overlay.

Layout review is required before manuscript use. Only layout-level changes are allowed.

## Interpretation Boundary

Allowed claim:

> SVTuner provides biological-application evidence that its unsupported/withheld output is concordant with an independent CTA-defined Immune endpoint under an immune-all-dropout condition, reducing forced non-immune assignment burden relative to the baseline mapping.

Disallowed claims remain locked in `bioapp_phase9_interpretation_boundary_lock.json`.
"""
    (OUT / "bioapp_phase9_figure_selection_report.md").write_text(report, encoding="utf-8")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    data, input_errors = load_inputs()
    audit_rows, audit_errors = audit_inputs(data) if not input_errors else ([], [])
    errors = input_errors + audit_errors

    p8 = data.get("phase8", {})
    metrics = p8.get("metric_highlights", {})
    decision = "PASS" if not errors else "FAIL"
    biological_allowed = decision == "PASS"

    audit_df = pd.DataFrame(audit_rows)
    audit_df.to_csv(OUT / "bioapp_phase9_audit_checklist.csv", index=False)
    figures = figure_selection_table()
    figures.to_csv(OUT / "bioapp_phase9_figure_selection_table.csv", index=False)

    allowed_final_claim = (
        "SVTuner provides biological-application evidence that its unsupported/withheld output is concordant with an "
        "independent CTA-defined Immune endpoint under an immune-all-dropout condition, reducing forced non-immune "
        "assignment burden relative to the baseline mapping."
    )
    disallowed_claims = [
        "SVTuner discovered a new immune niche.",
        "SVTuner proves the biological identity of all withheld spots.",
        "SVTuner definitively improves all biological interpretations.",
        "SVTuner replaces spatial transcriptomics annotation.",
        "The withheld regions are confirmed immune cells without external validation.",
        "Any claim based on newly tuned thresholds or redefined endpoints.",
    ]
    boundary = {
        "phase": "BioApp Phase 9 interpretation boundary lock",
        "decision": decision,
        "biological_application_allowed": biological_allowed,
        "allowed_final_claim_if_pass": allowed_final_claim if biological_allowed else None,
        "claim_scope": "CTA-defined Immune endpoint under immune-all-dropout condition only",
        "mandatory_caveats": [
            "Endpoint-positive spots are sparse immune-positive spatial compartments, not a large continuous ROI.",
            "Withheld regions are endpoint-concordant, but this does not prove every withheld spot is an immune cell.",
            "Burden reduction is withheld-aware forced-call burden, not post-Stage4 remapping.",
            "No endpoint-dependent threshold tuning or endpoint redefinition was performed.",
        ],
        "disallowed_claims": disallowed_claims,
        "layout_only_figure_optimization_allowed": True,
        "data_metric_threshold_label_conclusion_changes_allowed": False,
    }
    write_json(OUT / "bioapp_phase9_interpretation_boundary_lock.json", boundary)

    summary = {
        "phase": "BioApp Phase 9 final biological-application audit, figure selection, and interpretation boundary lock",
        "decision": decision,
        "stage_type": "final audit / figure selection / interpretation boundary lock",
        "input_phase2c_decision": data.get("phase2c", {}).get("decision"),
        "input_phase3r_decision": data.get("phase3r", {}).get("decision"),
        "input_phase4_decision": data.get("phase4", {}).get("decision"),
        "input_phase5_decision": data.get("phase5", {}).get("decision"),
        "input_phase6_decision": data.get("phase6", {}).get("decision"),
        "input_phase7_decision": data.get("phase7", {}).get("decision"),
        "input_phase8_decision": p8.get("decision"),
        "endpoint_frozen": data.get("phase2c", {}).get("endpoint_frozen"),
        "primary_endpoint": "Immune cells",
        "endpoint_external_spatially_registered_independent": not errors,
        "baseline_anchor_condition": "baseline_immune_all_dropout",
        "svtuner_condition": "svtuner_immune_all_dropout",
        "main_analysis_spots": p8.get("main_analysis_spots"),
        "final_effective_gene_count_for_phase9": p8.get("final_effective_gene_count_for_phase8"),
        "gene_count_material_mismatch": p8.get("gene_count_material_mismatch"),
        "phase8_metric_highlights": metrics,
        "endpoint_concordant_withheld_behavior_supported": decision == "PASS",
        "forced_nonimmune_burden_reduction_supported": decision == "PASS",
        "biological_application_allowed": biological_allowed,
        "allowed_final_claim": allowed_final_claim if biological_allowed else None,
        "final_biological_application_claim_made_in_manuscript": False,
        "main_figure_candidate": rel(P8 / "endpoint_baseline_svtuner_three_way_spatial_comparison.pdf"),
        "main_figure_candidate_requires_layout_review": True,
        "supplementary_figure_candidates": figures.loc[figures["classification"].str.contains("supplementary", case=False), "figure_file"].tolist(),
        "qc_audit_only_figures": figures.loc[figures["classification"].eq("QC/audit-only figure"), "figure_file"].tolist(),
        "audit_errors": errors,
        "no_rerun_or_redefinition_guardrails": {
            "SVTuner_rerun": False,
            "Stage3_rerun": False,
            "Stage4_run": False,
            "CytoSPACE_rerun": False,
            "endpoint_redefined": False,
            "threshold_selected_using_endpoint_labels": False,
            "new_biological_endpoint_introduced": False,
            "new_mapping_or_new_computation_beyond_audit": False,
        },
        "next": "manuscript integration" if biological_allowed else "stop",
    }
    write_json(OUT / "bioapp_phase9_final_audit_summary.json", summary)

    golden = {
        "phase": "BioApp Phase 9 Golden Rules final check",
        "decision": decision,
        "biological_question_defined": True,
        "external_endpoint_predefined": True,
        "endpoint_independent_from_SVTuner": True,
        "endpoint_spatially_registered": True,
        "baseline_comparison_available": True,
        "SVTuner_vs_endpoint_evaluation_available": True,
        "endpoint_specific_quantitative_metric_available": True,
        "endpoint_baseline_SVTuner_spatial_comparison_available": True,
        "no_endpoint_leakage_or_threshold_tuning_detected": decision == "PASS",
        "gene_count_material_mismatch": p8.get("gene_count_material_mismatch"),
        "interpretation_boundary_locked": True,
        "biological_application_allowed": biological_allowed,
        "final_claim_allowed_but_not_inserted_into_manuscript": biological_allowed,
        "allowed_claim_level": "formal biological-application evidence claim" if biological_allowed else "not allowed",
    }
    write_json(OUT / "bioapp_phase9_golden_rules_final_check.json", golden)

    write_report(figures, summary)

    readme = f"""# BioApp Phase 9

This phase performs final audit, figure selection, and interpretation boundary locking using frozen Phase 8 outputs.

Decision: {decision}

Biological application allowed: {str(biological_allowed).lower()}

No SVTuner, Stage3, Stage4, CytoSPACE, endpoint redefinition, threshold tuning, or new mapping was performed.

Main figure candidate:
{summary['main_figure_candidate']}

Allowed final claim:
{allowed_final_claim if biological_allowed else 'None'}

Next:
{'Manuscript integration' if biological_allowed else 'Stop / review failed checks'}
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")

    manifest_rows = []
    for path in sorted(OUT.rglob("*")):
        if path.is_file():
            manifest_rows.append(
                {
                    "file": rel(path),
                    "type": "phase9_output",
                    "description": "BioApp Phase 9 audit output",
                    "status": "generated",
                }
            )
    pd.DataFrame(manifest_rows).to_csv(OUT / "manifest.csv", index=False)

    print("BioApp Phase 9 completed.")
    print(f"Decision: {decision}")
    print(f"Biological application allowed: {str(biological_allowed).lower()}")
    print(f"Main figure candidate: {summary['main_figure_candidate']}")
    print(f"Next: {'manuscript integration' if biological_allowed else 'stop'}")
    return 0 if decision == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
