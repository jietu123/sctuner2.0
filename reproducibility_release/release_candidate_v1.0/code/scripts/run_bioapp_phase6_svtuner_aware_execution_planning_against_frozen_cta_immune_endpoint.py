#!/usr/bin/env python3
"""BioApp Phase 6: SVTuner-aware execution planning only.

This phase creates planning artifacts. It never runs SVTuner, Stage3, Stage4,
CytoSPACE, prevention analysis, or SVTuner-vs-endpoint metrics.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase6_svtuner_aware_execution_planning_against_frozen_cta_immune_endpoint"
P2C = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping"
P3R = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase3r_review_resolution_for_cta_endpoint_baseline_feasibility"
P4 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase4_formal_cytospace_baseline_execution_against_frozen_cta_immune_endpoint"
P5 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase5_endpoint_specific_evaluation_of_cytospace_baselines_against_frozen_cta_immune_endpoint"

PHASE = "BioApp Phase 6 SVTuner-aware execution planning against frozen CTA Immune endpoint"
STAGE_TYPE = "SVTuner-aware execution planning / biological application candidate preparation"
PRIMARY_ENDPOINT = "Immune cells"
PRIMARY_ENDPOINT_NOTE = "sparse immune-positive spatial compartments, not a large continuous ROI"
ALLOWED_CLAIMS = [
    "SVTuner-aware execution was planned against the frozen CTA Immune endpoint.",
    "The primary SVTuner execution condition and output contract were frozen.",
    "Future SVTuner-vs-endpoint metrics were aligned to Phase 5 baseline metrics.",
]
DISALLOWED_CLAIMS = [
    "SVTuner improves endpoint recovery.",
    "SVTuner prevents contradicted interpretation.",
    "Biological application completed.",
    "Biological discovery made.",
]


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return str(path).replace("\\", "/")


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def write_json(path: Path, data: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, ensure_ascii=False)
        fh.write("\n")


def rows_cols(path: Path) -> tuple[int | None, int | None, bool]:
    if not path.exists() or not path.is_file():
        return None, None, False
    try:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
            reader = csv.reader(fh)
            header = next(reader)
            n = sum(1 for _ in reader)
        return n, len(header), True
    except Exception:
        return None, None, False


def validate_inputs() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any], list[str]]:
    errors: list[str] = []
    p2c = read_json(P2C / "bioapp_phase2c_endpoint_freeze_summary.json")
    p3r = read_json(P3R / "bioapp_phase3r_review_resolution_summary.json")
    p4 = read_json(P4 / "bioapp_phase4_formal_cytospace_baseline_execution_summary.json")
    p5 = read_json(P5 / "bioapp_phase5_endpoint_specific_baseline_evaluation_summary.json")
    checks = [
        (p2c.get("decision") == "PASS", "Phase 2C decision is not PASS"),
        (p2c.get("endpoint_frozen") is True, "Phase 2C endpoint_frozen is not true"),
        (p2c.get("primary_endpoint") == PRIMARY_ENDPOINT, "Phase 2C primary endpoint is not Immune cells"),
        (p3r.get("decision") == "PASS", "Phase 3R decision is not PASS"),
        (p3r.get("full_2248_coverage_achieved") is True, "Phase 3R full coverage missing"),
        (p3r.get("analysis_universe_type") == "full_frozen_endpoint_spot_universe", "Phase 3R universe invalid"),
        (int(p3r.get("analysis_universe_spots", -1)) == 2248, "Phase 3R spot count invalid"),
        (int(p3r.get("n_gene_overlap_exact", -1)) == 2000, "Phase 3R gene overlap must be 2000"),
        (p4.get("decision") == "PASS", "Phase 4 decision is not PASS"),
        (p4.get("all_required_baseline_runs_completed") is True, "Phase 4 baseline runs incomplete"),
        (p4.get("ready_for_phase5_endpoint_evaluation") is True, "Phase 4 not ready for Phase 5"),
        (p4.get("SVTuner_run") is False, "Phase 4 SVTuner_run must be false"),
        (p4.get("Stage4_run") is False, "Phase 4 Stage4_run must be false"),
        (p5.get("decision") == "PASS", "Phase 5 decision is not PASS"),
        (p5.get("baseline_vs_endpoint_metric_computed") is True, "Phase 5 metric not computed"),
        (p5.get("ready_for_phase6_svtuner_execution") is True, "Phase 5 not ready for Phase 6"),
        (p5.get("biological_application_allowed") is False, "Phase 5 biological_application_allowed must be false"),
        (p5.get("SVTuner_run") is False, "Phase 5 SVTuner_run must be false"),
        (p5.get("Stage4_run") is False, "Phase 5 Stage4_run must be false"),
        (p5.get("prevention_analysis_run") is False, "Phase 5 prevention_analysis_run must be false"),
        (p5.get("endpoint_redefined") is False, "Phase 5 endpoint_redefined must be false"),
    ]
    for ok, msg in checks:
        if not ok:
            errors.append(msg)
    return p2c, p3r, p4, p5, errors


def input_manifest() -> tuple[list[dict[str, Any]], bool]:
    items = [
        (P3R / "formal_gene_intersection.txt", "formal_gene_intersection", True, "Phase 3R"),
        (P3R / "st_expression_full_2248_data.csv.gz", "ST_expression_input", True, "Phase 3R"),
        (P3R / "st_coordinates_full_2248.csv", "ST_coordinates", True, "Phase 3R"),
        (P3R / "cytospace_input_manifest.json", "harmonized_input_manifest", True, "Phase 3R"),
        (P3R / "immune_all_dropout_design.json", "immune_all_dropout_reference_design", True, "Phase 3R"),
        (P4 / "baseline_immune_all_dropout" / "reference_metadata_used.csv", "immune_all_dropout_reference_metadata", True, "Phase 4"),
        (P4 / "baseline_immune_all_dropout" / "inputs" / "reference_expression_used.csv.gz", "immune_all_dropout_reference_expression", True, "Phase 4"),
        (P4 / "baseline_immune_all_dropout" / "cytospace_output", "baseline_anchor_output_dir", True, "Phase 4"),
        (P5 / "baseline_endpoint_metric_summary.csv", "phase5_baseline_metric_table", True, "Phase 5"),
        (P2C / "spot_level_endpoint_freeze.csv", "endpoint_freeze_table", True, "Phase 2C"),
        (P2C / "primary_endpoint_positive_spots.txt", "primary_endpoint_positive_spots", True, "Phase 2C"),
        (P2C / "primary_endpoint_negative_spots.txt", "primary_endpoint_negative_spots", True, "Phase 2C"),
        (P2C / "primary_endpoint_ambiguous_spots.txt", "primary_endpoint_ambiguous_spots", True, "Phase 2C"),
        (P2C / "primary_endpoint_excluded_spots.txt", "primary_endpoint_excluded_spots", True, "Phase 2C"),
    ]
    rows = []
    complete = True
    for path, role, required, source in items:
        n_rows, n_cols, readable = rows_cols(path) if path.suffix.lower() in {".csv", ".txt"} else (None, None, path.exists())
        exists = path.exists()
        if required and not exists:
            complete = False
        rows.append(
            {
                "file_path": rel(path),
                "role": role,
                "exists": exists,
                "required_for_phase7": required,
                "source_phase": source,
                "read_success": readable,
                "n_rows": n_rows,
                "n_cols": n_cols,
                "notes": "",
            }
        )
    return rows, complete


def entrypoint_inventory() -> tuple[pd.DataFrame, bool, bool]:
    candidates = [
        ROOT / "src" / "svtuner" / "run.py",
        ROOT / "src" / "svtuner" / "cli.py",
        ROOT / "src" / "stages" / "stage3_type_plugin.py",
        ROOT / "src" / "stages" / "stage3b_st_unsupported.py",
        ROOT / "src" / "stages" / "stage4_cytospace.py",
        ROOT / "scripts" / "run_project_mainline.py",
    ]
    rows = []
    for path in candidates:
        text = path.read_text(encoding="utf-8", errors="ignore") if path.exists() else ""
        lower = text.lower()
        contains_stage3 = "stage3" in lower
        contains_stage4 = "stage4" in lower
        contains_cytospace = "cytospace" in lower
        rows.append(
            {
                "file_path": rel(path),
                "candidate_entrypoint_type": "cli_or_stage_module" if path.exists() else "missing",
                "contains_stage3": contains_stage3,
                "contains_stage4": contains_stage4,
                "contains_cytospace_call": contains_cytospace,
                "baseline_only_mode_possible": "filter_mode" in lower and "none" in lower,
                "svtuner_mode_possible": path.exists() and ("stage3" in lower or "unsupported" in lower or "svtuner" in lower),
                "notes": "candidate only; not executed in Phase 6",
            }
        )
    df = pd.DataFrame(rows)
    found = bool(df["svtuner_mode_possible"].any())
    explicit_main_entrypoint = (ROOT / "src" / "svtuner" / "run.py").exists()
    explicit_stage3_module = (ROOT / "src" / "stages" / "stage3_type_plugin.py").exists()
    explicit_stage4_module = (ROOT / "src" / "stages" / "stage4_cytospace.py").exists()
    ambiguous = not (explicit_main_entrypoint and explicit_stage3_module and explicit_stage4_module)
    return df, found, ambiguous


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    p2c, p3r, p4, p5, errors = validate_inputs()

    condition_selection = {
        "primary_svtuner_condition": "baseline_immune_all_dropout",
        "selection_rule": "predefined immune-reference incompleteness scenario",
        "selected_based_on_SVTuner_output": False,
        "selected_based_on_future_mapping_output": False,
        "selected_based_on_endpoint_redefinition": False,
        "baseline_anchor_phase": "Phase 5",
        "baseline_anchor_run": "baseline_immune_all_dropout",
    }
    write_json(OUT_DIR / "primary_svtuner_condition_selection.json", condition_selection)

    condition_plan = {
        "required_future_svtuner_runs": ["svtuner_immune_all_dropout"],
        "optional_future_svtuner_runs": ["svtuner_full_reference_context", "svtuner_nonimmune_all_dropout_control_context"],
        "deferred_sensitivity_runs": ["random_all_cell_size_matched_sensitivity"],
        "main_condition": "immune-all-dropout",
        "main_condition_reason": "frozen CTA Immune endpoint tests immune reference incompleteness",
    }
    write_json(OUT_DIR / "svtuner_condition_plan.json", condition_plan)

    manifest_rows, manifest_complete = input_manifest()
    write_json(OUT_DIR / "svtuner_future_input_manifest.json", {"items": manifest_rows})
    pd.DataFrame(manifest_rows).to_csv(OUT_DIR / "svtuner_future_input_manifest.csv", index=False)

    output_fields = [
        "barcode",
        "primary_endpoint_status",
        "baseline_anchor_run",
        "svtuner_condition",
        "reference_supported_status",
        "reference_unrepresented_score",
        "withheld_score",
        "withheld_binary",
        "supported_assignment_score",
        "svtuner_immune_score",
        "svtuner_nonimmune_score",
        "svtuner_dominant_label",
        "svtuner_dominant_fraction",
        "baseline_immune_score_from_phase5",
        "baseline_nonimmune_score_from_phase5",
        "baseline_dominant_label_from_phase5",
        "baseline_dominant_fraction_from_phase5",
        "score_parse_status",
        "svtuner_output_available",
    ]
    output_contract = {
        "required_spot_level_fields": output_fields,
        "primary_condition": "svtuner_immune_all_dropout",
        "baseline_anchor": "baseline_immune_all_dropout",
        "withheld_binary_threshold_note": "Threshold must be fixed in Phase 7 before evaluation and cannot be chosen using endpoint performance.",
        "immune_all_dropout_note": "For immune-all-dropout, immune labels may be unavailable; primary endpoint alignment should use withheld/reference-unrepresented enrichment.",
    }
    write_json(OUT_DIR / "svtuner_output_contract.json", output_contract)
    (OUT_DIR / "svtuner_output_contract.md").write_text(
        "# SVTuner Output Contract\n\n"
        "Phase 7 must produce a spot-level table with these fields:\n\n"
        + "\n".join(f"- `{f}`" for f in output_fields)
        + "\n\nBaseline fields must be imported from Phase 5 without recomputation or redefinition.\n",
        encoding="utf-8",
    )

    metric_plan = {
        "endpoint_analysis_set": {
            "positive": "CTA Immune-positive spots",
            "negative": "CTA Immune-negative spots",
            "ambiguous": "excluded from main metrics",
            "excluded": "excluded from main metrics",
            "main_analysis_set": "positive + negative",
        },
        "baseline_anchor": "baseline_immune_all_dropout",
        "planned_metrics_not_computed_in_phase6": [
            "withheld enrichment in CTA Immune-positive spots",
            "withheld AUROC / AUPRC",
            "endpoint-positive nonimmune assignment reduction",
            "contradiction prevention rate",
            "endpoint-concordant interpretation score",
        ],
        "phase6_computation": "none",
    }
    write_json(OUT_DIR / "svtuner_metric_alignment_plan.json", metric_plan)
    (OUT_DIR / "svtuner_metric_alignment_plan.md").write_text(
        "# SVTuner Metric Alignment Plan\n\n"
        "Future SVTuner-vs-endpoint metrics must use the same endpoint analysis set as Phase 5.\n\n"
        "Baseline anchor: `baseline_immune_all_dropout`.\n\n"
        "Phase 6 does not compute any SVTuner metric.\n",
        encoding="utf-8",
    )

    threshold_policy = {
        "primary_metric_uses_continuous_score": True,
        "binary_withheld_threshold_required": False,
        "binary_threshold_policy": "predefined only if later needed",
        "allowed_strategies": [
            "use predefined default SVTuner threshold from method configuration",
            "use threshold fixed from prior method-validation experiments",
            "use unsupervised fixed quantile based only on unsupported score distribution, not endpoint labels",
            "use no binary threshold and report continuous withheld_score as primary",
        ],
        "forbidden_strategies": [
            "maximize CTA endpoint AUROC",
            "tune threshold based on Phase 5 baseline results",
            "tune threshold based on visual appearance",
        ],
        "threshold_selected_using_endpoint_labels": False,
    }
    write_json(OUT_DIR / "svtuner_threshold_policy.json", threshold_policy)

    entry_df, entry_found, entry_ambiguous = entrypoint_inventory()
    entry_df.to_csv(OUT_DIR / "svtuner_entrypoint_inventory.csv", index=False)

    dry_plan = f"""SVTuner command dry plan

Future run name:
svtuner_immune_all_dropout

Input ST expression:
{rel(P3R / 'st_expression_full_2248_data.csv.gz')}

Input coordinates:
{rel(P3R / 'st_coordinates_full_2248.csv')}

Immune-all-dropout reference expression / metadata:
{rel(P4 / 'baseline_immune_all_dropout' / 'inputs' / 'reference_expression_used.csv.gz')}
{rel(P4 / 'baseline_immune_all_dropout' / 'reference_metadata_used.csv')}

Formal gene intersection:
{rel(P3R / 'formal_gene_intersection.txt')}

Expected output directory:
visualizations/bioapp_experiment/bioapp_phase7_svtuner_aware_execution_against_frozen_cta_immune_endpoint/svtuner_immune_all_dropout

Expected spot-level output:
svtuner_immune_all_dropout_spot_level_output.csv

Expected command template:
python <svtuner_entrypoint> --condition svtuner_immune_all_dropout --input-manifest <svtuner_future_input_manifest.json> --output <phase7_output_dir>

Environment notes:
Use the project SVTuner/CytoSPACE environment. Phase 7 must explicitly disable endpoint-dependent threshold tuning.

No SVTuner command was executed in Phase 6.
"""
    (OUT_DIR / "svtuner_command_dry_plan.txt").write_text(dry_plan, encoding="utf-8")

    phase7_gate = {
        "phase7_name": "BioApp Phase 7 鈥?SVTuner-aware execution against frozen CTA Immune endpoint",
        "phase7_allowed_only_if_phase6_pass": True,
        "phase7_requires_manual_authorization": True,
        "allowed_in_phase7": [
            "run SVTuner-aware execution for immune-all-dropout condition",
            "generate raw SVTuner outputs",
            "generate run logs",
            "generate output inventory",
        ],
        "not_allowed_in_phase7": [
            "compute final prevention analysis unless explicitly assigned to a later evaluation phase",
            "claim biological application completed",
            "modify frozen endpoint",
            "choose endpoint-dependent thresholds",
        ],
        "phase8_expected": "SVTuner-vs-endpoint evaluation and baseline/SVTuner comparison",
    }
    write_json(OUT_DIR / "phase7_execution_gate.json", phase7_gate)

    anchor = {
        "baseline_anchor_condition": "baseline_immune_all_dropout",
        "phase5_metric_anchor": (p5.get("baseline_metric_highlights") or {}).get("baseline_immune_all_dropout", {}),
    }
    write_json(OUT_DIR / "phase5_baseline_anchor_snapshot.json", anchor)

    generated = {
        "svtuner_input_manifest_generated": True,
        "svtuner_output_contract_generated": True,
        "svtuner_metric_alignment_plan_generated": True,
        "svtuner_threshold_policy_generated": True,
        "svtuner_command_dry_plan_generated": True,
    }
    planning_complete = manifest_complete and all(generated.values())
    if errors:
        decision = "FAIL"
        reasons = errors
        ready = False
    elif entry_ambiguous:
        decision = "REVIEW_REQUIRED"
        reasons = ["SVTuner entrypoint ambiguous"]
        ready = False
    elif planning_complete:
        decision = "PASS"
        reasons = ["SVTuner-aware execution planning complete; Phase 7 requires manual authorization"]
        ready = True
    else:
        decision = "REVIEW_REQUIRED"
        reasons = ["input manifest or planning artifacts incomplete"]
        ready = False

    summary = {
        "phase": PHASE,
        "decision": decision,
        "decision_reasons": reasons,
        "stage_type": STAGE_TYPE,
        "input_phase2c_decision": p2c.get("decision"),
        "input_phase3r_decision": p3r.get("decision"),
        "input_phase4_decision": p4.get("decision"),
        "input_phase5_decision": p5.get("decision"),
        "endpoint_frozen": True,
        "primary_endpoint": PRIMARY_ENDPOINT,
        "primary_endpoint_spatial_note": PRIMARY_ENDPOINT_NOTE,
        "phase5_baseline_evaluation_available": True,
        "baseline_anchor_condition": "baseline_immune_all_dropout",
        "primary_svtuner_condition": "svtuner_immune_all_dropout",
        "primary_condition_selected_based_on_SVTuner_output": False,
        "analysis_universe_type": "full_frozen_endpoint_spot_universe",
        "analysis_universe_spots": 2248,
        "main_analysis_spots_positive_negative": int(p5.get("main_analysis_spots", 1888)),
        "n_gene_overlap": 2000,
        **generated,
        "svtuner_entrypoint_found": entry_found,
        "svtuner_entrypoint_ambiguous": entry_ambiguous,
        "SVTuner_run": False,
        "SVTuner_Stage3_run": False,
        "SVTuner_Stage4_run": False,
        "Stage4_run": False,
        "CytoSPACE_rerun": False,
        "prevention_analysis_run": False,
        "SVTuner_vs_endpoint_metric_computed": False,
        "endpoint_redefined": False,
        "expression_markers_used_to_define_endpoint": False,
        "mapping_outputs_used_to_define_endpoint": False,
        "threshold_selected_using_endpoint_labels": False,
        "ready_for_phase7_svtuner_execution": ready,
        "phase7_requires_manual_authorization": True,
        "biological_application_allowed": False,
        "allowed_claims": ALLOWED_CLAIMS,
        "disallowed_claims": DISALLOWED_CLAIMS,
        "output_dir": rel(OUT_DIR),
    }
    write_json(OUT_DIR / "bioapp_phase6_svtuner_execution_planning_summary.json", summary)
    golden = {
        "stage_type": STAGE_TYPE,
        "biological_question_defined": True,
        "external_endpoint_predefined": True,
        "endpoint_independent_from_SVTuner": True,
        "endpoint_spatially_registered": True,
        "baseline_comparison_available": True,
        "endpoint_specific_improvement_defined": True,
        "endpoint_specific_quantitative_metric_available": True,
        "endpoint_baseline_SVTuner_spatial_comparison_available": False,
        "interpretation_boundary_defined": True,
        "biological_application_allowed": False,
        "allowed_claim_level": "SVTuner-aware execution planning only",
        "decision": decision,
    }
    write_json(OUT_DIR / "bioapp_phase6_golden_rules_v2_1_check.json", golden)

    next_text = (
        "BioApp Phase 7 - SVTuner-aware execution against frozen CTA Immune endpoint\nManual authorization required before execution."
        if decision == "PASS"
        else "Review SVTuner entrypoint, output contract, or metric alignment plan before SVTuner execution.\nDo not run SVTuner / Stage4 / prevention analysis."
        if decision == "REVIEW_REQUIRED"
        else "Stop. Fix planning boundary violation or invalid upstream inputs before retrying Phase 6."
    )
    decision_txt = f"""BioApp Phase 6 - SVTuner-aware execution planning against frozen CTA Immune endpoint

Decision: {decision}

Input endpoint:
CTA-defined Immune cells

Endpoint status:
Frozen at Phase 2C

Primary endpoint note:
Sparse immune-positive spatial compartments, not a large continuous ROI.

Baseline anchor:
baseline_immune_all_dropout

Primary future SVTuner condition:
svtuner_immune_all_dropout

Analysis universe:
2248 frozen endpoint spots
2000 exact-overlap genes
main analysis set = positive + negative = {summary['main_analysis_spots_positive_negative']}

Planning outputs:
svtuner input manifest generated = {generated['svtuner_input_manifest_generated']}
svtuner output contract generated = {generated['svtuner_output_contract_generated']}
metric alignment plan generated = {generated['svtuner_metric_alignment_plan_generated']}
threshold policy generated = {generated['svtuner_threshold_policy_generated']}
command dry plan generated = {generated['svtuner_command_dry_plan_generated']}
entrypoint found = {entry_found}
entrypoint ambiguous = {entry_ambiguous}

Boundary checks:
SVTuner run: false
SVTuner Stage3 run: false
SVTuner Stage4 run: false
Stage4 run: false
CytoSPACE rerun: false
Prevention analysis run: false
SVTuner-vs-endpoint metric computed: false
Endpoint redefined: false
Threshold selected using endpoint labels: false

Allowed claims:
{chr(10).join('- ' + x for x in ALLOWED_CLAIMS)}

Disallowed claims:
{chr(10).join('- ' + x for x in DISALLOWED_CLAIMS)}

Next:
{next_text}
"""
    (OUT_DIR / "decision.txt").write_text(decision_txt, encoding="utf-8")
    readme = """# BioApp Phase 6

Purpose: SVTuner-aware execution planning against the frozen CTA Immune endpoint.

Inputs are Phase 2C, Phase 3R, Phase 4, and Phase 5 outputs. The frozen CTA endpoint is not redefined.

Phase 5 `baseline_immune_all_dropout` is the baseline anchor for future comparison.

This phase only plans future SVTuner execution. It does not run SVTuner, Stage4, CytoSPACE, prevention analysis, or SVTuner-vs-endpoint metrics.
"""
    (OUT_DIR / "README.md").write_text(readme, encoding="utf-8")
    manifest = []
    for path in sorted(OUT_DIR.iterdir()):
        if path.is_file():
            manifest.append(
                {
                    "file": rel(path),
                    "type": path.suffix.lstrip(".") or "text",
                    "description": "BioApp Phase 6 planning artifact",
                    "created_by_phase": PHASE,
                    "status": "generated",
                    "notes": "planning only; no SVTuner execution",
                }
            )
    pd.DataFrame(manifest).to_csv(OUT_DIR / "manifest.csv", index=False)

    print("BioApp Phase 6 completed.")
    print()
    print("Decision:")
    print(decision)
    print()
    print("Input endpoint:")
    print("CTA-defined Immune cells")
    print()
    print("Endpoint frozen:")
    print("true")
    print()
    print("Baseline anchor:")
    print("baseline_immune_all_dropout")
    print()
    print("Primary future SVTuner condition:")
    print("svtuner_immune_all_dropout")
    print()
    print("Analysis universe:")
    print("2248 frozen endpoint spots")
    print("2000 exact-overlap genes")
    print(f"main analysis spots = {summary['main_analysis_spots_positive_negative']}")
    print()
    print("Planning outputs:")
    print(f"svtuner input manifest generated = {generated['svtuner_input_manifest_generated']}")
    print(f"svtuner output contract generated = {generated['svtuner_output_contract_generated']}")
    print(f"metric alignment plan generated = {generated['svtuner_metric_alignment_plan_generated']}")
    print(f"threshold policy generated = {generated['svtuner_threshold_policy_generated']}")
    print(f"command dry plan generated = {generated['svtuner_command_dry_plan_generated']}")
    print(f"entrypoint found = {entry_found}")
    print(f"entrypoint ambiguous = {entry_ambiguous}")
    print()
    print("Boundary checks:")
    print("SVTuner run: false")
    print("SVTuner Stage3 run: false")
    print("SVTuner Stage4 run: false")
    print("Stage4 run: false")
    print("CytoSPACE rerun: false")
    print("Prevention analysis run: false")
    print("SVTuner-vs-endpoint metric computed: false")
    print("Endpoint redefined: false")
    print("Threshold selected using endpoint labels: false")
    print()
    print("Ready for Phase 7 SVTuner-aware execution:")
    print(str(ready).lower())
    print()
    print("Next:")
    print(next_text)
    return 0 if decision in {"PASS", "REVIEW_REQUIRED"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
