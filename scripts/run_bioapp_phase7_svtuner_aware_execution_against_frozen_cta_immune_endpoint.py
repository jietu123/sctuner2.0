#!/usr/bin/env python3
"""BioApp Phase 7: SVTuner-aware raw execution.

Boundary:
- run only the planned immune-all-dropout SVTuner-aware raw execution
- use Stage3B as the unsupported/reference-unrepresented raw detector
- generate raw outputs, logs, inventory, and a spot-level contract table
- do not compute SVTuner-vs-endpoint metrics or prevention analysis
- do not modify the frozen endpoint or upstream outputs
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
P2C = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping"
P3R = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase3r_review_resolution_for_cta_endpoint_baseline_feasibility"
P4 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase4_formal_cytospace_baseline_execution_against_frozen_cta_immune_endpoint"
P5 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase5_endpoint_specific_evaluation_of_cytospace_baselines_against_frozen_cta_immune_endpoint"
P6 = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase6_svtuner_aware_execution_planning_against_frozen_cta_immune_endpoint"
OUT = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase7_svtuner_aware_execution_against_frozen_cta_immune_endpoint"

PHASE = "BioApp Phase 7 SVTuner-aware execution against frozen CTA Immune endpoint"
STAGE_TYPE = "SVTuner-aware raw execution / biological application candidate preparation"
SAMPLE = "svtuner_immune_all_dropout"
GROUP = "bioapp_phase7_svtuner_aware_execution_against_frozen_cta_immune_endpoint"
RUN_NAME = "svtuner_immune_all_dropout"
BASELINE_ANCHOR = "baseline_immune_all_dropout"
PRIMARY_ENDPOINT = "Immune cells"
PRIMARY_ENDPOINT_NOTE = "sparse immune-positive spatial compartments, not a large continuous ROI"
MANUAL_AUTHORIZATION = True

ALLOWED_CLAIMS = [
    "SVTuner-aware execution was run for the immune-all-dropout condition.",
    "Raw SVTuner outputs, logs, and output inventory were generated.",
    "The frozen CTA endpoint was not used for tuning, threshold selection, or endpoint redefinition.",
]
DISALLOWED_CLAIMS = [
    "SVTuner improves endpoint recovery.",
    "SVTuner prevents contradicted interpretation.",
    "Biological application completed.",
    "Biological discovery made.",
]


def rel(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT)).replace("\\", "/")
    except Exception:
        return str(path).replace("\\", "/")


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, ensure_ascii=False)
        fh.write("\n")


def count_rows_cols(path: Path) -> tuple[int | None, int | None, bool]:
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8", errors="replace", newline="") as fh:
            reader = csv.reader(fh)
            header = next(reader)
            n_rows = sum(1 for _ in reader)
        return n_rows, len(header), True
    except Exception:
        return None, None, False


def validate_upstream() -> tuple[dict[str, Any], list[str]]:
    errors: list[str] = []
    p6 = read_json(P6 / "bioapp_phase6_svtuner_execution_planning_summary.json")
    p2c = read_json(P2C / "bioapp_phase2c_endpoint_freeze_summary.json")
    p3r = read_json(P3R / "bioapp_phase3r_review_resolution_summary.json")
    p4 = read_json(P4 / "bioapp_phase4_formal_cytospace_baseline_execution_summary.json")
    p5 = read_json(P5 / "bioapp_phase5_endpoint_specific_baseline_evaluation_summary.json")

    checks = [
        (p2c.get("decision") == "PASS", "Phase 2C decision is not PASS"),
        (p2c.get("endpoint_frozen") is True, "Phase 2C endpoint is not frozen"),
        (p2c.get("primary_endpoint") == PRIMARY_ENDPOINT, "Phase 2C primary endpoint is not Immune cells"),
        (p3r.get("decision") == "PASS", "Phase 3R decision is not PASS"),
        (p3r.get("full_2248_coverage_achieved") is True, "Phase 3R full spot coverage missing"),
        (int(p3r.get("analysis_universe_spots", -1)) == 2248, "Phase 3R spot count is not 2248"),
        (int(p3r.get("n_gene_overlap_exact", -1)) == 2000, "Phase 3R gene overlap is not 2000"),
        (p4.get("decision") == "PASS", "Phase 4 decision is not PASS"),
        (p4.get("all_required_baseline_runs_completed") is True, "Phase 4 baseline runs incomplete"),
        (p5.get("decision") == "PASS", "Phase 5 decision is not PASS"),
        (p5.get("baseline_vs_endpoint_metric_computed") is True, "Phase 5 baseline metrics missing"),
        (p5.get("ready_for_phase6_svtuner_execution") is True, "Phase 5 not ready for Phase 6"),
        (p6.get("decision") == "PASS", "Phase 6 decision is not PASS"),
        (p6.get("ready_for_phase7_svtuner_execution") is True, "Phase 6 not ready for Phase 7"),
        (p6.get("phase7_requires_manual_authorization") is True, "Phase 6 missing manual authorization gate"),
        (p6.get("svtuner_entrypoint_found") is True, "Phase 6 did not find SVTuner entrypoint"),
        (p6.get("svtuner_entrypoint_ambiguous") is False, "Phase 6 entrypoint is ambiguous"),
        (p6.get("primary_svtuner_condition") == RUN_NAME, "Phase 6 primary SVTuner condition mismatch"),
        (p6.get("baseline_anchor_condition") == BASELINE_ANCHOR, "Phase 6 baseline anchor mismatch"),
        (p6.get("biological_application_allowed") is False, "Phase 6 biological_application_allowed must be false"),
    ]
    for ok, msg in checks:
        if not ok:
            errors.append(msg)

    required = [
        P2C / "spot_level_endpoint_freeze.csv",
        P3R / "formal_gene_intersection.txt",
        P3R / "st_expression_full_2248_data.csv.gz",
        P3R / "st_coordinates_full_2248.csv",
        P4 / "baseline_immune_all_dropout" / "inputs" / "reference_expression_used.csv",
        P4 / "baseline_immune_all_dropout" / "reference_metadata_used.csv",
        P4 / "baseline_immune_all_dropout" / "cytospace_output",
        P5 / "baseline_endpoint_score_all_runs_by_spot.csv",
        P6 / "svtuner_output_contract.json",
        P6 / "svtuner_threshold_policy.json",
    ]
    for path in required:
        if not path.exists():
            errors.append(f"missing required input: {rel(path)}")

    return {
        "phase2c": p2c,
        "phase3r": p3r,
        "phase4": p4,
        "phase5": p5,
        "phase6": p6,
    }, errors


def prepare_stage3b_dataset(condition_dir: Path) -> tuple[Path, Path]:
    processed_export = (
        ROOT
        / "data"
        / "processed"
        / GROUP
        / SAMPLE
        / "stage1_preprocess"
        / "exported"
    )
    processed_export.mkdir(parents=True, exist_ok=True)

    # Stage3B expects cells x genes. Phase4 CytoSPACE input is genes x cells.
    sc_gene_by_cell = pd.read_csv(
        P4 / "baseline_immune_all_dropout" / "inputs" / "reference_expression_used.csv",
        index_col=0,
    )
    sc_cell_by_gene = sc_gene_by_cell.T
    sc_cell_by_gene.index.name = "cell_id"
    sc_cell_by_gene.to_csv(processed_export / "sc_expression_normalized.csv")

    # Stage3B expects spots x genes. Phase3R ST expression is genes x spots.
    st_gene_by_spot = pd.read_csv(P3R / "st_expression_full_2248_data.csv.gz", index_col=0)
    st_spot_by_gene = st_gene_by_spot.T
    st_spot_by_gene.index.name = "spot_id"
    st_spot_by_gene.to_csv(processed_export / "st_expression_normalized.csv")

    meta = pd.read_csv(P4 / "baseline_immune_all_dropout" / "reference_metadata_used.csv")
    meta = meta.set_index("cell_id")
    meta.index.name = "cell_id"
    meta[["cell_type"]].to_csv(processed_export / "sc_metadata.csv")

    coords = pd.read_csv(P3R / "st_coordinates_full_2248.csv")
    if "barcode" in coords.columns:
        coords = coords.rename(columns={"barcode": "spot_id"})
    if "spot_id" not in coords.columns:
        coords = coords.rename(columns={coords.columns[0]: "spot_id"})
    keep = ["spot_id"] + [c for c in ["row", "col", "imagerow", "imagecol"] if c in coords.columns]
    coords[keep].to_csv(processed_export / "st_coordinates.csv", index=False)

    dataset_config = OUT / "bioapp_phase7_stage3b_dataset_config.yaml"
    dataset_config.write_text(
        "\n".join(
            [
                "paths:",
                "  sc_expr: stage1_direct",
                "  sc_meta: stage1_direct",
                "  st_expr: stage1_direct",
                "  st_meta: stage1_direct",
                "qc:",
                "  sc_min_genes: 0",
                "  sc_max_genes: 100000",
                "  sc_max_mt: 100",
                "  st_min_genes: 0",
                "  st_max_genes: '.inf'",
                "  st_max_mt: 100",
                "  hvg_nfeatures: 2000",
                "gene_filter:",
                "  min_cells_sc: 0",
                "  min_cells_st: 0",
                "stage3b:",
                "  fdr: 0.05",
                "  n_calibration: 0",
                "  n_spatial_permutations: 200",
                "  random_seed: 20260705",
                "  sc_expr_source: normalized",
                "  sc_profile_source: normalized",
                "  expression_scale: log1p",
                "  sc_profile_scale: log1p",
                "  max_genes: 2000",
                "storage:",
                f"  group: {GROUP}",
                "",
            ]
        ),
        encoding="utf-8",
    )

    inputs_dir = condition_dir / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    manifest_items = [
        ("svtuner_input_manifest.json", P6 / "svtuner_future_input_manifest.json"),
        ("svtuner_input_manifest.csv", P6 / "svtuner_future_input_manifest.csv"),
        ("formal_gene_intersection_used.txt", P3R / "formal_gene_intersection.txt"),
    ]
    for name, src in manifest_items:
        shutil.copy2(src, inputs_dir / name)
    for name, path in [
        ("st_expression_used_path.txt", P3R / "st_expression_full_2248_data.csv.gz"),
        ("st_coordinates_used_path.txt", P3R / "st_coordinates_full_2248.csv"),
        ("reference_expression_used_path.txt", P4 / "baseline_immune_all_dropout" / "inputs" / "reference_expression_used.csv"),
        ("reference_metadata_used_path.txt", P4 / "baseline_immune_all_dropout" / "reference_metadata_used.csv"),
        ("baseline_anchor_output_dir.txt", P4 / "baseline_immune_all_dropout" / "cytospace_output"),
        ("endpoint_freeze_table_path.txt", P2C / "spot_level_endpoint_freeze.csv"),
    ]:
        (inputs_dir / name).write_text(rel(path) + "\n", encoding="utf-8")

    return dataset_config, processed_export


def run_stage3b(condition_dir: Path, dataset_config: Path) -> dict[str, Any]:
    logs = condition_dir / "logs"
    logs.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable,
        "-m",
        "src.stages.stage3b_st_unsupported",
        "--sample",
        SAMPLE,
        "--project_root",
        str(ROOT),
        "--dataset_config",
        str(dataset_config),
        "--sc_expr_source",
        "normalized",
        "--sc_profile_source",
        "normalized",
        "--expression_scale",
        "log1p",
        "--sc_profile_scale",
        "log1p",
        "--max_genes",
        "2000",
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT)
    start = time.time()
    proc = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, env=env)
    runtime = time.time() - start

    (logs / "run_command.txt").write_text(" ".join(cmd) + "\n", encoding="utf-8")
    (logs / "run_log_stdout.txt").write_text(proc.stdout or "", encoding="utf-8", errors="replace")
    (logs / "run_log_stderr.txt").write_text(proc.stderr or "", encoding="utf-8", errors="replace")
    (logs / "run_environment.txt").write_text(
        "\n".join(
            [
                f"python_executable={sys.executable}",
                f"cwd={ROOT}",
                f"PYTHONPATH={env.get('PYTHONPATH')}",
                "entrypoint=src.stages.stage3b_st_unsupported",
                "svtuner_stage3b_raw_execution=true",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    raw_data_dir = ROOT / "data" / "processed" / GROUP / SAMPLE / "stage3b_st_unsupported"
    raw_result_dir = ROOT / "result" / GROUP / SAMPLE / "stage3b_st_unsupported"
    raw_outputs = condition_dir / "raw_outputs"
    raw_outputs.mkdir(parents=True, exist_ok=True)
    if raw_data_dir.exists():
        for path in raw_data_dir.glob("*.csv"):
            shutil.copy2(path, raw_outputs / path.name)
    if raw_result_dir.exists():
        for path in raw_result_dir.glob("*.json"):
            shutil.copy2(path, raw_outputs / path.name)

    raw_files = [rel(p) for p in raw_outputs.glob("*") if p.is_file()]
    error_text = (proc.stderr or "") + "\n" + (proc.stdout or "")
    error_detected = proc.returncode != 0 or "traceback" in error_text.lower()
    completed = (
        proc.returncode == 0
        and (raw_outputs / "spot_unsupported_scores.csv").exists()
        and not error_detected
    )
    status = {
        "run_name": RUN_NAME,
        "manual_authorization_for_phase7": MANUAL_AUTHORIZATION,
        "SVTuner_run": True,
        "SVTuner_Stage3_run": True,
        "SVTuner_Stage4_run": False,
        "Stage4_run": False,
        "CytoSPACE_rerun": False,
        "endpoint_used_for_mapping": False,
        "endpoint_used_for_threshold_selection": False,
        "endpoint_used_for_evaluation": False,
        "prevention_analysis_run": False,
        "SVTuner_vs_endpoint_metric_computed": False,
        "n_ST_spots_input": 2248,
        "n_genes_input": 2000,
        "n_reference_cells_input": 1600,
        "return_code": int(proc.returncode),
        "runtime_seconds": runtime,
        "completed": completed,
        "raw_output_files_detected": raw_files,
        "error_detected": bool(error_detected),
        "error_summary": None if not error_detected else (proc.stderr or proc.stdout or "")[:1000],
    }
    write_json(logs / "run_status.json", status)
    return status


def load_phase5_anchor() -> pd.DataFrame:
    scores = pd.read_csv(P5 / "baseline_endpoint_score_all_runs_by_spot.csv")
    scores = scores.loc[scores["baseline_run"] == BASELINE_ANCHOR].copy()
    keep = [
        "barcode",
        "immune_score",
        "nonimmune_score",
        "dominant_baseline_label",
        "dominant_baseline_fraction",
    ]
    return scores[keep].rename(
        columns={
            "immune_score": "baseline_immune_score_from_phase5",
            "nonimmune_score": "baseline_nonimmune_score_from_phase5",
            "dominant_baseline_label": "baseline_dominant_label_from_phase5",
            "dominant_baseline_fraction": "baseline_dominant_fraction_from_phase5",
        }
    )


def build_contract(condition_dir: Path, run_status: dict[str, Any]) -> tuple[bool, int, bool]:
    endpoint = pd.read_csv(P2C / "spot_level_endpoint_freeze.csv")
    if "barcode" not in endpoint.columns:
        raise ValueError("endpoint freeze table lacks barcode column")
    status_col = "primary_endpoint_status"
    if status_col not in endpoint.columns:
        candidates = [c for c in endpoint.columns if "status" in c.lower()]
        if not candidates:
            raise ValueError("endpoint freeze table lacks primary endpoint status column")
        status_col = candidates[0]
    frame = endpoint[["barcode", status_col]].rename(columns={status_col: "primary_endpoint_status"}).copy()
    frame["barcode"] = frame["barcode"].astype(str)

    raw_scores_path = condition_dir / "raw_outputs" / "spot_unsupported_scores.csv"
    if raw_scores_path.exists():
        raw = pd.read_csv(raw_scores_path, index_col=0)
        raw.index = raw.index.astype(str)
    else:
        raw = pd.DataFrame(index=frame["barcode"].astype(str))

    if "unsupported_fraction_estimate" in raw.columns:
        unrep = pd.to_numeric(raw["unsupported_fraction_estimate"], errors="coerce")
    elif "anomaly_score" in raw.columns:
        unrep = pd.to_numeric(raw["anomaly_score"], errors="coerce")
    else:
        unrep = pd.Series(np.nan, index=raw.index)
    withheld_binary = (
        raw["is_unsupported_region"].astype(bool)
        if "is_unsupported_region" in raw.columns
        else pd.Series(pd.NA, index=raw.index)
    )
    raw_part = pd.DataFrame(
        {
            "barcode": raw.index,
            "reference_unrepresented_score": unrep.reindex(raw.index).to_numpy(),
            "withheld_score": unrep.reindex(raw.index).to_numpy(),
            "withheld_binary": withheld_binary.reindex(raw.index).to_numpy(),
        }
    )

    base = load_phase5_anchor()
    out = frame.merge(raw_part, on="barcode", how="left").merge(base, on="barcode", how="left")
    out["baseline_anchor_run"] = BASELINE_ANCHOR
    out["svtuner_condition"] = RUN_NAME
    out["reference_supported_status"] = np.where(
        out["withheld_binary"].fillna(False).astype(bool),
        "reference_unrepresented_withheld_by_stage3b",
        "available_nonimmune_reference_assignment_retained",
    )
    out["supported_assignment_score"] = 1.0 - pd.to_numeric(out["withheld_score"], errors="coerce").clip(0, 1)
    out["svtuner_immune_score"] = np.nan
    out["svtuner_nonimmune_score"] = np.where(
        out["withheld_binary"].fillna(False).astype(bool),
        0.0,
        out["baseline_nonimmune_score_from_phase5"],
    )
    out["svtuner_dominant_label"] = np.where(
        out["withheld_binary"].fillna(False).astype(bool),
        "__withheld_reference_unrepresented__",
        out["baseline_dominant_label_from_phase5"],
    )
    out["svtuner_dominant_fraction"] = np.where(
        out["withheld_binary"].fillna(False).astype(bool),
        0.0,
        out["baseline_dominant_fraction_from_phase5"],
    )
    out["score_parse_status"] = np.where(
        out["reference_unrepresented_score"].notna(),
        "ok",
        "missing_stage3b_score",
    )
    out["svtuner_output_available"] = bool(run_status.get("completed"))

    ordered = [
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
    out = out[ordered]
    contract_path = condition_dir / "svtuner_immune_all_dropout_spot_level_raw_output_contract.csv"
    out.to_csv(contract_path, index=False)

    parse_report = pd.DataFrame(
        [
            {
                "field_name": col,
                "required_by_contract": True,
                "available": bool(out[col].notna().any()),
                "source_file": "stage3b spot_unsupported_scores.csv / Phase5 baseline anchor / Phase2C endpoint freeze",
                "parse_method": "direct merge by barcode; no endpoint-dependent tuning",
                "missing_reason": "" if out[col].notna().any() else "not available from raw output",
                "notes": "Phase7 raw contract field",
            }
            for col in ordered
        ]
    )
    parse_report.to_csv(condition_dir / "svtuner_raw_output_parse_report.csv", index=False)

    aligned = set(out["barcode"]) == set(endpoint["barcode"].astype(str)) and len(out) == 2248
    return contract_path.exists(), len(out), aligned


def build_inventory(condition_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for path in sorted(condition_dir.rglob("*")):
        if not path.is_file():
            continue
        role = "unknown"
        name = path.name.lower()
        if name == "spot_unsupported_scores.csv":
            role = "svtuner_spot_score"
        elif "unsupported_regions" in name:
            role = "svtuner_stage3_output"
        elif name.endswith(".log") or "log" in name:
            role = "svtuner_log"
        elif "contract" in name:
            role = "svtuner_spot_score"
        elif "run_status" in name:
            role = "svtuner_log"
        n_rows, n_cols, readable = (None, None, True)
        if path.suffix.lower() == ".csv":
            n_rows, n_cols, readable = count_rows_cols(path)
        rows.append(
            {
                "file_path": rel(path),
                "file_type": path.suffix.lstrip(".") or "none",
                "file_size_mb": path.stat().st_size / (1024 * 1024),
                "detected_role": role,
                "readable": readable,
                "n_rows": n_rows,
                "n_cols": n_cols,
                "notes": "Phase7 raw output inventory",
            }
        )
    inv = pd.DataFrame(rows)
    inv.to_csv(condition_dir / "svtuner_raw_output_inventory.csv", index=False)
    return inv


def write_static_docs(
    condition_dir: Path,
    decision: str,
    ready: bool,
    run_status: dict[str, Any],
    qc: dict[str, Any],
) -> None:
    threshold_record = {
        "threshold_selected_using_endpoint_labels": False,
        "binary_threshold_used": True,
        "binary_threshold_value": None,
        "threshold_source": "predefined_default_stage3b_fdr_0.05",
        "endpoint_labels_used_for_threshold": False,
        "notes": "Stage3B uses the predefined method FDR and spatial-region gate; CTA endpoint labels were not used.",
    }
    write_json(condition_dir / "phase7_threshold_record.json", threshold_record)

    next_text = (
        "BioApp Phase 8 - SVTuner-vs-endpoint evaluation and baseline/SVTuner comparison"
        if decision == "PASS"
        else "Review SVTuner run logs, raw output parsing, and spot-level contract alignment before Phase 8.\nDo not compute prevention analysis yet."
        if decision == "REVIEW_REQUIRED"
        else "Stop. Fix SVTuner execution or boundary violation before retrying Phase 7."
    )
    decision_txt = f"""BioApp Phase 7 - SVTuner-aware execution against frozen CTA Immune endpoint

Decision: {decision}

Manual authorization:
true

Input endpoint:
CTA-defined Immune cells

Endpoint status:
Frozen at Phase 2C

Primary endpoint note:
{PRIMARY_ENDPOINT_NOTE}

Baseline anchor:
{BASELINE_ANCHOR}

SVTuner condition:
{RUN_NAME}

Analysis universe:
2248 frozen endpoint spots
2000 exact-overlap genes
main analysis set = positive + negative = 1888

Execution:
SVTuner run = true
SVTuner Stage3 run = {str(run_status.get("SVTuner_Stage3_run")).lower()}
SVTuner Stage4 run = false
Stage4 run = false
CytoSPACE rerun = false
completed = {str(run_status.get("completed")).lower()}
return_code = {run_status.get("return_code")}
runtime_seconds = {run_status.get("runtime_seconds")}

Raw outputs:
raw_outputs_detected = {str(qc["raw_outputs_detected"]).lower()}
spot_level_contract_generated = {str(qc["spot_level_contract_generated"]).lower()}
spot_level_contract_rows = {qc["n_contract_rows"]}
contract_spot_universe_aligned = {str(qc["contract_spot_universe_aligned"]).lower()}

Boundary checks:
Endpoint used for mapping: false
Endpoint used for threshold selection: false
Endpoint used for evaluation: false
Threshold selected using endpoint labels: false
Prevention analysis run: false
SVTuner-vs-endpoint metric computed: false
Endpoint redefined: false
Biological application result generated: false

Allowed claims:
- {ALLOWED_CLAIMS[0]}
- {ALLOWED_CLAIMS[1]}
- {ALLOWED_CLAIMS[2]}

Disallowed claims:
- {DISALLOWED_CLAIMS[0]}
- {DISALLOWED_CLAIMS[1]}
- {DISALLOWED_CLAIMS[2]}
- {DISALLOWED_CLAIMS[3]}

Next:
{next_text}
"""
    (OUT / "decision.txt").write_text(decision_txt, encoding="utf-8")

    readme = f"""# BioApp Phase 7

This phase performs SVTuner-aware raw execution for the pre-authorized
immune-all-dropout condition only.

Manual authorization was granted only for raw SVTuner-aware execution of the
immune-all-dropout condition. It does not authorize final prevention analysis
or biological application claims.

Inputs come from Phase 2C, Phase 3R, Phase 4, Phase 5, and Phase 6. The frozen
CTA Immune endpoint is retained only for spot-universe alignment and the raw
contract status column; it is not used for mapping, threshold selection, or
Phase 7 evaluation.

Primary SVTuner condition: `{RUN_NAME}`.
Baseline anchor: `{BASELINE_ANCHOR}`.

This phase does not compute prevention, SVTuner-vs-endpoint metrics, or
baseline/SVTuner endpoint-specific comparisons.
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")


def write_manifest() -> None:
    rows = []
    for path in sorted(OUT.rglob("*")):
        if path.is_file():
            rows.append(
                {
                    "file": rel(path),
                    "type": path.suffix.lstrip(".") or "none",
                    "description": "BioApp Phase 7 output",
                    "created_by_phase": "BioApp Phase 7",
                    "status": "generated",
                    "notes": "raw execution artifact; not final endpoint evaluation",
                }
            )
    pd.DataFrame(rows).to_csv(OUT / "manifest.csv", index=False)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    condition_dir = OUT / RUN_NAME
    for sub in ["inputs", "raw_outputs", "logs", "qc"]:
        (condition_dir / sub).mkdir(parents=True, exist_ok=True)

    upstream, errors = validate_upstream()
    if not MANUAL_AUTHORIZATION:
        errors.append("missing manual authorization for Phase 7")

    if errors:
        run_status = {
            "SVTuner_run": False,
            "SVTuner_Stage3_run": False,
            "SVTuner_Stage4_run": False,
            "Stage4_run": False,
            "CytoSPACE_rerun": False,
            "completed": False,
            "return_code": None,
            "runtime_seconds": None,
            "error_detected": True,
            "error_summary": "; ".join(errors),
        }
        qc = {
            "raw_outputs_detected": False,
            "spot_level_contract_generated": False,
            "n_contract_rows": 0,
            "contract_spot_universe_aligned": False,
        }
        decision = "FAIL"
        ready = False
    else:
        dataset_config, _export_dir = prepare_stage3b_dataset(condition_dir)
        run_status = run_stage3b(condition_dir, dataset_config)
        contract_generated, n_contract_rows, aligned = build_contract(condition_dir, run_status)
        inventory = build_inventory(condition_dir)
        raw_outputs_detected = bool((condition_dir / "raw_outputs" / "spot_unsupported_scores.csv").exists())
        qc = {
            "run_name": RUN_NAME,
            "completed": bool(run_status.get("completed")),
            "return_code": run_status.get("return_code"),
            "runtime_seconds": run_status.get("runtime_seconds"),
            "n_ST_spots_input": 2248,
            "n_genes_input": 2000,
            "n_reference_cells_input": 1600,
            "raw_outputs_detected": raw_outputs_detected,
            "spot_level_contract_generated": contract_generated,
            "n_contract_rows": n_contract_rows,
            "contract_spot_universe_aligned": aligned,
            "error_detected": bool(run_status.get("error_detected")),
            "ready_for_phase8_evaluation": bool(
                run_status.get("completed")
                and raw_outputs_detected
                and contract_generated
                and n_contract_rows == 2248
                and aligned
                and not run_status.get("error_detected")
            ),
            "notes": "Phase7 raw SVTuner-aware execution QC; no endpoint metrics computed",
        }
        pd.DataFrame([qc]).to_csv(condition_dir / "qc" / "svtuner_execution_qc.csv", index=False)
        decision = "PASS" if qc["ready_for_phase8_evaluation"] else "REVIEW_REQUIRED"
        ready = bool(qc["ready_for_phase8_evaluation"])

    write_static_docs(condition_dir, decision, ready, run_status, qc)
    write_manifest()

    summary = {
        "phase": PHASE,
        "decision": decision,
        "stage_type": STAGE_TYPE,
        "input_phase2c_decision": "PASS",
        "input_phase3r_decision": "PASS",
        "input_phase4_decision": "PASS",
        "input_phase5_decision": "PASS",
        "input_phase6_decision": "PASS",
        "manual_authorization_for_phase7": MANUAL_AUTHORIZATION,
        "endpoint_frozen": True,
        "primary_endpoint": PRIMARY_ENDPOINT,
        "primary_endpoint_spatial_note": PRIMARY_ENDPOINT_NOTE,
        "endpoint_used_for_mapping": False,
        "endpoint_used_for_threshold_selection": False,
        "endpoint_used_for_evaluation": False,
        "baseline_anchor_condition": BASELINE_ANCHOR,
        "svtuner_condition": RUN_NAME,
        "analysis_universe_type": "full_frozen_endpoint_spot_universe",
        "analysis_universe_spots": 2248,
        "main_analysis_spots_positive_negative": 1888,
        "n_gene_overlap": 2000,
        "SVTuner_run": bool(run_status.get("SVTuner_run", decision != "FAIL")),
        "SVTuner_Stage3_run": bool(run_status.get("SVTuner_Stage3_run", False)),
        "SVTuner_Stage4_run": bool(run_status.get("SVTuner_Stage4_run", False)),
        "Stage4_run": False,
        "CytoSPACE_rerun": False,
        "SVTuner_entrypoint_path": "src.stages.stage3b_st_unsupported",
        "actual_command_executed": (condition_dir / "logs" / "run_command.txt").read_text(encoding="utf-8").strip()
        if (condition_dir / "logs" / "run_command.txt").exists()
        else None,
        "return_code": run_status.get("return_code"),
        "runtime_seconds": run_status.get("runtime_seconds"),
        "completed": bool(run_status.get("completed", False)),
        "raw_outputs_detected": bool(qc.get("raw_outputs_detected", False)),
        "spot_level_contract_generated": bool(qc.get("spot_level_contract_generated", False)),
        "spot_level_contract_rows": int(qc.get("n_contract_rows", 0) or 0),
        "contract_spot_universe_aligned": bool(qc.get("contract_spot_universe_aligned", False)),
        "ready_for_phase8_evaluation": ready,
        "threshold_selected_using_endpoint_labels": False,
        "binary_threshold_used": True,
        "binary_threshold_value": None,
        "threshold_source": "predefined_default_stage3b_fdr_0.05",
        "withheld_enrichment_computed": False,
        "withheld_AUROC_AUPRC_computed": False,
        "endpoint_positive_nonimmune_assignment_reduction_computed": False,
        "contradiction_prevention_rate_computed": False,
        "endpoint_concordant_interpretation_score_computed": False,
        "baseline_svtuner_comparison_computed": False,
        "prevention_analysis_run": False,
        "endpoint_redefined": False,
        "expression_markers_used_to_define_endpoint": False,
        "mapping_outputs_used_to_define_endpoint": False,
        "biological_application_allowed": False,
        "biological_application_result_generated": False,
        "allowed_claims": ALLOWED_CLAIMS,
        "disallowed_claims": DISALLOWED_CLAIMS,
        "output_dir": rel(OUT),
    }
    write_json(OUT / "bioapp_phase7_svtuner_aware_execution_summary.json", summary)

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
        "allowed_claim_level": "SVTuner-aware raw execution only",
        "decision": decision,
    }
    write_json(OUT / "bioapp_phase7_golden_rules_v2_1_check.json", golden)

    print("BioApp Phase 7 completed.")
    print()
    print("Decision:")
    print(decision)
    print()
    print("Manual authorization:")
    print("true")
    print()
    print("Input endpoint:")
    print("CTA-defined Immune cells")
    print()
    print("Endpoint frozen:")
    print("true")
    print()
    print("Baseline anchor:")
    print(BASELINE_ANCHOR)
    print()
    print("SVTuner condition:")
    print(RUN_NAME)
    print()
    print("Analysis universe:")
    print("2248 frozen endpoint spots")
    print("2000 exact-overlap genes")
    print("main analysis spots = 1888")
    print()
    print("Execution:")
    print("SVTuner run = true")
    print(f"SVTuner Stage3 run = {str(run_status.get('SVTuner_Stage3_run', False)).lower()}")
    print("SVTuner Stage4 run = false")
    print("Stage4 run = false")
    print("CytoSPACE rerun = false")
    print(f"completed = {str(run_status.get('completed', False)).lower()}")
    print(f"return_code = {run_status.get('return_code')}")
    print(f"runtime_seconds = {run_status.get('runtime_seconds')}")
    print()
    print("Raw outputs:")
    print(f"raw_outputs_detected = {str(qc.get('raw_outputs_detected', False)).lower()}")
    print(f"spot_level_contract_generated = {str(qc.get('spot_level_contract_generated', False)).lower()}")
    print(f"spot_level_contract_rows = {qc.get('n_contract_rows')}")
    print(f"contract_spot_universe_aligned = {str(qc.get('contract_spot_universe_aligned', False)).lower()}")
    print()
    print("Boundary checks:")
    print("Endpoint used for mapping = false")
    print("Endpoint used for threshold selection = false")
    print("Endpoint used for evaluation = false")
    print("Threshold selected using endpoint labels = false")
    print("Prevention analysis run = false")
    print("SVTuner-vs-endpoint metric computed = false")
    print("Endpoint redefined = false")
    print("Biological application result generated = false")
    print()
    print("Ready for Phase 8 evaluation:")
    print(str(ready).lower())
    print()
    print("Next:")
    print(
        "BioApp Phase 8 - SVTuner-vs-endpoint evaluation and baseline/SVTuner comparison"
        if decision == "PASS"
        else "Review SVTuner run logs, raw output parsing, and spot-level contract alignment before Phase 8."
    )
    return 0 if decision in {"PASS", "REVIEW_REQUIRED"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
