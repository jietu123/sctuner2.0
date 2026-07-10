#!/usr/bin/env python3
"""BioApp Phase 4: formal CytoSPACE baseline execution.

Boundary:
- run CytoSPACE baseline only
- do not run SVTuner, Stage3, Stage4
- do not compute contradiction/prevention or endpoint-specific metrics
- do not redefine endpoint
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

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PHASE3R_DIR = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase3r_review_resolution_for_cta_endpoint_baseline_feasibility"
PHASE2C_DIR = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping"
OUT_DIR = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase4_formal_cytospace_baseline_execution_against_frozen_cta_immune_endpoint"
CYTOSPACE_PY = Path("E:/ANACONDA/envs/cytospace_v1.1.0_py310/python.exe")

PHASE = "BioApp Phase 4 formal CytoSPACE baseline execution against frozen CTA Immune endpoint"
STAGE_TYPE = "formal CytoSPACE baseline execution / biological application candidate preparation"
PRIMARY_ENDPOINT = "Immune cells"
PRIMARY_ENDPOINT_NOTE = "sparse immune-positive spatial compartments, not a large continuous ROI"
LABEL_COLUMN = "cell_type"
IMMUNE_LABELS = [
    "B cells",
    "CD8 T cells",
    "Monocytes and Macrophages",
    "NK cells",
    "CD4 T cells",
    "Plasma cells",
    "T-cells",
]
ALLOWED_CLAIMS = [
    "Formal CytoSPACE baseline runs were executed against the frozen CTA Immune endpoint input universe.",
    "Full-reference, immune-all-dropout, and nonimmune-all-dropout control baseline outputs were generated or audited.",
    "Baseline outputs are ready for a later endpoint-specific evaluation phase if all required runs completed.",
]
DISALLOWED_CLAIMS = [
    "Baseline creates false niche calls.",
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


def ensure_out() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)


def count_rows_cols(path: Path) -> tuple[int | None, int | None, bool]:
    try:
        if path.suffix == ".gz":
            opener = gzip.open
            mode = "rt"
        else:
            opener = open
            mode = "r"
        with opener(path, mode, encoding="utf-8", errors="replace", newline="") as fh:
            reader = csv.reader(fh)
            header = next(reader)
            n_rows = sum(1 for _ in reader)
        return n_rows, len(header), True
    except Exception:
        return None, None, False


def validate_inputs() -> tuple[dict[str, Any], dict[str, Any], list[str]]:
    errors: list[str] = []
    p3r_path = PHASE3R_DIR / "bioapp_phase3r_review_resolution_summary.json"
    p2c_path = PHASE2C_DIR / "bioapp_phase2c_endpoint_freeze_summary.json"
    if not p3r_path.exists():
        errors.append(f"missing Phase 3R summary: {rel(p3r_path)}")
        return {}, {}, errors
    if not p2c_path.exists():
        errors.append(f"missing Phase 2C summary: {rel(p2c_path)}")
        return {}, {}, errors
    p3r = read_json(p3r_path)
    p2c = read_json(p2c_path)

    checks = [
        (p3r.get("decision") == "PASS", "Phase 3R decision is not PASS"),
        (p3r.get("ready_for_phase4") is True, "Phase 3R ready_for_phase4 is not true"),
        (p3r.get("biological_application_allowed") is False, "Phase 3R biological_application_allowed must be false"),
        (p3r.get("full_2248_coverage_achieved") is True, "Phase 3R full_2248_coverage_achieved is not true"),
        (p3r.get("analysis_universe_type") == "full_frozen_endpoint_spot_universe", "Phase 3R analysis universe is not full"),
        (int(p3r.get("analysis_universe_spots", -1)) == 2248, "Phase 3R analysis universe spots is not 2248"),
        (int(p3r.get("n_gene_overlap_exact", 0)) >= 500, "Phase 3R exact gene overlap is <500"),
        (p3r.get("cytospace_input_manifest_ready") is True, "Phase 3R CytoSPACE manifest not ready"),
        (p2c.get("decision") == "PASS", "Phase 2C decision is not PASS"),
        (p2c.get("endpoint_frozen") is True, "Phase 2C endpoint_frozen is not true"),
        (p2c.get("primary_endpoint") == PRIMARY_ENDPOINT, "Phase 2C primary endpoint is not Immune cells"),
    ]
    for ok, msg in checks:
        if not ok:
            errors.append(msg)

    required = [
        "cytospace_input_manifest.json",
        "control_strategy_resolution.json",
        "immune_all_dropout_design.json",
        "nonimmune_all_dropout_control_design.json",
        "random_all_cell_size_matched_control_dry_design.json",
        "st_expression_full_2248_data.csv.gz",
        "st_expression_full_2248_counts.csv.gz",
        "st_coordinates_full_2248.csv",
        "formal_gene_intersection.txt",
    ]
    for name in required:
        if not (PHASE3R_DIR / name).exists():
            errors.append(f"missing Phase 3R required file: {name}")
    return p3r, p2c, errors


def write_runner_helper() -> Path:
    helper = OUT_DIR / "phase4_run_single_cytospace.py"
    helper.write_text(
        r'''
import argparse
import json
import sys
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True)
args = parser.parse_args()
cfg = json.loads(Path(args.config).read_text(encoding="utf-8"))

sys.path.insert(0, cfg["cytospace_package_path"])
# datatable.fread can mis-handle non-ASCII Windows paths in this workspace.
# Force the official CytoSPACE reader to use its pandas fallback.
import cytospace.common.common as cytospace_common
cytospace_common.dt = None
from cytospace.cytospace import main_cytospace

main_cytospace(
    scRNA_path=cfg["scRNA_path"],
    cell_type_path=cfg["cell_type_path"],
    n_cells_per_spot_path=None,
    st_cell_type_path=None,
    cell_type_fraction_estimation_path=cfg["cell_type_fraction_estimation_path"],
    spaceranger_path=None,
    st_path=cfg["st_path"],
    coordinates_path=cfg["coordinates_path"],
    output_folder=cfg["output_folder"],
    output_prefix="",
    mean_cell_numbers=cfg.get("mean_cell_numbers", 5),
    downsample_off=True,
    scRNA_max_transcripts_per_cell=1500,
    solver_method=cfg.get("solver_method", "lap_CSPR"),
    distance_metric="Pearson_correlation",
    sampling_method="duplicates",
    single_cell=False,
    number_of_selected_spots=10000,
    sampling_sub_spots=False,
    number_of_selected_sub_spots=10000,
    number_of_processors=cfg.get("number_of_processors", 1),
    seed=cfg.get("seed", 20260705),
    plot_off=True,
    geometry="honeycomb",
    max_num_cells_plot=50000,
    num_column=3,
)
'''.strip()
        + "\n",
        encoding="utf-8",
    )
    return helper


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, list[str], pd.DataFrame]:
    ref_expr = pd.read_csv(PHASE3R_DIR / "cytospace_input_manifest.json")
    raise RuntimeError("unreachable")


def prepare_common_input() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    sc_expr_path = ROOT / read_json(PHASE3R_DIR / "cytospace_input_manifest.json")["reference_expression"]
    sc_meta_path = ROOT / read_json(PHASE3R_DIR / "cytospace_input_manifest.json")["reference_metadata"]
    genes = [g.strip() for g in (PHASE3R_DIR / "formal_gene_intersection.txt").read_text(encoding="utf-8").splitlines() if g.strip()]
    st_expr = pd.read_csv(PHASE3R_DIR / "st_expression_full_2248_data.csv.gz", index_col=0)
    st_expr = st_expr.loc[genes]
    st_coords_full = pd.read_csv(PHASE3R_DIR / "st_coordinates_full_2248.csv")
    if {"barcode", "row", "col"}.issubset(st_coords_full.columns):
        st_coords = st_coords_full[["barcode", "row", "col"]].copy()
    elif {"barcode", "imagerow", "imagecol"}.issubset(st_coords_full.columns):
        st_coords = st_coords_full[["barcode", "imagerow", "imagecol"]].copy()
        st_coords = st_coords.rename(columns={"imagerow": "row", "imagecol": "col"})
    else:
        raise ValueError("ST coordinates must contain barcode,row,col or barcode,imagerow,imagecol")
    st_coords = st_coords.set_index("barcode")
    st_coords = st_coords.loc[st_expr.columns]
    sc_expr_cells_by_genes = pd.read_csv(sc_expr_path)
    sc_meta = pd.read_csv(sc_meta_path)
    if "cell_id" not in sc_expr_cells_by_genes.columns or "cell_id" not in sc_meta.columns:
        raise ValueError("reference expression and metadata must contain cell_id")
    missing = [g for g in genes if g not in sc_expr_cells_by_genes.columns]
    if missing:
        raise ValueError(f"formal genes missing from scRNA expression: {len(missing)}")
    return sc_expr_cells_by_genes, sc_meta, st_expr, genes, st_coords


def subset_reference(
    run_name: str,
    sc_expr_cells_by_genes: pd.DataFrame,
    sc_meta: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if run_name == "baseline_full_reference":
        meta = sc_meta.copy()
    elif run_name == "baseline_immune_all_dropout":
        meta = sc_meta.loc[~sc_meta[LABEL_COLUMN].isin(IMMUNE_LABELS)].copy()
    elif run_name == "baseline_nonimmune_all_dropout_control":
        meta = sc_meta.loc[sc_meta[LABEL_COLUMN].isin(IMMUNE_LABELS)].copy()
    else:
        raise ValueError(f"unknown run: {run_name}")
    expr = sc_expr_cells_by_genes.loc[sc_expr_cells_by_genes["cell_id"].isin(meta["cell_id"])].copy()
    expr = expr.set_index("cell_id").loc[meta["cell_id"]]
    return expr, meta


def write_run_inputs(
    run_dir: Path,
    run_name: str,
    sc_expr_cells_by_genes: pd.DataFrame,
    sc_meta: pd.DataFrame,
    st_expr: pd.DataFrame,
    st_coords: pd.DataFrame,
    genes: list[str],
) -> dict[str, Any]:
    inputs_dir = run_dir / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    ref_expr_subset, ref_meta_subset = subset_reference(run_name, sc_expr_cells_by_genes, sc_meta)
    ref_expr_gene_by_cell = ref_expr_subset[genes].transpose()
    ref_expr_path = inputs_dir / "reference_expression_used.csv"
    ref_expr_gz = inputs_dir / "reference_expression_used.csv.gz"
    ref_expr_gene_by_cell.to_csv(ref_expr_path)
    with ref_expr_path.open("rb") as src, gzip.open(ref_expr_gz, "wb") as dst:
        shutil.copyfileobj(src, dst)

    ref_meta_used = ref_meta_subset[["cell_id", LABEL_COLUMN]].copy()
    ref_meta_used.to_csv(run_dir / "reference_metadata_used.csv", index=False)
    ref_meta_used.set_index("cell_id").to_csv(inputs_dir / "cell_type_labels.csv")
    (run_dir / "reference_cell_list_used.txt").write_text("\n".join(ref_meta_used["cell_id"].astype(str)) + "\n", encoding="utf-8")
    counts = ref_meta_used[LABEL_COLUMN].value_counts().rename_axis("cell_type").reset_index(name="n_cells")
    counts["fraction"] = counts["n_cells"] / len(ref_meta_used)
    counts.to_csv(run_dir / "reference_label_counts_used.csv", index=False)
    fractions = pd.DataFrame([counts.set_index("cell_type")["fraction"]])
    fractions.index = ["fraction"]
    fractions.to_csv(inputs_dir / "cell_type_fractions.csv")

    st_path = inputs_dir / "st_expression_used.csv"
    st_expr.to_csv(st_path)
    coords_path = inputs_dir / "st_coordinates_used.csv"
    st_coords.to_csv(coords_path)
    (run_dir / "st_spot_list_used.txt").write_text("\n".join(st_expr.columns.astype(str)) + "\n", encoding="utf-8")
    (run_dir / "gene_list_used.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")

    return {
        "scRNA_path": str(ref_expr_path.resolve()),
        "cell_type_path": str((inputs_dir / "cell_type_labels.csv").resolve()),
        "cell_type_fraction_estimation_path": str((inputs_dir / "cell_type_fractions.csv").resolve()),
        "st_path": str(st_path.resolve()),
        "coordinates_path": str(coords_path.resolve()),
        "reference_expression_used": rel(ref_expr_gz),
        "reference_metadata_used": rel(run_dir / "reference_metadata_used.csv"),
        "n_reference_cells_input": int(len(ref_meta_used)),
        "n_reference_labels_input": int(ref_meta_used[LABEL_COLUMN].nunique()),
        "reference_labels_input": sorted(ref_meta_used[LABEL_COLUMN].astype(str).unique().tolist()),
        "n_ST_spots_input": int(st_expr.shape[1]),
        "n_genes_input": int(len(genes)),
    }


def detect_output_role(path: Path) -> str:
    name = path.name.lower()
    if name == "assigned_locations.csv":
        return "cytospace_assignment"
    if "assigned_expression" in str(path).lower():
        return "cytospace_assigned_expression"
    if "log" in name:
        return "cytospace_log"
    if "cell" in name and "spot" in name:
        return "cytospace_cell_to_spot"
    if "fraction" in name or "composition" in name:
        return "cytospace_spot_composition"
    return "unknown"


def inventory_outputs(run_name: str, run_dir: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    output_dir = run_dir / "cytospace_output"
    if not output_dir.exists():
        return rows
    for path in sorted(output_dir.rglob("*")):
        if not path.is_file():
            continue
        n_rows, n_cols, readable = (None, None, False)
        if path.suffix.lower() in {".csv", ".txt", ".tsv"}:
            n_rows, n_cols, readable = count_rows_cols(path)
        rows.append(
            {
                "run_name": run_name,
                "file_path": rel(path),
                "file_type": path.suffix.lstrip("."),
                "file_size_mb": round(path.stat().st_size / (1024 * 1024), 6),
                "detected_role": detect_output_role(path),
                "readable": readable,
                "n_rows": n_rows,
                "n_cols": n_cols,
                "notes": "",
            }
        )
    return rows


def run_cytospace(run_name: str, common: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str], pd.DataFrame]) -> dict[str, Any]:
    sc_expr, sc_meta, st_expr, genes, st_coords = common
    run_dir = OUT_DIR / run_name
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    output_dir = run_dir / "cytospace_output"
    output_dir.mkdir(parents=True, exist_ok=True)
    input_info = write_run_inputs(run_dir, run_name, sc_expr, sc_meta, st_expr, st_coords, genes)
    cfg = {
        **input_info,
        "run_name": run_name,
        "cytospace_package_path": str((ROOT / "external" / "cytospace").resolve()),
        "output_folder": str(output_dir.resolve()),
        "solver_method": "lap_CSPR",
        "mean_cell_numbers": 5,
        "number_of_processors": 1,
        "seed": 20260705,
    }
    cfg_path = run_dir / "cytospace_run_config.json"
    write_json(cfg_path, cfg)
    helper = write_runner_helper()
    command = [str(CYTOSPACE_PY), str(helper), "--config", str(cfg_path)]
    (run_dir / "run_command.txt").write_text(" ".join(command) + "\n", encoding="utf-8")
    env = os.environ.copy()
    env["NUMBA_DISABLE_JIT"] = "1"
    env["NUMBA_CACHE_DIR"] = str((Path(os.environ.get("TEMP", str(OUT_DIR))) / "numba_cache_svtuner").resolve())
    env["PYTHONPATH"] = str((ROOT / "external" / "cytospace").resolve())
    env["NUMEXPR_MAX_THREADS"] = "1"
    (run_dir / "run_environment.txt").write_text(
        "\n".join(
            [
                f"python={CYTOSPACE_PY}",
                "NUMBA_DISABLE_JIT=1",
                f"NUMBA_CACHE_DIR={env['NUMBA_CACHE_DIR']}",
                f"PYTHONPATH={env['PYTHONPATH']}",
                "entrypoint=external/cytospace/cytospace/cytospace.py:main_cytospace",
                "solver_method=lap_CSPR",
                "downsample_off=true",
                "plot_off=true",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    input_manifest = {
        "run_name": run_name,
        "endpoint_used_for_mapping": False,
        "endpoint_used_for_evaluation": False,
        "st_expression_input_type": "log1p_CPM_from_counts_in_Phase3R",
        **cfg,
    }
    write_json(run_dir / "input_manifest.json", input_manifest)
    pd.DataFrame([input_manifest]).to_csv(run_dir / "input_manifest.csv", index=False)

    start = time.perf_counter()
    proc = subprocess.run(command, cwd=str(ROOT), capture_output=True, text=True, env=env, timeout=1800)
    runtime = time.perf_counter() - start
    (run_dir / "run_log_stdout.txt").write_text(proc.stdout, encoding="utf-8", errors="replace")
    (run_dir / "run_log_stderr.txt").write_text(proc.stderr, encoding="utf-8", errors="replace")

    output_files = [rel(p) for p in output_dir.rglob("*") if p.is_file()]
    error_text = (proc.stderr or "") + "\n" + (proc.stdout or "")
    error_detected = proc.returncode != 0 or any(tok in error_text.lower() for tok in ["traceback", "error:", "exception"])
    completed = proc.returncode == 0 and (output_dir / "assigned_locations.csv").exists()
    status = {
        "run_name": run_name,
        "CytoSPACE_run": True,
        "SVTuner_run": False,
        "Stage4_run": False,
        "endpoint_used_for_mapping": False,
        "endpoint_used_for_evaluation": False,
        "n_ST_spots_input": input_info["n_ST_spots_input"],
        "n_genes_input": input_info["n_genes_input"],
        "n_reference_cells_input": input_info["n_reference_cells_input"],
        "n_reference_labels_input": input_info["n_reference_labels_input"],
        "reference_labels_input": input_info["reference_labels_input"],
        "return_code": proc.returncode,
        "runtime_seconds": round(runtime, 3),
        "completed": completed,
        "output_files_detected": output_files,
        "error_detected": error_detected,
        "error_summary": None if not error_detected else "non-zero return code or error text detected in logs",
    }
    write_json(run_dir / "run_status.json", status)
    return status


def prepare_random_sensitivity(sc_meta: pd.DataFrame) -> dict[str, Any]:
    seed = 20260705
    n_remove = 2414
    rng = pd.Series(sc_meta["cell_id"].astype(str)).sample(n=n_remove, random_state=seed)
    path = OUT_DIR / "random_all_cell_size_matched_control_cell_list_seed20260705.txt"
    path.write_text("\n".join(rng.tolist()) + "\n", encoding="utf-8")
    design = {
        "design": "random_all_cell_size_matched_control_phase4_prepared",
        "seed": seed,
        "n_cells_to_remove": n_remove,
        "n_cells_remaining": int(len(sc_meta) - n_remove),
        "cell_list_path": rel(path),
        "mapping_run": False,
        "phase_to_execute": "later optional sensitivity phase after manual authorization",
    }
    write_json(OUT_DIR / "random_all_cell_size_matched_control_design_phase4_prepared.json", design)
    return {
        "random_all_cell_size_matched_sensitivity_prepared": True,
        "random_sensitivity_mapping_run": False,
        "random_sensitivity_cell_list_path": rel(path),
    }


def write_phase_outputs(summary: dict[str, Any], statuses: dict[str, dict[str, Any]], inventory: list[dict[str, Any]], qc_rows: list[dict[str, Any]]) -> None:
    write_json(OUT_DIR / "bioapp_phase4_formal_cytospace_baseline_execution_summary.json", summary)
    golden = {
        "stage_type": STAGE_TYPE,
        "biological_question_defined": True,
        "external_endpoint_predefined": True,
        "endpoint_independent_from_SVTuner": True,
        "endpoint_spatially_registered": True,
        "baseline_comparison_available": summary["all_required_baseline_runs_completed"],
        "endpoint_specific_improvement_defined": True,
        "endpoint_specific_quantitative_metric_available": False,
        "endpoint_baseline_SVTuner_spatial_comparison_available": False,
        "interpretation_boundary_defined": True,
        "biological_application_allowed": False,
        "allowed_claim_level": "formal baseline execution only",
        "decision": summary["decision"],
    }
    write_json(OUT_DIR / "bioapp_phase4_golden_rules_v2_1_check.json", golden)
    pd.DataFrame(inventory).to_csv(OUT_DIR / "baseline_output_inventory.csv", index=False)
    pd.DataFrame(qc_rows).to_csv(OUT_DIR / "baseline_execution_qc.csv", index=False)

    next_text = (
        "BioApp Phase 5 - endpoint-specific evaluation of CytoSPACE baseline outputs against frozen CTA Immune endpoint"
        if summary["decision"] == "PASS"
        else "Review CytoSPACE run logs, output inventory, and input compatibility before endpoint-specific evaluation.\nDo not run SVTuner/SVTuner Stage4/prevention analysis."
        if summary["decision"] == "REVIEW_REQUIRED"
        else "Stop. Fix CytoSPACE baseline execution or boundary violation before retrying Phase 4."
    )
    lines = [
        "BioApp Phase 4 - formal CytoSPACE baseline execution against frozen CTA Immune endpoint",
        "",
        f"Decision: {summary['decision']}",
        "",
        "Input endpoint:",
        "CTA-defined Immune cells",
        "",
        "Endpoint status:",
        "Frozen at Phase 2C",
        "",
        "Primary endpoint note:",
        "Sparse immune-positive spatial compartments, not a large continuous ROI.",
        "",
        "Analysis universe:",
        "full_frozen_endpoint_spot_universe",
        "spots = 2248",
        "genes = 2000",
        "",
        "Baseline runs:",
    ]
    for name, status in statuses.items():
        lines.append(f"{name}: completed = {status.get('completed')}")
    lines += [
        "",
        "Control strategy:",
        "true non-immune size-matched control possible = false",
        "primary control = nonimmune-all-dropout control",
        "supplementary control = random-all-cell size-matched dropout dry/sensitivity design",
        "random sensitivity mapping run = false",
        "",
        "Boundary checks:",
        "CytoSPACE run: true",
        "SVTuner run: false",
        "Stage4 run: false",
        "SVTuner Stage3 run: false",
        "SVTuner Stage4 run: false",
        "Contradiction analysis run: false",
        "Prevention analysis run: false",
        "Endpoint redefined: false",
        "Baseline-vs-endpoint metric computed: false",
        "SVTuner-vs-endpoint metric computed: false",
        "",
        "Allowed claims:",
        *[f"- {x}" for x in ALLOWED_CLAIMS],
        "",
        "Disallowed claims:",
        *[f"- {x}" for x in DISALLOWED_CLAIMS],
        "",
        "Next:",
        next_text,
    ]
    (OUT_DIR / "decision.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    readme = f"""# BioApp Phase 4

Purpose: execute formal CytoSPACE baseline runs against the frozen CTA-defined Immune cells endpoint input universe.

Inputs come from Phase 2C and Phase 3R. The endpoint is only an anchor for later evaluation and was not used for mapping.

ST expression uses the full 2248 frozen endpoint spot universe. The CytoSPACE ST data matrix is `log1p_CPM_from_counts_in_Phase3R`; it was not used to define or select the endpoint.

Exact-overlap genes used: 2000.

Baseline runs:
- baseline_full_reference
- baseline_immune_all_dropout
- baseline_nonimmune_all_dropout_control

Control limitation: true non-immune size-matched control is impossible because non-immune pool size is 1600 and immune dropout size is 2414. The primary available label-level clean control is nonimmune-all-dropout.

This phase does not compute endpoint metrics and does not run SVTuner, Stage3, or Stage4.

Decision: `{summary['decision']}`
"""
    (OUT_DIR / "README.md").write_text(readme, encoding="utf-8")
    manifest_rows = []
    for path in sorted(OUT_DIR.rglob("*")):
        if path.is_file():
            manifest_rows.append(
                {
                    "file": rel(path),
                    "type": path.suffix.lstrip(".") or "text",
                    "description": "BioApp Phase 4 artifact",
                    "created_by_phase": PHASE,
                    "status": "generated",
                    "notes": "baseline execution artifact; no endpoint metric",
                }
            )
    pd.DataFrame(manifest_rows).to_csv(OUT_DIR / "manifest.csv", index=False)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    p3r, p2c, input_errors = validate_inputs()
    statuses: dict[str, dict[str, Any]] = {}
    inventory: list[dict[str, Any]] = []
    qc_rows: list[dict[str, Any]] = []
    boundary_violation = False
    decision_reasons: list[str] = []
    manual_authorization = True

    if input_errors:
        decision = "FAIL"
        decision_reasons = input_errors
        common = None
        random_info = {
            "random_all_cell_size_matched_sensitivity_prepared": False,
            "random_sensitivity_mapping_run": False,
        }
    else:
        common = prepare_common_input()
        _, sc_meta, _, _, _ = common
        random_info = prepare_random_sensitivity(sc_meta)
        for run_name in [
            "baseline_full_reference",
            "baseline_immune_all_dropout",
            "baseline_nonimmune_all_dropout_control",
        ]:
            status = run_cytospace(run_name, common)
            statuses[run_name] = status
            inv = inventory_outputs(run_name, OUT_DIR / run_name)
            inventory.extend(inv)
            output_detected = len(inv) > 0 and any(x["detected_role"] == "cytospace_assignment" for x in inv)
            ready = (
                status["completed"]
                and output_detected
                and status["n_ST_spots_input"] == 2248
                and status["n_genes_input"] == 2000
                and not boundary_violation
            )
            qc_rows.append(
                {
                    "run_name": run_name,
                    "completed": status["completed"],
                    "return_code": status["return_code"],
                    "runtime_seconds": status["runtime_seconds"],
                    "n_ST_spots_input": status["n_ST_spots_input"],
                    "n_genes_input": status["n_genes_input"],
                    "n_reference_cells_input": status["n_reference_cells_input"],
                    "n_reference_labels_input": status["n_reference_labels_input"],
                    "output_detected": output_detected,
                    "error_detected": status["error_detected"],
                    "ready_for_phase5_endpoint_evaluation": ready,
                    "notes": "" if ready else "run incomplete, missing assignment output, or universe mismatch",
                }
            )
        all_completed = all(s.get("completed") for s in statuses.values()) and all(r["ready_for_phase5_endpoint_evaluation"] for r in qc_rows)
        decision = "PASS" if all_completed else "REVIEW_REQUIRED"
        decision_reasons = (
            ["all required baseline runs completed and raw outputs detected"]
            if all_completed
            else ["one or more baseline runs failed, output missing, or QC not ready"]
        )

    baseline_runs_summary = {}
    for name in [
        "baseline_full_reference",
        "baseline_immune_all_dropout",
        "baseline_nonimmune_all_dropout_control",
    ]:
        s = statuses.get(name, {})
        baseline_runs_summary[name] = {
            "attempted": bool(s),
            "completed": s.get("completed", False),
            "n_reference_cells_input": s.get("n_reference_cells_input"),
            "output_detected": any(i["run_name"] == name and i["detected_role"] == "cytospace_assignment" for i in inventory),
            "ready_for_phase5_endpoint_evaluation": any(
                q["run_name"] == name and q["ready_for_phase5_endpoint_evaluation"] for q in qc_rows
            ),
        }
    all_required_completed = decision == "PASS"
    summary = {
        "phase": PHASE,
        "decision": decision,
        "decision_reasons": decision_reasons,
        "stage_type": STAGE_TYPE,
        "input_phase2c_decision": p2c.get("decision") if p2c else None,
        "input_phase3r_decision": p3r.get("decision") if p3r else None,
        "manual_authorization_for_phase4": manual_authorization,
        "endpoint_frozen": True,
        "primary_endpoint": PRIMARY_ENDPOINT,
        "primary_endpoint_spatial_note": PRIMARY_ENDPOINT_NOTE,
        "endpoint_used_for_mapping": False,
        "endpoint_used_for_evaluation": False,
        "analysis_universe_type": "full_frozen_endpoint_spot_universe",
        "analysis_universe_spots": 2248,
        "n_gene_overlap_used": 2000,
        "ST_expression_input_type": "log1p_CPM_from_counts_in_Phase3R",
        "ST_expression_input_path": rel(PHASE3R_DIR / "st_expression_full_2248_data.csv.gz"),
        "ST_expression_note": "The data layer was generated from counts in Phase 3R and was not used to define the endpoint.",
        "reference_label_column": LABEL_COLUMN,
        "total_reference_cells": 4014,
        "immune_reference_cells": 2414,
        "nonimmune_reference_cells": 1600,
        "control_strategy": {
            "true_nonimmune_size_matched_control_possible": False,
            "primary_control": "nonimmune-all-dropout control",
            "primary_control_size_matched": False,
            "supplementary_control": "random-all-cell size-matched dropout dry/sensitivity design",
            "random_sensitivity_mapping_run": False,
        },
        **random_info,
        "baseline_runs": baseline_runs_summary,
        "all_required_baseline_runs_completed": all_required_completed,
        "ready_for_phase5_endpoint_evaluation": all_required_completed,
        "CytoSPACE_run": bool(statuses),
        "SVTuner_run": False,
        "Stage4_run": False,
        "SVTuner_Stage3_run": False,
        "SVTuner_Stage4_run": False,
        "contradiction_analysis_run": False,
        "prevention_analysis_run": False,
        "endpoint_redefined": False,
        "expression_markers_used_to_define_endpoint": False,
        "mapping_outputs_used_to_define_endpoint": False,
        "baseline_vs_endpoint_metric_computed": False,
        "SVTuner_vs_endpoint_metric_computed": False,
        "boundary_violation": boundary_violation,
        "biological_application_allowed": False,
        "cytospace_entrypoint_found": CYTOSPACE_PY.exists(),
        "cytospace_entrypoint_path": "external/cytospace/cytospace/cytospace.py:main_cytospace",
        "cytospace_environment": str(CYTOSPACE_PY),
        "command_template": "python phase4_run_single_cytospace.py --config <run_config>",
        "actual_command_executed": {k: (OUT_DIR / k / "run_command.txt").read_text(encoding="utf-8").strip() if (OUT_DIR / k / "run_command.txt").exists() else None for k in baseline_runs_summary},
        "allowed_claims": ALLOWED_CLAIMS,
        "disallowed_claims": DISALLOWED_CLAIMS,
        "output_dir": rel(OUT_DIR),
    }
    write_phase_outputs(summary, statuses, inventory, qc_rows)

    print("BioApp Phase 4 completed.")
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
    print("Analysis universe:")
    print("2248 frozen endpoint spots")
    print("2000 exact-overlap genes")
    print()
    print("Baseline runs:")
    for name, s in baseline_runs_summary.items():
        print(f"{name}: attempted={s['attempted']}, completed={s['completed']}, output_detected={s['output_detected']}")
    print()
    print("Control strategy:")
    print("true_nonimmune_size_matched_control_possible = false")
    print("primary_control = nonimmune-all-dropout control")
    print("supplementary_control = random-all-cell size-matched dropout dry/sensitivity design")
    print("random_sensitivity_mapping_run = false")
    print()
    print("Boundary checks:")
    print("CytoSPACE run: true")
    print("SVTuner run: false")
    print("Stage4 run: false")
    print("SVTuner Stage3 run: false")
    print("SVTuner Stage4 run: false")
    print("Contradiction analysis run: false")
    print("Prevention analysis run: false")
    print("Endpoint redefined: false")
    print("Baseline-vs-endpoint metric computed: false")
    print("SVTuner-vs-endpoint metric computed: false")
    print()
    print("Ready for Phase 5 endpoint evaluation:")
    print(str(all_required_completed).lower())
    print()
    print("Next:")
    if decision == "PASS":
        print("BioApp Phase 5 - endpoint-specific evaluation of CytoSPACE baseline outputs against frozen CTA Immune endpoint")
    elif decision == "REVIEW_REQUIRED":
        print("Review CytoSPACE run logs, output inventory, and input compatibility before endpoint-specific evaluation.")
        print("Do not run SVTuner/SVTuner Stage4/prevention analysis.")
    else:
        print("Stop. Fix CytoSPACE baseline execution or boundary violation before retrying Phase 4.")
    return 0 if decision in {"PASS", "REVIEW_REQUIRED"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
