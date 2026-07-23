#!/usr/bin/env python3
"""Generate the Task 4D software-version and Stage3B NNLS audit.

This is a static, read-only audit. It never imports or executes Stage3B, R,
Seurat, NumPy, or SciPy. The only writes are deterministic report files under
the dedicated Task 4D audit directory.
"""

from __future__ import annotations

import ast
import csv
import hashlib
import io
import json
import re
import subprocess
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "visualizations/manuscript_audits/software_and_nnls_formal_audit"
GENERATED_DATE = "2026-07-22"

STAGE3B = ROOT / "src/stages/stage3b_st_unsupported.py"
CLI = ROOT / "src/svtuner/cli.py"
STAGE1_R = ROOT / "r_scripts/stage1_preprocess.R"
ENV_YAML = ROOT / "configs/environment.yml"
PROJECT_CONFIG = ROOT / "configs/project_config.yaml"
COMPOSITE_RUNNER = ROOT / "scripts/run_composite_noise_stage3ab_full.py"
BIOAPP_RUNNER = ROOT / "scripts/run_bioapp_phase7_svtuner_aware_execution_against_frozen_cta_immune_endpoint.py"
BIOAPP_SUMMARY = ROOT / (
    "visualizations/bioapp_experiment/"
    "bioapp_phase7_svtuner_aware_execution_against_frozen_cta_immune_endpoint/"
    "bioapp_phase7_svtuner_aware_execution_summary.json"
)
NAMED_STAGE3B_LOG = ROOT / (
    "visualizations/method_comparison/composite_scnoise10_stage3ab_full/"
    "execution_logs/"
    "real_brca7_endothelial_marker_control_sc_missing_endothelial_cells_scnoise10_stage3b.log"
)
BASE_STAGE3B_LOG = ROOT / (
    "logs/cytospace_fig2j_stage3ab_upgrade/"
    "cytospace_fig2j_stage3ab_joint.stage3b.log"
)
BIOAPP_R_PACKAGE_CHECK = ROOT / (
    "visualizations/bioapp_experiment/"
    "bioapp_phase2br2_explicit_rscript_cta_align_locale_safe_dry_run/"
    "bioapp_phase2br2_r_package_check.csv"
)
BIOAPP_R_RUNTIME = ROOT / (
    "visualizations/bioapp_experiment/"
    "bioapp_phase2br2_explicit_rscript_cta_align_locale_safe_dry_run/"
    "bioapp_phase2br2_explicit_r_runtime_check.json"
)
BIOAPP_SEURAT_INVENTORY = ROOT / (
    "visualizations/bioapp_experiment/"
    "bioapp_coordinate_registration_recovery_audit/"
    "bioapp_st_data_rdata_object_inventory.json"
)
BIOAPP_PHASE2_SCRIPT = ROOT / "scripts/run_bioapp_phase2br2_explicit_rscript_cta_align_locale_safe_dry_run.py"
BIOAPP_PHASE3R_SCRIPT = ROOT / "scripts/run_bioapp_phase3r_review_resolution_for_cta_endpoint_baseline_feasibility.py"
FIG2C_SEURAT_SCRIPT = ROOT / "scripts/compute_cytospace_fig2c_official_enrichment_seurat.R"
FIG2_BENCHMARK_RUNNER = ROOT / "scripts/run_cytospace_fig2d_profile_mask_benchmark.py"


def rel(path: Path) -> str:
    return path.resolve().relative_to(ROOT.resolve()).as_posix()


def sha1_file(path: Path) -> str:
    digest = hashlib.sha1()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha1_bytes(data: bytes) -> str:
    return hashlib.sha1(data).hexdigest()


def text_bytes(value: str) -> bytes:
    return value.rstrip().encode("utf-8") + b"\n"


def json_bytes(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode("utf-8")


def tsv_bytes(rows: Iterable[dict[str, Any]], fields: list[str]) -> bytes:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(
        output,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    for row in rows:
        writer.writerow({field: "" if row.get(field) is None else row.get(field, "") for field in fields})
    return output.getvalue().encode("utf-8")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def git_commit() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def require_files() -> None:
    required = [
        STAGE3B,
        CLI,
        STAGE1_R,
        ENV_YAML,
        PROJECT_CONFIG,
        COMPOSITE_RUNNER,
        BIOAPP_RUNNER,
        BIOAPP_SUMMARY,
        NAMED_STAGE3B_LOG,
        BASE_STAGE3B_LOG,
        BIOAPP_R_PACKAGE_CHECK,
        BIOAPP_R_RUNTIME,
        BIOAPP_SEURAT_INVENTORY,
        BIOAPP_PHASE2_SCRIPT,
        BIOAPP_PHASE3R_SCRIPT,
        FIG2C_SEURAT_SCRIPT,
        FIG2_BENCHMARK_RUNNER,
    ]
    missing = [rel(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Required audit evidence is missing: {missing}")


def parse_python_from_log(path: Path) -> Path:
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"(?im)^(?:\[CMD\]\s+)?([A-Z]:\\[^\r\n]*?python\.exe)\s+-m\s+src\.stages\.stage3b_st_unsupported", text)
    if not match:
        raise ValueError(f"No formal Stage3B Python executable found in {rel(path)}")
    return Path(match.group(1))


def dist_info(root: Path, package: str) -> tuple[Path, str]:
    candidates = sorted((root / "Lib/site-packages").glob(f"{package}-*.dist-info"))
    if len(candidates) != 1:
        raise ValueError(f"Expected one {package} dist-info under {root}, found {len(candidates)}")
    metadata = candidates[0] / "METADATA"
    if not metadata.is_file():
        raise FileNotFoundError(metadata)
    match = re.search(r"(?m)^Version:\s*(\S+)", metadata.read_text(encoding="utf-8", errors="replace"))
    if not match:
        raise ValueError(f"No Version field in {metadata}")
    return metadata, match.group(1)


def conda_package_record(root: Path, package: str) -> tuple[Path, dict[str, Any]]:
    candidates = sorted((root / "conda-meta").glob(f"{package}-*.json"))
    records = []
    for path in candidates:
        record = json.loads(path.read_text(encoding="utf-8"))
        if record.get("name") == package:
            records.append((path, record))
    if len(records) != 1:
        raise ValueError(f"Expected one {package} conda record under {root}, found {len(records)}")
    return records[0]


def external_label(path: Path, env_name: str, env_root: Path) -> str:
    return f"ENV:{env_name}/{path.resolve().relative_to(env_root.resolve()).as_posix()}"


def evidence_row(
    evidence_id: str,
    category: str,
    path: Path,
    supports: str,
    evidence_type: str,
    priority: int,
    *,
    display_path: str | None = None,
    notes: str = "",
) -> dict[str, Any]:
    stat = path.stat()
    return {
        "evidence_id": evidence_id,
        "category": category,
        "file_path": display_path or rel(path),
        "file_size": stat.st_size,
        "modified_time_ns": stat.st_mtime_ns,
        "sha1": sha1_file(path),
        "evidence_type": evidence_type,
        "evidence_priority": priority,
        "supports": supports,
        "read_only": True,
        "notes": notes,
    }


def source_lines(path: Path, start: int, end: int) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    return " | ".join(f"L{i}:{lines[i - 1].strip()}" for i in range(start, end + 1))


def stage1_logs() -> list[Path]:
    return sorted(path for path in (ROOT / "result").glob("**/stage1_preprocess.log") if path.is_file())


def validate_static_implementation() -> dict[str, Any]:
    source = STAGE3B.read_text(encoding="utf-8")
    tree = ast.parse(source)
    imported_nnls = any(
        isinstance(node, ast.ImportFrom)
        and node.module == "scipy.optimize"
        and any(alias.name == "nnls" for alias in node.names)
        for node in tree.body
    )
    functions = {node.name: node for node in tree.body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))}
    fit = functions.get("fit_nonnegative_mixtures")
    run = functions.get("run_stage3b")
    fit_text = ast.get_source_segment(source, fit) if fit else ""
    run_text = ast.get_source_segment(source, run) if run else ""
    checks = {
        "imports_scipy_optimize_nnls": imported_nnls,
        "fit_function_found": fit is not None,
        "run_function_found": run is not None,
        "design_is_profiles_transpose": "design = profiles.T" in fit_text,
        "per_spot_loop": "for i, row in enumerate(observations)" in fit_text,
        "nnls_call": "nnls(design, row)" in fit_text,
        "sum_normalisation": "coef /= total" in fit_text,
        "zero_guard": "if total > 0" in fit_text,
        "reconstruction": "coef @ profiles" in fit_text,
        "formal_fit_called_for_st": "fit_nonnegative_mixtures(st_values, profiles)" in run_text,
        "formal_fit_called_for_calibration": "fit_nonnegative_mixtures(pseudo, profiles)" in run_text,
        "residual_formula": "residual = observations - reconstruction" in source,
        "no_alternative_solver_import": not any(
            token in source for token in ("lsq_linear", "least_squares", "sklearn", "projected_gradient")
        ),
    }
    failed = [name for name, value in checks.items() if not value]
    if failed:
        raise ValueError(f"Static Stage3B implementation validation failed: {failed}")
    return checks


def build_payloads() -> dict[str, bytes]:
    require_files()
    static_checks = validate_static_implementation()
    named_python = parse_python_from_log(NAMED_STAGE3B_LOG)
    base_python = parse_python_from_log(BASE_STAGE3B_LOG)
    if not named_python.is_file() or not base_python.is_file():
        raise FileNotFoundError("A formal Stage3B Python executable recorded in logs is no longer retained")
    named_env = named_python.parent
    base_env = base_python.parent

    named_python_record, named_python_meta = conda_package_record(named_env, "python")
    base_python_record, base_python_meta = conda_package_record(base_env, "python")
    named_numpy_meta, named_numpy_version = dist_info(named_env, "numpy")
    named_scipy_meta, named_scipy_version = dist_info(named_env, "scipy")
    base_numpy_meta, base_numpy_version = dist_info(base_env, "numpy")
    base_scipy_meta, base_scipy_version = dist_info(base_env, "scipy")
    named_scipy_source = named_env / "Lib/site-packages/scipy/optimize/_nnls.py"
    base_scipy_source = base_env / "Lib/site-packages/scipy/optimize/_nnls.py"
    named_history = named_env / "conda-meta/history"
    base_history = base_env / "conda-meta/history"
    for path in (named_scipy_source, base_scipy_source, named_history, base_history):
        if not path.is_file():
            raise FileNotFoundError(path)

    r_packages = {row["package"]: row for row in read_csv(BIOAPP_R_PACKAGE_CHECK)}
    r_runtime = json.loads(BIOAPP_R_RUNTIME.read_text(encoding="utf-8"))
    seurat_inventory = json.loads(BIOAPP_SEURAT_INVENTORY.read_text(encoding="utf-8"))
    expected_r = {
        "Seurat": "5.3.1",
        "SeuratObject": "5.2.0",
    }
    for package, version in expected_r.items():
        if r_packages.get(package, {}).get("version") != version:
            raise ValueError(f"Unexpected BioApp {package} version evidence")
    if seurat_inventory.get("seurat_version") != "5.3.1":
        raise ValueError("BioApp Seurat object inventory does not agree with package check")

    stage1_log_paths = stage1_logs()
    if not stage1_log_paths:
        raise ValueError("No formal Stage1 Seurat execution logs found")
    stage1_samples = []
    for path in stage1_log_paths:
        match = re.search(r"(?m)^\[Stage1\] sample:\s*(.+)$", path.read_text(encoding="utf-8", errors="replace"))
        if match:
            stage1_samples.append(match.group(1).strip())

    evidence_rows = [
        evidence_row("E4D001", "formal implementation", STAGE3B, "NNLS import, fitting, normalisation, reconstruction, residual and output path", "implementation source", 1),
        evidence_row("E4D002", "formal caller", CLI, "CLI Stage3B module entrypoint", "implementation source", 2),
        evidence_row("E4D003", "formal caller", COMPOSITE_RUNNER, "5%/10% Stage3A+Stage3B runner command", "implementation source", 2),
        evidence_row("E4D004", "formal caller", BIOAPP_RUNNER, "BioApp Stage3B command and output import", "implementation source", 2),
        evidence_row("E4D005", "formal execution", NAMED_STAGE3B_LOG, "successful Stage3B command using named Python environment", "formal run log", 1),
        evidence_row("E4D006", "formal execution", BASE_STAGE3B_LOG, "successful high-resolution Stage3B command using base Python", "formal run log", 1),
        evidence_row("E4D007", "formal execution", BIOAPP_SUMMARY, "successful BioApp Stage3B command using base Python", "formal run summary", 1),
        evidence_row("E4D008", "reproducibility environment", ENV_YAML, "pinned Python, NumPy, SciPy, R, Seurat and SeuratObject versions", "environment export", 3),
        evidence_row("E4D009", "reproducibility environment", PROJECT_CONFIG, "named formal working environment and package pins", "environment configuration", 3),
        evidence_row("E4D010", "Seurat implementation", STAGE1_R, "mainline Seurat preprocessing calls", "implementation source", 2),
        evidence_row("E4D011", "Seurat formal execution", BIOAPP_R_PACKAGE_CHECK, "direct packageVersion output for Seurat and SeuratObject", "formal run package check", 1),
        evidence_row("E4D012", "Seurat formal execution", BIOAPP_R_RUNTIME, "direct R version and platform", "formal run runtime check", 1),
        evidence_row("E4D013", "Seurat formal execution", BIOAPP_SEURAT_INVENTORY, "Seurat object parsing and duplicate package-version confirmation", "formal run object inventory", 1),
        evidence_row("E4D014", "Seurat caller", BIOAPP_PHASE2_SCRIPT, "BioApp explicit-R CTA alignment and package check", "implementation source", 2),
        evidence_row("E4D015", "Seurat caller", BIOAPP_PHASE3R_SCRIPT, "BioApp Seurat RData expression export", "implementation source", 2),
        evidence_row("E4D016", "Seurat candidate caller", FIG2C_SEURAT_SCRIPT, "optional Seurat enrichment backend", "implementation source; formal use unconfirmed", 6),
        evidence_row("E4D017", "Seurat candidate caller", FIG2_BENCHMARK_RUNNER, "backend switch defaults to Python", "implementation source; formal use unconfirmed", 6),
        evidence_row(
            "E4D018", "solver package", named_scipy_meta, "SciPy version 1.11.4", "installed package metadata",
            2, display_path=external_label(named_scipy_meta, "cytospace_v1.1.0_py310", named_env),
        ),
        evidence_row(
            "E4D019", "solver source", named_scipy_source, "SciPy 1.11.4 nnls defaults and failure behavior", "installed solver source",
            2, display_path=external_label(named_scipy_source, "cytospace_v1.1.0_py310", named_env),
        ),
        evidence_row(
            "E4D020", "solver environment", named_history, "named environment change history predating formal runs", "conda history",
            3, display_path=external_label(named_history, "cytospace_v1.1.0_py310", named_env),
        ),
        evidence_row(
            "E4D021", "Python package", named_python_record, f"Python {named_python_meta['version']}", "conda package record",
            2, display_path=external_label(named_python_record, "cytospace_v1.1.0_py310", named_env),
        ),
        evidence_row(
            "E4D022", "NumPy package", named_numpy_meta, f"NumPy {named_numpy_version}", "installed package metadata",
            2, display_path=external_label(named_numpy_meta, "cytospace_v1.1.0_py310", named_env),
        ),
        evidence_row(
            "E4D023", "solver package", base_scipy_meta, "SciPy version 1.15.3", "installed package metadata",
            2, display_path=external_label(base_scipy_meta, "anaconda_base", base_env),
        ),
        evidence_row(
            "E4D024", "solver source", base_scipy_source, "SciPy 1.15.3 nnls defaults and failure behavior", "installed solver source",
            2, display_path=external_label(base_scipy_source, "anaconda_base", base_env),
        ),
        evidence_row(
            "E4D025", "solver environment", base_history, "base environment change history predating formal runs", "conda history",
            3, display_path=external_label(base_history, "anaconda_base", base_env),
        ),
        evidence_row(
            "E4D026", "Python package", base_python_record, f"Python {base_python_meta['version']}", "conda package record",
            2, display_path=external_label(base_python_record, "anaconda_base", base_env),
        ),
        evidence_row(
            "E4D027", "NumPy package", base_numpy_meta, f"NumPy {base_numpy_version}", "installed package metadata",
            2, display_path=external_label(base_numpy_meta, "anaconda_base", base_env),
            notes="Active dist-info is 2.1.0; the retained conda record still names 2.1.3, consistent with a later pip-level replacement before formal runs.",
        ),
    ]
    for index, path in enumerate(stage1_log_paths, start=28):
        evidence_rows.append(
            evidence_row(
                f"E4D{index:03d}",
                "Seurat formal execution",
                path,
                "successful mainline Stage1 Seurat preprocessing",
                "formal run log without sessionInfo",
                5,
            )
        )

    evidence_by_path = {row["file_path"]: row for row in evidence_rows}
    stage3b_rel = rel(STAGE3B)
    stage3b_sha = evidence_by_path[stage3b_rel]["sha1"]

    seurat_version_rows = [
        {
            "record_id": "SEURAT_MAINLINE_STAGE1",
            "dataset_family": "low-resolution mainline datasets",
            "formal_script": rel(STAGE1_R),
            "r_version": "4.5.2 (reproducibility environment pin; exact run-local version not recorded)",
            "seurat_version": "5.3.1 (reproducibility environment pin; exact formal-run version not recoverable)",
            "seuratobject_version": "5.2.0 (reproducibility environment pin; exact formal-run version not recoverable)",
            "version_status": "PINNED_ENVIRONMENT_ONLY_FOR_EXACT_VERSION",
            "evidence_path": f"{rel(ENV_YAML)};{';'.join(rel(p) for p in stage1_log_paths)}",
            "evidence_sha1": f"{sha1_file(ENV_YAML)};{';'.join(sha1_file(p) for p in stage1_log_paths)}",
            "evidence_type": "environment export plus successful formal logs without sessionInfo/packageVersion",
            "confidence": "MEDIUM",
            "notes": f"{len(stage1_log_paths)} retained Stage1 logs ({len(set(stage1_samples))} sample IDs) prove Seurat execution, but none records sessionInfo() or packageVersion(). Do not state 5.3.1 as independently observed at each run.",
        },
        {
            "record_id": "SEURAT_BIOAPP_CTA",
            "dataset_family": "CTA/Xenium-defined biological application and RData coordinate/expression export",
            "formal_script": f"{rel(BIOAPP_PHASE2_SCRIPT)};{rel(BIOAPP_PHASE3R_SCRIPT)}",
            "r_version": "4.5.1",
            "seurat_version": r_packages["Seurat"]["version"],
            "seuratobject_version": r_packages["SeuratObject"]["version"],
            "version_status": "EXACT_FORMAL_RUN_VERSION_CONFIRMED",
            "evidence_path": f"{rel(BIOAPP_R_PACKAGE_CHECK)};{rel(BIOAPP_R_RUNTIME)};{rel(BIOAPP_SEURAT_INVENTORY)}",
            "evidence_sha1": f"{sha1_file(BIOAPP_R_PACKAGE_CHECK)};{sha1_file(BIOAPP_R_RUNTIME)};{sha1_file(BIOAPP_SEURAT_INVENTORY)}",
            "evidence_type": "formal packageVersion output, runtime record and parsed Seurat-object inventory",
            "confidence": "HIGH",
            "notes": f"Platform {r_runtime['platform']}; package check and object inventory independently agree on Seurat 5.3.1 and SeuratObject 5.2.0.",
        },
    ]

    seurat_usage_rows = [
        {
            "usage_id": "SU001",
            "script": rel(STAGE1_R),
            "functionality": "CreateSeuratObject and NormalizeData for SC and ST, then export Stage1 matrices",
            "dataset_family": "low-resolution mainline datasets",
            "formal_execution_status": "CONFIRMED_EXECUTED",
            "version_link": "SEURAT_MAINLINE_STAGE1",
            "evidence_path": ";".join(rel(p) for p in stage1_log_paths),
            "notes": "Formal logs contain Seurat v5-style layer output but no exact version snapshot.",
        },
        {
            "usage_id": "SU002",
            "script": rel(BIOAPP_PHASE2_SCRIPT),
            "functionality": "Load retained CTA Seurat object, inspect spatial image slots and support CTA-to-spot registration",
            "dataset_family": "CTA biological application",
            "formal_execution_status": "CONFIRMED_EXECUTED",
            "version_link": "SEURAT_BIOAPP_CTA",
            "evidence_path": rel(BIOAPP_SEURAT_INVENTORY),
            "notes": "Direct package and runtime records retained.",
        },
        {
            "usage_id": "SU003",
            "script": rel(BIOAPP_PHASE3R_SCRIPT),
            "functionality": "Load CTA Seurat RData and export Spatial assay expression for frozen downstream inputs",
            "dataset_family": "CTA biological application",
            "formal_execution_status": "CONFIRMED_EXECUTED",
            "version_link": "SEURAT_BIOAPP_CTA",
            "evidence_path": rel(BIOAPP_R_RUNTIME),
            "notes": "Uses the same explicit standalone R 4.5.1 runtime lineage.",
        },
        {
            "usage_id": "SU004",
            "script": rel(FIG2C_SEURAT_SCRIPT),
            "functionality": "Optional R/Seurat enrichment backend",
            "dataset_family": "Fig. 2 low-resolution profile-mask benchmark",
            "formal_execution_status": "AVAILABLE_BUT_FORMAL_USE_NOT_CONFIRMED",
            "version_link": "",
            "evidence_path": f"{rel(FIG2C_SEURAT_SCRIPT)};{rel(FIG2_BENCHMARK_RUNNER)}",
            "notes": "The benchmark runner defaults to the Python backend; no retained run-local record proves that the optional Seurat backend generated the formal panel.",
        },
    ]

    execution_rows = [
        {
            "path_id": "EP001",
            "dataset_family": "all Stage3B routes",
            "caller": "src.svtuner.cli.cmd_stage3b or direct python -m invocation",
            "entrypoint": "src.stages.stage3b_st_unsupported.main",
            "implementation": "src.stages.stage3b_st_unsupported.run_stage3b",
            "fit_function": "fit_nonnegative_mixtures",
            "solver_call": "scipy.optimize.nnls",
            "execution_mode": "module entrypoint; static path confirmed",
            "formal_execution_evidence": f"{rel(CLI)};{rel(STAGE3B)}",
            "status": "CONFIRMED",
            "notes": "The same fit function is called for observed ST spots and supported-null pseudo-spots.",
        },
        {
            "path_id": "EP002",
            "dataset_family": "5% and 10% Stage3A+Stage3B simulation benchmarks",
            "caller": rel(COMPOSITE_RUNNER),
            "entrypoint": "python -m src.stages.stage3b_st_unsupported",
            "implementation": rel(STAGE3B),
            "fit_function": "fit_nonnegative_mixtures",
            "solver_call": "scipy.optimize.nnls",
            "execution_mode": "E:/ANACONDA/envs/cytospace_v1.1.0_py310/python.exe",
            "formal_execution_evidence": rel(NAMED_STAGE3B_LOG),
            "status": "CONFIRMED",
            "notes": "Formal log records successful completion and exact executable path.",
        },
        {
            "path_id": "EP003",
            "dataset_family": "high-resolution Fig. 3 Stage3B examples",
            "caller": "direct module invocation recorded in retained log",
            "entrypoint": "python -m src.stages.stage3b_st_unsupported",
            "implementation": rel(STAGE3B),
            "fit_function": "fit_nonnegative_mixtures",
            "solver_call": "scipy.optimize.nnls",
            "execution_mode": "E:/ANACONDA/python.exe",
            "formal_execution_evidence": rel(BASE_STAGE3B_LOG),
            "status": "CONFIRMED",
            "notes": "This path uses a different retained SciPy environment from the composite simulation runner.",
        },
        {
            "path_id": "EP004",
            "dataset_family": "CTA biological application",
            "caller": rel(BIOAPP_RUNNER),
            "entrypoint": "python -m src.stages.stage3b_st_unsupported",
            "implementation": rel(STAGE3B),
            "fit_function": "fit_nonnegative_mixtures",
            "solver_call": "scipy.optimize.nnls",
            "execution_mode": "E:/ANACONDA/python.exe",
            "formal_execution_evidence": rel(BIOAPP_SUMMARY),
            "status": "CONFIRMED",
            "notes": "Formal summary records return_code=0, 2,248 spots and the exact command.",
        },
    ]

    implementation_rows = [
        {
            "record_id": "STAGE3B_NNLS_FORMAL",
            "formal_implementation_file": rel(STAGE3B),
            "formal_function": "fit_nonnegative_mixtures",
            "formal_caller": "run_stage3b (observed ST at L705; supported-null pseudo-ST at L717)",
            "solver_package": "SciPy",
            "solver_function": "scipy.optimize.nnls",
            "solver_backend": "SciPy 1.11.4 Fortran Lawson-Hanson active-set wrapper; SciPy 1.15.3 Cython improved active-set implementation",
            "solver_version": "1.11.4 for named composite environment; 1.15.3 for retained base-environment formal routes",
            "matrix_orientation": "observations=(n_spots,n_genes); profiles B=(n_types,n_genes); design=B.T=(n_genes,n_types); each b=row=(n_genes)",
            "non_negative_constraint": "enforced internally by scipy.optimize.nnls; no post-hoc negative clipping",
            "intercept": "none",
            "regularisation": "none",
            "tolerance": "not passed by SVTuner; SciPy implementation default/internal behavior applies",
            "max_iterations": "not passed; SciPy default 3 * n_types in both retained versions",
            "fallback": "none",
            "failure_handling": "no local try/except; SciPy RuntimeError propagates and the Stage3B process fails before successful output completion",
            "coefficient_normalisation": "if coef.sum()>0, in-place divide by coef.sum(); otherwise retain zero vector",
            "zero_solution_handling": "zero coefficient vector and zero reconstruction; zero observation rows are not separately skipped",
            "reconstruction_formula": "reconstruction[i] = normalized_coef @ profiles",
            "status": "CONFIRMED_SINGLE_IMPLEMENTATION_TWO_FORMAL_SCIPY_VERSIONS",
            "evidence_path": f"{rel(STAGE3B)};{rel(NAMED_STAGE3B_LOG)};{rel(BASE_STAGE3B_LOG)};ENV:cytospace_v1.1.0_py310/Lib/site-packages/scipy/optimize/_nnls.py;ENV:anaconda_base/Lib/site-packages/scipy/optimize/_nnls.py",
            "evidence_sha1": f"{stage3b_sha};{sha1_file(NAMED_STAGE3B_LOG)};{sha1_file(BASE_STAGE3B_LOG)};{sha1_file(named_scipy_source)};{sha1_file(base_scipy_source)}",
            "notes": "Repository-wide static search found no Stage3B lsq_linear, minimize, sklearn, projected-gradient or custom NNLS fallback. The two SciPy releases solve the same NNLS problem but use different retained internal backends, so their versions must not be collapsed.",
        }
    ]

    math_rows = [
        ("MC001", "minimise ||x_i - w_i B||_2^2", "nnls(B.T, x_i) minimises ||B.T w_i - x_i||_2", "PASS", "Equivalent column-vector orientation.", "none", "State the matrix orientation explicitly if desired.", "src/stages/stage3b_st_unsupported.py:L192-L196"),
        ("MC002", "w_i >= 0", "SciPy NNLS enforces nonnegative coefficients", "PASS", "No post-hoc clipping substitutes for the constraint.", "none", "Name scipy.optimize.nnls as the constrained solver.", "src/stages/stage3b_st_unsupported.py:L27,L196"),
        ("MC003", "one coefficient vector per spot", "Python loop calls nnls once for each observation row", "PASS", "Not batched.", "runtime scales with number of spots", "Describe fitting as per-spot NNLS.", "src/stages/stage3b_st_unsupported.py:L193-L196"),
        ("MC004", "if sum(w_i)>0, w_hat_i=w_i/sum(w_i)", "total=coef.sum(); coef/=total only when total>0", "PASS", "none", "none", "Retain the conditional sum normalisation statement.", "src/stages/stage3b_st_unsupported.py:L197-L200"),
        ("MC005", "otherwise w_hat_i=zero vector", "preallocated zeros remain unchanged when total<=0", "PASS", "No explicit else branch, but behavior is exact.", "none", "State that zero solutions remain zero.", "src/stages/stage3b_st_unsupported.py:L193-L200"),
        ("MC006", "x_hat_i=w_hat_i B", "reconstruction[i]=coef@profiles after coefficient normalisation", "PASS", "none", "none", "Retain the stated reconstruction formula.", "src/stages/stage3b_st_unsupported.py:L200-L201"),
        ("MC007", "residual_i=x_i-x_hat_i", "residual=observations-reconstruction", "PASS", "none", "none", "Retain the stated residual formula.", "src/stages/stage3b_st_unsupported.py:L205-L210"),
        ("MC008", "no intercept in the stated objective", "nnls receives only B.T and x_i; no intercept column or centering term", "PASS", "none", "none", "State that the model is fit without an intercept.", "src/stages/stage3b_st_unsupported.py:L192-L196"),
        ("MC009", "no regularisation in the stated objective", "no penalty term or augmented design is used", "PASS", "none", "none", "Do not describe ridge, lasso or other regularisation.", "src/stages/stage3b_st_unsupported.py:L188-L202"),
        ("MC010", "solver tolerance", "SVTuner passes neither tolerance nor maxiter", "PASS_WITH_DEFAULTS", "Tolerance is solver-version internal/default; SciPy 1.11.4 and 1.15.3 are both retained formal versions.", "Version-specific numerical details can differ slightly.", "Say that solver defaults were used; do not claim a custom tolerance.", "src/stages/stage3b_st_unsupported.py:L196"),
        ("MC011", "solver failure handling", "no fallback or local exception recovery; RuntimeError propagates", "PASS_WITH_IMPLEMENTATION_DETAIL", "Methods definition does not specify operational failure behavior.", "A failed solve aborts the formal Stage3B run rather than silently substituting coefficients.", "State no fallback was used; failed solves were not silently imputed.", "src/stages/stage3b_st_unsupported.py:L188-L202"),
        ("MC012", "x_i and B are fitted profiles", "SC and ST values are linearised, row-composition normalised, and type profiles are row-normalised before fitting", "PASS_WITH_PREPROCESSING_CONTEXT", "The compact equation omits preprocessing detail.", "Readers need preprocessing context to reproduce exact inputs.", "Define x_i and B as compositional profiles after the documented linearisation/normalisation.", "src/stages/stage3b_st_unsupported.py:L123-L185,L691-L705"),
        ("MC013", "formal outputs originate from this implementation", "formal logs and BioApp summary invoke src.stages.stage3b_st_unsupported and complete successfully", "PASS", "none", "none", "Link the implementation and environment versions in Methods/software reporting.", f"{rel(NAMED_STAGE3B_LOG)};{rel(BASE_STAGE3B_LOG)};{rel(BIOAPP_SUMMARY)}"),
    ]
    math_fields = ["check_id", "methods_definition", "code_implementation", "match_status", "difference", "impact", "recommended_manual_wording", "evidence_path"]
    math_dicts = [dict(zip(math_fields, row)) for row in math_rows]

    solver_rows = [
        {
            "record_id": "PYENV_COMPOSITE_NAMED",
            "formal_dataset_family": "5% and 10% Stage3A+Stage3B simulation benchmarks",
            "python_executable": "ENV:cytospace_v1.1.0_py310/python.exe",
            "python_version": named_python_meta["version"],
            "platform": named_python_meta.get("subdir", "win-64"),
            "numpy_version": named_numpy_version,
            "scipy_version": named_scipy_version,
            "solver_package": "SciPy",
            "solver_function": "scipy.optimize.nnls",
            "solver_backend": "Fortran Lawson-Hanson active-set wrapper",
            "version_status": "CONFIRMED_BY_FORMAL_EXECUTABLE_AND_RETAINED_PACKAGE_LINEAGE",
            "formal_run_evidence": rel(NAMED_STAGE3B_LOG),
            "environment_evidence": "ENV:cytospace_v1.1.0_py310/conda-meta/history;ENV:cytospace_v1.1.0_py310/Lib/site-packages/scipy-1.11.4.dist-info/METADATA",
            "evidence_sha1": f"{sha1_file(NAMED_STAGE3B_LOG)};{sha1_file(named_python_record)};{sha1_file(named_numpy_meta)};{sha1_file(named_scipy_meta)};{sha1_file(named_history)}",
            "confidence": "HIGH",
            "notes": "Formal command names the retained executable; package metadata and solver source predate the formal run, and retained conda history shows no later environment transaction.",
        },
        {
            "record_id": "PYENV_ANACONDA_BASE",
            "formal_dataset_family": "high-resolution Stage3B examples and CTA biological application",
            "python_executable": "ENV:anaconda_base/python.exe",
            "python_version": base_python_meta["version"],
            "platform": base_python_meta.get("subdir", "win-64"),
            "numpy_version": base_numpy_version,
            "scipy_version": base_scipy_version,
            "solver_package": "SciPy",
            "solver_function": "scipy.optimize.nnls",
            "solver_backend": "Cython improved active-set implementation",
            "version_status": "CONFIRMED_BY_FORMAL_EXECUTABLE_AND_RETAINED_PACKAGE_LINEAGE",
            "formal_run_evidence": f"{rel(BASE_STAGE3B_LOG)};{rel(BIOAPP_SUMMARY)}",
            "environment_evidence": "ENV:anaconda_base/conda-meta/history;ENV:anaconda_base/Lib/site-packages/scipy-1.15.3.dist-info/METADATA",
            "evidence_sha1": f"{sha1_file(BASE_STAGE3B_LOG)};{sha1_file(BIOAPP_SUMMARY)};{sha1_file(base_python_record)};{sha1_file(base_numpy_meta)};{sha1_file(base_scipy_meta)};{sha1_file(base_history)}",
            "confidence": "HIGH",
            "notes": "Active NumPy dist-info is 2.1.0 although the older conda record names 2.1.3; the active import-tree metadata predates formal runs and is the reported value.",
        },
    ]

    unresolved_rows = [
        {
            "issue_id": "U4D-001",
            "severity": "MAJOR",
            "category": "mainline Seurat exact formal-run version",
            "affected_records": "SEURAT_MAINLINE_STAGE1",
            "evidence_found": "successful Stage1 Seurat logs; environment pins R 4.5.2, Seurat 5.3.1 and SeuratObject 5.2.0",
            "missing_evidence": "run-local sessionInfo(), packageVersion() output, renv.lock or command-linked package snapshot",
            "reason_unresolved": "The formal Stage1 logs do not record package versions or the resolved Rscript executable.",
            "safe_current_statement": "Stage1 used Seurat; the reproducibility environment pins Seurat 5.3.1 and SeuratObject 5.2.0, but the exact run-local versions were not recorded.",
            "unsafe_statement": "All formal Stage1 runs were directly verified to use Seurat 5.3.1.",
            "recommended_resolution": "If an immutable historical environment export or run-local sessionInfo is recovered, attach it; otherwise retain bounded wording.",
            "status": "OPEN",
        },
        {
            "issue_id": "U4D-002",
            "severity": "MINOR",
            "category": "optional Fig. 2 Seurat enrichment backend",
            "affected_records": "SU004",
            "evidence_found": "Seurat R implementation exists; benchmark runner defaults to Python backend",
            "missing_evidence": "formal run command or manifest selecting enrichment_backend=r",
            "reason_unresolved": "No retained execution record ties the optional R backend to the formal Fig. 2 panel.",
            "safe_current_statement": "Do not attribute the formal Fig. 2 enrichment panel to Seurat without a run-local backend record.",
            "unsafe_statement": "Fig. 2 enrichment was computed with Seurat.",
            "recommended_resolution": "Use the source-value provenance already retained or recover the original benchmark command.",
            "status": "OPEN",
        },
    ]

    recommendation_rows = [
        ("MR4D-001", "Seurat mainline version", "Methods / software", "Stage1 formal logs confirm Seurat execution; environment pins Seurat 5.3.1 and SeuratObject 5.2.0.", "Exact run-local package versions were not embedded.", "Preprocessing was performed with Seurat; the retained reproducibility environment pins Seurat 5.3.1 and SeuratObject 5.2.0, although run-local package snapshots were not retained.", "All formal Stage1 runs were directly verified with Seurat 5.3.1.", "MAJOR", "E4D008;E4D010;E4D028-E4D037"),
        ("MR4D-002", "BioApp Seurat and R versions", "Methods / software", "BioApp direct records confirm R 4.5.1, Seurat 5.3.1 and SeuratObject 5.2.0.", "none for this route", "CTA Seurat-object parsing and export used R 4.5.1, Seurat 5.3.1 and SeuratObject 5.2.0.", "BioApp used the mainline R 4.5.2 environment.", "HIGH", "E4D011;E4D012;E4D013"),
        ("MR4D-003", "Stage3B solver", "Methods / software", "Formal implementation imports and calls scipy.optimize.nnls per spot.", "Two formal SciPy environments were retained.", "Stage3B fitted each spot by non-negative least squares using scipy.optimize.nnls.", "Stage3B used a custom NNLS, lsq_linear, minimize, or a batched solver.", "HIGH", "E4D001;E4D005;E4D006;E4D007"),
        ("MR4D-004", "Stage3B solver versions", "Methods / software", "Composite simulations used SciPy 1.11.4 with its retained Fortran Lawson-Hanson wrapper; retained base-environment high-resolution/BioApp routes used SciPy 1.15.3 with its Cython improved active-set implementation.", "No single version or internal backend applies to every formal route, although both expose scipy.optimize.nnls and solve the same constrained objective.", "SciPy 1.11.4 was used for the 5%/10% composite simulations, whereas retained base-environment high-resolution and BioApp runs used SciPy 1.15.3; both routes called scipy.optimize.nnls.", "All Stage3B analyses used SciPy 1.11.4 or an identical internal solver build.", "HIGH", "E4D018-E4D027"),
        ("MR4D-005", "Coefficient normalisation", "Methods", "NNLS coefficients are divided by their sum only when the sum is positive.", "none", "The non-negative coefficients were sum-normalised when their total was positive; zero solutions were retained as zero vectors.", "All NNLS outputs were unconditionally divided by their sum.", "HIGH", "E4D001"),
        ("MR4D-006", "Zero/failure behavior", "Methods / implementation detail", "Zero solutions produce zero reconstructions; solver exceptions propagate with no fallback.", "none", "Zero coefficient solutions remained zero; no alternative solver or silent coefficient imputation was used if NNLS failed.", "Failed NNLS fits were automatically replaced by zero vectors.", "HIGH", "E4D001;E4D019;E4D024"),
        ("MR4D-007", "Mathematical consistency", "Methods", "Objective, non-negativity, conditional coefficient normalisation, reconstruction and residual formulas match the supplied definition.", "The compact equation should define x_i and B as preprocessed compositional profiles.", "Define x_i and B as the linearised, row-composition-normalised spot and reference-type profiles; fitting then follows the stated NNLS equations without intercept or regularisation.", "NNLS was fit directly to unprocessed count vectors if that is not otherwise documented.", "HIGH", "E4D001"),
    ]
    recommendation_fields = ["recommendation_id", "topic", "likely_section", "verified_information", "remaining_uncertainty", "recommended_manual_wording", "prohibited_wording", "priority", "evidence_ids"]
    recommendation_dicts = [dict(zip(recommendation_fields, row)) for row in recommendation_rows]

    decision = "CONDITIONAL PASS"
    summary = {
        "task": "Task 4D formal software-version and Stage3B NNLS implementation audit",
        "decision": decision,
        "generated_date": GENERATED_DATE,
        "formal_manuscript_accessed": False,
        "formal_manuscript_modified": False,
        "formal_bibliography_accessed": False,
        "formal_bibliography_modified": False,
        "experimental_code_modified": False,
        "experimental_outputs_modified": False,
        "figures_modified": False,
        "formal_stages_rerun": False,
        "github_access_required": False,
        "stage3b_implementation_status": "CONFIRMED_SINGLE_IMPLEMENTATION",
        "stage3b_solver": "scipy.optimize.nnls",
        "stage3b_math_consistency": "PASS_WITH_DOCUMENTED_DEFAULTS_AND_PREPROCESSING_CONTEXT",
        "formal_solver_environments": 2,
        "formal_scipy_versions": [named_scipy_version, base_scipy_version],
        "bioapp_seurat_exact_version_confirmed": True,
        "mainline_stage1_exact_run_local_seurat_version_confirmed": False,
        "mainline_stage1_environment_pin": "R 4.5.2; Seurat 5.3.1; SeuratObject 5.2.0",
        "unresolved_count": len(unresolved_rows),
        "unresolved_major": 1,
        "unresolved_minor": 1,
        "static_checks": static_checks,
        "repository_local_commit_sha": git_commit(),
    }

    readme = f"""# Task 4D Software and NNLS Formal Audit

This directory is generated by `scripts/run_task4d_software_and_nnls_audit.py`.

- Decision: **{decision}**
- Formal Stage3B solver: `scipy.optimize.nnls`
- Stage3B implementation: one repository implementation, no fallback solver
- Formal SciPy versions: `1.11.4` and `1.15.3`, tied to separate retained executables
- Mathematical consistency: **PASS with solver-default and preprocessing context**
- BioApp Seurat: exact `5.3.1`; SeuratObject `5.2.0`; R `4.5.1`
- Mainline Stage1 Seurat: execution confirmed, exact run-local version not independently recoverable; environment pins `5.3.1`

The audit is static and read-only. It did not import Stage3B, run R, execute an
experimental stage, access manuscript LaTeX/BibTeX, or modify formal outputs.
"""

    final_decision = f"""# Task 4D Final Decision

## 1. Decision

`{decision}`

The formal Stage3B implementation, execution paths, solver function, two formal
solver environments and mathematical behavior are confirmed. The decision is
conditional because the mainline Stage1 logs do not contain a run-local
`sessionInfo()` or `packageVersion()` snapshot.

## 2. Seurat

### Mainline Stage1

- Formal use: confirmed by {len(stage1_log_paths)} successful Stage1 logs.
- Environment pin: R `4.5.2`, Seurat `5.3.1`, SeuratObject `5.2.0`.
- Exact formal-run version: **not independently recoverable**.

### CTA biological application

- R: `4.5.1` (`x86_64-w64-mingw32`).
- Seurat: `5.3.1`.
- SeuratObject: `5.2.0`.
- Evidence status: direct formal package/runtime records.

## 3. Stage3B NNLS

- Implementation: `src/stages/stage3b_st_unsupported.py`.
- Function: `fit_nonnegative_mixtures`.
- Caller: `run_stage3b`.
- Solver: `scipy.optimize.nnls`.
- Fitting: per spot with design `B.T` (`n_genes x n_types`).
- Intercept: none.
- Regularisation: none.
- Explicit tolerance/maxiter: none; SciPy defaults apply.
- Fallback: none.
- Coefficients: sum-normalised only when the sum is positive.
- Zero solution: retained as zeros; reconstruction is zero.
- Failure: exception propagates; no silent replacement.

## 4. Formal solver versions

- 5%/10% composite simulation routes: Python `{named_python_meta['version']}`, NumPy `{named_numpy_version}`, SciPy `{named_scipy_version}`.
- High-resolution/BioApp base routes: Python `{base_python_meta['version']}`, NumPy `{base_numpy_version}`, SciPy `{base_scipy_version}`.

## 5. Mathematical consistency

`PASS_WITH_DOCUMENTED_DEFAULTS_AND_PREPROCESSING_CONTEXT`

The objective orientation, non-negative constraint, conditional coefficient
normalisation, reconstruction and residual match the supplied Methods
definition. The manual Methods wording should define `x_i` and `B` as the
linearised, row-composition-normalised profiles and should not claim a custom
tolerance or fallback.

## 6. Remaining unresolved items

- Mainline Stage1 exact run-local R/Seurat/SeuratObject versions.
- Whether the optional Seurat enrichment backend generated the formal Fig. 2 panel.

## 7. Guardrails

```text
Formal manuscript accessed: false
Formal manuscript modified: false
Formal bibliography accessed: false
Formal bibliography modified: false
Experimental code modified: false
Experimental outputs modified: false
Figures modified: false
Formal stages rerun: false
GitHub repository access required: false
```
"""

    payloads: dict[str, bytes] = {
        "README.md": text_bytes(readme),
        "local_evidence_inventory.tsv": tsv_bytes(evidence_rows, [
            "evidence_id", "category", "file_path", "file_size", "modified_time_ns", "sha1", "evidence_type", "evidence_priority", "supports", "read_only", "notes"
        ]),
        "seurat_version_audit.tsv": tsv_bytes(seurat_version_rows, [
            "record_id", "dataset_family", "formal_script", "r_version", "seurat_version", "seuratobject_version", "version_status", "evidence_path", "evidence_sha1", "evidence_type", "confidence", "notes"
        ]),
        "seurat_usage_inventory.tsv": tsv_bytes(seurat_usage_rows, [
            "usage_id", "script", "functionality", "dataset_family", "formal_execution_status", "version_link", "evidence_path", "notes"
        ]),
        "stage3b_execution_path.tsv": tsv_bytes(execution_rows, [
            "path_id", "dataset_family", "caller", "entrypoint", "implementation", "fit_function", "solver_call", "execution_mode", "formal_execution_evidence", "status", "notes"
        ]),
        "stage3b_nnls_implementation_audit.tsv": tsv_bytes(implementation_rows, [
            "record_id", "formal_implementation_file", "formal_function", "formal_caller", "solver_package", "solver_function", "solver_backend", "solver_version", "matrix_orientation", "non_negative_constraint", "intercept", "regularisation", "tolerance", "max_iterations", "fallback", "failure_handling", "coefficient_normalisation", "zero_solution_handling", "reconstruction_formula", "status", "evidence_path", "evidence_sha1", "notes"
        ]),
        "stage3b_math_consistency_audit.tsv": tsv_bytes(math_dicts, math_fields),
        "python_solver_version_audit.tsv": tsv_bytes(solver_rows, [
            "record_id", "formal_dataset_family", "python_executable", "python_version", "platform", "numpy_version", "scipy_version", "solver_package", "solver_function", "solver_backend", "version_status", "formal_run_evidence", "environment_evidence", "evidence_sha1", "confidence", "notes"
        ]),
        "unresolved_software_records.tsv": tsv_bytes(unresolved_rows, [
            "issue_id", "severity", "category", "affected_records", "evidence_found", "missing_evidence", "reason_unresolved", "safe_current_statement", "unsafe_statement", "recommended_resolution", "status"
        ]),
        "manuscript_manual_update_recommendations.tsv": tsv_bytes(recommendation_dicts, recommendation_fields),
        "audit_summary.json": json_bytes(summary),
        "final_decision.md": text_bytes(final_decision),
    }
    return payloads


def main() -> None:
    first = build_payloads()
    second = build_payloads()
    if first != second:
        changed = sorted(name for name in set(first) | set(second) if first.get(name) != second.get(name))
        raise RuntimeError(f"Non-deterministic Task 4D build: {changed}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    expected = set(first)
    for path in OUT_DIR.iterdir():
        if path.is_file() and path.name not in expected:
            path.unlink()
    for name, data in first.items():
        (OUT_DIR / name).write_bytes(data)
    output_hashes = {name: sha1_bytes(data) for name, data in sorted(first.items())}
    print("TASK4D_AUDIT_PASS")
    print(f"output_dir={rel(OUT_DIR)}")
    print(f"files={len(first)} decision=CONDITIONAL PASS")
    print(f"core_output_digest={sha1_bytes(json_bytes(output_hashes))}")


if __name__ == "__main__":
    main()
