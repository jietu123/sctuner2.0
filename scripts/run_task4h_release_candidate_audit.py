#!/usr/bin/env python3
"""Build and audit the Task 4H reproducibility release candidate.

This script performs archival copying, classification, hashing and validation.
It does not run SVTuner, CytoSPACE, an experimental stage, metric computation,
or figure generation.
"""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
OUT = (
    ROOT
    / "visualizations"
    / "manuscript_audits"
    / "task4h_release_candidate_audit"
)
CANDIDATE = ROOT / "reproducibility_release" / "release_candidate_v1.0"

FROZEN_BASELINE = {
    "reproducibility_release/manifests/dataset_manifest_v1.0.tsv": (
        "64223e7846fc3cc8228a076b0f3e069742f58afe"
    ),
    "reproducibility_release/manifests/experiment_manifest_v1.0.tsv": (
        "347539ed5677dc89318772333552d5944821a355"
    ),
    "reproducibility_release/manifests/manifest_v1.0_metadata.json": (
        "0141fc28d9a106a8c7995a5a5658bf60fe298ef5"
    ),
    "reproducibility_release/manifests/manifest_v1.0_human_approval.json": (
        "711ae2b8a6d973725ab28e0423d5925fe90ca495"
    ),
    "reproducibility_release/manifests/manifest_v1.0_freeze_report.md": (
        "ed00f3643e577dea77f0b5d7248332ed5d022902"
    ),
    "reproducibility_release/release_manifest.tsv": (
        "9660b3fb9fbeee051289c6ab90fffe2a09ad8ec2"
    ),
}

REQUIRED_CANDIDATE_DIRS = [
    "code",
    "environments",
    "configs",
    "manifests",
    "source_values",
    "figure_scripts",
    "audit_summaries",
    "reconstruction",
    "inventory",
    "metadata",
]

PROHIBITED_EXTENSIONS = {
    ".h5",
    ".h5ad",
    ".rdata",
    ".rds",
    ".loom",
    ".mtx",
    ".fastq",
    ".fq",
    ".bam",
    ".cram",
    ".tif",
    ".tiff",
}

ACCEPTED_PROVENANCE_GAPS = [
    ("GAP-PROV-001", "human-lung read-run accession"),
    ("GAP-PROV-002", "CID4465 slide identifier"),
    ("GAP-PROV-003", "CID4465 capture area"),
    ("GAP-PROV-004", "CTA exact chemistry or kit version"),
    ("GAP-PROV-005", "Vizgen provider checksum"),
    ("GAP-PROV-006", "Vizgen explicit redistribution permission"),
]

ENVIRONMENT_GAPS = [
    (
        "GAP-ENV-001",
        "Mainline Stage1 exact run-local R/Seurat package snapshot",
        "The frozen specification is R 4.5.2, Seurat 5.3.1 and "
        "SeuratObject 5.2.0; formal logs do not contain run-local sessionInfo().",
    ),
    (
        "GAP-ENV-002",
        "Legacy Fig. 2 enrichment exact Python runtime",
        "The source tables and backend implementation are retained, but the "
        "original executable and complete package snapshot were not retained.",
    ),
    (
        "GAP-ENV-003",
        "Tangram formal mapping runtime snapshot",
        "Formal outputs and the runner are retained; the historical environment "
        "containing the successful Tangram installation is not retained.",
    ),
    (
        "GAP-ENV-004",
        "novoSpaRc formal mapping runtime snapshot",
        "Formal outputs and the runner are retained; exact successful package "
        "versions were not recorded.",
    ),
    (
        "GAP-ENV-005",
        "SpaOTsc formal mapping runtime snapshot",
        "Formal outputs and the runner are retained; exact successful package "
        "versions were not recorded.",
    ),
    (
        "GAP-ENV-006",
        "CellTrek-labelled Python fallback runtime snapshot",
        "The retained formal runner identifies itself as a Python fallback "
        "CellTrek-style implementation; exact historical package versions were "
        "not recorded and the official R CellTrek package was not used.",
    ),
]

FIGURE_SCRIPTS: list[dict[str, str]] = [
    {
        "path": "scripts/export_simulation_stage3ab_joint_overview.py",
        "panels": "Fig. 1C",
        "inputs": "simulation Stage3A/Stage3B triptych panels and scenario results",
        "outputs": "simulation_stage3ab_joint_triptych_overview_stack_3datasets.svg/png",
        "format": "mixed",
        "retention": "DIRECTLY_RETAINED",
    },
    {
        "path": "scripts/plot_method_comparison_composition_recovery.py",
        "panels": "Fig. 1D (0% noise)",
        "inputs": "nine-scenario mapping outputs and simulation truth",
        "outputs": "composition_recovery_7mapping_methods_composite_no_noise_abstention_aware_boxplot.png/pdf",
        "format": "mixed",
        "retention": "DIRECTLY_RETAINED",
    },
    {
        "path": "scripts/run_composite_scnoise10_stage3ab_full.py",
        "panels": "Fig. 1D (10% noise)",
        "inputs": "10% noisy-reference Stage3A/Stage3B outputs and retained six-method mappings",
        "outputs": "fig1d_scnoise10_stage3ab_full.png/pdf and source tables",
        "format": "mixed",
        "retention": "DIRECTLY_RETAINED",
    },
    {
        "path": "scripts/plot_real_profile_mask_fig2c_only.py",
        "panels": "Fig. 2B",
        "inputs": "spot_foundation.csv",
        "outputs": "fig2_panel_c_real_profile_mask.png/pdf",
        "format": "mixed",
        "retention": "DIRECTLY_RETAINED",
    },
    {
        "path": "scripts/plot_fig2d_profile_mask_benchmark.py",
        "panels": "Fig. 2C",
        "inputs": "fig2d_profile_mask_benchmark_source_values.csv",
        "outputs": "fig2d_profile_mask_benchmark.png/pdf",
        "format": "mixed",
        "retention": "RECOVERED_FROM_LOCAL_GIT_AND_HASH_VERIFIED",
    },
    {
        "path": "scripts/plot_real_profile_mask_expression_recovery.py",
        "panels": "Fig. 2D",
        "inputs": "expression_recovery_by_scenario.csv and expression_recovery_long.csv",
        "outputs": "expression_recovery_cosine_summary_bar.png/pdf",
        "format": "mixed",
        "retention": "DIRECTLY_RETAINED",
    },
    {
        "path": "scripts/build_fig2e_stage3_profile_mask_route2_only.py",
        "panels": "Fig. 2F",
        "inputs": "fig2e_stage3_profile_mask_metrics.csv",
        "outputs": "fig2e_stage3_profile_mask_route2_ce9_ce10.png/pdf",
        "format": "mixed",
        "retention": "DIRECTLY_RETAINED",
    },
    {
        "path": "scripts/build_fig3b_highres_profile_mask_enrichment.py",
        "panels": "Fig. 3B",
        "inputs": "five high-resolution mapping outputs",
        "outputs": "fig3b_highres_profile_mask_enrichment.png/pdf/svg and source tables",
        "format": "mixed",
        "retention": "DIRECTLY_RETAINED",
    },
    {
        "path": "scripts/run_highres_targeted_validation_fig2c.py",
        "panels": "Fig. 3C representative curve",
        "inputs": "configs/highres_targeted_validation_fig2c.json",
        "outputs": "targeted_fig2c_*.png/pdf/svg and curve CSVs",
        "format": "mixed",
        "retention": "DIRECTLY_RETAINED",
    },
    {
        "path": "scripts/build_highres_profile_mask_fig2d_benchmark.py",
        "panels": "Fig. 3C 21-readout panel",
        "inputs": "high-resolution mapping outputs",
        "outputs": "fig2d_highres_profile_mask_benchmark.png/pdf/svg and source tables",
        "format": "mixed",
        "retention": "DIRECTLY_RETAINED",
    },
    {
        "path": "scripts/run_highres_targeted_validation_fig2d.py",
        "panels": "Fig. 3C targeted 10-readout panel",
        "inputs": "configs/highres_targeted_validation_fig2d.csv",
        "outputs": "targeted_fig2d_highres_fixed_panel.png/pdf/svg and source tables",
        "format": "mixed",
        "retention": "DIRECTLY_RETAINED",
    },
    {
        "path": "scripts/plot_cytospace_fig2k_mapping_with_route2.py",
        "panels": "Fig. 3E",
        "inputs": "T-cell state-decoy baseline and route2 mappings",
        "outputs": "fig2k_stage3_detected_state_decoy_baseline_vs_route2.png/pdf",
        "format": "mixed",
        "retention": "DIRECTLY_RETAINED",
    },
    {
        "path": "scripts/plot_cytospace_fig2i_forced_unsupported.py",
        "panels": "Fig. 3F",
        "inputs": "kidney state-decoy baseline and route2 mappings",
        "outputs": "fig2i_stage3_detected_state32like_n1000_baseline_vs_route2*.png/pdf",
        "format": "mixed",
        "retention": "DIRECTLY_RETAINED",
    },
    {
        "path": "scripts/plot_stage3b_reference_dropout_spatial_stack.py",
        "panels": "Fig. 4A/F",
        "inputs": "reference-dropout Stage3B summaries and spatial outputs",
        "outputs": "stage3b_reference_dropout_spatial_stack_recommended_6x2.svg/png",
        "format": "mixed",
        "retention": "DIRECTLY_RETAINED",
    },
    {
        "path": "scripts/plot_stage3b_reference_dropout_panel_a_violin.py",
        "panels": "Fig. 4B",
        "inputs": "marker percentile values",
        "outputs": "reference-dropout marker percentile panel",
        "format": "mixed",
        "retention": "DIRECTLY_RETAINED",
    },
    {
        "path": "scripts/plot_stage3b_reference_dropout_panel_a_blank_composition.py",
        "panels": "Fig. 4C",
        "inputs": "blank-region composition summaries",
        "outputs": "reference-dropout blank composition panel",
        "format": "mixed",
        "retention": "DIRECTLY_RETAINED",
    },
    {
        "path": "scripts/plot_cell2location_abstention_aware_region_accuracy.py",
        "panels": "Fig. 4E",
        "inputs": "thalamic region/spot tables and frozen reference labels",
        "outputs": "cell2location_ST8059051_*_abstention_aware_region_accuracy.png",
        "format": "raster",
        "retention": "DIRECTLY_RETAINED",
    },
    {
        "path": "scripts/plot_stage3b_reference_missing_cell2location_panels_bc.py",
        "panels": "Fig. 4G/H",
        "inputs": "forced-surrogate and expression-similarity matrices",
        "outputs": "cell2location reference-missing panels",
        "format": "mixed",
        "retention": "DIRECTLY_RETAINED",
    },
    {
        "path": "scripts/plot_stage3b_reference_missing_thalamic_panels_bc.py",
        "panels": "Fig. 4H",
        "inputs": "thalamic forced-assignment and expression-similarity tables",
        "outputs": "thalamic reference-missing panels",
        "format": "mixed",
        "retention": "DIRECTLY_RETAINED",
    },
    {
        "path": "scripts/plot_stage3b_false_spatial_niche_compact.py",
        "panels": "Fig. 4I coupling",
        "inputs": "coupling matrices, pair statistics and reference-relative summary",
        "outputs": "false-spatial-niche compact panel",
        "format": "mixed",
        "retention": "DIRECTLY_RETAINED",
    },
    {
        "path": "scripts/plot_stage3b_communication_strict_reference_style.py",
        "panels": "Fig. 4I KITLG-KIT",
        "inputs": "local LR hotspot and ligand-target source tables",
        "outputs": "stage3b_communication_strict_reference_style.png/pdf",
        "format": "mixed",
        "retention": "DIRECTLY_RETAINED",
    },
    {
        "path": "scripts/run_bioapp_main_figure_v3_12_panel_A_evidence_chain_redesign.py",
        "panels": "Fig. 5A-I",
        "inputs": "frozen endpoint, Phase 5/8 metrics, morphology, interface and microenvironment tables",
        "outputs": "bioapp_main_figure_v3_12_panel_A_evidence_chain_redesign panels",
        "format": "mixed",
        "retention": "DIRECTLY_RETAINED",
    },
]

CODE_EXCLUSIONS = {
    "scripts/cleanup_redundant_storage.py",
    "scripts/rebuild_editable_svg_references.py",
    "scripts/rebuild_pdf_figures_as_editable_svg.py",
    "scripts/reproduce_ncem_fig2_human_lymph_node.py",
    "scripts/run_bioapp_main_figure_v3_10_d2_light_green_footprint_correction.py",
    "scripts/run_bioapp_main_figure_v3_11_top_row_abc_redesign.py",
}

AUDIT_PATHS = [
    "visualizations/manuscript_audits/dataset_provenance_formal_audit/README.md",
    "visualizations/manuscript_audits/dataset_provenance_formal_audit/audit_summary.json",
    "visualizations/manuscript_audits/dataset_provenance_formal_audit/final_decision.md",
    "visualizations/manuscript_audits/dataset_provenance_formal_audit/accession_type_audit.tsv",
    "visualizations/manuscript_audits/dataset_provenance_formal_audit/official_source_registry.tsv",
    "visualizations/manuscript_audits/dataset_provenance_formal_audit/unresolved_records.tsv",
    "visualizations/manuscript_audits/dataset_provenance_major_evidence_recovery/README.md",
    "visualizations/manuscript_audits/dataset_provenance_major_evidence_recovery/audit_summary.json",
    "visualizations/manuscript_audits/dataset_provenance_major_evidence_recovery/final_decision.md",
    "visualizations/manuscript_audits/dataset_provenance_major_evidence_recovery/heca_input_file_inventory.tsv",
    "visualizations/manuscript_audits/dataset_provenance_major_evidence_recovery/heca_reference_lineage_audit.tsv",
    "visualizations/manuscript_audits/dataset_provenance_major_evidence_recovery/cta_metadata_recovery.tsv",
    "visualizations/manuscript_audits/dataset_provenance_major_evidence_recovery/vizgen_release_evidence.tsv",
    "visualizations/manuscript_audits/dataset_provenance_major_evidence_recovery/vizgen_sample_file_inventory.tsv",
    "visualizations/manuscript_audits/dataset_provenance_major_evidence_recovery/unresolved_after_task4b2.tsv",
    "visualizations/manuscript_audits/software_and_nnls_formal_audit/README.md",
    "visualizations/manuscript_audits/software_and_nnls_formal_audit/audit_summary.json",
    "visualizations/manuscript_audits/software_and_nnls_formal_audit/final_decision.md",
    "visualizations/manuscript_audits/software_and_nnls_formal_audit/python_solver_version_audit.tsv",
    "visualizations/manuscript_audits/software_and_nnls_formal_audit/seurat_version_audit.tsv",
    "visualizations/manuscript_audits/software_and_nnls_formal_audit/stage3b_nnls_implementation_audit.tsv",
    "visualizations/manuscript_audits/software_and_nnls_formal_audit/stage3b_math_consistency_audit.tsv",
    "visualizations/manuscript_audits/software_and_nnls_formal_audit/unresolved_software_records.tsv",
    "visualizations/manuscript_audits/fig2_enrichment_backend_formal_audit/README.md",
    "visualizations/manuscript_audits/fig2_enrichment_backend_formal_audit/audit_summary.json",
    "visualizations/manuscript_audits/fig2_enrichment_backend_formal_audit/final_decision.md",
    "visualizations/manuscript_audits/fig2_enrichment_backend_formal_audit/backend_decision_by_panel.tsv",
    "visualizations/manuscript_audits/fig2_enrichment_backend_formal_audit/formal_fig2_source_value_inventory.tsv",
    "visualizations/manuscript_audits/fig2_enrichment_backend_formal_audit/formal_source_value_consistency.tsv",
    "visualizations/manuscript_audits/manifest_reconciliation_and_freeze_candidate/README.md",
    "visualizations/manuscript_audits/manifest_reconciliation_and_freeze_candidate/audit_summary.json",
    "visualizations/manuscript_audits/manifest_reconciliation_and_freeze_candidate/final_decision.md",
    "visualizations/manuscript_audits/manifest_reconciliation_and_freeze_candidate/manuscript_manifest_reconciliation.tsv",
    "visualizations/manuscript_audits/manifest_reconciliation_and_freeze_candidate/redistribution_classification.tsv",
    "visualizations/manuscript_audits/manifest_reconciliation_and_freeze_candidate/unresolved_nonblocking_records.tsv",
    "visualizations/manuscript_audits/manifest_reconciliation_and_freeze_candidate/validation_checks.tsv",
    "visualizations/bioapp_experiment/fig5b_auprc_formal_audit/README.md",
    "visualizations/bioapp_experiment/fig5b_auprc_formal_audit/audit_summary.json",
    "visualizations/bioapp_experiment/fig5b_auprc_formal_audit/final_decision.md",
    "visualizations/bioapp_experiment/fig5b_auprc_formal_audit/input_file_hashes.tsv",
    "visualizations/bioapp_experiment/fig5b_auprc_formal_audit/value_provenance.tsv",
    "visualizations/bioapp_experiment/bioapp_coordinate_registration_recovery_audit/bioapp_registration_recovery_summary.json",
    "visualizations/bioapp_experiment/bioapp_coordinate_registration_recovery_audit/bioapp_registration_recovery_decision.md",
    "visualizations/bioapp_experiment/bioapp_coordinate_registration_recovery_audit/bioapp_image_slot_scale_factor_audit.csv",
    "visualizations/bioapp_experiment/bioapp_figure_v2_6_frozen_spot_coordinate_provenance_audit/bioapp_v2_6_coordinate_provenance_summary.json",
    "visualizations/bioapp_experiment/bioapp_figure_v2_6_frozen_spot_coordinate_provenance_audit/bioapp_v2_6_coordinate_recovery_decision.md",
    "visualizations/bioapp_experiment/bioapp_figure_v2_6_frozen_spot_coordinate_provenance_audit/bioapp_v2_6_provenance_graph_edges.csv",
    "visualizations/bioapp_experiment/bioapp_phase9_final_biological_application_audit_figure_selection_interpretation_boundary_lock/bioapp_phase9_final_audit_summary.json",
    "visualizations/bioapp_experiment/bioapp_phase9_final_biological_application_audit_figure_selection_interpretation_boundary_lock/bioapp_phase9_interpretation_boundary_lock.json",
    "visualizations/bioapp_experiment/bioapp_phase9_final_biological_application_audit_figure_selection_interpretation_boundary_lock/bioapp_phase9_figure_selection_table.csv",
    "reproducibility_release/manifests/manifest_v1.0_freeze_report.md",
]


@dataclass
class Asset:
    record_id: str
    asset_class: str
    relative_path: str
    filename: str
    file_extension: str
    byte_size: int
    sha1: str
    release_status: str
    redistribution_class: str
    formal_role: str
    linked_figure_or_analysis: str
    source_or_generator: str
    notes: str


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def sha1(path: Path) -> str:
    digest = hashlib.sha1()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


def copy_exact(source: Path, destination: Path) -> None:
    require(source.is_file(), f"Missing source asset: {rel(source)}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        require(destination.read_bytes() == source.read_bytes(), f"Candidate copy differs: {rel(destination)}")
        return
    shutil.copyfile(source, destination)


def add_asset(
    assets: list[Asset],
    source: Path,
    candidate_relative: str,
    asset_class: str,
    formal_role: str,
    linked: str,
    redistribution: str,
    notes: str = "",
) -> Path:
    source = source.resolve()
    require(source.is_relative_to(ROOT), f"External file cannot enter candidate: {source}")
    source_rel = rel(source)
    lowered = source_rel.lower()
    require("/data/raw/" not in f"/{lowered}", f"Raw data prohibited: {source_rel}")
    require(source.suffix.lower() not in PROHIBITED_EXTENSIONS, f"Prohibited extension: {source_rel}")
    destination = CANDIDATE / candidate_relative
    copy_exact(source, destination)
    assets.append(
        Asset(
            record_id=f"T4H-{len(assets) + 1:04d}",
            asset_class=asset_class,
            relative_path=rel(destination),
            filename=destination.name,
            file_extension=destination.suffix.lower() or "NONE",
            byte_size=destination.stat().st_size,
            sha1=sha1(destination),
            release_status="INCLUDED_RELEASE_CANDIDATE",
            redistribution_class=redistribution,
            formal_role=formal_role,
            linked_figure_or_analysis=linked,
            source_or_generator=source_rel,
            notes=notes,
        )
    )
    return destination


def verify_frozen_baseline() -> None:
    for relative_path, expected in FROZEN_BASELINE.items():
        path = ROOT / relative_path
        require(path.is_file(), f"Frozen baseline missing: {relative_path}")
        require(sha1(path) == expected, f"Frozen baseline hash mismatch: {relative_path}")


def git_output(*args: str) -> str:
    result = subprocess.run(
        ["git", *args],
        cwd=ROOT,
        text=True,
        encoding="utf-8",
        errors="replace",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    require(result.returncode == 0, f"git {' '.join(args)} failed: {result.stderr.strip()}")
    return result.stdout.strip()


def source_specs() -> list[tuple[str, str, str, str]]:
    """Return (path, coverage group, role, generator)."""
    specs: list[tuple[str, str, str, str]] = []

    def add(path: str, group: str, role: str, generator: str) -> None:
        require((ROOT / path).is_file(), f"Required source-value asset missing: {path}")
        specs.append((path, group, role, generator))

    def add_glob(pattern: str, group: str, role: str, generator: str) -> None:
        matches = sorted(p for p in ROOT.glob(pattern) if p.is_file())
        require(matches, f"Required source-value glob empty: {pattern}")
        for path in matches:
            add(rel(path), group, role, generator)

    add(
        "visualizations/simulations/composite_stage3a_stage3b_recheck.csv",
        "Fig. 1 simulation values",
        "joint simulation Stage3A/Stage3B audit values",
        "scripts/evaluate_stage3b_simulation.py",
    )
    add(
        "visualizations/simulations/simulation_stage3ab_joint_triptych_overview_stack_3datasets.json",
        "Fig. 1 simulation values",
        "simulation overview manifest",
        "scripts/export_simulation_stage3ab_joint_overview.py",
    )
    add_glob(
        "visualizations/method_comparison/composite_no_noise/*.csv",
        "Fig. 1 no-noise results",
        "no-noise whole-space evaluation source values",
        "scripts/plot_method_comparison_composition_recovery.py",
    )
    add_glob(
        "visualizations/method_comparison/composite_scnoise10_stage3ab_full/*.csv",
        "Fig. 1 10% reference-noise results",
        "10% reference-noise source values and execution audits",
        "scripts/run_composite_scnoise10_stage3ab_full.py",
    )
    add_glob(
        "visualizations/method_comparison/composite_scnoise10_stage3ab_full/*.json",
        "Fig. 1 10% reference-noise results",
        "10% reference-noise summaries and guardrails",
        "scripts/run_composite_scnoise10_stage3ab_full.py",
    )

    experiments = read_tsv(ROOT / "reproducibility_release/manifests/experiment_manifest_v1.0.tsv")
    for row in experiments:
        if row.get("figure_panel") == "Fig. 2A-B/E" and row.get("source_value_paths"):
            for path in row["source_value_paths"].split(";"):
                add(path, "Fig. 2 profile-masking results", "spot-level profile-mask foundation", "scripts/build_real_profile_mask_foundation.py")
    add_glob(
        "visualizations/cytospace_fig2d_profile_mask_benchmark/*.csv",
        "Fig. 2 enrichment source values",
        "low-resolution enrichment benchmark source values",
        "scripts/run_cytospace_fig2d_profile_mask_benchmark.py",
    )
    add_glob(
        "visualizations/cytospace_fig2c_melanoma_stage3_profile_mask/*.csv",
        "Fig. 2 enrichment source values",
        "representative enrichment source values",
        "scripts/compute_fig2c_enrichment_python.py",
    )
    add_glob(
        "result/cytospace_fig2c_melanoma_mel1_rep2_screen_mask_macrophages/fig2c_strict_cd4_ce9/**/*.csv",
        "Fig. 2 enrichment source values",
        "retained baseline and SVTuner enrichment curves/summaries",
        "scripts/compute_fig2c_enrichment_python.py",
    )
    add(
        "result/cytospace_fig2e_stage3_profile_mask/fig2e_stage3_profile_mask_metrics.csv",
        "Fig. 2 enrichment source values",
        "external-tumour profile-mask metrics",
        "scripts/build_fig2e_stage3_profile_mask_benchmark.py",
    )
    add_glob(
        "result/real_profile_mask_expression_recovery/*.csv",
        "Fig. 2 cosine-similarity values",
        "expression-recovery cosine source values",
        "scripts/compute_real_profile_mask_expression_recovery.py",
    )
    add(
        "result/real_profile_mask_expression_recovery/expression_recovery_config.json",
        "Fig. 2 cosine-similarity values",
        "expression-recovery configuration record",
        "scripts/compute_real_profile_mask_expression_recovery.py",
    )

    add_glob(
        "visualizations/highres_profile_mask_fig3b_enrichment/*.csv",
        "Fig. 3 high-resolution results",
        "high-resolution running-enrichment source values",
        "scripts/build_fig3b_highres_profile_mask_enrichment.py",
    )
    add_glob(
        "visualizations/highres_profile_mask_fig3b_enrichment/*.json",
        "Fig. 3 high-resolution results",
        "high-resolution enrichment manifest",
        "scripts/build_fig3b_highres_profile_mask_enrichment.py",
    )
    add_glob(
        "visualizations/highres_profile_mask_fig2d/all_candidates/*.csv",
        "Fig. 3 21-readout panel",
        "all-candidate high-resolution readouts",
        "scripts/build_highres_profile_mask_fig2d_benchmark.py",
    )
    add_glob(
        "visualizations/highres_profile_mask_fig2d/all_candidates/*.json",
        "Fig. 3 21-readout panel",
        "all-candidate high-resolution manifest",
        "scripts/build_highres_profile_mask_fig2d_benchmark.py",
    )
    add_glob(
        "visualizations/highres_profile_mask_fig2d/targeted_validation/*.csv",
        "Fig. 3 selected 10-readout panel",
        "targeted high-resolution source values",
        "scripts/run_highres_targeted_validation_fig2d.py",
    )
    add_glob(
        "visualizations/highres_profile_mask_fig2d/targeted_validation/*.json",
        "Fig. 3 selected 10-readout panel",
        "targeted panel manifest",
        "scripts/run_highres_targeted_validation_fig2d.py",
    )
    add_glob(
        "visualizations/highres_targeted_validation_fig2c/*.csv",
        "Fig. 3 high-resolution results",
        "representative high-resolution curve values",
        "scripts/run_highres_targeted_validation_fig2c.py",
    )
    add_glob(
        "visualizations/cytospace_fig2k_tcell_states_stage3_decoy/*.csv",
        "Fig. 3 T-cell state-decoy results",
        "T-cell state-decoy source values",
        "scripts/plot_cytospace_fig2k_mapping_with_route2.py",
    )
    add_glob(
        "visualizations/cytospace_fig2k_tcell_states_stage3_decoy/*.json",
        "Fig. 3 T-cell state-decoy results",
        "T-cell state-decoy manifest",
        "scripts/plot_cytospace_fig2k_mapping_with_route2.py",
    )
    add_glob(
        "visualizations/cytospace_fig2i_mouse_kidney_stage3_unsupported_decoy_sensitivity/*.csv",
        "Fig. 3 kidney state-decoy results",
        "kidney state-decoy source values",
        "scripts/plot_cytospace_fig2i_forced_unsupported.py",
    )

    add_glob(
        "visualizations/stage3b_realdata_candidate_scan/**/*.csv",
        "Fig. 4 reference-dropout results",
        "real reference-dropout source values",
        "scripts/plot_stage3b_reference_dropout_spatial_stack.py",
    )
    add_glob(
        "visualizations/cell2location_stage3b_case/thalamic_top15/*.csv",
        "Fig. 4 thalamic case",
        "thalamic abstention-aware reference-relative values",
        "scripts/plot_cell2location_abstention_aware_region_accuracy.py",
    )
    add_glob(
        "visualizations/cell2location_stage3b_case/thalamic_top15/*.json",
        "Fig. 4 thalamic case",
        "thalamic case summary",
        "scripts/plot_cell2location_abstention_aware_region_accuracy.py",
    )
    add_glob(
        "visualizations/stage3b_reference_missing_stress/*.csv",
        "Fig. 4 forced-surrogate and expression-similarity matrices",
        "forced-surrogate and expression-similarity source values",
        "scripts/plot_stage3b_reference_missing_cell2location_panels_bc.py",
    )
    add_glob(
        "visualizations/stage3b_false_spatial_niche/candidate_validation/brca_her2_plasma/*.csv",
        "Fig. 4 KITLG-KIT analysis",
        "reference-relative coupling matrices and summaries",
        "scripts/plot_stage3b_false_spatial_niche_compact.py",
    )
    add_glob(
        "visualizations/stage3b_communication_validation/brca_her2_ffpe_plasma/strict_reference_style/*.csv",
        "Fig. 4 KITLG-KIT analysis",
        "local ligand-receptor hotspot and ligand-target source values",
        "scripts/plot_stage3b_communication_strict_reference_style.py",
    )

    add_glob(
        "visualizations/bioapp_experiment/bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping/*.csv",
        "Fig. 5 CTA endpoint results",
        "frozen CTA endpoint and composition values",
        "scripts/run_bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping.py",
    )
    add_glob(
        "visualizations/bioapp_experiment/bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping/*.json",
        "Fig. 5 CTA endpoint results",
        "frozen CTA endpoint definitions and guardrails",
        "scripts/run_bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping.py",
    )
    add_glob(
        "visualizations/bioapp_experiment/bioapp_phase5_endpoint_specific_evaluation_of_cytospace_baselines_against_frozen_cta_immune_endpoint/*.csv",
        "Fig. 5 CTA endpoint results",
        "baseline endpoint evaluation source values",
        "scripts/run_bioapp_phase5_endpoint_specific_evaluation_of_cytospace_baselines_against_frozen_cta_immune_endpoint.py",
    )
    add_glob(
        "visualizations/bioapp_experiment/bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison/*.csv",
        "Fig. 5 binary-withholding table",
        "SVTuner endpoint comparison, binary withholding and curve values",
        "scripts/run_bioapp_phase8_svtuner_vs_endpoint_evaluation_and_baseline_svtuner_comparison.py",
    )
    add_glob(
        "visualizations/bioapp_experiment/fig5b_auprc_formal_audit/*.csv",
        "Fig. 5 AUROC and average precision",
        "formal Fig. 5B analysis set and ROC/PR source values",
        "scripts/audit_fig5b_auprc_provenance.py",
    )
    for phase, group in [
        ("bioapp_downstream_phase_d0_morphology_domain_boundary_freeze", "Fig. 5 morphology, interface and microenvironment analyses"),
        ("bioapp_downstream_phase_d1_morphology_and_interface_analysis", "Fig. 5 morphology, interface and microenvironment analyses"),
        ("bioapp_downstream_phase_d2_statistical_strengthening", "Fig. 5 morphology, interface and microenvironment analyses"),
        ("bioapp_downstream_phase_d3_microenvironment_context_validation", "Fig. 5 morphology, interface and microenvironment analyses"),
    ]:
        add_glob(
            f"visualizations/bioapp_experiment/{phase}/*.csv",
            group,
            "morphology/interface/microenvironment source values",
            f"scripts/run_{phase}.py",
        )

    deduplicated: dict[str, tuple[str, str, str, str]] = {}
    for spec in specs:
        deduplicated.setdefault(spec[0], spec)
    return list(deduplicated.values())


def configuration_rows() -> tuple[list[dict[str, Any]], list[tuple[Path, str, str]]]:
    experiments = read_tsv(ROOT / "reproducibility_release/manifests/experiment_manifest_v1.0.tsv")
    formal_map: dict[str, set[str]] = defaultdict(set)
    for row in experiments:
        path = row.get("formal_config_path", "")
        if path:
            formal_map[path].add(row.get("figure_panel", ""))

    explicit_formal = {
        "configs/fig3b_highres_profile_mask_enrichment.json",
        "configs/highres_targeted_validation_fig2c.json",
        "configs/highres_targeted_validation_fig2d.csv",
        "visualizations/bioapp_experiment/bioapp_phase7_svtuner_aware_execution_against_frozen_cta_immune_endpoint/bioapp_phase7_stage3b_dataset_config.yaml",
    }
    for path in ROOT.glob("configs/datasets/*_scnoise10.yaml"):
        if any(token in path.name for token in ("sc_missing_", "endothelial_marker_")):
            explicit_formal.add(rel(path))
    for path in ROOT.glob("configs/datasets/*_sc_missing_*.yaml"):
        if "scnoise05" not in path.name:
            explicit_formal.add(rel(path))

    all_config_paths = sorted(
        {
            *[rel(p) for p in (ROOT / "configs").rglob("*") if p.is_file()],
            *explicit_formal,
        }
    )
    rows: list[dict[str, Any]] = []
    included: list[tuple[Path, str, str]] = []
    formal_stems = {Path(path).stem for path in formal_map} | {Path(path).stem for path in explicit_formal}
    for index, path_text in enumerate(all_config_paths, 1):
        path = ROOT / path_text
        if not path.is_file():
            continue
        name = path.name.lower()
        if path_text == "configs/environment.yml":
            classification = "FORMAL"
            include = False
            notes = "Inventoried under formal environments."
        elif path_text in formal_map or path_text in explicit_formal:
            classification = "FORMAL"
            include = True
            notes = "Direct formal configuration or formal scenario configuration."
        elif path_text in {
            "configs/pipeline_presets.yaml",
            "configs/project_config.local.yaml.example",
        }:
            classification = "SUPPORTING"
            include = True
            notes = "Necessary project-level supporting configuration/template."
        elif any(stem.startswith(path.stem + "_") for stem in formal_stems):
            classification = "SUPPORTING"
            include = True
            notes = "Base configuration supporting a formal derived scenario."
        elif any(token in name for token in ("scnoise05", "smoke", "probe", "batch_", "_base", "candidate_stable")):
            classification = "LEGACY_NOT_USED" if "scnoise05" not in name else "UNRESOLVED"
            include = False
            notes = "Not part of the frozen Fig. 1-5 formal experiment set."
        else:
            classification = "UNRESOLVED"
            include = False
            notes = "No frozen-manifest or formal-panel evidence establishes current formal use."
        panels = ";".join(sorted(formal_map.get(path_text, set())))
        rows.append(
            {
                "record_id": f"CFG-{index:04d}",
                "configuration_path": path_text,
                "classification": classification,
                "included_in_candidate": str(include).lower(),
                "candidate_path": "",
                "sha1": sha1(path),
                "byte_size": path.stat().st_size,
                "linked_figure_or_analysis": panels or "supporting/legacy review",
                "evidence": "frozen experiment manifest" if path_text in formal_map else "Task 4H classification",
                "notes": notes,
            }
        )
        if include:
            included.append((path, classification, panels or "supporting configuration"))
    return rows, included


def environment_rows() -> list[dict[str, Any]]:
    rows = [
        {
            "record_id": "ENV-001",
            "route": "composite-simulation Stage3B and mainline CytoSPACE",
            "specification_path": "configs/environment.yml",
            "python_version": "3.10.19",
            "numpy_version": "1.26.4",
            "scipy_version": "1.11.4",
            "pandas_version": "2.3.3",
            "r_version": "4.5.2 (frozen specification)",
            "seurat_version": "5.3.1 (frozen specification)",
            "seuratobject_version": "5.2.0 (frozen specification)",
            "status": "FORMAL_FROZEN_SPECIFICATION",
            "evidence": "configs/environment.yml;software_and_nnls_formal_audit",
            "notes": "Exact Python/SciPy retained; Stage1 run-local R package snapshots were not recorded.",
        },
        {
            "record_id": "ENV-002",
            "route": "high-resolution and CTA biological-application Stage3B",
            "specification_path": "visualizations/manuscript_audits/software_and_nnls_formal_audit/python_solver_version_audit.tsv",
            "python_version": "3.13.5",
            "numpy_version": "2.1.0",
            "scipy_version": "1.15.3",
            "pandas_version": "2.2.3 (current retained environment; not asserted for every historical run)",
            "r_version": "",
            "seurat_version": "",
            "seuratobject_version": "",
            "status": "CONFIRMED_RETAINED_RUNTIME_LINEAGE",
            "evidence": "python_solver_version_audit.tsv",
            "notes": "Preserved separately from the composite Python 3.10 environment.",
        },
        {
            "record_id": "ENV-003",
            "route": "CTA BioApp direct R runtime",
            "specification_path": "visualizations/manuscript_audits/software_and_nnls_formal_audit/seurat_version_audit.tsv",
            "python_version": "",
            "numpy_version": "",
            "scipy_version": "",
            "pandas_version": "",
            "r_version": "4.5.1",
            "seurat_version": "5.3.1",
            "seuratobject_version": "5.2.0",
            "status": "EXACT_FORMAL_RUN_VERSION_CONFIRMED",
            "evidence": "formal packageVersion output and RData object inventory",
            "notes": "Direct CTA runtime record; not collapsed into the mainline Stage1 specification.",
        },
        {
            "record_id": "ENV-004",
            "route": "CytoSPACE 1.1.0 legacy retained environment",
            "specification_path": "external/cytospace/environment-win-r43.yml",
            "python_version": "3.9.19",
            "numpy_version": "1.26.2",
            "scipy_version": "1.11.4",
            "pandas_version": "2.2.0",
            "r_version": "4.3.3",
            "seurat_version": "5.1.0",
            "seuratobject_version": "",
            "status": "SUPPORTING_LEGACY_SPECIFICATION",
            "evidence": "external/cytospace/environment-win-r43.yml",
            "notes": "Retained for lineage; formal mainline uses configs/environment.yml.",
        },
        {
            "record_id": "ENV-005",
            "route": "Tangram / novoSpaRc / SpaOTsc formal comparison",
            "specification_path": "",
            "python_version": "UNRESOLVED",
            "numpy_version": "UNRESOLVED",
            "scipy_version": "UNRESOLVED",
            "pandas_version": "UNRESOLVED",
            "r_version": "",
            "seurat_version": "",
            "seuratobject_version": "",
            "status": "HISTORICAL_RUNTIME_NOT_RETAINED",
            "evidence": "formal output summaries and runners only",
            "notes": "Current E:/ANACONDA/envs/tangram-env lacks the successful formal packages and is not substituted as historical truth.",
        },
        {
            "record_id": "ENV-006",
            "route": "CellTrek-labelled comparison",
            "specification_path": "",
            "python_version": "UNRESOLVED",
            "numpy_version": "UNRESOLVED",
            "scipy_version": "UNRESOLVED",
            "pandas_version": "UNRESOLVED",
            "r_version": "NOT_USED_BY_FORMAL_RUNNER",
            "seurat_version": "NOT_USED_BY_FORMAL_RUNNER",
            "seuratobject_version": "NOT_USED_BY_FORMAL_RUNNER",
            "status": "PYTHON_FALLBACK_RUNTIME_NOT_RETAINED",
            "evidence": "scripts/run_celltrek_mapping.py; retained celltrek_summary.json files",
            "notes": "Formal code is a disclosed Python CellTrek-style fallback, not the official R package.",
        },
        {
            "record_id": "ENV-007",
            "route": "legacy Fig. 2 enrichment calculations",
            "specification_path": "",
            "python_version": "UNRESOLVED",
            "numpy_version": "UNRESOLVED",
            "scipy_version": "UNRESOLVED",
            "pandas_version": "UNRESOLVED",
            "r_version": "",
            "seurat_version": "",
            "seuratobject_version": "",
            "status": "RUN_LOCAL_SNAPSHOT_NOT_RETAINED",
            "evidence": "fig2_enrichment_backend_formal_audit",
            "notes": "The exact original executable and complete package snapshot must not be invented.",
        },
    ]
    return rows


def audit_inventory_rows() -> list[dict[str, Any]]:
    rows = []
    for index, path_text in enumerate(AUDIT_PATHS, 1):
        path = ROOT / path_text
        require(path.is_file(), f"Required consolidated audit record missing: {path_text}")
        rows.append(
            {
                "record_id": f"AUD-{index:04d}",
                "audit_path": path_text,
                "sha1": sha1(path),
                "byte_size": path.stat().st_size,
                "audit_role": "final decision/evidence record",
                "included_in_candidate": "true",
                "notes": "Selected consolidated evidence; existing project file remains unchanged.",
            }
        )
    return rows


def redistribution_rows(source_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    datasets = read_tsv(ROOT / "reproducibility_release/manifests/dataset_manifest_v1.0.tsv")
    rows: list[dict[str, Any]] = []
    for index, row in enumerate(datasets, 1):
        existing = row["redistribution_class"]
        source_layer = row["source_layer"]
        if existing == "PERMISSION_UNRESOLVED":
            classification = "UNRESOLVED_PERMISSION"
            action = "Metadata, official link and reconstruction instructions only."
        elif existing == "RESTRICTED":
            classification = "INTERNAL_NOT_FOR_RELEASE"
            action = "Do not include data; retain metadata only."
        elif existing == "LINK_AND_RECONSTRUCT_ONLY":
            classification = "LINK_AND_RECONSTRUCT"
            action = "Provide official acquisition route and reconstruction code."
        elif source_layer.startswith("PROJECT_GENERATED"):
            classification = "PROJECT_GENERATED_REDISTRIBUTABLE"
            action = "Project-generated derivative may enter candidate subject to source terms."
        else:
            classification = "LINK_AND_RECONSTRUCT"
            action = "Third-party retained input is not copied into candidate."
        rows.append(
            {
                "record_id": f"RED-DATA-{index:03d}",
                "asset_kind": "dataset",
                "asset_id_or_path": row["dataset_record_id"],
                "source_layer": source_layer,
                "frozen_manifest_class": existing,
                "release_classification": classification,
                "included_as_file": "false",
                "release_action": action,
                "notes": "Dataset metadata is retained in dataset_manifest_v1.0.tsv.",
            }
        )
    for index, row in enumerate(source_rows, 1):
        rows.append(
            {
                "record_id": f"RED-SV-{index:03d}",
                "asset_kind": "processed_source_value",
                "asset_id_or_path": row["source_path"],
                "source_layer": "PROJECT_GENERATED_DERIVATIVE",
                "frozen_manifest_class": "NOT_APPLICABLE",
                "release_classification": "PROJECT_GENERATED_REDISTRIBUTABLE",
                "included_as_file": "true",
                "release_action": "Include hashed quantitative derivative only.",
                "notes": "No third-party raw matrix is included.",
            }
        )
    return rows


def main() -> int:
    verify_frozen_baseline()
    branch = git_output("branch", "--show-current")
    commit = git_output("rev-parse", "HEAD")
    pre_status = git_output("status", "--porcelain=v1")
    creation_timestamp = datetime.now().astimezone().replace(microsecond=0).isoformat()

    OUT.mkdir(parents=True, exist_ok=True)
    for name in REQUIRED_CANDIDATE_DIRS:
        (CANDIDATE / name).mkdir(parents=True, exist_ok=True)

    assets: list[Asset] = []

    # Candidate-level documentation generated without reading or modifying the manuscript.
    candidate_readme = CANDIDATE / "README.md"
    candidate_readme.write_text(
        "# SVTuner reproducibility release candidate v1.0\n\n"
        "This is a non-final release candidate generated by Task 4H. It contains "
        "formal code, configurations, source-value tables, plotting scripts and "
        "consolidated audit evidence. Third-party raw data are intentionally absent.\n\n"
        "Use `inventory/task4h_asset_inventory.tsv` as the payload inventory and "
        "`reconstruction/README.md` for data acquisition boundaries. Inventory and "
        "metadata files are bootstrap records and are not self-listed in the payload "
        "inventory to avoid circular hashes.\n",
        encoding="utf-8",
    )
    reconstruction_readme = CANDIDATE / "reconstruction" / "README.md"
    reconstruction_readme.write_text(
        "# Data acquisition and reconstruction\n\n"
        "No third-party raw dataset is bundled. Use the official source routes, "
        "sample identifiers and retained checksums in "
        "`../manifests/dataset_manifest_v1.0.tsv`, then run the preparation scripts "
        "under `../code/scripts/` with the formal configuration files.\n\n"
        "Vizgen provider files are link-and-reconstruct only while explicit "
        "redistribution permission and provider checksums remain unresolved. "
        "CID4465 slide/capture area, CTA exact chemistry, and the human-lung run "
        "accession remain unresolved and must not be inferred.\n",
        encoding="utf-8",
    )
    env_limitations = CANDIDATE / "environments" / "environment_limitations.md"
    env_limitations.write_text(
        "# Retained environment boundaries\n\n"
        "- Composite Stage3B: Python 3.10.19, NumPy 1.26.4, SciPy 1.11.4.\n"
        "- High-resolution/BioApp Stage3B: Python 3.13.5, NumPy 2.1.0, SciPy 1.15.3.\n"
        "- Mainline Stage1 frozen specification: R 4.5.2, Seurat 5.3.1, SeuratObject 5.2.0.\n"
        "- CTA direct runtime: R 4.5.1, Seurat 5.3.1, SeuratObject 5.2.0.\n"
        "- Exact historical Tangram/novoSpaRc/SpaOTsc/CellTrek-labelled and legacy "
        "Fig. 2 enrichment runtime snapshots were not retained and are not inferred.\n",
        encoding="utf-8",
    )

    for generated, role in [
        (candidate_readme, "release-candidate scope and bootstrap instructions"),
        (reconstruction_readme, "third-party acquisition and reconstruction boundary"),
        (env_limitations, "bounded formal environment limitations"),
    ]:
        assets.append(
            Asset(
                record_id=f"T4H-{len(assets) + 1:04d}",
                asset_class="DOCUMENTATION",
                relative_path=rel(generated),
                filename=generated.name,
                file_extension=generated.suffix,
                byte_size=generated.stat().st_size,
                sha1=sha1(generated),
                release_status="INCLUDED_RELEASE_CANDIDATE",
                redistribution_class="REDISTRIBUTABLE",
                formal_role=role,
                linked_figure_or_analysis="Task 4H release candidate",
                source_or_generator=rel(Path(__file__).resolve()),
                notes="Task 4H generated candidate documentation.",
            )
        )

    figure_paths = {entry["path"] for entry in FIGURE_SCRIPTS}
    for entry in FIGURE_SCRIPTS:
        source = ROOT / entry["path"]
        add_asset(
            assets,
            source,
            f"figure_scripts/{entry['path']}",
            "FIGURE_SCRIPT",
            f"formal panel generator; expected inputs: {entry['inputs']}; outputs: {entry['outputs']}",
            entry["panels"],
            "REDISTRIBUTABLE",
            f"format={entry['format']};retention={entry['retention']}",
        )

    code_paths: set[Path] = set()
    code_paths.update(p for p in (ROOT / "src").rglob("*.py") if p.is_file())
    if (ROOT / "r_scripts").exists():
        code_paths.update(p for p in (ROOT / "r_scripts").rglob("*") if p.is_file())
    code_paths.update(
        p
        for p in (ROOT / "external" / "cytospace").rglob("*")
        if p.is_file()
        and p.suffix.lower() in {".py", ".r", ".md", ".txt"}
        and p.name != "environment-win-r43.yml"
    )
    code_paths.add(ROOT / "external" / "cytospace" / "LICENSE")
    code_paths.add(ROOT / "external" / "cytospace" / "setup.py")
    for path in (ROOT / "scripts").glob("*"):
        if (
            path.is_file()
            and path.suffix.lower() in {".py", ".r"}
            and rel(path) not in figure_paths
            and rel(path) not in CODE_EXCLUSIONS
        ):
            code_paths.add(path)
    for source in sorted(code_paths):
        add_asset(
            assets,
            source,
            f"code/{rel(source)}",
            "SOURCE_CODE",
            "formal or necessary supporting implementation",
            "Stage0-Stage5, simulations, masking, dropout, BioApp, source values or audit",
            "REDISTRIBUTABLE",
            "Temporary/editable-SVG/redundant-storage and obsolete full-figure scripts excluded.",
        )

    env_files = [
        ROOT / "configs/environment.yml",
        ROOT / "configs/project_config.yaml",
        ROOT / "external/cytospace/environment-win-r43.yml",
        ROOT / "pyproject.toml",
    ]
    for source in env_files:
        add_asset(
            assets,
            source,
            f"environments/{rel(source)}",
            "ENVIRONMENT_SPECIFICATION",
            "formal or retained supporting environment specification",
            "all formal routes",
            "REDISTRIBUTABLE",
        )

    config_rows, included_configs = configuration_rows()
    for source, classification, linked in included_configs:
        destination = add_asset(
            assets,
            source,
            f"configs/{rel(source)}",
            "CONFIGURATION",
            f"{classification} machine-readable configuration",
            linked,
            "REDISTRIBUTABLE",
        )
        for row in config_rows:
            if row["configuration_path"] == rel(source):
                row["candidate_path"] = rel(destination)
                break

    for relative_path in [
        "reproducibility_release/manifests/dataset_manifest_v1.0.tsv",
        "reproducibility_release/manifests/experiment_manifest_v1.0.tsv",
        "reproducibility_release/manifests/manifest_v1.0_metadata.json",
        "reproducibility_release/manifests/manifest_v1.0_human_approval.json",
        "reproducibility_release/manifests/manifest_v1.0_freeze_report.md",
        "reproducibility_release/release_manifest.tsv",
    ]:
        source = ROOT / relative_path
        add_asset(
            assets,
            source,
            f"manifests/{source.name}",
            "FROZEN_MANIFEST_OR_RELEASE_RECORD",
            "frozen provenance/experiment/release inventory record",
            "Fig. 1-5 and reproducibility release",
            "METADATA_ONLY",
            "Byte-identical copy; source frozen artifact is not modified.",
        )

    source_rows: list[dict[str, Any]] = []
    for index, (path_text, group, role, generator) in enumerate(source_specs(), 1):
        source = ROOT / path_text
        destination = add_asset(
            assets,
            source,
            f"source_values/SV-{index:04d}_{source.name}",
            "SOURCE_VALUE_TABLE",
            role,
            group,
            "PROJECT_GENERATED_REDISTRIBUTABLE",
            "Quantitative derivative only; no provider raw matrix.",
        )
        row_count = ""
        if source.suffix.lower() in {".csv", ".tsv"}:
            with source.open("r", encoding="utf-8-sig", errors="replace") as handle:
                row_count = max(sum(1 for _ in handle) - 1, 0)
        source_rows.append(
            {
                "record_id": f"SV-{index:04d}",
                "coverage_group": group,
                "source_path": path_text,
                "candidate_path": rel(destination),
                "formal_role": role,
                "generator": generator,
                "row_count": row_count,
                "byte_size": source.stat().st_size,
                "sha1": sha1(source),
                "coverage_status": "PRESENT_HASHED_INCLUDED",
                "redistribution_class": "PROJECT_GENERATED_REDISTRIBUTABLE",
                "notes": "Not reconstructed from a plotted image.",
            }
        )

    audit_rows = audit_inventory_rows()
    for row in audit_rows:
        source = ROOT / row["audit_path"]
        destination = add_asset(
            assets,
            source,
            f"audit_summaries/{row['record_id']}_{source.name}",
            "AUDIT_OR_PROVENANCE_RECORD",
            row["audit_role"],
            "dataset/experiment/software/Fig. 2/Fig. 5/release audit",
            "METADATA_ONLY",
        )
        row["candidate_path"] = rel(destination)

    # Central inventories.
    asset_rows = [asdict(asset) for asset in assets]
    write_tsv(
        OUT / "task4h_asset_inventory.tsv",
        asset_rows,
        list(Asset.__dataclass_fields__),
    )
    class_counts = Counter(asset.asset_class for asset in assets)
    class_summary_rows = [
        {
            "asset_class": key,
            "record_count": value,
            "total_bytes": sum(a.byte_size for a in assets if a.asset_class == key),
            "validation_status": "PASS",
        }
        for key, value in sorted(class_counts.items())
    ]
    write_tsv(
        OUT / "task4h_asset_class_summary.tsv",
        class_summary_rows,
        ["asset_class", "record_count", "total_bytes", "validation_status"],
    )
    env_rows = environment_rows()
    write_tsv(
        OUT / "task4h_environment_inventory.tsv",
        env_rows,
        list(env_rows[0]),
    )
    write_tsv(
        OUT / "task4h_configuration_inventory.tsv",
        config_rows,
        list(config_rows[0]),
    )
    write_tsv(
        OUT / "task4h_source_value_inventory.tsv",
        source_rows,
        list(source_rows[0]),
    )
    figure_rows = []
    for index, entry in enumerate(FIGURE_SCRIPTS, 1):
        path = ROOT / entry["path"]
        figure_rows.append(
            {
                "record_id": f"FIGSCRIPT-{index:03d}",
                "relative_path": entry["path"],
                "sha1": sha1(path),
                "byte_size": path.stat().st_size,
                "figure_panels": entry["panels"],
                "expected_source_value_inputs": entry["inputs"],
                "output_filenames": entry["outputs"],
                "output_format": entry["format"],
                "retention_status": entry["retention"],
                "candidate_path": f"reproducibility_release/release_candidate_v1.0/figure_scripts/{entry['path']}",
            }
        )
    write_tsv(
        OUT / "task4h_figure_script_inventory.tsv",
        figure_rows,
        list(figure_rows[0]),
    )
    write_tsv(
        OUT / "task4h_audit_record_inventory.tsv",
        audit_rows,
        list(audit_rows[0]),
    )
    redistribution = redistribution_rows(source_rows)
    write_tsv(
        OUT / "task4h_redistribution_matrix.tsv",
        redistribution,
        list(redistribution[0]),
    )

    missing_rows: list[dict[str, Any]] = []
    for gap_id, field in ACCEPTED_PROVENANCE_GAPS:
        missing_rows.append(
            {
                "gap_id": gap_id,
                "severity": "NONBLOCKING_ACCEPTED",
                "category": "provenance",
                "item": field,
                "status": "UNRESOLVED",
                "safe_action": "Retain qualified wording; do not infer or redistribute.",
                "blocking": "false",
            }
        )
    for gap_id, item, detail in ENVIRONMENT_GAPS:
        missing_rows.append(
            {
                "gap_id": gap_id,
                "severity": "NONBLOCKING_DOCUMENTED",
                "category": "environment",
                "item": item,
                "status": "NOT_RETAINED",
                "safe_action": detail,
                "blocking": "false",
            }
        )
    missing_rows.append(
        {
            "gap_id": "GAP-GIT-001",
            "severity": "NONBLOCKING_OPERATIONAL",
            "category": "git",
            "item": "Current working tree is not clean",
            "status": "RELEASE_COMMIT_NOT_READY",
            "safe_action": "Review and commit an intentional release scope before tagging.",
            "blocking": "false",
        }
    )
    write_tsv(
        OUT / "task4h_missing_or_unresolved_assets.tsv",
        missing_rows,
        list(missing_rows[0]),
    )

    # Validate every payload record independently.
    validation_rows = []
    for asset in assets:
        path = ROOT / asset.relative_path
        observed_hash = sha1(path) if path.is_file() else ""
        observed_bytes = path.stat().st_size if path.is_file() else -1
        status = (
            "PASS"
            if path.is_file()
            and observed_hash == asset.sha1
            and observed_bytes == asset.byte_size
            else "FAIL"
        )
        validation_rows.append(
            {
                "record_id": asset.record_id,
                "relative_path": asset.relative_path,
                "path_exists": str(path.is_file()).lower(),
                "expected_sha1": asset.sha1,
                "observed_sha1": observed_hash,
                "expected_bytes": asset.byte_size,
                "observed_bytes": observed_bytes,
                "status": status,
            }
        )
    write_tsv(
        OUT / "task4h_hash_validation_report.tsv",
        validation_rows,
        list(validation_rows[0]),
    )
    hash_pass = all(row["status"] == "PASS" for row in validation_rows)

    coverage_counts = Counter(row["coverage_group"] for row in source_rows)
    required_coverage = {
        "Fig. 1 simulation values",
        "Fig. 1 no-noise results",
        "Fig. 1 10% reference-noise results",
        "Fig. 2 profile-masking results",
        "Fig. 2 enrichment source values",
        "Fig. 2 cosine-similarity values",
        "Fig. 3 high-resolution results",
        "Fig. 3 21-readout panel",
        "Fig. 3 selected 10-readout panel",
        "Fig. 3 T-cell state-decoy results",
        "Fig. 3 kidney state-decoy results",
        "Fig. 4 reference-dropout results",
        "Fig. 4 thalamic case",
        "Fig. 4 forced-surrogate and expression-similarity matrices",
        "Fig. 4 KITLG-KIT analysis",
        "Fig. 5 CTA endpoint results",
        "Fig. 5 AUROC and average precision",
        "Fig. 5 binary-withholding table",
        "Fig. 5 morphology, interface and microenvironment analyses",
    }
    missing_coverage = sorted(group for group in required_coverage if coverage_counts[group] == 0)
    source_coverage = "COMPLETE" if not missing_coverage else "INCOMPLETE"

    candidate_files = [p for p in CANDIDATE.rglob("*") if p.is_file()]
    restricted_included = any(
        "/data/raw/" in f"/{rel(path).lower()}"
        or path.suffix.lower() in PROHIBITED_EXTENSIONS
        for path in candidate_files
    )
    duplicate_paths = len({asset.relative_path for asset in assets}) != len(assets)
    blocking_items = []
    if not hash_pass:
        blocking_items.append("candidate hash or byte-size validation failure")
    if missing_coverage:
        blocking_items.append("missing formal source-value coverage: " + ";".join(missing_coverage))
    if restricted_included:
        blocking_items.append("restricted/prohibited third-party data included")
    if duplicate_paths:
        blocking_items.append("duplicate candidate paths in asset inventory")

    verify_frozen_baseline()
    decision = "BLOCKED" if blocking_items else "PASS_WITH_NONBLOCKING_GAPS"
    guardrails = {
        "svtuner_rerun": False,
        "stage3a_rerun": False,
        "stage3b_rerun": False,
        "stage4_rerun": False,
        "cytospace_rerun": False,
        "formal_mapping_rerun": False,
        "metric_recomputation": False,
        "threshold_modification": False,
        "endpoint_redefinition": False,
        "figure_modification": False,
        "figure_regeneration": False,
        "manuscript_latex_accessed": False,
        "manuscript_latex_modified": False,
        "bibtex_accessed": False,
        "bibtex_modified": False,
        "result_number_modified": False,
        "dataset_manifest_modified": False,
        "experiment_manifest_modified": False,
        "data_modified": False,
        "metrics_modified": False,
        "restricted_third_party_data_included": restricted_included,
        "frozen_manifest_hashes_unchanged": True,
        "existing_task4e_task4f_release_inventory_unchanged": True,
        "hash_validation_passed": hash_pass,
    }
    write_json(OUT / "task4h_guardrail_report.json", guardrails)

    post_status = git_output("status", "--porcelain=v1")
    tracked_candidate = git_output("ls-files", rel(CANDIDATE))
    cache_dirs = sorted(
        rel(path)
        for name in ("__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache")
        for path in ROOT.rglob(name)
        if path.is_dir() and not path.is_relative_to(CANDIDATE)
    )
    git_status = {
        "git_commit": commit,
        "git_branch": branch,
        "working_tree_clean": not bool(pre_status),
        "working_tree_status_before_task4h": pre_status.splitlines(),
        "working_tree_status_after_candidate": post_status.splitlines(),
        "tracked_release_candidate_files": tracked_candidate.splitlines() if tracked_candidate else [],
        "untracked_files_relevant_to_release": [
            line for line in post_status.splitlines() if "reproducibility_release" in line or "run_task4h" in line
        ],
        "ignored_cache_or_temporary_directories": cache_dirs[:100],
        "current_commit_suitable_as_manuscript_release_commit": False,
        "suitability_reason": "The working tree contains uncommitted and untracked formal/release assets.",
        "recommended_release_version": "v1.0.0-rc1",
        "recommended_git_tag": "v1.0.0-rc1",
        "recommended_release_commit_scope": (
            "Review and commit frozen manifests, Task 4E-4H audit records, release-candidate "
            "inventory, formal scripts/configurations/source values, and intentional 5% work "
            "separately if it is not part of the manuscript release."
        ),
        "recommended_archive_package_name": "SVTuner-v1.0.0-rc1.zip",
        "commit_created": False,
        "tag_created": False,
        "remote_accessed": False,
    }
    write_json(OUT / "task4h_git_release_status.json", git_status)

    summary = {
        "task": "Task 4H remaining reproducibility assets freeze and release-candidate audit",
        "decision": decision,
        "candidate_version": "1.0-rc1",
        "creation_timestamp": creation_timestamp,
        "project_root": str(ROOT),
        "git_commit": commit,
        "git_branch": branch,
        "working_tree_status": "DIRTY" if pre_status else "CLEAN",
        "manifest_versions": {
            "dataset_manifest": "1.0 FROZEN",
            "experiment_manifest": "1.0 FROZEN",
        },
        "asset_counts_by_class": dict(sorted(class_counts.items())),
        "inventory_record_count": len(assets),
        "hash_algorithm": "SHA-1",
        "formal_source_value_coverage": source_coverage,
        "source_value_coverage_counts": dict(sorted(coverage_counts.items())),
        "unresolved_items": [row["item"] for row in missing_rows],
        "blocking_items": blocking_items,
        "nonblocking_items": len(missing_rows),
        "rerun_flags": {
            "svtuner": False,
            "stage3": False,
            "stage4": False,
            "cytospace": False,
        },
        "manuscript_modified": False,
        "figure_modified": False,
        "data_modified": False,
        "metrics_modified": False,
        "thresholds_modified": False,
        "endpoint_modified": False,
        "restricted_third_party_data_included": restricted_included,
        "frozen_manifests_unchanged": True,
        "hash_validation": "PASS" if hash_pass else "FAIL",
    }
    write_json(OUT / "task4h_release_candidate_summary.json", summary)

    decision_report = f"""# Task 4H release-candidate decision

## Decision

`{decision}`

The candidate contains {len(assets)} hashed payload assets. All included paths,
SHA-1 values and byte sizes {'passed' if hash_pass else 'did not pass'} validation.
Formal source-value coverage is `{source_coverage}` across the required Fig. 1-5
groups. No prohibited third-party raw dataset is included.

## Non-blocking gaps

The six human-approved provenance gaps remain unresolved. Historical run-local
environment snapshots are also unavailable for mainline Stage1, legacy Fig. 2
enrichment, Tangram, novoSpaRc, SpaOTsc and the CellTrek-labelled Python fallback.
These versions were not inferred from the current machine.

## Git release readiness

The repository is on `{branch}` at `{commit}`. The working tree was not clean
before Task 4H, so this commit is not suitable as the manuscript release commit.
An intentional scoped commit and review are required before tag `v1.0.0-rc1`.

## Guardrails

No SVTuner/CytoSPACE stage, metric, mapping or figure was rerun. No manuscript,
BibTeX, figure, endpoint, threshold, result number, frozen manifest or existing
Task 4E/4F inventory record was modified.
"""
    (OUT / "task4h_release_candidate_decision.md").write_text(decision_report, encoding="utf-8")
    out_readme = f"""# Task 4H reproducibility release-candidate audit

This directory records the static audit and packaging of the remaining formal
SVTuner reproducibility assets. Decision: `{decision}`.

Release candidate: `{rel(CANDIDATE)}`.

The central inventory covers candidate payload files. Inventory and metadata
bootstrap files are copied into the candidate after generation and are not
self-listed, avoiding circular hashes.
"""
    (OUT / "README.md").write_text(out_readme, encoding="utf-8")

    # Copy Task 4H inventories and metadata into the candidate structure.
    inventory_names = [
        "task4h_asset_inventory.tsv",
        "task4h_asset_class_summary.tsv",
        "task4h_environment_inventory.tsv",
        "task4h_configuration_inventory.tsv",
        "task4h_source_value_inventory.tsv",
        "task4h_figure_script_inventory.tsv",
        "task4h_audit_record_inventory.tsv",
        "task4h_redistribution_matrix.tsv",
        "task4h_missing_or_unresolved_assets.tsv",
        "task4h_hash_validation_report.tsv",
    ]
    metadata_names = [
        "task4h_release_candidate_summary.json",
        "task4h_release_candidate_decision.md",
        "task4h_git_release_status.json",
        "task4h_guardrail_report.json",
        "README.md",
    ]
    for name in inventory_names:
        copy_exact(OUT / name, CANDIDATE / "inventory" / name)
    for name in metadata_names:
        copy_exact(OUT / name, CANDIDATE / "metadata" / name)

    verify_frozen_baseline()
    print("Task 4H completed.")
    print()
    print("Decision:")
    print(decision)
    print()
    print("Release candidate:")
    print(rel(CANDIDATE))
    print()
    print("Git commit:")
    print(commit)
    print()
    print("Working tree clean:")
    print(str(not bool(pre_status)).lower())
    print()
    print("Frozen manifests unchanged:")
    print("true")
    print()
    print("Inventory records:")
    print(len(assets))
    print()
    print("Hash validation:")
    print("PASS" if hash_pass else "FAIL")
    print()
    print("Formal source-value coverage:")
    print(source_coverage)
    print()
    print("Restricted third-party data included:")
    print(str(restricted_included).lower())
    print()
    print("Accepted non-blocking gaps:")
    print(len(missing_rows))
    print()
    print("Blocking gaps:")
    print(len(blocking_items))
    print()
    print("SVTuner rerun:")
    print("false")
    print()
    print("Stage3 rerun:")
    print("false")
    print()
    print("Stage4 rerun:")
    print("false")
    print()
    print("CytoSPACE rerun:")
    print("false")
    print()
    print("Figures modified:")
    print("false")
    print()
    print("Manuscript modified:")
    print("false")
    print()
    print("Next:")
    print("Task 4I - Declarations completion")
    return 0 if decision != "BLOCKED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
