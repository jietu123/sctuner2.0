#!/usr/bin/env python
"""BioApp Phase 3 reference audit and baseline feasibility.

This stage audits reference metadata, immune-related labels, reference dropout
feasibility, and formal baseline input readiness for the frozen Phase 2C
CTA-defined Immune cells endpoint. It must not run CytoSPACE, SVTuner, Stage4,
baseline mapping, contradiction analysis, or prevention analysis.
"""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PHASE2C_DIR = (
    ROOT
    / "visualizations"
    / "bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping"
)
PHASE2BR2_DIR = (
    ROOT
    / "visualizations"
    / "bioapp_phase2br2_explicit_rscript_cta_align_locale_safe_dry_run"
)
OUT_DIR = (
    ROOT
    / "visualizations"
    / "bioapp_phase3_reference_audit_and_formal_baseline_feasibility_for_cta_endpoint"
)

PHASE_NAME = "BioApp Phase 3 reference audit and formal baseline feasibility for CTA-defined endpoint"
STAGE_TYPE = "reference audit / formal baseline feasibility / biological application candidate preparation"

PREFERRED_REFERENCE_DIRS = [
    ROOT
    / "data"
    / "processed"
    / "cytospace_fig2d_tme"
    / "cytospace_fig2d_tme_brca_her2_ffpe"
    / "stage1_preprocess"
    / "exported",
    ROOT
    / "data"
    / "processed"
    / "biological_application"
    / "cytospace_fig2d_tme_brca_her2_ffpe_profile_mask_t_cells_sc_missing_plasma_cells"
    / "stage1_preprocess"
    / "exported",
    ROOT
    / "data"
    / "processed"
    / "low_resolution_experiments"
    / "human_breast_cancer_real"
    / "stage1_preprocess"
    / "exported",
]

PHASE2C_REQUIRED = {
    "summary": "bioapp_phase2c_endpoint_freeze_summary.json",
    "golden_rules": "bioapp_phase2c_golden_rules_v2_1_check.json",
    "endpoint_definition": "endpoint_definition.json",
    "endpoint_freeze": "spot_level_endpoint_freeze.csv",
    "cta_composition": "spot_level_CTA_composition.csv",
    "primary_positive": "primary_endpoint_positive_spots.txt",
    "primary_negative": "primary_endpoint_negative_spots.txt",
    "primary_ambiguous": "primary_endpoint_ambiguous_spots.txt",
    "primary_excluded": "primary_endpoint_excluded_spots.txt",
}

LABEL_COLUMNS = [
    "cell_type",
    "celltype",
    "celltypes",
    "CellType",
    "cellType",
    "annotation",
    "Annotation",
    "label",
    "labels",
    "type",
    "Type",
    "cluster",
    "seurat_clusters",
    "major_cell_type",
    "minor_cell_type",
    "subtype",
]

IMMUNE_KEYWORDS = [
    "immune",
    "T cell",
    "T_cells",
    "T-cell",
    "CD4",
    "CD8",
    "Treg",
    "NK",
    "B cell",
    "B_cells",
    "B-cell",
    "Plasma",
    "Macrophage",
    "Monocyte",
    "Myeloid",
    "Dendritic",
    "DC",
    "Mast",
    "Neutrophil",
    "Lymphocyte",
    "Leukocyte",
]

ALLOWED_CLAIMS = [
    "Reference label vocabulary was audited against the frozen CTA Immune cells endpoint.",
    "A formal baseline feasibility plan was established without running CytoSPACE or SVTuner.",
    "Endpoint-specific metrics were predefined for later baseline/SVTuner comparison.",
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
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def delimiter_for(path: Path) -> str:
    if path.suffix.lower() == ".tsv" or path.name.lower().endswith(".tsv.gz"):
        return "\t"
    if path.suffix.lower() == ".txt":
        return "\t"
    return ","


def count_rows_fast(path: Path) -> int | None:
    try:
        with path.open("rb") as handle:
            return max(sum(1 for _ in handle) - 1, 0)
    except Exception:
        return None


def read_header(path: Path) -> list[str]:
    try:
        return pd.read_csv(path, sep=delimiter_for(path), nrows=0).columns.astype(str).tolist()
    except Exception:
        return []


def validate_phase2c() -> tuple[dict[str, Any], pd.DataFrame]:
    missing = [name for name in PHASE2C_REQUIRED.values() if not (PHASE2C_DIR / name).exists()]
    if missing:
        raise RuntimeError(f"missing frozen endpoint files from Phase 2C: {', '.join(missing)}")
    summary = json.loads((PHASE2C_DIR / PHASE2C_REQUIRED["summary"]).read_text(encoding="utf-8"))
    if (
        summary.get("decision") != "PASS"
        or summary.get("endpoint_frozen") is not True
        or summary.get("primary_endpoint") != "Immune cells"
        or summary.get("ready_for_phase3") is not True
    ):
        raise RuntimeError("Phase 2C not ready for Phase 3")
    if summary.get("CytoSPACE_run") or summary.get("SVTuner_run") or summary.get("Stage4_run"):
        raise RuntimeError("boundary violation detected in Phase 2C summary")
    if summary.get("expression_markers_used_to_define_endpoint") or summary.get("mapping_outputs_used_to_define_endpoint"):
        raise RuntimeError("endpoint redefinition boundary violation detected in Phase 2C summary")
    endpoint = pd.read_csv(PHASE2C_DIR / PHASE2C_REQUIRED["endpoint_freeze"])
    return summary, endpoint


def fail_outputs(reason: str) -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary = {
        "phase": PHASE_NAME,
        "decision": "FAIL",
        "reason": reason,
        "stage_type": STAGE_TYPE,
        "input_phase2c_decision": None,
        "endpoint_frozen": False,
        "primary_endpoint": None,
        "reference_found": False,
        "reference_metadata_found": False,
        "reference_expression_found": False,
        "immune_reference_available": False,
        "primary_dropout_feasible": False,
        "size_matched_control_feasible": False,
        "formal_baseline_feasibility": "failed",
        "gene_overlap_checked": False,
        "n_gene_overlap": None,
        "CytoSPACE_run": False,
        "SVTuner_run": False,
        "Stage4_run": False,
        "formal_baseline_mapping_run": False,
        "contradiction_analysis_run": False,
        "prevention_analysis_run": False,
        "endpoint_redefined": False,
        "expression_markers_used_to_define_endpoint": False,
        "mapping_outputs_used_to_define_endpoint": False,
        "endpoint_specific_metric_plan_generated": False,
        "ready_for_phase4": False,
        "allowed_claims": ALLOWED_CLAIMS,
        "disallowed_claims": DISALLOWED_CLAIMS,
    }
    write_json(OUT_DIR / "bioapp_phase3_reference_audit_summary.json", summary)
    write_json(
        OUT_DIR / "bioapp_phase3_golden_rules_v2_1_check.json",
        {
            "stage_type": STAGE_TYPE,
            "biological_question_defined": True,
            "external_endpoint_predefined": True,
            "endpoint_independent_from_SVTuner": True,
            "endpoint_spatially_registered": False,
            "baseline_comparison_available": False,
            "endpoint_specific_improvement_defined": False,
            "endpoint_specific_quantitative_metric_available": False,
            "endpoint_baseline_SVTuner_spatial_comparison_available": False,
            "interpretation_boundary_defined": True,
            "biological_application_allowed": False,
            "allowed_claim_level": "reference audit and baseline feasibility only",
            "decision": "FAIL",
            "reason": reason,
        },
    )
    write_text(
        OUT_DIR / "decision.txt",
        "\n".join(
            [
                "BioApp Phase 3 - reference audit and formal baseline feasibility for CTA-defined endpoint",
                "",
                "Decision: FAIL",
                "",
                f"Reason: {reason}",
                "",
                "Next:",
                "Stop. Fix reference/input audit or boundary violation before retrying Phase 3.",
                "",
            ]
        ),
    )
    print("BioApp Phase 3 completed.")
    print("\nDecision:\nFAIL")
    print("\nNext:\nStop. Fix reference/input audit or boundary violation before retrying Phase 3.")
    return 1


def detect_role(path: Path) -> tuple[str, bool]:
    low = path.as_posix().lower()
    name = path.name.lower()
    if "sc_metadata" in name or "celllabels" in name or "annotation" in low or "metadata" in low:
        return "scRNA_metadata", True
    if "sc_expression" in name or "scrna_gep" in name or ("sc" in low and "expression" in low):
        return "scRNA_expression", True
    if "st_expression" in name or "stdata" in name or "spatial" in low and "matrix" in low:
        return "ST_expression", False
    if "st_coordinates" in name or "coordinates" in name or "tissue_positions" in name:
        return "ST_coordinates", False
    if "cytospace" in low or "stage1_preprocess" in low:
        return "CytoSPACE_input", False
    return "unknown", False


def inventory_files() -> pd.DataFrame:
    candidate_paths: set[Path] = set()
    for path in PHASE2C_DIR.glob("*"):
        if path.is_file():
            candidate_paths.add(path)
    for path in PHASE2BR2_DIR.glob("*"):
        if path.is_file():
            candidate_paths.add(path)
    for base in PREFERRED_REFERENCE_DIRS:
        for name in [
            "sc_metadata.csv",
            "sc_expression_normalized.csv",
            "st_expression_normalized.csv",
            "st_coordinates.csv",
        ]:
            candidate_paths.add(base / name)
    processed = ROOT / "data" / "processed"
    if processed.exists():
        for meta in processed.glob("**/stage1_preprocess/exported/sc_metadata.csv"):
            base = meta.parent
            candidate_paths.add(meta)
            candidate_paths.add(base / "sc_expression_normalized.csv")
            candidate_paths.add(base / "st_expression_normalized.csv")
            candidate_paths.add(base / "st_coordinates.csv")
    for raw_root in [
        ROOT / "data" / "raw" / "新建文件夹",
        ROOT / "data" / "raw" / "bioapp_phase2_cta_minimal_sources",
        ROOT / "data" / "raw" / "cytospace_fig2c_melanoma" / "official_example" / "CytoSPACE_example_melanoma",
    ]:
        if raw_root.exists():
            for path in raw_root.rglob("*"):
                if path.is_file():
                    candidate_paths.add(path)
    rows: list[dict[str, Any]] = []
    for path in sorted(candidate_paths):
            if not path.exists() or not path.is_file():
                continue
            role, candidate = detect_role(path)
            read_success = False
            n_rows = None
            n_cols = None
            notes = ""
            try:
                if path.suffix.lower() in {".csv", ".tsv", ".txt"} and path.stat().st_size <= 150 * 1024 * 1024:
                    header = read_header(path)
                    n_cols = len(header)
                    n_rows = None
                    notes = "header read only for Phase 3 inventory"
                    read_success = bool(header)
                else:
                    read_success = path.exists()
                    notes = "binary or large file inventory only"
            except Exception as exc:  # noqa: BLE001
                notes = f"read failed: {exc}"
            rows.append(
                {
                    "file_path": rel(path),
                    "file_type": path.suffix.lower().lstrip(".") or "file",
                    "file_size_mb": round(path.stat().st_size / (1024 * 1024), 3),
                    "detected_role": role,
                    "read_success": read_success,
                    "n_rows": n_rows,
                    "n_cols": n_cols,
                    "candidate_reference": candidate,
                    "notes": notes,
                }
            )
    return pd.DataFrame(rows).sort_values(["candidate_reference", "detected_role", "file_path"], ascending=[False, True, True])


def choose_reference_dir() -> Path | None:
    for base in PREFERRED_REFERENCE_DIRS:
        if (base / "sc_metadata.csv").exists() and (base / "sc_expression_normalized.csv").exists():
            return base
    candidates = sorted((ROOT / "data").rglob("stage1_preprocess/exported/sc_metadata.csv"))
    for meta in candidates:
        base = meta.parent
        if (base / "sc_expression_normalized.csv").exists():
            return base
    return None


def immune_matches(label: str) -> tuple[bool, list[str], str]:
    low = str(label).lower()
    matched = [kw for kw in IMMUNE_KEYWORDS if kw.lower() in low]
    is_immune = bool(matched)
    if not is_immune:
        return False, [], "Non_immune"
    if any(x in low for x in ["cd4", "cd8", "t cell", "t-cell", "t_cells", "t-cells", "treg", "nk"]):
        return True, matched, "T_NK"
    if any(x in low for x in ["b cell", "b-cell", "b_cells", "plasma"]):
        return True, matched, "B_plasma"
    if any(x in low for x in ["macrophage", "monocyte", "myeloid", "dendritic", "dc", "neutrophil", "mast"]):
        return True, matched, "Myeloid"
    return True, matched, "Other_immune"


def audit_labels(metadata_files: list[Path]) -> tuple[pd.DataFrame, pd.DataFrame, str | None, pd.DataFrame | None, Path | None]:
    label_rows: list[dict[str, Any]] = []
    immune_rows: list[dict[str, Any]] = []
    recommended_col = None
    recommended_df = None
    recommended_file = None
    for meta_path in metadata_files:
        try:
            df = pd.read_csv(meta_path, sep=delimiter_for(meta_path))
        except Exception:
            continue
        candidate_cols = [col for col in LABEL_COLUMNS if col in df.columns]
        for col in candidate_cols:
            counts = df[col].fillna("NA").astype(str).value_counts()
            n_cells = int(len(df))
            n_unique = int(counts.shape[0])
            for label, count in counts.items():
                row = {
                    "source_file": rel(meta_path),
                    "label_column": col,
                    "n_cells": n_cells,
                    "n_unique_labels": n_unique,
                    "label": label,
                    "label_count": int(count),
                    "label_fraction": float(count / n_cells) if n_cells else 0.0,
                }
                label_rows.append(row)
                matched, keywords, group = immune_matches(label)
                immune_rows.append(
                    {
                        "source_file": rel(meta_path),
                        "label_column": col,
                        "label": label,
                        "label_count": int(count),
                        "label_fraction": float(count / n_cells) if n_cells else 0.0,
                        "immune_keyword_matched": matched,
                        "immune_group_suggested": group if matched else "Non_immune",
                        "matched_keywords": ";".join(keywords),
                    }
                )
        if recommended_df is None and "cell_type" in candidate_cols:
            recommended_col = "cell_type"
            recommended_df = df
            recommended_file = meta_path
    return pd.DataFrame(label_rows), pd.DataFrame(immune_rows), recommended_col, recommended_df, recommended_file


def inspect_st_rdata_genes() -> dict[str, Any]:
    rscript = Path("E:/R/R-4.5.1/bin/x64/Rscript.exe")
    rdata = PHASE2BR2_DIR / "_ascii_runtime" / "ST_data.RData"
    if not rscript.exists() or not rdata.exists():
        return {"read_success": False, "reason": "Rscript or ST_data.RData not found"}
    script = OUT_DIR / "phase3_st_rdata_gene_inventory.R"
    out_json = OUT_DIR / "phase3_st_rdata_gene_inventory.json"
    script.write_text(
        r"""
suppressPackageStartupMessages({
  library(Seurat)
  library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
rdata <- args[[1]]
out_json <- args[[2]]
env <- new.env()
objects_loaded <- load(rdata, envir = env)
seurat_name <- NULL
for (nm in objects_loaded) {
  obj0 <- get(nm, envir = env)
  if (inherits(obj0, "Seurat") && is.null(seurat_name)) seurat_name <- nm
}
if (is.null(seurat_name)) stop("No Seurat object found")
obj <- get(seurat_name, envir = env)
assay_name <- DefaultAssay(obj)
mat <- NULL
source <- NULL
try({
  mat <- GetAssayData(obj, assay = assay_name, layer = "data")
  source <- "layer:data"
}, silent = TRUE)
if (is.null(mat)) {
  try({
    mat <- GetAssayData(obj, assay = assay_name, layer = "counts")
    source <- "layer:counts"
  }, silent = TRUE)
}
if (is.null(mat)) {
  try({
    mat <- GetAssayData(obj, assay = assay_name, slot = "data")
    source <- "slot:data"
  }, silent = TRUE)
}
if (is.null(mat)) stop("Could not read assay matrix")
payload <- list(
  read_success = TRUE,
  seurat_object_name = seurat_name,
  assay_name = assay_name,
  matrix_source = source,
  n_ST_genes = nrow(mat),
  n_ST_spots = ncol(mat),
  ST_genes = rownames(mat),
  ST_spots = colnames(mat)
)
writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = FALSE), out_json)
""".strip(),
        encoding="utf-8",
    )
    proc = subprocess.run(
        [str(rscript), str(script), str(rdata), str(out_json)],
        cwd=str(ROOT),
        text=True,
        capture_output=True,
        timeout=20,
    )
    if proc.returncode != 0 or not out_json.exists():
        return {
            "read_success": False,
            "reason": "R ST expression inventory failed",
            "return_code": proc.returncode,
            "stderr": proc.stderr[-2000:],
        }
    return json.loads(out_json.read_text(encoding="utf-8"))


def expression_header(path: Path) -> tuple[list[str], int | None]:
    if not path.exists():
        return [], None
    header = read_header(path)
    rows = count_rows_fast(path) if path.stat().st_size <= 5 * 1024 * 1024 else None
    genes = header[1:] if len(header) > 1 else []
    return genes, rows


def build_dropout_design(ref_df: pd.DataFrame, label_col: str, immune_audit: pd.DataFrame) -> dict[str, Any]:
    counts = ref_df[label_col].fillna("NA").astype(str).value_counts()
    label_flags = {}
    for label in counts.index:
        matched, keywords, group = immune_matches(str(label))
        label_flags[str(label)] = {
            "is_immune": matched,
            "matched_keywords": keywords,
            "immune_group_suggested": group if matched else "Non_immune",
        }
    immune_labels = [label for label, payload in label_flags.items() if payload["is_immune"]]
    nonimmune_labels = [label for label, payload in label_flags.items() if not payload["is_immune"]]
    immune_removed = int(counts.loc[immune_labels].sum()) if immune_labels else 0
    total = int(counts.sum())
    remaining = total - immune_removed
    remaining_nonimmune_labels = len(nonimmune_labels)
    primary_feasible = immune_removed >= 100 and remaining >= 1000 and remaining_nonimmune_labels >= 2

    nonimmune_counts = counts.loc[nonimmune_labels].sort_values(ascending=False) if nonimmune_labels else pd.Series(dtype=int)
    selected_labels: list[str] = []
    selected_count = 0
    for label, count in nonimmune_counts.items():
        if selected_count < immune_removed:
            selected_labels.append(str(label))
            selected_count += int(count)
    if selected_count > immune_removed and len(selected_labels) > 1:
        without_last = selected_count - int(nonimmune_counts.loc[selected_labels[-1]])
        if abs(without_last - immune_removed) < abs(selected_count - immune_removed):
            selected_count = without_last
            selected_labels = selected_labels[:-1]

    if selected_count == 0 and not nonimmune_counts.empty:
        selected_labels = [str(nonimmune_counts.index[0])]
        selected_count = int(nonimmune_counts.iloc[0])

    control_remaining = total - selected_count
    ratio = float(selected_count / immune_removed) if immune_removed else 0.0
    size_feasible = bool(0.8 <= ratio <= 1.25 and control_remaining >= 1000)
    control_note = ""
    if not size_feasible:
        control_note = (
            "Label-level non-immune control cannot size-match the immune-all dropout within the 0.8-1.25 ratio; "
            "manual decision or a later cell-level control design is required."
        )

    return {
        "primary_dropout_design": {
            "name": "reference_immune_all_dropout",
            "definition": "remove all immune-related labels from the scRNA reference",
            "n_cells_removed": immune_removed,
            "n_cells_remaining": remaining,
            "removed_label_list": immune_labels,
            "remaining_label_list": nonimmune_labels,
            "removed_fraction": float(immune_removed / total) if total else 0.0,
            "primary_dropout_feasible": bool(primary_feasible),
        },
        "size_matched_control_dropout_design": {
            "name": "reference_immune_size_matched_nonimmune_dropout_control",
            "definition": "remove non-immune labels approximating the immune removed cell count",
            "candidate_nonimmune_labels": nonimmune_labels,
            "selected_control_labels_or_cells": selected_labels,
            "control_removed_cell_count": int(selected_count),
            "control_remaining_cell_count": int(control_remaining),
            "size_match_ratio": ratio,
            "size_matched_control_feasible": size_feasible,
            "notes": control_note,
        },
    }


def metric_plan() -> list[dict[str, Any]]:
    return [
        {
            "metric_name": "External-endpoint contradiction rate",
            "definition": "Among CTA Immune-positive spots, quantify the fraction or mass assigned by baseline to non-immune compartments.",
            "required_inputs": ["frozen CTA Immune-positive spots", "baseline spot-level immune/non-immune assignment scores"],
            "phase_to_compute": "Phase 4/5 after formal baseline output audit",
            "allowed_claim": "Baseline-vs-endpoint contradiction can be quantified after baseline outputs exist.",
            "disallowed_claim": "Do not claim baseline creates false niche calls in Phase 3.",
        },
        {
            "metric_name": "Endpoint-positive enrichment score",
            "definition": "Compare immune-related assignment score in CTA Immune-positive spots versus CTA Immune-negative spots.",
            "required_inputs": ["frozen endpoint labels", "mapping-derived immune assignment score"],
            "phase_to_compute": "Phase 4/5 after formal baseline output audit",
            "allowed_claim": "Endpoint enrichment metric is predefined before mapping comparison.",
            "disallowed_claim": "Do not report enrichment values in Phase 3.",
        },
        {
            "metric_name": "Endpoint-specific AUROC/AUPRC",
            "definition": "Use mapping-derived immune assignment score to classify CTA Immune-positive versus CTA Immune-negative spots.",
            "required_inputs": ["frozen endpoint labels", "mapping-derived immune assignment score"],
            "phase_to_compute": "Phase 4/5 after formal baseline output audit",
            "allowed_claim": "Classification metric is predefined before formal comparison.",
            "disallowed_claim": "Do not report AUROC/AUPRC in Phase 3.",
        },
        {
            "metric_name": "Contradiction prevention rate",
            "definition": "Among baseline contradicted CTA Immune-positive spots, quantify the fraction converted by SVTuner to reference-unrepresented or withheld.",
            "required_inputs": ["baseline contradiction candidates", "SVTuner Stage3/Stage4 outputs", "structured output audit"],
            "phase_to_compute": "Later SVTuner phase only",
            "allowed_claim": "Prevention metric is defined as a future SVTuner evaluation endpoint.",
            "disallowed_claim": "Do not claim prevention before SVTuner outputs and audit exist.",
        },
    ]


def write_manifest() -> None:
    rows = []
    for path in sorted(OUT_DIR.glob("*")):
        if path.is_file():
            rows.append(
                {
                    "file": rel(path),
                    "type": path.suffix.lstrip(".") or "file",
                    "description": "BioApp Phase 3 generated output",
                    "created_by_phase": "BioApp Phase 3",
                }
            )
    pd.DataFrame(rows).to_csv(OUT_DIR / "manifest.csv", index=False)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    try:
        phase2c_summary, endpoint_df = validate_phase2c()
    except Exception as exc:  # noqa: BLE001
        return fail_outputs(str(exc))

    endpoint_counts = endpoint_df["primary_endpoint_status"].value_counts().to_dict()
    inventory = inventory_files()
    inventory.to_csv(OUT_DIR / "reference_file_inventory.csv", index=False)

    ref_dir = choose_reference_dir()
    if ref_dir is None:
        return fail_outputs("no scRNA reference found")

    metadata_files = [ref_dir / "sc_metadata.csv"]
    label_audit, immune_audit, recommended_col, ref_df, recommended_file = audit_labels(metadata_files)
    label_audit.to_csv(OUT_DIR / "reference_label_audit.csv", index=False)
    immune_audit.to_csv(OUT_DIR / "immune_reference_label_audit.csv", index=False)

    if ref_df is None or recommended_col is None or recommended_file is None:
        return fail_outputs("no reference metadata found")

    total_reference_cells = int(len(ref_df))
    n_reference_labels = int(ref_df[recommended_col].nunique())
    immune_rows = immune_audit[immune_audit["immune_keyword_matched"] == True]  # noqa: E712
    immune_reference_cells = int(immune_rows["label_count"].sum())
    immune_reference_labels = immune_rows["label"].astype(str).tolist()
    immune_reference_available = bool(len(immune_reference_labels) >= 1 and immune_reference_cells >= 100)
    immune_reference_fraction = float(immune_reference_cells / total_reference_cells) if total_reference_cells else 0.0

    dropout = build_dropout_design(ref_df, recommended_col, immune_audit)
    primary_dropout_feasible = bool(dropout["primary_dropout_design"]["primary_dropout_feasible"])
    size_matched_control_feasible = bool(dropout["size_matched_control_dropout_design"]["size_matched_control_feasible"])
    write_json(OUT_DIR / "reference_dropout_design.json", dropout)

    sc_expr = ref_dir / "sc_expression_normalized.csv"
    st_expr_existing = ref_dir / "st_expression_normalized.csv"
    st_coords_existing = ref_dir / "st_coordinates.csv"
    sc_genes, n_sc_cells = expression_header(sc_expr)
    existing_st_genes, n_existing_st_spots = expression_header(st_expr_existing)
    st_rdata_inventory = inspect_st_rdata_genes()
    st_genes = st_rdata_inventory.get("ST_genes", []) if st_rdata_inventory.get("read_success") else []
    st_spots = st_rdata_inventory.get("ST_spots", []) if st_rdata_inventory.get("read_success") else []
    phase2c_barcodes = set(endpoint_df["barcode"].astype(str))
    st_rdata_barcode_overlap = len(phase2c_barcodes & set(map(str, st_spots))) if st_spots else 0
    existing_st_barcode_overlap = None
    if st_expr_existing.exists():
        try:
            existing_st_ids = pd.read_csv(st_expr_existing, usecols=[0]).iloc[:, 0].astype(str)
            existing_st_barcode_overlap = int(len(phase2c_barcodes & set(existing_st_ids)))
        except Exception:
            existing_st_barcode_overlap = None

    gene_overlap = sorted(set(map(str, st_genes)) & set(map(str, sc_genes))) if st_genes and sc_genes else []
    n_gene_overlap = len(gene_overlap)
    gene_overlap_checked = bool(st_genes and sc_genes)
    gene_overlap_df = pd.DataFrame(
        [
            {
                "sc_expression_file": rel(sc_expr),
                "st_expression_source": "Phase 2B-R2 ST_data.RData read-only assay inventory",
                "n_ST_spots": len(st_spots) if st_spots else None,
                "n_ST_genes": len(st_genes) if st_genes else None,
                "n_sc_cells": n_sc_cells,
                "n_sc_genes": len(sc_genes),
                "n_gene_overlap": n_gene_overlap if gene_overlap_checked else None,
                "gene_overlap_fraction_ST": float(n_gene_overlap / len(st_genes)) if st_genes else None,
                "gene_overlap_fraction_sc": float(n_gene_overlap / len(sc_genes)) if sc_genes else None,
                "phase2c_spot_overlap_with_ST_RData": st_rdata_barcode_overlap,
                "phase2c_spot_overlap_with_existing_exported_ST_expression": existing_st_barcode_overlap,
            }
        ]
    )
    if gene_overlap_checked:
        gene_overlap_df.to_csv(OUT_DIR / "gene_overlap_audit.csv", index=False)

    cytospace_entry = ROOT / "external" / "cytospace"
    baseline_items = [
        {
            "item": "Frozen Phase 2C endpoint",
            "status": "pass",
            "evidence_file": rel(PHASE2C_DIR / "spot_level_endpoint_freeze.csv"),
            "value": f"{len(endpoint_df)} spots",
            "required_for_phase4": True,
            "notes": "Primary endpoint is CTA-defined Immune cells.",
        },
        {
            "item": "scRNA metadata / labels",
            "status": "pass" if recommended_file.exists() else "fail",
            "evidence_file": rel(recommended_file),
            "value": f"{total_reference_cells} cells, {n_reference_labels} labels",
            "required_for_phase4": True,
            "notes": f"Recommended label column: {recommended_col}",
        },
        {
            "item": "scRNA expression matrix",
            "status": "pass" if sc_expr.exists() else "fail",
            "evidence_file": rel(sc_expr),
            "value": f"{n_sc_cells} cells, {len(sc_genes)} genes" if sc_genes else "not readable",
            "required_for_phase4": True,
            "notes": "CSV header read only; no mapping executed.",
        },
        {
            "item": "ST expression matrix for frozen endpoint",
            "status": "pass" if st_rdata_inventory.get("read_success") and st_rdata_barcode_overlap == len(endpoint_df) else "incomplete",
            "evidence_file": rel(PHASE2BR2_DIR / "_ascii_runtime" / "ST_data.RData"),
            "value": f"{len(st_spots)} spots, {len(st_genes)} genes" if st_genes else st_rdata_inventory.get("reason", "not readable"),
            "required_for_phase4": True,
            "notes": "Read-only RData assay inventory; no ST expression export or mapping executed.",
        },
        {
            "item": "Existing exported ST expression barcode compatibility",
            "status": "incomplete" if existing_st_barcode_overlap != len(endpoint_df) else "pass",
            "evidence_file": rel(st_expr_existing),
            "value": f"{existing_st_barcode_overlap} / {len(endpoint_df)} frozen endpoint spots overlap",
            "required_for_phase4": False,
            "notes": "Existing exported ST expression is not treated as sufficient unless it covers all frozen endpoint spots.",
        },
        {
            "item": "Gene overlap",
            "status": "pass" if gene_overlap_checked and n_gene_overlap >= 500 else "incomplete",
            "evidence_file": rel(OUT_DIR / "gene_overlap_audit.csv") if gene_overlap_checked else "",
            "value": n_gene_overlap if gene_overlap_checked else "not checked",
            "required_for_phase4": True,
            "notes": "Computed from Phase 2B-R2 ST_data.RData gene names and recommended scRNA expression header.",
        },
        {
            "item": "CytoSPACE backend/entry",
            "status": "pass" if cytospace_entry.exists() else "incomplete",
            "evidence_file": rel(cytospace_entry),
            "value": "available" if cytospace_entry.exists() else "not found",
            "required_for_phase4": True,
            "notes": "Entry exists; no CytoSPACE command executed.",
        },
        {
            "item": "Primary immune dropout design",
            "status": "pass" if primary_dropout_feasible else "incomplete",
            "evidence_file": rel(OUT_DIR / "reference_dropout_design.json"),
            "value": dropout["primary_dropout_design"]["n_cells_removed"],
            "required_for_phase4": True,
            "notes": "Design only; no reference files generated.",
        },
        {
            "item": "Size-matched non-immune control design",
            "status": "pass" if size_matched_control_feasible else "manual_review_required",
            "evidence_file": rel(OUT_DIR / "reference_dropout_design.json"),
            "value": dropout["size_matched_control_dropout_design"]["size_match_ratio"],
            "required_for_phase4": False,
            "notes": dropout["size_matched_control_dropout_design"]["notes"] or "Label-level control feasible.",
        },
    ]
    baseline_df = pd.DataFrame(baseline_items)
    baseline_df.to_csv(OUT_DIR / "formal_baseline_feasibility_audit.csv", index=False)

    formal_sufficient = (
        immune_reference_available
        and primary_dropout_feasible
        and sc_expr.exists()
        and st_rdata_inventory.get("read_success")
        and st_rdata_barcode_overlap == len(endpoint_df)
        and gene_overlap_checked
        and n_gene_overlap >= 500
        and cytospace_entry.exists()
    )
    formal_baseline_feasibility = "sufficient" if formal_sufficient else "incomplete"
    endpoint_specific_metric_plan_generated = True
    write_json(OUT_DIR / "endpoint_specific_metric_plan.json", {"metrics": metric_plan()})

    dry_plan = f"""BioApp Phase 3 baseline command dry plan

This file is a dry plan only. No CytoSPACE command was executed in Phase 3.

Frozen endpoint:
- CTA-defined Immune cells
- Endpoint table: {rel(PHASE2C_DIR / 'spot_level_endpoint_freeze.csv')}

Recommended reference input:
- scRNA metadata: {rel(recommended_file)}
- scRNA expression: {rel(sc_expr)}
- reference label column: {recommended_col}

ST input for future Phase 4:
- source RData: {rel(PHASE2BR2_DIR / '_ascii_runtime' / 'ST_data.RData')}
- read-only audit: n_ST_spots={len(st_spots) if st_spots else 'unknown'}, n_ST_genes={len(st_genes) if st_genes else 'unknown'}
- Phase 4 must export or provide a CytoSPACE-ready ST expression table covering all frozen endpoint spots before running baseline.

Planned branches:
1. full reference baseline would use the recommended scRNA reference without immune dropout.
2. immune-dropout reference baseline would use reference_immune_all_dropout.
3. size-matched control baseline would require manual decision because label-level non-immune size matching is {'feasible' if size_matched_control_feasible else 'not feasible'}.

Command template only:
python <cytospace_wrapper> --st-expression <future_ST_expression_csv> --st-coordinates <future_ST_coordinates_csv> --sc-expression <sc_expression_csv> --sc-metadata <sc_metadata_csv> --label-column {recommended_col} --output <future_output_dir>

Expected output files:
- assigned_locations.csv
- assigned_locations_with_labels.csv
- spot-level immune assignment summaries

No CytoSPACE command was executed in Phase 3.
"""
    write_text(OUT_DIR / "baseline_command_dry_plan.txt", dry_plan)

    if formal_sufficient and size_matched_control_feasible:
        decision = "PASS"
        ready_for_phase4 = True
        reason = ""
        next_text = "BioApp Phase 4 - formal CytoSPACE baseline execution against frozen CTA Immune endpoint"
    elif formal_sufficient:
        decision = "REVIEW_REQUIRED"
        ready_for_phase4 = False
        reason = "size-matched control requires manual decision before Phase 4"
        next_text = "Review reference label vocabulary, dropout design, and baseline input completeness before Phase 4.\nDo not run CytoSPACE/SVTuner/Stage4."
    else:
        decision = "REVIEW_REQUIRED"
        ready_for_phase4 = False
        reason = "formal baseline inputs are incomplete or need manual preparation before Phase 4"
        next_text = "Review reference label vocabulary, dropout design, and baseline input completeness before Phase 4.\nDo not run CytoSPACE/SVTuner/Stage4."

    if not immune_reference_available or ref_df is None:
        decision = "FAIL"
        ready_for_phase4 = False
        reason = "no usable immune-related scRNA reference labels found"
        formal_baseline_feasibility = "failed"
        next_text = "Stop. Fix reference/input audit or boundary violation before retrying Phase 3."

    summary = {
        "phase": PHASE_NAME,
        "decision": decision,
        "reason": reason,
        "stage_type": STAGE_TYPE,
        "input_phase2c_decision": phase2c_summary.get("decision"),
        "endpoint_frozen": True,
        "primary_endpoint": "Immune cells",
        "primary_endpoint_n_positive": int(endpoint_counts.get("positive", 0)),
        "primary_endpoint_n_negative": int(endpoint_counts.get("negative", 0)),
        "primary_endpoint_n_ambiguous": int(endpoint_counts.get("ambiguous", 0)),
        "primary_endpoint_n_excluded": int(endpoint_counts.get("excluded", 0)),
        "primary_endpoint_spatial_note": "sparse immune-positive spatial compartments, not a large continuous ROI",
        "reference_found": bool(ref_dir is not None),
        "reference_metadata_found": bool(recommended_file is not None and recommended_file.exists()),
        "reference_expression_found": bool(sc_expr.exists()),
        "recommended_reference_dir": rel(ref_dir),
        "recommended_reference_label_column": recommended_col,
        "total_reference_cells": total_reference_cells,
        "n_reference_labels": n_reference_labels,
        "immune_reference_available": immune_reference_available,
        "immune_reference_cells": immune_reference_cells,
        "immune_reference_fraction": immune_reference_fraction,
        "immune_reference_labels": immune_reference_labels,
        "primary_dropout_design": "reference_immune_all_dropout",
        "primary_dropout_feasible": primary_dropout_feasible,
        "size_matched_control_feasible": size_matched_control_feasible,
        "formal_baseline_feasibility": formal_baseline_feasibility,
        "gene_overlap_checked": gene_overlap_checked,
        "n_gene_overlap": n_gene_overlap if gene_overlap_checked else None,
        "CytoSPACE_run": False,
        "SVTuner_run": False,
        "Stage4_run": False,
        "formal_baseline_mapping_run": False,
        "contradiction_analysis_run": False,
        "prevention_analysis_run": False,
        "endpoint_redefined": False,
        "expression_markers_used_to_define_endpoint": False,
        "mapping_outputs_used_to_define_endpoint": False,
        "endpoint_specific_metric_plan_generated": endpoint_specific_metric_plan_generated,
        "ready_for_phase4": ready_for_phase4,
        "allowed_claims": ALLOWED_CLAIMS,
        "disallowed_claims": DISALLOWED_CLAIMS,
    }
    write_json(OUT_DIR / "bioapp_phase3_reference_audit_summary.json", summary)

    golden = {
        "stage_type": STAGE_TYPE,
        "biological_question_defined": True,
        "external_endpoint_predefined": True,
        "endpoint_independent_from_SVTuner": True,
        "endpoint_spatially_registered": True,
        "baseline_comparison_available": False,
        "endpoint_specific_improvement_defined": True,
        "endpoint_specific_quantitative_metric_available": False,
        "endpoint_baseline_SVTuner_spatial_comparison_available": False,
        "interpretation_boundary_defined": True,
        "biological_application_allowed": False,
        "allowed_claim_level": "reference audit and baseline feasibility only",
        "decision": decision,
    }
    write_json(OUT_DIR / "bioapp_phase3_golden_rules_v2_1_check.json", golden)

    decision_text = f"""BioApp Phase 3 - reference audit and formal baseline feasibility for CTA-defined endpoint

Decision: {decision}

Input endpoint:
CTA-defined Immune cells

Endpoint status:
Frozen at Phase 2C

Primary endpoint note:
Sparse immune-positive spatial compartments, not a large continuous ROI.

Reference audit:
reference_found = {str(ref_dir is not None).lower()}
reference_metadata_found = {str(recommended_file is not None and recommended_file.exists()).lower()}
recommended_label_column = {recommended_col}
immune_reference_available = {str(immune_reference_available).lower()}

Dropout design:
primary_dropout_design = reference_immune_all_dropout
primary_dropout_feasible = {str(primary_dropout_feasible).lower()}
size_matched_control_feasible = {str(size_matched_control_feasible).lower()}

Formal baseline feasibility:
{formal_baseline_feasibility}

Boundary checks:
CytoSPACE run: false
SVTuner run: false
Stage4 run: false
Formal baseline mapping run: false
Contradiction analysis run: false
Prevention analysis run: false
Endpoint redefined: false

Allowed claims:
{chr(10).join('- ' + x for x in ALLOWED_CLAIMS)}

Disallowed claims:
{chr(10).join('- ' + x for x in DISALLOWED_CLAIMS)}

Next:
{next_text}
"""
    write_text(OUT_DIR / "decision.txt", decision_text)

    readme = f"""# BioApp Phase 3 Reference Audit

Purpose: audit whether the frozen Phase 2C CTA-defined Immune cells endpoint has enough reference and input support for later formal baseline-vs-endpoint comparison.

Phase 2C input: `{rel(PHASE2C_DIR)}`.

Reference audit method: scan candidate files, read recommended scRNA metadata, audit label vocabulary, and match immune-related labels using predefined keywords.

Immune label keywords: {", ".join(IMMUNE_KEYWORDS)}.

Dropout design: primary design removes all immune-related labels; size-matched non-immune control is only designed, not executed.

Formal baseline feasibility: checks ST expression availability, scRNA expression, metadata, gene overlap, barcode compatibility, and CytoSPACE entry availability without running CytoSPACE.

Endpoint-specific metrics are only planned in `endpoint_specific_metric_plan.json`; no metric values are computed in Phase 3.

Decision: {decision}
Reason: {reason}

Allowed claims:
{chr(10).join('- ' + x for x in ALLOWED_CLAIMS)}

Disallowed claims:
{chr(10).join('- ' + x for x in DISALLOWED_CLAIMS)}
"""
    write_text(OUT_DIR / "README.md", readme)

    write_manifest()

    print("BioApp Phase 3 completed.")
    print("\nDecision:")
    print(decision)
    print("\nInput endpoint:")
    print("CTA-defined Immune cells")
    print("\nEndpoint frozen:")
    print("true")
    print("\nPrimary endpoint note:")
    print("sparse immune-positive spatial compartments, not a large continuous ROI")
    print("\nReference audit:")
    print(f"reference_found = {str(ref_dir is not None).lower()}")
    print(f"reference_metadata_found = {str(recommended_file is not None and recommended_file.exists()).lower()}")
    print(f"recommended_label_column = {recommended_col}")
    print(f"total_reference_cells = {total_reference_cells}")
    print(f"n_reference_labels = {n_reference_labels}")
    print(f"immune_reference_available = {str(immune_reference_available).lower()}")
    print(f"immune_reference_cells = {immune_reference_cells}")
    print(f"immune_reference_labels = {immune_reference_labels}")
    print("\nDropout design:")
    print(f"primary_dropout_feasible = {str(primary_dropout_feasible).lower()}")
    print(f"size_matched_control_feasible = {str(size_matched_control_feasible).lower()}")
    print("\nFormal baseline feasibility:")
    print(formal_baseline_feasibility)
    print("\nBoundary checks:")
    print("CytoSPACE run: false")
    print("SVTuner run: false")
    print("Stage4 run: false")
    print("Formal baseline mapping run: false")
    print("Contradiction analysis run: false")
    print("Prevention analysis run: false")
    print("Endpoint redefined: false")
    print("\nNext:")
    print(next_text)
    return 0 if decision in {"PASS", "REVIEW_REQUIRED"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
