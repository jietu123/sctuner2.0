#!/usr/bin/env python3
"""BioApp Phase 3R: resolve CTA endpoint baseline feasibility review items.

This phase is deliberately non-destructive:
- no CytoSPACE execution
- no SVTuner execution
- no Stage4 execution
- no formal contradiction/prevention computation
- no endpoint redefinition
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase3r_review_resolution_for_cta_endpoint_baseline_feasibility"
PHASE2C_DIR = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase2c_endpoint_freeze_with_validated_cta_to_spot_mapping"
PHASE3_DIR = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase3_reference_audit_and_formal_baseline_feasibility_for_cta_endpoint"
PHASE2BR2_ASCII = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase2br2_explicit_rscript_cta_align_locale_safe_dry_run" / "_ascii_runtime" / "ST_data.RData"
REF_DIR = ROOT / "data" / "processed" / "cytospace_fig2d_tme" / "cytospace_fig2d_tme_brca_her2_ffpe" / "stage1_preprocess" / "exported"
RSCRIPT = Path("E:/R/R-4.5.1/bin/x64/Rscript.exe")

PHASE = "BioApp Phase 3R review resolution for CTA endpoint baseline feasibility"
STAGE_TYPE = "baseline input harmonization / control-design review / biological application candidate preparation"
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
    "Phase 3 review items were resolved or explicitly identified.",
    "ST/reference input harmonization was audited for the frozen CTA Immune endpoint.",
    "A CytoSPACE input manifest was prepared without running CytoSPACE.",
    "Control strategy was reviewed without generating mapping results.",
]
DISALLOWED_CLAIMS = [
    "Baseline creates false niche calls.",
    "SVTuner improves endpoint recovery.",
    "SVTuner prevents contradicted interpretation.",
    "Biological application completed.",
    "Biological discovery made.",
]


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def write_json(path: Path, data: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, ensure_ascii=False)
        fh.write("\n")


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return str(path).replace("\\", "/")


def ensure_out() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)


def read_csv_header(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    mode = "rt" if path.suffix == ".gz" else "r"
    with opener(path, mode, encoding="utf-8", newline="") as fh:
        return next(csv.reader(fh))


def count_data_rows(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    mode = "rt" if path.suffix == ".gz" else "r"
    with opener(path, mode, encoding="utf-8", errors="replace") as fh:
        return max(sum(1 for _ in fh) - 1, 0)


def validate_inputs() -> tuple[dict[str, Any], dict[str, Any], list[str]]:
    errors: list[str] = []
    p2c_path = PHASE2C_DIR / "bioapp_phase2c_endpoint_freeze_summary.json"
    p3_path = PHASE3_DIR / "bioapp_phase3_reference_audit_summary.json"
    if not p2c_path.exists():
        errors.append(f"missing Phase 2C summary: {rel(p2c_path)}")
        return {}, {}, errors
    if not p3_path.exists():
        errors.append(f"missing Phase 3 summary: {rel(p3_path)}")
        return {}, {}, errors

    p2c = read_json(p2c_path)
    p3 = read_json(p3_path)
    checks = [
        (p2c.get("decision") == "PASS", "Phase 2C decision is not PASS"),
        (p2c.get("endpoint_frozen") is True, "Phase 2C endpoint_frozen is not true"),
        (p2c.get("primary_endpoint") == PRIMARY_ENDPOINT, "Phase 2C primary endpoint is not Immune cells"),
        (p2c.get("ready_for_phase3") is True, "Phase 2C ready_for_phase3 is not true"),
        (int(p2c.get("ST_spots_total", -1)) == 2248, "Phase 2C ST_spots_total is not 2248"),
        (p3.get("decision") == "REVIEW_REQUIRED", "Phase 3 decision is not REVIEW_REQUIRED"),
        (p3.get("endpoint_frozen") is True, "Phase 3 endpoint_frozen is not true"),
        (p3.get("primary_endpoint") == PRIMARY_ENDPOINT, "Phase 3 primary endpoint is not Immune cells"),
        (p3.get("reference_found") is True, "Phase 3 reference_found is not true"),
        (p3.get("reference_metadata_found") is True, "Phase 3 reference_metadata_found is not true"),
        (p3.get("reference_expression_found") is True, "Phase 3 reference_expression_found is not true"),
        (p3.get("immune_reference_available") is True, "Phase 3 immune_reference_available is not true"),
        (p3.get("primary_dropout_feasible") is True, "Phase 3 primary_dropout_feasible is not true"),
        (p3.get("ready_for_phase4") is False, "Phase 3 ready_for_phase4 should be false"),
    ]
    for ok, msg in checks:
        if not ok:
            errors.append(msg)
    return p2c, p3, errors


def write_r_export_script() -> Path:
    r_script = OUT_DIR / "phase3r_export_st_expression_from_rdata.R"
    r_code = r'''
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: Rscript phase3r_export_st_expression_from_rdata.R <ST_data.RData> <out_dir>")
}
rdata_path <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

write_lines <- function(x, path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(as.character(x), con = con)
}

env <- new.env(parent = emptyenv())
loaded <- load(rdata_path, envir = env)
obj_name <- NULL
for (nm in loaded) {
  candidate <- get(nm, envir = env)
  if (inherits(candidate, "Seurat")) {
    obj_name <- nm
    obj <- candidate
    break
  }
}
if (is.null(obj_name)) {
  stop("No Seurat object found in ST_data.RData")
}

assay_name <- if ("Spatial" %in% names(obj@assays)) "Spatial" else DefaultAssay(obj)

get_mat <- function(layer_name) {
  mat <- NULL
  try({
    mat <- GetAssayData(obj, assay = assay_name, layer = layer_name)
  }, silent = TRUE)
  if (is.null(mat)) {
    try({
      mat <- GetAssayData(obj, assay = assay_name, slot = layer_name)
    }, silent = TRUE)
  }
  if (is.null(mat)) {
    try({
      mat <- LayerData(obj[[assay_name]], layer = layer_name)
    }, silent = TRUE)
  }
  mat
}

write_mat <- function(mat, genes_by_spots_path, spots_by_genes_path) {
  dense <- as.matrix(mat)
  gb <- as.data.frame(dense, check.names = FALSE)
  gb <- cbind(gene = rownames(dense), gb)
  fwrite(gb, file = genes_by_spots_path, sep = ",", quote = TRUE, compress = "gzip")
  sb <- as.data.frame(t(dense), check.names = FALSE)
  sb <- cbind(spot_id = rownames(sb), sb)
  fwrite(sb, file = spots_by_genes_path, sep = ",", quote = TRUE, compress = "gzip")
}

counts <- get_mat("counts")
data <- get_mat("data")
if (is.null(counts)) {
  stop("No counts layer/slot found in Seurat object")
}

write_mat(
  counts,
  file.path(out_dir, "st_expression_full_2248_counts.csv.gz"),
  file.path(out_dir, "st_expression_full_2248_counts_spots_by_genes.csv.gz")
)

data_exported <- FALSE
data_source <- "missing"
if (!is.null(data) && length(dim(data)) == 2 && nrow(data) > 0 && ncol(data) > 0) {
  write_mat(
    data,
    file.path(out_dir, "st_expression_full_2248_data.csv.gz"),
    file.path(out_dir, "st_expression_full_2248_data_spots_by_genes.csv.gz")
  )
  data_exported <- TRUE
  data_source <- "Seurat data layer"
} else {
  # The raw ST object may contain only counts. For Phase 3R input preparation,
  # export a transparent log-normalized matrix derived from counts so the
  # downstream manifest has a complete expression candidate. This is not used
  # to define or modify the endpoint.
  norm <- log1p(t(t(counts) / Matrix::colSums(counts)) * 10000)
  write_mat(
    norm,
    file.path(out_dir, "st_expression_full_2248_data.csv.gz"),
    file.path(out_dir, "st_expression_full_2248_data_spots_by_genes.csv.gz")
  )
  data_exported <- TRUE
  data_source <- "log1p_CPM_from_counts_in_Phase3R"
}

spot_barcodes <- colnames(counts)
gene_names <- rownames(counts)
write_lines(gene_names, file.path(out_dir, "st_gene_list_full_2248.txt"))
write_lines(spot_barcodes, file.path(out_dir, "st_barcode_list_full_2248.txt"))

meta <- obj@meta.data
meta$barcode <- rownames(meta)
meta <- meta[, c("barcode", setdiff(colnames(meta), "barcode")), drop = FALSE]
fwrite(meta, file = file.path(out_dir, "st_spot_metadata_full_2248.csv"), sep = ",", quote = TRUE)

coord_df <- NULL
image_name <- NA_character_
if (length(names(obj@images)) > 0) {
  image_name <- if ("BCSA2TumB1" %in% names(obj@images)) "BCSA2TumB1" else names(obj@images)[[1]]
  coords <- GetTissueCoordinates(obj, image = image_name)
  coord_df <- as.data.frame(coords, check.names = FALSE)
  coord_df$barcode <- rownames(coord_df)
  coord_df <- coord_df[, c("barcode", setdiff(colnames(coord_df), "barcode")), drop = FALSE]
  fwrite(coord_df, file = file.path(out_dir, "st_coordinates_full_2248.csv"), sep = ",", quote = TRUE)
}

summary_lines <- c(
  paste0("object_name=", obj_name),
  paste0("assay_name=", assay_name),
  paste0("image_name=", image_name),
  paste0("n_spots=", length(spot_barcodes)),
  paste0("n_genes=", length(gene_names)),
  paste0("counts_exported=true"),
  paste0("data_exported=", tolower(as.character(data_exported))),
  paste0("data_source=", data_source),
  paste0("metadata_rows=", nrow(meta)),
  paste0("coordinates_rows=", ifelse(is.null(coord_df), NA, nrow(coord_df)))
)
write_lines(summary_lines, file.path(out_dir, "st_r_export_summary.txt"))
'''
    r_script.write_text(r_code.strip() + "\n", encoding="utf-8")
    return r_script


def run_r_export() -> dict[str, Any]:
    result: dict[str, Any] = {
        "attempted": True,
        "success": False,
        "return_code": None,
        "stdout_path": rel(OUT_DIR / "st_r_export_stdout.txt"),
        "stderr_path": rel(OUT_DIR / "st_r_export_stderr.txt"),
        "error": None,
    }
    if not RSCRIPT.exists():
        result["error"] = f"explicit Rscript not found: {RSCRIPT}"
        return result
    if not PHASE2BR2_ASCII.exists():
        result["error"] = f"ST_data.RData not found: {rel(PHASE2BR2_ASCII)}"
        return result
    r_script = write_r_export_script()
    ascii_out = Path(tempfile.mkdtemp(prefix="svtuner_phase3r_r_export_"))
    result["ascii_temp_output_dir"] = str(ascii_out)
    cmd = [str(RSCRIPT), str(r_script), str(PHASE2BR2_ASCII), str(ascii_out)]
    proc = subprocess.run(cmd, cwd=str(ROOT), capture_output=True, text=True, timeout=900)
    (OUT_DIR / "st_r_export_stdout.txt").write_text(proc.stdout, encoding="utf-8", errors="replace")
    (OUT_DIR / "st_r_export_stderr.txt").write_text(proc.stderr, encoding="utf-8", errors="replace")
    result["return_code"] = proc.returncode
    if proc.returncode == 0:
        for exported in ascii_out.iterdir():
            if exported.is_file():
                shutil.copy2(exported, OUT_DIR / exported.name)
    result["success"] = proc.returncode == 0 and (OUT_DIR / "st_expression_full_2248_counts.csv.gz").exists()
    if not result["success"]:
        result["error"] = f"R export failed with return code {proc.returncode}"
    return result


def audit_st_coverage(endpoint_df: pd.DataFrame) -> dict[str, Any]:
    barcode_file = OUT_DIR / "st_barcode_list_full_2248.txt"
    coord_file = OUT_DIR / "st_coordinates_full_2248.csv"
    if barcode_file.exists():
        st_barcodes = [x.strip() for x in barcode_file.read_text(encoding="utf-8").splitlines() if x.strip()]
    else:
        st_barcodes = []
    st_barcode_set = set(st_barcodes)
    if coord_file.exists():
        coords = pd.read_csv(coord_file)
        coord_barcodes = set(coords["barcode"].astype(str))
    else:
        coords = pd.DataFrame()
        coord_barcodes = set()

    endpoint_barcodes = endpoint_df["barcode"].astype(str).tolist()
    rows = []
    for barcode in endpoint_barcodes:
        in_expr = barcode in st_barcode_set
        in_coord = barcode in coord_barcodes
        rows.append(
            {
                "barcode": barcode,
                "in_phase2c_endpoint": True,
                "in_ST_expression": in_expr,
                "in_ST_coordinates": in_coord,
                "endpoint_status_primary": endpoint_df.loc[endpoint_df["barcode"].astype(str) == barcode, "primary_endpoint_status"].iloc[0],
                "coverage_status": "covered" if in_expr and in_coord else "missing",
                "notes": "" if in_expr and in_coord else "missing from ST expression or coordinates",
            }
        )
    pd.DataFrame(rows).to_csv(OUT_DIR / "st_endpoint_spot_coverage_audit.csv", index=False)

    return {
        "n_phase2c_endpoint_spots": len(endpoint_barcodes),
        "n_ST_expression_spots": len(st_barcodes),
        "n_ST_coordinate_spots": len(coord_barcodes),
        "n_overlap_endpoint_ST_expression": len(set(endpoint_barcodes) & st_barcode_set),
        "n_endpoint_missing_from_ST_expression": len(set(endpoint_barcodes) - st_barcode_set),
        "n_ST_expression_not_in_endpoint": len(st_barcode_set - set(endpoint_barcodes)),
        "n_overlap_endpoint_ST_coordinates": len(set(endpoint_barcodes) & coord_barcodes),
        "n_endpoint_missing_from_ST_coordinates": len(set(endpoint_barcodes) - coord_barcodes),
        "full_2248_coverage_achieved": len(endpoint_barcodes) == 2248
        and set(endpoint_barcodes) == st_barcode_set
        and set(endpoint_barcodes) == coord_barcodes,
    }


def audit_reference() -> dict[str, Any]:
    expr_path = REF_DIR / "sc_expression_normalized.csv"
    meta_path = REF_DIR / "sc_metadata.csv"
    ref: dict[str, Any] = {
        "reference_expression_found": expr_path.exists(),
        "reference_metadata_found": meta_path.exists(),
        "reference_label_column": LABEL_COLUMN,
    }
    if not expr_path.exists() or not meta_path.exists():
        return ref

    meta = pd.read_csv(meta_path)
    if LABEL_COLUMN not in meta.columns:
        ref["reference_label_column_found"] = False
        return ref
    ref["reference_label_column_found"] = True
    ref["total_reference_cells"] = int(len(meta))
    ref["n_reference_labels"] = int(meta[LABEL_COLUMN].nunique())
    ref["n_reference_metadata_duplicate_cells"] = int(meta["cell_id"].duplicated().sum()) if "cell_id" in meta.columns else None

    counts = (
        meta[LABEL_COLUMN]
        .fillna("NA")
        .astype(str)
        .value_counts()
        .rename_axis("cell_type")
        .reset_index(name="n_cells")
    )
    counts["fraction"] = counts["n_cells"] / len(meta)
    counts["immune_label"] = counts["cell_type"].isin(IMMUNE_LABELS)
    counts["immune_group_suggested"] = counts["immune_label"].map({True: "immune", False: "nonimmune"})
    counts.to_csv(OUT_DIR / "reference_cell_type_counts.csv", index=False)

    immune_cells = int(counts.loc[counts["immune_label"], "n_cells"].sum())
    nonimmune_cells = int(counts.loc[~counts["immune_label"], "n_cells"].sum())
    ref["immune_reference_cells"] = immune_cells
    ref["nonimmune_reference_cells"] = nonimmune_cells

    meta.to_csv(OUT_DIR / "reference_metadata_audit.csv", index=False)
    pd.DataFrame(
        [
            {
                "file": rel(expr_path),
                "exists": True,
                "n_rows_cells": count_data_rows(expr_path),
                "n_columns_including_cell_id": len(read_csv_header(expr_path)),
                "first_column": read_csv_header(expr_path)[0],
            }
        ]
    ).to_csv(OUT_DIR / "reference_expression_audit.csv", index=False)
    cell_ids = meta["cell_id"].astype(str).tolist() if "cell_id" in meta.columns else []
    (OUT_DIR / "reference_cell_list.txt").write_text("\n".join(cell_ids) + ("\n" if cell_ids else ""), encoding="utf-8")
    sc_genes = read_csv_header(expr_path)[1:]
    (OUT_DIR / "reference_gene_list.txt").write_text("\n".join(sc_genes) + "\n", encoding="utf-8")
    return ref


def audit_gene_overlap() -> dict[str, Any]:
    st_genes_path = OUT_DIR / "st_gene_list_full_2248.txt"
    sc_genes_path = OUT_DIR / "reference_gene_list.txt"
    if not st_genes_path.exists() or not sc_genes_path.exists():
        return {
            "gene_overlap_checked": False,
            "n_ST_genes": None,
            "n_scRNA_genes": None,
            "n_gene_overlap_exact": None,
            "n_gene_overlap_case_insensitive": None,
            "formal_gene_overlap_sufficient": False,
        }

    st_genes = [g.strip() for g in st_genes_path.read_text(encoding="utf-8").splitlines() if g.strip()]
    sc_genes = [g.strip() for g in sc_genes_path.read_text(encoding="utf-8").splitlines() if g.strip()]
    st_set = set(st_genes)
    sc_set = set(sc_genes)
    exact = sorted(st_set & sc_set)
    st_upper = {g.upper(): g for g in st_set}
    sc_upper = {g.upper(): g for g in sc_set}
    ci_upper = sorted(set(st_upper) & set(sc_upper))
    ci_st = {st_upper[u] for u in ci_upper}
    universe = sorted(st_set | sc_set)
    rows = []
    for gene in universe:
        rows.append(
            {
                "gene": gene,
                "in_ST": gene in st_set,
                "in_scRNA": gene in sc_set,
                "exact_overlap": gene in st_set and gene in sc_set,
                "case_insensitive_overlap": gene.upper() in set(st_upper) and gene.upper() in set(sc_upper),
                "used_for_formal_overlap": gene in exact,
            }
        )
    pd.DataFrame(rows).to_csv(OUT_DIR / "st_sc_gene_overlap_audit.csv", index=False)
    (OUT_DIR / "formal_gene_intersection.txt").write_text("\n".join(exact) + ("\n" if exact else ""), encoding="utf-8")
    (OUT_DIR / "st_only_genes.txt").write_text("\n".join(sorted(st_set - sc_set)) + "\n", encoding="utf-8")
    (OUT_DIR / "scrna_only_genes.txt").write_text("\n".join(sorted(sc_set - st_set)) + "\n", encoding="utf-8")
    return {
        "gene_overlap_checked": True,
        "n_ST_genes": len(st_set),
        "n_scRNA_genes": len(sc_set),
        "n_gene_overlap_exact": len(exact),
        "n_gene_overlap_case_insensitive": len(ci_upper),
        "gene_overlap_fraction_ST": len(exact) / len(st_set) if st_set else 0,
        "gene_overlap_fraction_sc": len(exact) / len(sc_set) if sc_set else 0,
        "formal_gene_overlap_sufficient": len(exact) >= 500,
        "case_insensitive_requires_manual_review": len(exact) < 500 and len(ci_upper) >= 500,
    }


def write_manifest_and_dry_plan(coverage: dict[str, Any], gene: dict[str, Any], ref: dict[str, Any]) -> dict[str, bool]:
    manifest_ready = bool(
        coverage.get("full_2248_coverage_achieved")
        and gene.get("formal_gene_overlap_sufficient")
        and ref.get("reference_expression_found")
        and ref.get("reference_metadata_found")
        and ref.get("reference_label_column_found")
    )
    manifest = {
        "phase": PHASE,
        "dry_plan_only": True,
        "cytospace_run": False,
        "reference_dir": rel(REF_DIR),
        "reference_expression": rel(REF_DIR / "sc_expression_normalized.csv"),
        "reference_metadata": rel(REF_DIR / "sc_metadata.csv"),
        "reference_label_column": LABEL_COLUMN,
        "st_expression_counts_genes_by_spots": rel(OUT_DIR / "st_expression_full_2248_counts.csv.gz"),
        "st_expression_counts_spots_by_genes": rel(OUT_DIR / "st_expression_full_2248_counts_spots_by_genes.csv.gz"),
        "st_expression_data_genes_by_spots": rel(OUT_DIR / "st_expression_full_2248_data.csv.gz"),
        "st_expression_data_spots_by_genes": rel(OUT_DIR / "st_expression_full_2248_data_spots_by_genes.csv.gz"),
        "st_coordinates": rel(OUT_DIR / "st_coordinates_full_2248.csv"),
        "formal_gene_intersection": rel(OUT_DIR / "formal_gene_intersection.txt"),
        "analysis_spot_universe": "full_frozen_endpoint_spot_universe"
        if coverage.get("full_2248_coverage_achieved")
        else "intersected_spot_universe_requires_review",
        "analysis_spot_count": coverage.get("n_phase2c_endpoint_spots")
        if coverage.get("full_2248_coverage_achieved")
        else coverage.get("n_overlap_endpoint_ST_expression"),
        "n_gene_overlap_exact": gene.get("n_gene_overlap_exact"),
        "manifest_ready_for_phase4": manifest_ready,
        "note": "No CytoSPACE command was executed in Phase 3R. This is a dry plan only.",
    }
    write_json(OUT_DIR / "cytospace_input_manifest.json", manifest)
    pd.DataFrame([manifest]).to_csv(OUT_DIR / "cytospace_input_manifest.csv", index=False)

    plan = "\n".join(
        [
            "BioApp Phase 3R baseline command dry plan",
            "",
            "No CytoSPACE command was executed in Phase 3R.",
            "This is a dry plan only.",
            "",
            f"Reference expression: {manifest['reference_expression']}",
            f"Reference metadata: {manifest['reference_metadata']}",
            f"ST expression, spots x genes: {manifest['st_expression_data_spots_by_genes']}",
            f"ST coordinates: {manifest['st_coordinates']}",
            f"Formal gene intersection: {manifest['formal_gene_intersection']}",
            f"Manifest ready for Phase 4: {manifest_ready}",
            "",
            "Future Phase 4 may construct the actual CytoSPACE command from this manifest only after manual authorization.",
        ]
    )
    (OUT_DIR / "baseline_command_dry_plan_phase3r.txt").write_text(plan + "\n", encoding="utf-8")
    return {
        "cytospace_input_manifest_generated": True,
        "cytospace_input_manifest_ready": manifest_ready,
        "baseline_command_dry_plan_generated": True,
    }


def resolve_control_strategy(ref: dict[str, Any]) -> dict[str, Any]:
    immune_removed = int(ref.get("immune_reference_cells") or 0)
    nonimmune_total = int(ref.get("nonimmune_reference_cells") or 0)
    possible = nonimmune_total >= immune_removed and immune_removed > 0
    if possible:
        primary = "true_nonimmune_cell_level_size_matched_control"
        supplementary = None
        limitation = "true nonimmune size-matched control is possible"
    else:
        primary = "nonimmune-all-dropout control"
        supplementary = "random-all-cell size-matched dropout dry design"
        limitation = "true nonimmune size-matched control impossible because nonimmune pool is smaller than immune dropout size"
    evaluated = [
        {
            "strategy": "A",
            "name": "nonimmune-all-dropout control",
            "feasible": nonimmune_total > 0,
            "n_cells_to_remove": nonimmune_total,
            "size_matched_to_immune_all_dropout": nonimmune_total == immune_removed,
            "note": "label-level clean but not size-matched when nonimmune pool is smaller",
        },
        {
            "strategy": "B",
            "name": "random-all-cell size-matched dropout control",
            "feasible": (immune_removed > 0 and (immune_removed <= (immune_removed + nonimmune_total))),
            "n_cells_to_remove": immune_removed,
            "size_matched_to_immune_all_dropout": True,
            "note": "dry design only; would require later fixed seed and explicit cell list",
        },
        {
            "strategy": "C",
            "name": "downsampled immune dropout matched to nonimmune-all removal",
            "feasible": nonimmune_total > 0 and immune_removed >= nonimmune_total,
            "n_cells_to_remove": min(immune_removed, nonimmune_total),
            "size_matched_to_nonimmune_all_dropout": True,
            "note": "does not simulate complete immune reference absence",
        },
        {
            "strategy": "D",
            "name": "no strict size-matched control",
            "feasible": True,
            "n_cells_to_remove": None,
            "note": "explicit limitation if used",
        },
    ]
    resolution = {
        "immune_removed_cell_count": immune_removed,
        "nonimmune_total_cells": nonimmune_total,
        "true_nonimmune_size_matched_control_possible": possible,
        "evaluated_strategies": evaluated,
        "recommended_primary_control": primary,
        "recommended_supplementary_control": supplementary,
        "control_limitations": limitation,
        "requires_manual_approval_before_phase4": True,
        "control_strategy_resolved": immune_removed > 0 and nonimmune_total > 0,
    }
    write_json(OUT_DIR / "control_strategy_resolution.json", resolution)
    lines = [
        "Control strategy resolution",
        "",
        f"immune_removed_cell_count = {immune_removed}",
        f"nonimmune_total_cells = {nonimmune_total}",
        f"true_nonimmune_size_matched_control_possible = {str(possible).lower()}",
        f"recommended_primary_control = {primary}",
        f"recommended_supplementary_control = {supplementary}",
        f"control_limitations = {limitation}",
        "requires_manual_approval_before_phase4 = true",
    ]
    (OUT_DIR / "control_strategy_resolution.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    dropout_designs = {
        "immune_all_dropout_design.json": {
            "design": "reference_immune_all_dropout",
            "labels_to_remove": IMMUNE_LABELS,
            "n_cells_to_remove": immune_removed,
            "n_cells_remaining": nonimmune_total,
            "expected_reference_metadata_filter_rule": f"exclude rows where {LABEL_COLUMN} is in immune labels",
            "expected_reference_expression_filter_rule": "keep expression rows/cells retained by filtered metadata",
            "random_seed_if_later_used": None,
            "phase_to_execute": "future Phase 4 after manual authorization",
            "mapping_result": False,
        },
        "nonimmune_all_dropout_control_design.json": {
            "design": "nonimmune_all_dropout_control",
            "labels_to_remove": "all labels not in immune labels",
            "n_cells_to_remove": nonimmune_total,
            "n_cells_remaining": immune_removed,
            "expected_reference_metadata_filter_rule": f"exclude rows where {LABEL_COLUMN} is not in immune labels",
            "expected_reference_expression_filter_rule": "keep expression rows/cells retained by filtered metadata",
            "random_seed_if_later_used": None,
            "phase_to_execute": "future Phase 4 after manual authorization",
            "mapping_result": False,
        },
        "random_all_cell_size_matched_control_dry_design.json": {
            "design": "random_all_cell_size_matched_control_dry_design",
            "labels_to_remove": "random cells from full reference; not generated in Phase 3R",
            "n_cells_to_remove": immune_removed,
            "n_cells_remaining": immune_removed + nonimmune_total - immune_removed,
            "expected_reference_metadata_filter_rule": "future fixed-seed random removal from full reference",
            "expected_reference_expression_filter_rule": "keep expression rows/cells retained by future random cell list",
            "random_seed_if_later_used": 20260705,
            "phase_to_execute": "future sensitivity stage only after manual authorization",
            "mapping_result": False,
        },
    }
    for filename, design in dropout_designs.items():
        write_json(OUT_DIR / filename, design)
    return resolution


def build_decision(
    input_errors: list[str],
    r_export: dict[str, Any],
    coverage: dict[str, Any],
    gene: dict[str, Any],
    ref: dict[str, Any],
    manifest_info: dict[str, bool],
    control: dict[str, Any],
) -> tuple[str, list[str], bool]:
    reasons: list[str] = []
    boundary_violation = False
    if input_errors:
        return "FAIL", input_errors, boundary_violation
    if not r_export.get("success"):
        reasons.append("ST_data.RData expression export failed")
        return "FAIL", reasons, boundary_violation
    if not coverage.get("full_2248_coverage_achieved"):
        reasons.append("ST expression/coordinates do not cover all 2248 frozen endpoint spots")
    if not gene.get("gene_overlap_checked"):
        reasons.append("gene overlap could not be completed")
    elif not gene.get("formal_gene_overlap_sufficient"):
        reasons.append("n_gene_overlap_exact < 500")
    if gene.get("case_insensitive_requires_manual_review"):
        reasons.append("case-insensitive gene overlap suggests naming issue requiring manual review")
    if not manifest_info.get("cytospace_input_manifest_ready"):
        reasons.append("CytoSPACE input manifest incomplete")
    if not control.get("control_strategy_resolved"):
        reasons.append("control strategy unresolved")

    if reasons:
        return "REVIEW_REQUIRED", reasons, boundary_violation
    return "PASS", ["Phase 3R review items resolved; Phase 4 can be considered after manual authorization"], boundary_violation


def write_readme_and_decision(summary: dict[str, Any], reasons: list[str]) -> None:
    next_text = (
        "BioApp Phase 4 鈥?formal CytoSPACE baseline execution against frozen CTA Immune endpoint"
        if summary["decision"] == "PASS"
        else "Resolve remaining ST/reference input or control-strategy issues before Phase 4.\nDo not run CytoSPACE/SVTuner/Stage4."
        if summary["decision"] == "REVIEW_REQUIRED"
        else "Stop. Fix input harmonization or boundary violation before retrying Phase 3R."
    )
    decision_text = f"""BioApp Phase 3R 鈥?review resolution for CTA endpoint baseline feasibility

Decision: {summary['decision']}

Input endpoint:
CTA-defined Immune cells

Endpoint status:
Frozen at Phase 2C

Primary endpoint note:
{PRIMARY_ENDPOINT_NOTE}

ST expression coverage:
full_2248_coverage_achieved = {str(summary['full_2248_coverage_achieved']).lower()}
n_overlap_endpoint_ST_expression = {summary['n_overlap_endpoint_ST_expression']}
n_endpoint_missing_from_ST_expression = {summary['n_endpoint_missing_from_ST_expression']}

Gene overlap:
n_gene_overlap_exact = {summary['n_gene_overlap_exact']}
formal_gene_overlap_sufficient = {str(summary['formal_gene_overlap_sufficient']).lower()}

CytoSPACE input manifest:
generated = {str(summary['cytospace_input_manifest_generated']).lower()}
ready = {str(summary.get('cytospace_input_manifest_ready', False)).lower()}

Control strategy:
true_nonimmune_size_matched_control_possible = {str(summary['true_nonimmune_size_matched_control_possible']).lower()}
recommended_primary_control = {summary['recommended_primary_control']}
recommended_supplementary_control = {summary['recommended_supplementary_control']}
control_strategy_resolved = {str(summary['control_strategy_resolved']).lower()}

Boundary checks:
CytoSPACE run: false
SVTuner run: false
Stage4 run: false
Formal baseline mapping run: false
Contradiction analysis run: false
Prevention analysis run: false
Endpoint redefined: false

Decision reasons:
{os.linesep.join('- ' + r for r in reasons)}

Allowed claims:
{os.linesep.join('- ' + c for c in ALLOWED_CLAIMS)}

Disallowed claims:
{os.linesep.join('- ' + c for c in DISALLOWED_CLAIMS)}

Next:
{next_text}
"""
    (OUT_DIR / "decision.txt").write_text(decision_text, encoding="utf-8")

    readme = f"""# BioApp Phase 3R

Purpose: resolve Phase 3 review items for the frozen CTA-defined Immune cells endpoint without running mapping or generating biological application results.

## Inputs

- Phase 2C frozen endpoint: `{rel(PHASE2C_DIR)}`
- Phase 3 reference audit: `{rel(PHASE3_DIR)}`
- Explicit Rscript: `{RSCRIPT}`
- ST object source: `{rel(PHASE2BR2_ASCII)}`
- Reference directory: `{rel(REF_DIR)}`

## Main Checks

- ST expression full endpoint coverage: `{summary['full_2248_coverage_achieved']}`
- Gene overlap exact count: `{summary['n_gene_overlap_exact']}`
- CytoSPACE manifest generated: `{summary['cytospace_input_manifest_generated']}`
- Control strategy resolved: `{summary['control_strategy_resolved']}`

## Interpretation Boundary

This phase is input harmonization and control-design review only. It does not run CytoSPACE, SVTuner, Stage4, contradiction analysis, or prevention analysis.

## Decision

`{summary['decision']}`

## Next

{next_text}
"""
    (OUT_DIR / "README.md").write_text(readme, encoding="utf-8")


def write_manifest() -> None:
    rows = []
    for path in sorted(OUT_DIR.iterdir()):
        if path.is_file():
            rows.append(
                {
                    "file": rel(path),
                    "type": path.suffix.lstrip(".") or "text",
                    "description": "BioApp Phase 3R output artifact",
                    "created_by_phase": PHASE,
                    "status": "generated",
                    "notes": "dry plan / audit artifact; no mapping result",
                }
            )
    pd.DataFrame(rows).to_csv(OUT_DIR / "manifest.csv", index=False)


def main() -> int:
    ensure_out()
    p2c, p3, input_errors = validate_inputs()

    endpoint_df = pd.read_csv(PHASE2C_DIR / "spot_level_endpoint_freeze.csv") if not input_errors else pd.DataFrame()

    r_export = {"attempted": False, "success": False}
    coverage = {
        "n_phase2c_endpoint_spots": None,
        "n_ST_expression_spots": None,
        "n_overlap_endpoint_ST_expression": None,
        "n_endpoint_missing_from_ST_expression": None,
        "full_2248_coverage_achieved": False,
    }
    ref = {}
    gene = {}
    manifest_info = {
        "cytospace_input_manifest_generated": False,
        "cytospace_input_manifest_ready": False,
        "baseline_command_dry_plan_generated": False,
    }
    control = {
        "immune_removed_cell_count": None,
        "nonimmune_total_cells": None,
        "true_nonimmune_size_matched_control_possible": False,
        "control_strategy_resolved": False,
        "recommended_primary_control": None,
        "recommended_supplementary_control": None,
    }

    if not input_errors:
        r_export = run_r_export()
        if r_export.get("success"):
            coverage = audit_st_coverage(endpoint_df)
        ref = audit_reference()
        gene = audit_gene_overlap()
        manifest_info = write_manifest_and_dry_plan(coverage, gene, ref)
        control = resolve_control_strategy(ref)

    decision, reasons, boundary_violation = build_decision(input_errors, r_export, coverage, gene, ref, manifest_info, control)
    ready_for_phase4 = decision == "PASS"
    analysis_universe_type = (
        "full_frozen_endpoint_spot_universe"
        if coverage.get("full_2248_coverage_achieved")
        else "intersected_spot_universe_requires_review"
        if coverage.get("n_overlap_endpoint_ST_expression")
        else "unresolved"
    )
    analysis_universe_spots = (
        coverage.get("n_phase2c_endpoint_spots")
        if analysis_universe_type == "full_frozen_endpoint_spot_universe"
        else coverage.get("n_overlap_endpoint_ST_expression")
    )

    summary = {
        "phase": PHASE,
        "decision": decision,
        "decision_reasons": reasons,
        "stage_type": STAGE_TYPE,
        "input_phase2c_decision": p2c.get("decision") if p2c else None,
        "input_phase3_decision": p3.get("decision") if p3 else None,
        "endpoint_frozen": True if p2c.get("endpoint_frozen") is True else False,
        "primary_endpoint": PRIMARY_ENDPOINT,
        "primary_endpoint_spatial_note": PRIMARY_ENDPOINT_NOTE,
        "ST_spots_frozen_endpoint_total": 2248 if not endpoint_df.empty else None,
        "ST_expression_full_2248_export_attempted": r_export.get("attempted", False),
        "ST_expression_full_2248_export_success": r_export.get("success", False),
        "ST_r_export_return_code": r_export.get("return_code"),
        "ST_r_export_error": r_export.get("error"),
        "n_ST_expression_spots": coverage.get("n_ST_expression_spots"),
        "n_overlap_endpoint_ST_expression": coverage.get("n_overlap_endpoint_ST_expression"),
        "n_endpoint_missing_from_ST_expression": coverage.get("n_endpoint_missing_from_ST_expression"),
        "full_2248_coverage_achieved": coverage.get("full_2248_coverage_achieved", False),
        "reference_expression_found": ref.get("reference_expression_found", False),
        "reference_metadata_found": ref.get("reference_metadata_found", False),
        "reference_label_column": LABEL_COLUMN,
        "total_reference_cells": ref.get("total_reference_cells"),
        "n_reference_labels": ref.get("n_reference_labels"),
        "immune_reference_cells": ref.get("immune_reference_cells"),
        "nonimmune_reference_cells": ref.get("nonimmune_reference_cells"),
        "gene_overlap_checked": gene.get("gene_overlap_checked", False),
        "n_ST_genes": gene.get("n_ST_genes"),
        "n_scRNA_genes": gene.get("n_scRNA_genes"),
        "n_gene_overlap_exact": gene.get("n_gene_overlap_exact"),
        "n_gene_overlap_case_insensitive": gene.get("n_gene_overlap_case_insensitive"),
        "formal_gene_overlap_sufficient": gene.get("formal_gene_overlap_sufficient", False),
        "analysis_universe_type": analysis_universe_type,
        "analysis_universe_spots": analysis_universe_spots,
        "cytospace_input_manifest_generated": manifest_info.get("cytospace_input_manifest_generated", False),
        "cytospace_input_manifest_ready": manifest_info.get("cytospace_input_manifest_ready", False),
        "baseline_command_dry_plan_generated": manifest_info.get("baseline_command_dry_plan_generated", False),
        "primary_dropout_design": "reference_immune_all_dropout",
        "immune_removed_cell_count": control.get("immune_removed_cell_count"),
        "nonimmune_total_cells": control.get("nonimmune_total_cells"),
        "true_nonimmune_size_matched_control_possible": control.get("true_nonimmune_size_matched_control_possible", False),
        "control_strategy_resolved": control.get("control_strategy_resolved", False),
        "recommended_primary_control": control.get("recommended_primary_control"),
        "recommended_supplementary_control": control.get("recommended_supplementary_control"),
        "ready_for_phase4": ready_for_phase4,
        "biological_application_allowed": False,
        "CytoSPACE_run": False,
        "SVTuner_run": False,
        "Stage4_run": False,
        "formal_baseline_mapping_run": False,
        "contradiction_analysis_run": False,
        "prevention_analysis_run": False,
        "endpoint_redefined": False,
        "expression_markers_used_to_define_endpoint": False,
        "mapping_outputs_used_to_define_endpoint": False,
        "boundary_violation": boundary_violation,
        "allowed_claims": ALLOWED_CLAIMS,
        "disallowed_claims": DISALLOWED_CLAIMS,
        "output_dir": rel(OUT_DIR),
    }
    write_json(OUT_DIR / "bioapp_phase3r_review_resolution_summary.json", summary)

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
        "allowed_claim_level": "input harmonization and control-design review only",
        "decision": decision,
    }
    write_json(OUT_DIR / "bioapp_phase3r_golden_rules_v2_1_check.json", golden)
    write_readme_and_decision(summary, reasons)
    write_manifest()

    next_text = (
        "BioApp Phase 4 鈥?formal CytoSPACE baseline execution against frozen CTA Immune endpoint"
        if ready_for_phase4
        else "Resolve remaining ST/reference input or control-strategy issues before Phase 4."
    )
    print("BioApp Phase 3R completed.")
    print()
    print("Decision:")
    print(decision)
    print()
    print("Input endpoint:")
    print("CTA-defined Immune cells")
    print()
    print("Endpoint frozen:")
    print(str(summary["endpoint_frozen"]).lower())
    print()
    print("Primary endpoint note:")
    print(PRIMARY_ENDPOINT_NOTE)
    print()
    print("ST expression coverage:")
    print(f"full_2248_coverage_achieved = {str(summary['full_2248_coverage_achieved']).lower()}")
    print(f"n_overlap_endpoint_ST_expression = {summary['n_overlap_endpoint_ST_expression']}")
    print(f"n_endpoint_missing_from_ST_expression = {summary['n_endpoint_missing_from_ST_expression']}")
    print()
    print("Gene overlap:")
    print(f"n_ST_genes = {summary['n_ST_genes']}")
    print(f"n_scRNA_genes = {summary['n_scRNA_genes']}")
    print(f"n_gene_overlap_exact = {summary['n_gene_overlap_exact']}")
    print(f"formal_gene_overlap_sufficient = {str(summary['formal_gene_overlap_sufficient']).lower()}")
    print()
    print("CytoSPACE input manifest:")
    print(f"generated = {str(summary['cytospace_input_manifest_generated']).lower()}")
    print(f"ready_for_phase4 = {str(summary['ready_for_phase4']).lower()}")
    print()
    print("Control strategy:")
    print(f"true_nonimmune_size_matched_control_possible = {str(summary['true_nonimmune_size_matched_control_possible']).lower()}")
    print(f"recommended_primary_control = {summary['recommended_primary_control']}")
    print(f"recommended_supplementary_control = {summary['recommended_supplementary_control']}")
    print(f"control_strategy_resolved = {str(summary['control_strategy_resolved']).lower()}")
    print()
    print("Boundary checks:")
    print("CytoSPACE run: false")
    print("SVTuner run: false")
    print("Stage4 run: false")
    print("Formal baseline mapping run: false")
    print("Contradiction analysis run: false")
    print("Prevention analysis run: false")
    print("Endpoint redefined: false")
    print()
    print("Next:")
    print(next_text)
    return 0 if decision in {"PASS", "REVIEW_REQUIRED"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
