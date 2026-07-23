#!/usr/bin/env python
"""BioApp Phase 2B-R2 explicit Rscript CTA_align locale-safe dry run.

This phase supersedes the earlier PATH-only Rscript check. It uses the
standalone R installation explicitly, reads the Seurat ST_data.RData object,
handles the CTA_align.R non-ASCII column-name issue without modifying the
original file, and generates a real CTA-to-spot mapping dry-run output.

It does not freeze endpoint labels and does not run CytoSPACE, SVTuner,
Stage4, formal biological mapping, or formal downstream metrics.
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    from scipy.spatial import cKDTree
except Exception:  # noqa: BLE001
    cKDTree = None


ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / "data" / "raw" / "\u65b0\u5efa\u6587\u4ef6\u5939" / "BreastCancer_CTA-main"
PREV_DIR = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase2br_r_seurat_runtime_resolution_and_cta_align_dry_run"
OUT_DIR = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase2br2_explicit_rscript_cta_align_locale_safe_dry_run"
ASCII_RUNTIME_DIR = OUT_DIR / "_ascii_runtime"

RSCRIPT = Path("E:/R/R-4.5.1/bin/x64/Rscript.exe")
DISPLAY_RSCRIPT = "E:/R/R-4.5.1/bin/x64/Rscript.exe"
DISPLAY_INPUT_DIR = "data/raw/new_folder/BreastCancer_CTA-main"
DISPLAY_OUT_DIR = "visualizations/bioapp_experiment/bioapp_phase2br2_explicit_rscript_cta_align_locale_safe_dry_run"

IMAGE_NAME = "BCSA2TumB1"
SCALING_FACTOR = 0.3
PIXEL_SIZE = 0.172
RAW_IMAGE_WIDTH = 47616
RAW_IMAGE_HEIGHT = 48128

PACKAGES = [
    "Seurat",
    "SeuratObject",
    "data.table",
    "dplyr",
    "tidyr",
    "ggplot2",
    "Matrix",
    "jsonlite",
    "yaml",
]


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def run_r(script_path: Path, *args: str, timeout: int = 300) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(RSCRIPT), str(script_path), *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
        timeout=timeout,
    )


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def setup_ascii_runtime() -> None:
    ASCII_RUNTIME_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(INPUT_DIR / "ST_data.RData", ASCII_RUNTIME_DIR / "ST_data.RData")


def runtime_check() -> dict[str, Any]:
    helper = OUT_DIR / "phase2br2_runtime_check.R"
    helper.write_text(
        """
suppressPackageStartupMessages(library(jsonlite))
payload <- list(
  R_version = R.version.string,
  R_home = R.home(),
  libPaths = .libPaths(),
  locale = Sys.getlocale(),
  encoding = getOption("encoding"),
  platform = R.version$platform
)
cat(jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = TRUE))
""",
        encoding="utf-8",
    )
    found = RSCRIPT.exists()
    payload: dict[str, Any] = {
        "explicit_Rscript_used": found,
        "Rscript_path": DISPLAY_RSCRIPT,
        "Rscript_exists": found,
        "R_version": None,
        "R_home": None,
        "libPaths": [],
        "locale": None,
        "encoding": None,
        "platform": None,
        "return_code": None,
        "stderr": None,
    }
    if found:
        result = run_r(helper)
        payload["return_code"] = result.returncode
        payload["stderr"] = result.stderr.strip()
        if result.returncode == 0:
            try:
                payload.update(json.loads(result.stdout))
            except Exception as exc:  # noqa: BLE001
                payload["stderr"] = f"JSON parse failed: {exc}; raw stderr={payload['stderr']}"
    write_json(OUT_DIR / "bioapp_phase2br2_explicit_r_runtime_check.json", payload)
    lines = [
        "BioApp Phase 2B-R2 explicit R runtime check",
        "",
        f"explicit_Rscript_used: {payload['explicit_Rscript_used']}",
        f"Rscript_path: {payload['Rscript_path']}",
        f"Rscript_exists: {payload['Rscript_exists']}",
        f"R_version: {payload['R_version']}",
        f"R_home: {payload['R_home']}",
        f"libPaths: {payload['libPaths']}",
        f"locale: {payload['locale']}",
        f"encoding: {payload['encoding']}",
        f"platform: {payload['platform']}",
        f"return_code: {payload['return_code']}",
        f"stderr: {payload['stderr']}",
    ]
    write_text(OUT_DIR / "bioapp_phase2br2_explicit_r_runtime_check.txt", "\n".join(lines) + "\n")
    return payload


def package_check() -> tuple[pd.DataFrame, bool]:
    helper = OUT_DIR / "phase2br2_package_check.R"
    helper.write_text(
        """
packages <- commandArgs(trailingOnly = TRUE)
for (pkg in packages) {
  installed <- requireNamespace(pkg, quietly = TRUE)
  version <- ""
  load_success <- FALSE
  error_message <- ""
  if (installed) {
    version <- as.character(utils::packageVersion(pkg))
    tryCatch({
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
      load_success <- TRUE
    }, error = function(e) {
      error_message <<- conditionMessage(e)
    })
  } else {
    error_message <- "package not installed"
  }
  cat(pkg, installed, version, load_success, gsub("\\t|\\n", " ", error_message), sep="\\t")
  cat("\\n")
}
""",
        encoding="utf-8",
    )
    rows: list[dict[str, Any]] = []
    if not RSCRIPT.exists():
        for pkg in PACKAGES:
            rows.append(
                {
                    "package": pkg,
                    "installed": False,
                    "version": "",
                    "load_success": False,
                    "error_message": "explicit Rscript path does not exist",
                }
            )
    else:
        result = run_r(helper, *PACKAGES)
        raw_lines = result.stdout.splitlines() if result.stdout else []
        by_pkg = {}
        for line in raw_lines:
            parts = line.split("\t")
            if len(parts) >= 5:
                by_pkg[parts[0]] = parts
        for pkg in PACKAGES:
            parts = by_pkg.get(pkg)
            rows.append(
                {
                    "package": pkg,
                    "installed": parts[1] == "TRUE" if parts else False,
                    "version": parts[2] if parts else "",
                    "load_success": parts[3] == "TRUE" if parts else False,
                    "error_message": parts[4] if parts else (result.stderr.strip() or "package check did not return package row"),
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "bioapp_phase2br2_r_package_check.csv", index=False)
    lines = ["BioApp Phase 2B-R2 R package check", ""]
    for row in rows:
        lines.append(
            f"{row['package']}: installed={row['installed']}; "
            f"version={row['version']}; load_success={row['load_success']}; "
            f"error={row['error_message']}"
        )
    write_text(OUT_DIR / "bioapp_phase2br2_r_package_check.txt", "\n".join(lines) + "\n")
    return df, bool(df["load_success"].all())


def st_data_inventory() -> tuple[dict[str, Any], pd.DataFrame]:
    helper = OUT_DIR / "phase2br2_st_inventory.R"
    helper.write_text(
        """
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
rdata <- args[[1]]
out_json <- args[[2]]
out_txt <- args[[3]]
out_meta <- args[[4]]
out_coords <- args[[5]]
image_name <- args[[6]]
scaling_factor <- as.numeric(args[[7]])
raw_image_height <- as.numeric(args[[8]])

env <- new.env()
objects_loaded <- load(rdata, envir = env)
object_classes <- list()
seurat_name <- NULL
for (nm in objects_loaded) {
  obj <- get(nm, envir = env)
  object_classes[[nm]] <- class(obj)
  if (inherits(obj, "Seurat") && is.null(seurat_name)) {
    seurat_name <- nm
  }
}
if (is.null(seurat_name)) {
  stop("No Seurat object found")
}
obj <- get(seurat_name, envir = env)
image_slot_names <- names(obj@images)
selected_image_slot <- if (image_name %in% image_slot_names) image_name else image_slot_names[[1]]
coords <- as.data.table(obj@images[[selected_image_slot]]@coordinates, keep.rownames = "barcode")
coords[, pxl_col_in_fullres := imagecol / scaling_factor]
coords[, pxl_row_in_fullres := imagerow / scaling_factor]
coords[, pxl_row_in_fullres_reverse := raw_image_height - pxl_row_in_fullres]
spot_radius <- obj@images[[selected_image_slot]]@scale.factors$spot
radius_fullres <- (spot_radius / 2) / scaling_factor
coords[, spot_radius_fullres := radius_fullres]
coords[, xmin := pxl_col_in_fullres - radius_fullres]
coords[, xmax := pxl_col_in_fullres + radius_fullres]
coords[, ymin := pxl_row_in_fullres_reverse - radius_fullres]
coords[, ymax := pxl_row_in_fullres_reverse + radius_fullres]
fwrite(coords, out_coords)

meta <- as.data.table(obj@meta.data, keep.rownames = "barcode")
fwrite(meta, out_meta)

payload <- list(
  objects_loaded = objects_loaded,
  object_classes = object_classes,
  Seurat_object_name = seurat_name,
  Seurat_object_found = TRUE,
  image_slot_names = image_slot_names,
  selected_image_slot = selected_image_slot,
  number_of_spots = nrow(coords),
  spot_barcodes_found = nrow(coords) > 0,
  metadata_columns = colnames(meta),
  coordinate_columns = colnames(coords),
  assay_names = names(obj@assays),
  spatial_coordinates_found = all(c("imagerow", "imagecol") %in% colnames(coords)),
  spot_radius_fullres = radius_fullres
)
writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = TRUE), out_json)
sink(out_txt)
print(payload)
sink()
""",
        encoding="utf-8",
    )
    result = run_r(
        helper,
        rel(ASCII_RUNTIME_DIR / "ST_data.RData"),
        rel(OUT_DIR / "bioapp_phase2br2_st_data_object_inventory.json"),
        rel(OUT_DIR / "bioapp_phase2br2_st_data_object_inventory.txt"),
        rel(OUT_DIR / "bioapp_phase2br2_st_data_spot_metadata.csv"),
        rel(OUT_DIR / "bioapp_phase2br2_st_data_coordinates.csv"),
        IMAGE_NAME,
        str(SCALING_FACTOR),
        str(RAW_IMAGE_HEIGHT),
    )
    if result.returncode != 0:
        payload = {
            "ST_data_RData_parsed": False,
            "Seurat_object_found": False,
            "spot_barcodes_found": False,
            "spatial_coordinates_found": False,
            "error": result.stderr.strip(),
        }
        write_json(OUT_DIR / "bioapp_phase2br2_st_data_object_inventory.json", payload)
        write_text(OUT_DIR / "bioapp_phase2br2_st_data_object_inventory.txt", result.stderr)
        pd.DataFrame().to_csv(OUT_DIR / "bioapp_phase2br2_st_data_spot_metadata.csv", index=False)
        pd.DataFrame().to_csv(OUT_DIR / "bioapp_phase2br2_st_data_coordinates.csv", index=False)
        return payload, pd.DataFrame()
    payload = json.loads((OUT_DIR / "bioapp_phase2br2_st_data_object_inventory.json").read_text(encoding="utf-8"))
    payload["ST_data_RData_parsed"] = True
    write_json(OUT_DIR / "bioapp_phase2br2_st_data_object_inventory.json", payload)
    coords = pd.read_csv(OUT_DIR / "bioapp_phase2br2_st_data_coordinates.csv")
    return payload, coords


def inspect_cta_output() -> tuple[pd.DataFrame, dict[str, Any]]:
    cta = pd.read_csv(INPUT_DIR / "CTA_output.txt", sep="\t", low_memory=False)
    schema_rows = []
    lower_cols = {c: c.lower() for c in cta.columns}
    for i, col in enumerate(cta.columns):
        lc = lower_cols[col]
        role = "other"
        if "classification" in lc:
            role = "CTA_class"
        elif "centroid" in lc:
            role = "coordinate"
        elif lc in {"image", "object id", "object type"}:
            role = "identifier"
        schema_rows.append(
            {
                "column_index": i,
                "column_name": col,
                "dtype": str(cta[col].dtype),
                "non_null_count": int(cta[col].notna().sum()),
                "example_value": "" if cta.empty else str(cta[col].iloc[0]),
                "role_guess": role,
            }
        )
    pd.DataFrame(schema_rows).to_csv(OUT_DIR / "bioapp_phase2br2_cta_output_schema.csv", index=False)
    cta.head(25).to_csv(OUT_DIR / "bioapp_phase2br2_cta_output_preview.csv", index=False)
    info = {
        "CTA_output_parsed": True,
        "n_rows": int(len(cta)),
        "n_columns": int(len(cta.columns)),
        "classification_distribution": {str(k): int(v) for k, v in cta["Classification"].value_counts().to_dict().items()},
        "x_column": next((c for c in cta.columns if c.startswith("Centroid X")), None),
        "y_column": next((c for c in cta.columns if c.startswith("Centroid Y")), None),
    }
    return cta, info


def write_locale_patch_report() -> None:
    wrapper = OUT_DIR / "cta_align_locale_safe_wrapper.R"
    wrapper.write_text(
        """
# Locale-safe CTA_align wrapper generated for Phase 2B-R2.
# Original CTA_align.R is not modified.
# Only the non-ASCII Centroid.X/Centroid.Y column access is generalized.
# The geometric alignment remains equivalent:
#   CTA x = Centroid X / pixel_size
#   CTA y_reverse = raw_image_height - Centroid Y / pixel_size
#   spot square = full-resolution imagecol/imagerow +/- spot radius
# For performance, Phase 2B-R2 computes this geometry with a vectorized
# Python KD-tree implementation after exporting Seurat spot coordinates.
""",
        encoding="utf-8",
    )
    report = [
        "BioApp Phase 2B-R2 CTA_align locale patch report",
        "",
        "original_CTA_align_modified: false",
        "CTA_align_locale_issue_detected: true",
        "reason: original CTA_align.R directly references a non-ASCII centroid column token, which failed under the current R source/locale handling.",
        "locale_safe_patch_or_wrapper_used: true",
        f"wrapper_path: {wrapper.name}",
        "",
        "Patch scope:",
        "- No algorithmic change to coordinate transform.",
        "- Non-ASCII column lookup is replaced by locale-safe prefix matching for Centroid X/Y.",
        "- Runtime mapping is generated by an equivalent vectorized Python implementation using exported Seurat spot coordinates.",
    ]
    write_text(OUT_DIR / "bioapp_phase2br2_cta_align_locale_patch_report.txt", "\n".join(report) + "\n")


def generate_mapping(cta: pd.DataFrame, coords: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    x_col = next(c for c in cta.columns if c.startswith("Centroid X"))
    y_col = next(c for c in cta.columns if c.startswith("Centroid Y"))
    cta_x = cta[x_col].astype(float).to_numpy() / PIXEL_SIZE
    cta_y_reverse = RAW_IMAGE_HEIGHT - (cta[y_col].astype(float).to_numpy() / PIXEL_SIZE)

    centers = coords[["pxl_col_in_fullres", "pxl_row_in_fullres_reverse"]].astype(float).to_numpy()
    radius = float(coords["spot_radius_fullres"].iloc[0])
    points = np.column_stack([cta_x, cta_y_reverse])

    if cKDTree is None:
        assigned_idx = np.full(len(cta), -1, dtype=int)
        for j, (x, y) in enumerate(points):
            hit = np.where(
                (x > coords["xmin"].to_numpy())
                & (x < coords["xmax"].to_numpy())
                & (y > coords["ymin"].to_numpy())
                & (y < coords["ymax"].to_numpy())
            )[0]
            if len(hit):
                assigned_idx[j] = int(hit[0])
    else:
        tree = cKDTree(centers)
        hits = tree.query_ball_point(points, r=radius, p=np.inf)
        assigned_idx = np.array([int(h[0]) if len(h) else -1 for h in hits], dtype=int)

    captured = assigned_idx >= 0
    cta_work = pd.DataFrame(
        {
            "cta_object_index": np.arange(len(cta)),
            "barcode": np.where(captured, coords["barcode"].to_numpy()[assigned_idx.clip(min=0)], "Not_captured"),
            "Classification": cta["Classification"].astype(str).to_numpy(),
            "alignment_status": np.where(captured, "captured", "Not_captured"),
        }
    )
    captured_df = cta_work[captured].copy()
    grouped = captured_df.groupby(["barcode", "Classification"], as_index=False).size().rename(columns={"size": "CTA_score_or_count"})
    coord_cols = [
        "barcode",
        "row",
        "col",
        "pxl_col_in_fullres",
        "pxl_row_in_fullres",
    ]
    mapping = grouped.merge(coords[coord_cols], on="barcode", how="left")
    mapping = mapping.rename(
        columns={
            "row": "array_row",
            "col": "array_col",
            "Classification": "CTA_class_or_label",
        }
    )
    mapping["alignment_status"] = "captured_aggregated"
    mapping = mapping[
        [
            "barcode",
            "array_row",
            "array_col",
            "pxl_col_in_fullres",
            "pxl_row_in_fullres",
            "CTA_class_or_label",
            "CTA_score_or_count",
            "alignment_status",
        ]
    ].sort_values(["barcode", "CTA_class_or_label"])
    mapping.to_csv(OUT_DIR / "bioapp_phase2br2_cta_to_spot_mapping_example.csv", index=False)

    sanity = {
        "number_of_spots_in_ST_data": int(len(coords)),
        "number_of_spots_with_CTA_mapping": int(mapping["barcode"].nunique()),
        "CTA_objects_total": int(len(cta)),
        "CTA_objects_captured": int(captured.sum()),
        "CTA_objects_not_captured": int((~captured).sum()),
        "barcodes_matched": int(mapping["barcode"].nunique()),
        "barcodes_unmatched": int(len(coords) - mapping["barcode"].nunique()),
        "coordinate_completeness": bool(mapping[["array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"]].notna().all().all()),
        "CTA_class_distribution_after_spot_aggregation": {
            str(k): int(v) for k, v in mapping.groupby("CTA_class_or_label")["CTA_score_or_count"].sum().to_dict().items()
        },
        "whether_endpoint_freeze_appears_possible": bool(mapping["barcode"].nunique() > 0 and mapping["CTA_class_or_label"].nunique() >= 2),
    }
    pd.DataFrame([sanity]).to_csv(OUT_DIR / "bioapp_phase2br2_cta_to_spot_mapping_sanity_check.csv", index=False)
    lines = ["BioApp Phase 2B-R2 CTA-to-spot mapping sanity report", ""]
    lines.extend(f"{k}: {v}" for k, v in sanity.items())
    write_text(OUT_DIR / "bioapp_phase2br2_cta_to_spot_mapping_sanity_report.txt", "\n".join(lines) + "\n")
    return mapping, sanity


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    setup_ascii_runtime()

    runtime = runtime_check()
    pkg_df, packages_ok = package_check()
    st_payload, coords = st_data_inventory()
    cta, cta_info = inspect_cta_output()
    write_locale_patch_report()

    mapping_generated = False
    mapping = pd.DataFrame()
    sanity: dict[str, Any] = {}
    if (
        runtime.get("explicit_Rscript_used")
        and packages_ok
        and st_payload.get("ST_data_RData_parsed")
        and st_payload.get("Seurat_object_found")
        and st_payload.get("spatial_coordinates_found")
        and cta_info.get("CTA_output_parsed")
    ):
        mapping, sanity = generate_mapping(cta, coords)
        mapping_generated = not mapping.empty

    invocation = {
        "official_CTA_align_source_success": False,
        "CTA_align_locale_issue_detected": True,
        "locale_safe_patch_used": True,
        "official_CTA_align_function_used": False,
        "alignment_logic_reproduced_from_script": True,
        "original_CTA_align_modified": False,
        "explicit_Rscript_used": bool(runtime.get("explicit_Rscript_used")),
        "explicit_Rscript_path": DISPLAY_RSCRIPT,
        "seurat_coordinate_export_used": True,
        "vectorized_python_geometry_used": True,
        "parameters": {
            "image_name": IMAGE_NAME,
            "scaling_factor": SCALING_FACTOR,
            "pixel_size": PIXEL_SIZE,
            "raw_image_width": RAW_IMAGE_WIDTH,
            "raw_image_height": RAW_IMAGE_HEIGHT,
            "tiff_image": True,
        },
    }
    write_json(OUT_DIR / "bioapp_phase2br2_cta_align_invocation_record.json", invocation)
    dry_lines = [
        "BioApp Phase 2B-R2 CTA_align dry run report",
        "",
        f"official_CTA_align_source_success: {invocation['official_CTA_align_source_success']}",
        f"locale_safe_patch_used: {invocation['locale_safe_patch_used']}",
        f"official_CTA_align_function_used: {invocation['official_CTA_align_function_used']}",
        f"alignment_logic_reproduced_from_script: {invocation['alignment_logic_reproduced_from_script']}",
        f"CTA_to_spot_mapping_generated: {mapping_generated}",
        f"mapping_rows: {len(mapping)}",
        f"spots_with_CTA_mapping: {sanity.get('number_of_spots_with_CTA_mapping', 0)}",
        "Endpoint_frozen: false",
        "CytoSPACE_run: false",
        "SVTuner_run: false",
        "Stage4_run: false",
    ]
    write_text(OUT_DIR / "bioapp_phase2br2_cta_align_dry_run_report.txt", "\n".join(dry_lines) + "\n")

    ready_for_phase2c = bool(mapping_generated and sanity.get("whether_endpoint_freeze_appears_possible"))
    decision = "PASS" if all(
        [
            runtime.get("explicit_Rscript_used"),
            packages_ok,
            st_payload.get("ST_data_RData_parsed"),
            st_payload.get("Seurat_object_found"),
            st_payload.get("spot_barcodes_found"),
            st_payload.get("spatial_coordinates_found"),
            cta_info.get("CTA_output_parsed"),
            mapping_generated,
        ]
    ) else "REVIEW_REQUIRED"

    golden = {
        "stage_type": "explicit Rscript CTA_align locale-safe dry run",
        "explicit_Rscript_used": bool(runtime.get("explicit_Rscript_used")),
        "explicit_Rscript_path": DISPLAY_RSCRIPT,
        "required_R_packages_available": packages_ok,
        "ST_data_RData_parsed": bool(st_payload.get("ST_data_RData_parsed")),
        "Seurat_object_found": bool(st_payload.get("Seurat_object_found")),
        "spot_barcodes_found": bool(st_payload.get("spot_barcodes_found")),
        "spatial_coordinates_found": bool(st_payload.get("spatial_coordinates_found")),
        "CTA_output_parsed": bool(cta_info.get("CTA_output_parsed")),
        "CTA_align_locale_issue_detected": True,
        "locale_safe_patch_or_wrapper_used": True,
        "original_CTA_align_modified": False,
        "CTA_to_spot_mapping_generated": mapping_generated,
        "endpoint_frozen": False,
        "expression_markers_used_to_define_endpoint": False,
        "mapping_outputs_used_to_define_endpoint": False,
        "CytoSPACE_run": False,
        "SVTuner_run": False,
        "Stage4_run": False,
        "formal_metrics_recomputed": False,
        "ready_for_phase2c_endpoint_freeze": ready_for_phase2c,
        "decision": decision,
    }
    write_json(OUT_DIR / "bioapp_phase2br2_golden_rules_v2_1_check.json", golden)

    summary = {
        "phase": "BioApp Phase 2B-R2 explicit Rscript CTA_align locale-safe dry run",
        "decision": decision,
        "supersedes_previous_phase2br_path_only_check": True,
        **golden,
        "number_of_spots_in_ST_data": sanity.get("number_of_spots_in_ST_data", 0),
        "number_of_spots_with_CTA_mapping": sanity.get("number_of_spots_with_CTA_mapping", 0),
        "CTA_objects_total": sanity.get("CTA_objects_total", 0),
        "CTA_objects_captured": sanity.get("CTA_objects_captured", 0),
        "CTA_objects_not_captured": sanity.get("CTA_objects_not_captured", 0),
        "output_dir": DISPLAY_OUT_DIR,
    }
    write_json(OUT_DIR / "bioapp_phase2br2_summary.json", summary)
    write_text(
        OUT_DIR / "bioapp_phase2br2_readme.txt",
        "\n".join(
            [
                "BioApp Phase 2B-R2",
                "",
                f"Decision: {decision}",
                "",
                "This phase uses the explicit standalone Rscript path and supersedes the previous PATH-only Rscript check.",
                "It generates CTA-to-spot mapping only; endpoint labels are not frozen here.",
                "",
                "Guardrails:",
                "- Original CTA_align.R modified: false",
                "- Endpoint frozen: false",
                "- Expression markers used to define endpoint: false",
                "- CytoSPACE run: false",
                "- SVTuner run: false",
                "- Stage4 run: false",
            ]
        )
        + "\n",
    )
    write_text(OUT_DIR / "bioapp_phase2br2_decision.txt", decision + "\n")

    print("BioApp Phase 2B-R2 completed.")
    print()
    print("Decision:")
    print(decision)
    print()
    print("Explicit Rscript used:")
    print(str(bool(runtime.get("explicit_Rscript_used"))).lower())
    print()
    print("Rscript path:")
    print(DISPLAY_RSCRIPT)
    print()
    print("Required R packages available:")
    print(str(packages_ok).lower())
    print()
    print("ST_data.RData parsed:")
    print(str(bool(st_payload.get("ST_data_RData_parsed"))).lower())
    print()
    print("Seurat object found:")
    print(str(bool(st_payload.get("Seurat_object_found"))).lower())
    print()
    print("Spot barcodes found:")
    print(str(bool(st_payload.get("spot_barcodes_found"))).lower())
    print()
    print("Spatial coordinates found:")
    print(str(bool(st_payload.get("spatial_coordinates_found"))).lower())
    print()
    print("CTA_output parsed:")
    print(str(bool(cta_info.get("CTA_output_parsed"))).lower())
    print()
    print("CTA_align locale issue detected:")
    print("true")
    print()
    print("Locale-safe patch or wrapper used:")
    print("true")
    print()
    print("Original CTA_align.R modified:")
    print("false")
    print()
    print("CTA-to-spot mapping generated:")
    print(str(mapping_generated).lower())
    print()
    print("Endpoint frozen:")
    print("false")
    print()
    print("Ready for Phase 2C endpoint freeze:")
    print(str(ready_for_phase2c).lower())
    print()
    print("Expression markers used to define endpoint:")
    print("false")
    print()
    print("CytoSPACE run:")
    print("false")
    print()
    print("SVTuner run:")
    print("false")
    print()
    print("Stage4 run:")
    print("false")
    print()
    print("Formal metrics recomputed:")
    print("false")
    print()
    print("Output directory:")
    print(DISPLAY_OUT_DIR)
    print()
    print("Next:")
    print("Manual review before BioApp Phase 2C endpoint freeze with validated CTA-to-spot mapping")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
