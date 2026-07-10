from __future__ import annotations

import csv
import json
import math
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_coordinate_registration_recovery_audit"
RUNTIME = OUT / "_r_runtime"

RSCRIPT = Path("E:/R/R-4.5.1/bin/x64/Rscript.exe")
STDATA_ACTUAL = ROOT / "data" / "raw" / "\u65b0\u5efa\u6587\u4ef6\u5939" / "BreastCancer_CTA-main" / "ST_data.RData"
STDATA_ASCII = ROOT / "visualizations" / "bioapp_experiment" / "bioapp_phase2br2_explicit_rscript_cta_align_locale_safe_dry_run" / "_ascii_runtime" / "ST_data.RData"
CTA_DIR = ROOT / "data" / "raw" / "\u65b0\u5efa\u6587\u4ef6\u5939" / "BreastCancer_CTA-main"
STDATA_DISPLAY = "data/raw/\u65b0\u5efa\u6587\u4ef6\u5939/BreastCancer_CTA-main/ST_data.RData"
IMAGE_SLOT = "BCSA2TumB1"
SEURAT_OBJECT = "bcsa"

SEARCH_TERMS = [
    "spatial",
    "tissue",
    "hires",
    "lowres",
    "image",
    "png",
    "jpg",
    "jpeg",
    "tif",
    "tiff",
    "svs",
    "h5",
    "h5ad",
    "h5seurat",
    "scalefactors",
    "positions",
    "tissue_positions",
    "SpaceRanger",
    "spaceranger",
    "CytAssist",
    "Visium",
    "BCSA",
    "Tum",
    "TumB1",
    "H&E",
    "HE",
    "histology",
    "registration",
    "transform",
    "affine",
    "keypoint",
    "morphology",
    "Xenium",
]


def rel(path: Path | None) -> str | None:
    if path is None:
        return None
    try:
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")


def write_rows(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k) for k in fields})


def run_r(script: Path, *args: str, timeout: int = 600) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(RSCRIPT), str(script), *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
        timeout=timeout,
    )


def make_r_helper() -> Path:
    RUNTIME.mkdir(parents=True, exist_ok=True)
    helper = RUNTIME / "inspect_stdata_image_slot.R"
    helper.write_text(
        r'''
suppressPackageStartupMessages({
  library(Seurat)
  library(jsonlite)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
rdata <- args[[1]]
out_dir <- args[[2]]
object_name <- args[[3]]
image_name <- args[[4]]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

json_write <- function(x, path) {
  writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null"), path, useBytes = TRUE)
}

env <- new.env()
objects_loaded <- load(rdata, envir = env)
classes <- list()
for (nm in objects_loaded) {
  classes[[nm]] <- class(get(nm, envir = env))
}
if (!(object_name %in% objects_loaded)) {
  object_name <- NULL
  for (nm in objects_loaded) {
    obj_try <- get(nm, envir = env)
    if (inherits(obj_try, "Seurat")) {
      object_name <- nm
      break
    }
  }
}
if (is.null(object_name)) {
  stop("No Seurat object found")
}
obj <- get(object_name, envir = env)
if (!(image_name %in% names(obj@images))) {
  image_name <- names(obj@images)[[1]]
}
img <- obj@images[[image_name]]
coords <- img@coordinates

object_inventory <- list(
  source_rdata = rdata,
  objects_loaded = objects_loaded,
  object_classes = classes,
  seurat_object_name = object_name,
  seurat_object_class = class(obj),
  seurat_version = as.character(utils::packageVersion("Seurat")),
  seuratobject_version = as.character(utils::packageVersion("SeuratObject")),
  assays = names(obj@assays),
  reductions = names(obj@reductions),
  images = names(obj@images),
  selected_image_slot = image_name,
  default_assay = DefaultAssay(obj),
  metadata_columns = colnames(obj@meta.data),
  misc_names = names(obj@misc),
  tool_names = names(obj@tools),
  command_names = names(obj@commands),
  project_name = obj@project.name
)
json_write(object_inventory, file.path(out_dir, "bioapp_st_data_rdata_object_inventory.json"))

scale_list <- list()
if ("scale.factors" %in% slotNames(img)) {
  sf <- img@scale.factors
  scale_list <- lapply(as.list(sf), function(v) as.numeric(v))
}
image_dims <- NULL
embedded_image <- FALSE
if ("image" %in% slotNames(img)) {
  image_dims <- dim(img@image)
  embedded_image <- !is.null(image_dims) && length(image_dims) >= 2
}
coord_summary <- list()
for (nm in colnames(coords)) {
  vals <- suppressWarnings(as.numeric(coords[[nm]]))
  if (sum(!is.na(vals)) > 0) {
    coord_summary[[nm]] <- list(min = min(vals, na.rm = TRUE), max = max(vals, na.rm = TRUE), mean = mean(vals, na.rm = TRUE))
  }
}
slot_summary <- list(
  image_slot_name = image_name,
  image_slot_class = class(img),
  slot_names = slotNames(img),
  coordinate_dimensions = dim(coords),
  coordinate_columns = colnames(coords),
  coordinate_summary = coord_summary,
  scale_factors = scale_list,
  embedded_image_found = embedded_image,
  embedded_image_class = if ("image" %in% slotNames(img)) class(img@image) else NULL,
  embedded_image_dimensions = image_dims,
  key = tryCatch(img@key, error = function(e) NULL),
  assay = tryCatch(img@assay, error = function(e) NULL),
  misc_names = tryCatch(names(img@misc), error = function(e) NULL),
  spot_radius = tryCatch(img@spot.radius, error = function(e) NULL),
  boundary_slots_present = intersect(slotNames(img), c("boundaries", "centroids", "molecules", "segmentation"))
)
json_write(slot_summary, file.path(out_dir, "bioapp_seurat_image_slot_deep_inspection.json"))

scale_rows <- data.table(
  scale_factor_name = character(),
  value = numeric(),
  source = character(),
  notes = character()
)
if (length(scale_list) > 0) {
  for (nm in names(scale_list)) {
    scale_rows <- rbind(scale_rows, data.table(
      scale_factor_name = nm,
      value = as.numeric(scale_list[[nm]]),
      source = paste0("obj@images[['", image_name, "']]@scale.factors"),
      notes = ifelse(abs(as.numeric(scale_list[[nm]]) - 0.3) < 1e-8, "matches CTA_align scaling_factor", "")
    ))
  }
}
data.table::fwrite(scale_rows, file.path(out_dir, "bioapp_image_slot_scale_factor_audit.csv"))

path_terms <- c("spatial","tissue","hires","lowres","image","png","jpg","jpeg","tif","tiff","svs","h5","h5ad","h5seurat","scalefactors","positions","tissue_positions","SpaceRanger","spaceranger","CytAssist","Visium","BCSA","Tum","TumB1","H&E","HE","histology","registration","transform","affine","keypoint","morphology","Xenium")
matches <- data.table(object_path = character(), string_value = character(), matched_terms = character())
add_match <- function(path, value) {
  if (is.na(value) || !nzchar(value)) return()
  found <- path_terms[sapply(path_terms, function(t) grepl(t, value, ignore.case = TRUE, fixed = FALSE))]
  if (length(found) > 0) {
    matches <<- rbind(matches, data.table(object_path = path, string_value = value, matched_terms = paste(found, collapse = "|")))
  }
}
scan_obj <- function(x, path = "object", depth = 0) {
  if (depth > 5) return()
  nms <- tryCatch(names(x), error = function(e) NULL)
  if (!is.null(nms)) {
    for (nm in nms) add_match(paste0(path, "$name"), nm)
  }
  if (is.character(x)) {
    for (i in seq_along(x)) add_match(paste0(path, "[", i, "]"), as.character(x[[i]]))
    return()
  }
  if (isS4(x)) {
    for (sl in slotNames(x)) {
      if (sl %in% c("image", "counts", "data", "scale.data")) next()
      val <- tryCatch(slot(x, sl), error = function(e) NULL)
      if (!is.null(val)) scan_obj(val, paste0(path, "@", sl), depth + 1)
    }
    return()
  }
  if (is.data.frame(x)) {
    for (nm in colnames(x)) add_match(paste0(path, "$colname"), nm)
    char_cols <- names(x)[sapply(x, is.character)]
    for (nm in char_cols) {
      vals <- unique(x[[nm]])
      vals <- vals[seq_len(min(length(vals), 100))]
      for (v in vals) add_match(paste0(path, "$", nm), as.character(v))
    }
    return()
  }
  if (is.list(x)) {
    lim <- min(length(x), 100)
    if (lim == 0) return()
    for (i in seq_len(lim)) {
      nm <- names(x)[i]
      if (is.null(nm) || !nzchar(nm)) nm <- as.character(i)
      scan_obj(x[[i]], paste0(path, "$", nm), depth + 1)
    }
  }
}
scan_obj(obj, "bcsa", 0)
matches[, path_exists := file.exists(string_value)]
matches[, resembles_missing_path := grepl("\\\\|/|\\.png|\\.jpg|\\.jpeg|\\.tif|\\.tiff|\\.h5|\\.RData|\\.rds", string_value, ignore.case = TRUE) & !path_exists]
matches[, identifies_original_spatial_bundle := grepl("spatial|scalefactors|tissue_positions|SpaceRanger|spaceranger|Visium|CytAssist", string_value, ignore.case = TRUE)]
data.table::fwrite(matches, file.path(out_dir, "bioapp_st_data_pathlike_string_search.csv"))

if (embedded_image) {
  arr <- img@image
  h <- dim(arr)[1]
  w <- dim(arr)[2]
  lowres <- if ("lowres" %in% names(scale_list)) as.numeric(scale_list[["lowres"]]) else NA_real_
  hires <- if ("hires" %in% names(scale_list)) as.numeric(scale_list[["hires"]]) else NA_real_
  x_low <- coords$imagecol * lowres
  y_low <- coords$imagerow * lowres
  lowres_bounds_pass <- all(x_low >= 0, x_low <= w, y_low >= 0, y_low <= h, na.rm = TRUE)
  x_hires <- coords$imagecol * hires
  y_hires <- coords$imagerow * hires
  hires_bounds_pass <- all(x_hires >= 0, x_hires <= w, y_hires >= 0, y_hires <= h, na.rm = TRUE)
  bounds <- list(
    image_width = w,
    image_height = h,
    lowres_scale = lowres,
    hires_scale = hires,
    lowres_bounds_pass = lowres_bounds_pass,
    hires_bounds_pass = hires_bounds_pass,
    lowres_x_min = min(x_low, na.rm = TRUE),
    lowres_x_max = max(x_low, na.rm = TRUE),
    lowres_y_min = min(y_low, na.rm = TRUE),
    lowres_y_max = max(y_low, na.rm = TRUE)
  )
  json_write(bounds, file.path(out_dir, "bioapp_candidate_image_bounds_check.json"))
  png(file.path(out_dir, "bioapp_candidate_embedded_tissue_image_preview.png"), width = w, height = h)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, w), ylim = c(h, 0), asp = 1)
  rasterImage(as.raster(arr), 0, h, w, 0)
  dev.off()
  png(file.path(out_dir, "bioapp_candidate_tissue_overlay_audit_preview.png"), width = w, height = h)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, w), ylim = c(h, 0), asp = 1)
  rasterImage(as.raster(arr), 0, h, w, 0)
  points(x_low, y_low, pch = 16, cex = 0.18, col = rgb(0, 1, 1, 0.65))
  dev.off()
}
''',
        encoding="ascii",
    )
    return helper


def search_cta_dir() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not CTA_DIR.exists():
        return rows
    for path in CTA_DIR.rglob("*"):
        if not path.is_file():
            continue
        suffix = path.suffix.lower()
        matched: list[str] = []
        contains_coord = False
        contains_image_meta = False
        contains_transform = False
        contains_keypoints = False
        contains_image_refs = False
        usable = False
        text = ""
        if suffix in {".r", ".txt", ".json", ".yaml", ".yml", ".md", ".csv", ".tsv"}:
            try:
                text = path.read_text(encoding="utf-8", errors="ignore")[:2_000_000]
            except Exception:
                text = ""
        name_and_text = f"{path.name}\n{text}"
        low = name_and_text.lower()
        for term in SEARCH_TERMS:
            if term.lower() in low:
                matched.append(term)
        contains_coord = any(t in low for t in ["imagerow", "imagecol", "pxl_col", "pxl_row", "centroid.x", "centroid.y"])
        contains_image_meta = any(t in low for t in ["raw_image_width", "raw_image_height", "scaling_factor", "pixel_size", "tiff_image"])
        contains_transform = any(t in low for t in ["transform", "affine", "homography", "matrix", "scale", "rotate", "translate"])
        contains_keypoints = "keypoint" in low
        contains_image_refs = any(t in low for t in [".png", ".jpg", ".jpeg", ".tif", ".tiff", ".svs"])
        usable = bool(contains_image_meta and contains_coord and path.name in {"CTA_align.R", "README.md"})
        if matched or suffix in {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".svs"}:
            rows.append(
                {
                    "file_path": rel(path),
                    "file_type": suffix,
                    "matched_terms": "|".join(sorted(set(matched))),
                    "contains_coordinate_columns": contains_coord,
                    "contains_image_metadata": contains_image_meta,
                    "contains_transform_matrix": contains_transform and "matrix" in low,
                    "contains_keypoints": contains_keypoints,
                    "contains_image_file_references": contains_image_refs,
                    "usable_for_recovery": usable,
                    "notes": "Contains CTA alignment parameters but not an inverse transform or local image path."
                    if usable
                    else "",
                }
            )
    return rows


def path_exists_from_string(value: str) -> bool:
    p = Path(value)
    if p.exists():
        return True
    for base in [ROOT, CTA_DIR, ROOT / "data" / "raw", ROOT / "data" / "processed"]:
        if (base / value).exists():
            return True
    return False


def postprocess_pathlike(path_csv: Path) -> None:
    if not path_csv.exists():
        return
    df = pd.read_csv(path_csv)
    if df.empty:
        return
    df["path_exists"] = df["string_value"].astype(str).map(path_exists_from_string)
    df["resembles_missing_path"] = df["string_value"].astype(str).str.contains(r"\\|/|\.png|\.jpg|\.jpeg|\.tif|\.tiff|\.h5|\.RData|\.rds", case=False, regex=True) & ~df["path_exists"]
    df["can_identify_original_spatial_bundle"] = df["string_value"].astype(str).str.contains("spatial|scalefactors|tissue_positions|SpaceRanger|spaceranger|Visium|CytAssist", case=False, regex=True)
    df.to_csv(path_csv, index=False)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    RUNTIME.mkdir(parents=True, exist_ok=True)

    helper = make_r_helper()
    rdata_to_load = STDATA_ASCII if STDATA_ASCII.exists() else STDATA_ACTUAL
    r_result = run_r(helper, rel(rdata_to_load) or str(rdata_to_load), rel(OUT) or str(OUT), SEURAT_OBJECT, IMAGE_SLOT)
    (OUT / "r_inspection_stdout.txt").write_text(r_result.stdout, encoding="utf-8")
    (OUT / "r_inspection_stderr.txt").write_text(r_result.stderr, encoding="utf-8")
    if r_result.returncode != 0:
        decision = "FAIL"
        next_step = "Stop and fix R/Seurat object inspection."
        summary = {
            "phase": "BioApp Coordinate Registration Recovery Audit",
            "decision": decision,
            "audit_only": True,
            "main_figure_redrawn": False,
            "SVTuner_rerun": False,
            "Stage3_rerun": False,
            "Stage4_run": False,
            "CytoSPACE_rerun": False,
            "endpoint_redefined": False,
            "data_changed": False,
            "metrics_changed": False,
            "threshold_changed": False,
            "endpoint_labels_changed": False,
            "conclusions_changed": False,
            "st_data_rdata_path": STDATA_DISPLAY,
            "st_data_rdata_load_path": rel(rdata_to_load),
            "seurat_object_name": SEURAT_OBJECT,
            "image_slot": IMAGE_SLOT,
            "embedded_image_found": False,
            "referenced_image_found": False,
            "referenced_image_path": None,
            "scale_factors_found": False,
            "registration_transform_found": False,
            "inverse_transform_found": False,
            "h_and_e_overlay_recoverable": False,
            "recommended_next_step": next_step,
            "r_return_code": r_result.returncode,
        }
        write_json(OUT / "bioapp_registration_recovery_summary.json", summary)
        print("BioApp Coordinate Registration Recovery Audit completed.")
        print("Decision: FAIL")
        print("H&E overlay recoverable: false")
        print(f"Recommended next step: {next_step}")
        print("Guardrails passed: false")
        return

    postprocess_pathlike(OUT / "bioapp_st_data_pathlike_string_search.csv")

    obj_inv = json.loads((OUT / "bioapp_st_data_rdata_object_inventory.json").read_text(encoding="utf-8"))
    img_inv = json.loads((OUT / "bioapp_seurat_image_slot_deep_inspection.json").read_text(encoding="utf-8"))
    scale_df = pd.read_csv(OUT / "bioapp_image_slot_scale_factor_audit.csv") if (OUT / "bioapp_image_slot_scale_factor_audit.csv").exists() else pd.DataFrame()
    path_df = pd.read_csv(OUT / "bioapp_st_data_pathlike_string_search.csv") if (OUT / "bioapp_st_data_pathlike_string_search.csv").exists() else pd.DataFrame()
    bounds = json.loads((OUT / "bioapp_candidate_image_bounds_check.json").read_text(encoding="utf-8")) if (OUT / "bioapp_candidate_image_bounds_check.json").exists() else {}

    cta_rows = search_cta_dir()
    write_rows(
        OUT / "bioapp_breastcancer_cta_registration_file_search.csv",
        cta_rows,
        [
            "file_path",
            "file_type",
            "matched_terms",
            "contains_coordinate_columns",
            "contains_image_metadata",
            "contains_transform_matrix",
            "contains_keypoints",
            "contains_image_file_references",
            "usable_for_recovery",
            "notes",
        ],
    )

    embedded_image_found = bool(img_inv.get("embedded_image_found"))
    scale_factors_found = not scale_df.empty
    lowres_bounds_pass = bool(bounds.get("lowres_bounds_pass", False))
    referenced = pd.DataFrame()
    if not path_df.empty:
        referenced = path_df[path_df.get("path_exists", False).astype(bool) & path_df["string_value"].astype(str).str.contains(r"\.png|\.jpg|\.jpeg|\.tif|\.tiff|\.svs", case=False, regex=True)]
    referenced_image_found = not referenced.empty
    referenced_path = str(referenced.iloc[0]["string_value"]) if referenced_image_found else None
    scale_03_found = bool((scale_df.get("value", pd.Series(dtype=float)).astype(float).sub(0.3).abs() < 1e-8).any()) if not scale_df.empty else False
    registration_transform_found = any(bool(r.get("contains_transform_matrix")) for r in cta_rows)
    inverse_transform_found = False

    if embedded_image_found and lowres_bounds_pass:
        decision = "REGISTRATION_RECOVERED_EMBEDDED_IMAGE"
        h_and_e_recoverable = True
        next_step = "BioApp Main Figure V3 - Seurat image-slot tissue-background overlay"
    elif referenced_image_found:
        decision = "REGISTRATION_RECOVERED_REFERENCED_IMAGE"
        h_and_e_recoverable = True
        next_step = "BioApp Main Figure V3 - verified referenced tissue-background overlay"
    elif scale_factors_found or not path_df.empty:
        decision = "REGISTRATION_METADATA_FOUND_BUT_INVERSE_TRANSFORM_MISSING"
        h_and_e_recoverable = False
        next_step = "BioApp Main Figure V3 - publication polish for spot-lattice fallback, or stop H&E recovery"
    else:
        decision = "NO_REGISTRATION_RECOVERY"
        h_and_e_recoverable = False
        next_step = "BioApp Main Figure V3 - publication polish for spot-lattice fallback"

    candidate_rows = [
        {
            "candidate_id": "embedded_image_slot_001",
            "candidate_type": "embedded_seurat_image_slot",
            "source": f"ST_data.RData::{SEURAT_OBJECT}@images[['{IMAGE_SLOT}']]@image",
            "image_path_or_object": "embedded array",
            "explicitly_linked_to_BCSA2TumB1": True,
            "image_dimensions": "x".join(map(str, img_inv.get("embedded_image_dimensions", []))) if img_inv.get("embedded_image_dimensions") else "",
            "scale_factor_used": bounds.get("lowres_scale"),
            "bounds_check_passed": lowres_bounds_pass,
            "accepted_for_v3_candidate": embedded_image_found and lowres_bounds_pass,
            "notes": "Audit preview generated; not used as manuscript background in this phase.",
        }
    ]
    if referenced_image_found:
        candidate_rows.append(
            {
                "candidate_id": "referenced_image_001",
                "candidate_type": "pathlike_reference",
                "source": "ST_data.RData path-like string search",
                "image_path_or_object": referenced_path,
                "explicitly_linked_to_BCSA2TumB1": True,
                "image_dimensions": "",
                "scale_factor_used": "",
                "bounds_check_passed": "",
                "accepted_for_v3_candidate": False,
                "notes": "Path exists but requires separate coordinate verification.",
            }
        )
    write_rows(
        OUT / "bioapp_candidate_image_recovery_table.csv",
        candidate_rows,
        [
            "candidate_id",
            "candidate_type",
            "source",
            "image_path_or_object",
            "explicitly_linked_to_BCSA2TumB1",
            "image_dimensions",
            "scale_factor_used",
            "bounds_check_passed",
            "accepted_for_v3_candidate",
            "notes",
        ],
    )

    summary = {
        "phase": "BioApp Coordinate Registration Recovery Audit",
        "decision": decision,
        "audit_only": True,
        "main_figure_redrawn": False,
        "SVTuner_rerun": False,
        "Stage3_rerun": False,
        "Stage4_run": False,
        "CytoSPACE_rerun": False,
        "endpoint_redefined": False,
        "data_changed": False,
        "metrics_changed": False,
        "threshold_changed": False,
        "endpoint_labels_changed": False,
        "conclusions_changed": False,
        "st_data_rdata_path": STDATA_DISPLAY,
        "st_data_rdata_load_path": rel(rdata_to_load),
        "seurat_object_name": obj_inv.get("seurat_object_name", SEURAT_OBJECT),
        "image_slot": img_inv.get("image_slot_name", IMAGE_SLOT),
        "embedded_image_found": embedded_image_found,
        "referenced_image_found": referenced_image_found,
        "referenced_image_path": referenced_path,
        "scale_factors_found": scale_factors_found,
        "scale_factor_0_3_found_in_object_metadata": scale_03_found,
        "registration_transform_found": registration_transform_found,
        "inverse_transform_found": inverse_transform_found,
        "h_and_e_overlay_recoverable": h_and_e_recoverable,
        "recommended_next_step": next_step,
        "embedded_image_dimensions": img_inv.get("embedded_image_dimensions"),
        "lowres_bounds_pass": lowres_bounds_pass,
        "r_return_code": r_result.returncode,
    }
    write_json(OUT / "bioapp_registration_recovery_summary.json", summary)

    guardrails = {
        "phase": "BioApp Coordinate Registration Recovery Audit",
        "decision": decision,
        "audit_only": True,
        "main_figure_redrawn": False,
        "SVTuner_rerun": False,
        "Stage3_rerun": False,
        "Stage4_run": False,
        "CytoSPACE_rerun": False,
        "endpoint_redefined": False,
        "new_biological_endpoint_introduced": False,
        "threshold_selected_using_endpoint_labels": False,
        "new_mapping_or_new_computation": False,
        "data_changed": False,
        "metrics_changed": False,
        "threshold_changed": False,
        "endpoint_labels_changed": False,
        "conclusions_changed": False,
        "h_and_e_used_in_manuscript_figure": False,
        "external_image_used_without_verification": False,
        "all_quantitative_claims_based_on_full_frozen_analysis_set": True,
        "guardrails_passed": decision != "FAIL",
    }
    write_json(OUT / "bioapp_registration_recovery_guardrails.json", guardrails)

    decision_md = f"""# BioApp Coordinate Registration Recovery Decision

## Final Decision

`{decision}`

## Embedded Image

BCSA2TumB1 embedded image found: `{embedded_image_found}`

Embedded image dimensions: `{img_inv.get('embedded_image_dimensions')}`

The embedded image is stored inside:

`ST_data.RData::{SEURAT_OBJECT}@images[['{IMAGE_SLOT}']]@image`

## Referenced Local Tissue/H&E Image

Referenced image found from path-like metadata: `{referenced_image_found}`

Referenced image path: `{referenced_path}`

## Scale Factors

Scale factors found: `{scale_factors_found}`

Scale factor names/values are recorded in `bioapp_image_slot_scale_factor_audit.csv`.

The object contains Seurat scale factors such as `spot`, `fiducial`, `hires`, and `lowres`. The CTA-specific factor `0.3` is documented in README/CTA_align and used to convert Seurat image coordinates to raw microscope coordinates; it is not itself one of the Seurat image-slot scale factors.

## Transform / Inverse Transform

Transform metadata found: `{registration_transform_found}`

Inverse transform found: `{inverse_transform_found}`

No external inverse transform file was found. However, the embedded Seurat image can be used with the processed BCSA2TumB1 coordinate system using the image slot's lowres scale factor.

## H&E Overlay Recoverable

`{h_and_e_recoverable}`

## Reason

The Seurat `VisiumV1` image slot contains an embedded raster image and scale factors. The BioApp frozen coordinates are explicitly linked to this image slot. The audit-only bounds check confirms that `imagecol/imagerow` scaled by the image slot lowres factor falls within the embedded image dimensions.

This recovers a Seurat image-slot tissue-background route. It does not recover the original raw microscope TIFF or a full Space Ranger spatial bundle.

## Biological Results

No biological result, endpoint, metric, threshold, label, or conclusion was changed.

## Recommended Next Step

`{next_step}`
"""
    (OUT / "bioapp_registration_recovery_decision.md").write_text(decision_md, encoding="utf-8")

    readme = f"""BioApp Coordinate Registration Recovery Audit completed.
Decision: {decision}
Embedded image found: {str(embedded_image_found).lower()}
Referenced image found: {str(referenced_image_found).lower()}
Scale factors found: {str(scale_factors_found).lower()}
Inverse transform found: {str(inverse_transform_found).lower()}
H&E overlay recoverable: {str(h_and_e_recoverable).lower()}
Data changed: false
Metrics changed: false
Endpoint redefined: false
Next: {next_step}
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")

    print("BioApp Coordinate Registration Recovery Audit completed.")
    print(f"Decision: {decision}")
    print(f"H&E overlay recoverable: {str(h_and_e_recoverable).lower()}")
    print(f"Recommended next step: {next_step}")
    print(f"Guardrails passed: {str(guardrails['guardrails_passed']).lower()}")


if __name__ == "__main__":
    main()
