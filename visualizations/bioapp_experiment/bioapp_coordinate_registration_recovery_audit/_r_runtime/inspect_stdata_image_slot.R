
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
