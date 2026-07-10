
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
