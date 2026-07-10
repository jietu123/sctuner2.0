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