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
