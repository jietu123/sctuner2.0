# Stage 1: Preprocess scRNA and ST data (Seurat, mild/adjustable filtering)
# Usage (from project root):
#   Rscript r_scripts/stage1_preprocess.R --sample real_brca

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a)) a else b

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(jsonlite)
  library(yaml)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
sample_id <- if (length(args) >= 2 && args[1] == "--sample") args[2] else "real_brca"

# Resolve project root: script is under r_scripts/, so parent is project root
get_script_path <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1])))
  }
  if (!is.null(sys.frames()) && length(sys.frames()) > 0) {
    fr <- sys.frames()[[1]]
    if (!is.null(fr$ofile)) {
      return(normalizePath(fr$ofile))
    }
  }
  normalizePath(file.path("r_scripts", "stage1_preprocess.R"))
}

script_path <- get_script_path()
project_root <- normalizePath(file.path(dirname(script_path), ".."))

# Paths
input_dir <- file.path(project_root, "data", "raw", sample_id)
processed_dir <- file.path(project_root, "data", "processed", sample_id, "stage1_preprocess")
summary_dir <- file.path(project_root, "result", sample_id, "stage1_preprocess")
qc_dir <- file.path(summary_dir, "qc_plots")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# Defaults
default_cfg <- list(
  paths = list(
    sc_expr = sprintf("%s_scRNA_GEP.txt", sample_id),
    sc_meta = sprintf("%s_scRNA_celllabels.txt", sample_id),
    st_expr = sprintf("%s_STdata_GEP.txt", sample_id),
    st_meta = sprintf("%s_STdata_coordinates.txt", sample_id),
    svg_marker_whitelist = NULL
  ),
  qc = list(
    sc_min_genes = 200,
    sc_max_genes = 6000,
    sc_max_mt = 10,
    st_min_genes = 100,
    st_max_genes = Inf,
    st_max_mt = 20,
    hvg_nfeatures = 2000,
    mt_pattern = "^MT-"
  ),
  gene_filter = list(
    min_cells_sc = 0,
    min_cells_st = 0
  )
)

read_yaml_safe <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(list())
  tryCatch(yaml::read_yaml(path), error = function(e) list())
}

dataset_cfg_path <- file.path(project_root, "configs", "datasets", paste0(sample_id, ".yaml"))
project_cfg_path <- file.path(project_root, "configs", "project_config.yaml")

cfg <- default_cfg
cfg <- modifyList(cfg, read_yaml_safe(project_cfg_path))
cfg <- modifyList(cfg, read_yaml_safe(dataset_cfg_path))

resolve_path <- function(p) {
  if (is.null(p)) return(NULL)
  if (grepl("^([A-Za-z]:|/)", p)) return(normalizePath(p, winslash = "/"))
  normalizePath(file.path(input_dir, p), winslash = "/", mustWork = FALSE)
}

paths <- lapply(cfg$paths, resolve_path)

qc <- cfg$qc
gene_filter <- cfg$gene_filter

guess_prefixes <- unique(c(sample_id, sub("^.*_", "", sample_id)))

find_input <- function(key, suffix) {
  if (!is.null(paths[[key]]) && file.exists(paths[[key]])) {
    return(paths[[key]])
  }
  for (pfx in guess_prefixes) {
    cand <- file.path(input_dir, sprintf("%s_%s", pfx, suffix))
    if (file.exists(cand)) return(cand)
  }
  # fallback to first prefix even if missing (will error downstream)
  file.path(input_dir, sprintf("%s_%s", guess_prefixes[[1]], suffix))
}

read_expression <- function(path) {
  dt <- fread(path)
  genes <- dt[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- genes
  mode(mat) <- "numeric"
  Matrix(mat, sparse = TRUE)
}

load_whitelist <- function(path) {
  if (is.null(path) || !file.exists(path)) return(character(0))
  trimws(readLines(path))
}

filter_genes <- function(mat, min_cells, whitelist) {
  if (min_cells <= 0) return(mat)
  keep <- rowSums(mat > 0) >= min_cells
  if (length(whitelist) > 0) {
    keep[rownames(mat) %in% whitelist] <- TRUE
  }
  mat[keep, , drop = FALSE]
}

load_scrna <- function() {
  expr_path <- find_input("sc_expr", "scRNA_GEP.txt")
  meta_path <- find_input("sc_meta", "scRNA_celllabels.txt")
  mat_raw <- read_expression(expr_path)
  meta <- fread(meta_path)
  setnames(meta, c("cell_id", "cell_type"))
  setkey(meta, cell_id)
  common_cells <- intersect(colnames(mat_raw), meta$cell_id)
  mat_raw <- mat_raw[, common_cells, drop = FALSE]
  meta_use <- meta[common_cells]
  rownames(meta_use) <- meta_use$cell_id

  whitelist <- load_whitelist(paths$svg_marker_whitelist)
  mat_filtered <- filter_genes(mat_raw, gene_filter$min_cells_sc %||% 0, whitelist)

  stats <- list(
    n_cells_raw = ncol(mat_raw),
    n_genes_raw = nrow(mat_raw),
    n_cells_after_gene_filter = ncol(mat_filtered),
    n_genes_after_gene_filter = nrow(mat_filtered),
    whitelist_retained = sum(rownames(mat_filtered) %in% whitelist)
  )

  seurat_obj <- CreateSeuratObject(mat_filtered, meta.data = meta_use)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = qc$mt_pattern %||% "^MT-")
  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA >= qc$sc_min_genes &
      nFeature_RNA <= qc$sc_max_genes &
      percent.mt <= qc$sc_max_mt
  )
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = qc$hvg_nfeatures)
  list(obj = seurat_obj, stats = stats, whitelist = whitelist)
}

load_st <- function() {
  expr_path <- find_input("st_expr", "STdata_GEP.txt")
  coord_path <- find_input("st_meta", "STdata_coordinates.txt")
  mat_raw <- read_expression(expr_path)
  coords <- fread(coord_path)
  setnames(coords, c("spot_id", "row", "col"))
  setkey(coords, spot_id)
  common_spots <- intersect(colnames(mat_raw), coords$spot_id)
  mat_raw <- mat_raw[, common_spots, drop = FALSE]
  coords_use <- coords[common_spots]
  rownames(coords_use) <- coords_use$spot_id

  whitelist <- load_whitelist(paths$svg_marker_whitelist)
  mat_filtered <- filter_genes(mat_raw, gene_filter$min_cells_st %||% 0, whitelist)

  stats <- list(
    n_cells_raw = ncol(mat_raw),
    n_genes_raw = nrow(mat_raw),
    n_cells_after_gene_filter = ncol(mat_filtered),
    n_genes_after_gene_filter = nrow(mat_filtered),
    whitelist_retained = sum(rownames(mat_filtered) %in% whitelist)
  )

  seurat_obj <- CreateSeuratObject(mat_filtered, meta.data = coords_use)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = qc$mt_pattern %||% "^MT-")
  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA >= qc$st_min_genes &
      nFeature_RNA <= qc$st_max_genes &
      percent.mt <= qc$st_max_mt
  )
  seurat_obj <- NormalizeData(seurat_obj)
  list(obj = seurat_obj, stats = stats, whitelist = whitelist)
}

plot_hist <- function(df, col, title, path) {
  p <- ggplot(df, aes_string(x = col)) +
    geom_histogram(bins = 50, fill = "#3182bd", color = "white", alpha = 0.8) +
    theme_minimal() +
    ggtitle(title)
  ggsave(path, p, width = 6, height = 4, dpi = 150)
}

message("[Stage1] Loading scRNA...")
sc_res <- load_scrna()
sc <- sc_res$obj

message("[Stage1] Loading ST...")
st_res <- load_st()
st <- st_res$obj

message("[Stage1] Aligning genes...")
common_genes <- intersect(rownames(sc), rownames(st))
common_genes_path <- file.path(processed_dir, "common_genes.txt")
writeLines(common_genes, common_genes_path)
sc <- subset(sc, features = common_genes)
st <- subset(st, features = common_genes)

# HVG info
hvg_list <- VariableFeatures(sc)
hvg_path <- file.path(processed_dir, "hvg_genes.txt")
writeLines(hvg_list, hvg_path)

# QC stats
mt_stats <- function(obj) {
  pct <- obj$percent.mt
  list(mean = mean(pct), median = median(pct), max = max(pct))
}

summary_list <- list(
  sample = sample_id,
  input_dir = input_dir,
  processed_dir = processed_dir,
  summary_dir = summary_dir,
  qc_dir = qc_dir,
  qc_params = qc,
  gene_filter = gene_filter,
  paths = paths,
  stats = list(
    sc = c(sc_res$stats, n_cells_filtered = ncol(sc), n_genes_filtered = nrow(sc)),
    st = c(st_res$stats, n_cells_filtered = ncol(st), n_genes_filtered = nrow(st))
  ),
  genes = list(
    common = length(common_genes),
    common_genes_path = common_genes_path,
    hvg_requested = qc$hvg_nfeatures,
    hvg_actual = length(hvg_list),
    hvg_genes_path = hvg_path
  ),
  mt_pct = list(
    sc = mt_stats(sc),
    st = mt_stats(st)
  )
)

# QC plots (basic histograms)
try({
  sc_df <- sc@meta.data
  st_df <- st@meta.data
  plot_hist(sc_df, "nCount_RNA", "scRNA nCount_RNA", file.path(qc_dir, "sc_qc_nCount.png"))
  plot_hist(sc_df, "nFeature_RNA", "scRNA nFeature_RNA", file.path(qc_dir, "sc_qc_nFeature.png"))
  plot_hist(sc_df, "percent.mt", "scRNA percent.mt", file.path(qc_dir, "sc_qc_percent_mt.png"))
  plot_hist(st_df, "nCount_RNA", "ST nCount_RNA", file.path(qc_dir, "st_qc_nCount.png"))
  plot_hist(st_df, "nFeature_RNA", "ST nFeature_RNA", file.path(qc_dir, "st_qc_nFeature.png"))
}, silent = TRUE)

message("[Stage1] Saving outputs...")
saveRDS(sc, file.path(processed_dir, "sc_processed.rds"))
saveRDS(st, file.path(processed_dir, "st_processed.rds"))
write_json(summary_list, file.path(summary_dir, "stage1_summary.json"), auto_unbox = TRUE, pretty = TRUE)

message("[Stage1] Done.")
