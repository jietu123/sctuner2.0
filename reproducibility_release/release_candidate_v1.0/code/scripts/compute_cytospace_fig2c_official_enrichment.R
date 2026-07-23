#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(fgsea)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else "."
project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)

cytospace_dir <- if (length(args) >= 2) {
  normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
} else {
  file.path(project_root, "result", "cytospace_fig2c_official_slide1", "cytospace_output")
}
expr_dir <- file.path(cytospace_dir, "assigned_expression")
loc_path <- file.path(cytospace_dir, "assigned_locations.csv")
gene_set_path <- file.path(
  project_root,
  "data", "raw", "cytospace_fig2c_melanoma", "prepared", "gene_sets",
  "tcell_exhaustion_zheng_cell2017.txt"
)
out_dir <- if (length(args) >= 3) {
  args[[3]]
} else {
  file.path(project_root, "result", "cytospace_fig2c_official_slide1", "fig2c_enrichment")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(file.path(expr_dir, "matrix.mtx")))
stopifnot(file.exists(file.path(expr_dir, "genes.tsv")))
stopifnot(file.exists(file.path(expr_dir, "barcodes.tsv")))
stopifnot(file.exists(loc_path))
stopifnot(file.exists(gene_set_path))

message("[Fig2c official] reading assigned expression")
counts <- readMM(file.path(expr_dir, "matrix.mtx"))
genes_table <- read.delim(file.path(expr_dir, "genes.tsv"), header = FALSE, stringsAsFactors = FALSE)
barcodes <- readLines(file.path(expr_dir, "barcodes.tsv"), warn = FALSE)
genes <- make.unique(genes_table[[1]])
rownames(counts) <- genes
colnames(counts) <- barcodes

meta <- read.csv(loc_path, stringsAsFactors = FALSE, check.names = FALSE)
rownames(meta) <- meta$UniqueCID
meta <- meta[colnames(counts), , drop = FALSE]

exhaustion_genes <- unique(readLines(gene_set_path, warn = FALSE))
exhaustion_genes <- exhaustion_genes[nzchar(exhaustion_genes)]
exhaustion_genes <- intersect(exhaustion_genes, rownames(counts))
if (length(exhaustion_genes) < 10) {
  stop("Too few exhaustion genes found in assigned expression: ", length(exhaustion_genes))
}

melanoma <- meta[meta$CellType == "Melanoma cells", , drop = FALSE]
if (nrow(melanoma) < 5) {
  stop("Need at least five mapped melanoma cells; found ", nrow(melanoma))
}
tumor_xy <- as.matrix(melanoma[, c("X", "Y")])

running_enrichment <- function(stats, geneset) {
  stats <- sort(stats, decreasing = TRUE)
  stats <- stats[is.finite(stats)]
  hits <- names(stats) %in% geneset
  n_hits <- sum(hits)
  n_miss <- length(stats) - n_hits
  if (n_hits == 0 || n_miss == 0) {
    stop("Invalid ranked list for enrichment curve")
  }
  hit_weights <- abs(stats[hits])
  hit_norm <- sum(hit_weights)
  increments <- ifelse(hits, abs(stats) / hit_norm, -1 / n_miss)
  es <- cumsum(increments)
  data.frame(
    rank = seq_along(stats),
    gene = names(stats),
    stat = as.numeric(stats),
    hit = hits,
    running_es = as.numeric(es),
    stringsAsFactors = FALSE
  )
}

compute_for_type <- function(cell_type) {
  message("[Fig2c official] processing ", cell_type)
  target_meta <- meta[meta$CellType == cell_type, , drop = FALSE]
  if (nrow(target_meta) < 20) {
    stop("Too few mapped cells for ", cell_type, ": ", nrow(target_meta))
  }

  target_xy <- as.matrix(target_meta[, c("X", "Y")])
  dmat <- as.matrix(dist(rbind(target_xy, tumor_xy)))
  dmat <- dmat[seq_len(nrow(target_xy)), nrow(target_xy) + seq_len(nrow(tumor_xy)), drop = FALSE]
  mean_nearest5 <- apply(dmat, 1, function(x) mean(sort(x, partial = 5)[1:5]))
  cutoff <- median(mean_nearest5)
  distance_group <- ifelse(mean_nearest5 <= cutoff, "close", "far")
  names(distance_group) <- rownames(target_meta)

  # Equivalent to Seurat NormalizeData(normalization.method="LogNormalize",
  # scale.factor=10000) followed by FoldChange(slot="data", base=2,
  # pseudocount.use=1) for close versus far cells.
  sub_counts <- counts[, rownames(target_meta), drop = FALSE]
  lib_size <- Matrix::colSums(sub_counts)
  lib_size[lib_size == 0] <- 1
  norm_counts <- t(t(sub_counts) / lib_size) * 10000
  norm_data <- log1p(norm_counts)
  close_cells <- names(distance_group)[distance_group == "close"]
  far_cells <- names(distance_group)[distance_group == "far"]
  close_mean <- Matrix::rowMeans(expm1(norm_data[, close_cells, drop = FALSE]))
  far_mean <- Matrix::rowMeans(expm1(norm_data[, far_cells, drop = FALSE]))
  stats <- log2((close_mean + 1) / (far_mean + 1))
  names(stats) <- rownames(norm_data)
  stats <- stats[is.finite(stats)]
  stats <- sort(stats, decreasing = TRUE)

  set.seed(1)
  fg <- suppressWarnings(fgsea(
    pathways = list(Exhaustion = exhaustion_genes),
    stats = stats,
    nperm = 10000,
    minSize = 1,
    maxSize = 5000
  ))

  dist_out <- data.frame(
    UniqueCID = rownames(target_meta),
    OriginalCID = target_meta$OriginalCID,
    CellType = target_meta$CellType,
    SpotID = target_meta$SpotID,
    X = target_meta$X,
    Y = target_meta$Y,
    mean_distance_to_nearest5_melanoma = as.numeric(mean_nearest5),
    distance_group = distance_group,
    stringsAsFactors = FALSE
  )

  ranked_out <- data.frame(
    gene = names(stats),
    avg_log2FC_close_vs_far = as.numeric(stats),
    in_exhaustion_gene_set = names(stats) %in% exhaustion_genes,
    stringsAsFactors = FALSE
  )

  curve_out <- running_enrichment(stats, exhaustion_genes)
  curve_out$CellType <- cell_type

  prefix <- gsub("[^A-Za-z0-9]+", "_", tolower(cell_type))
  write.csv(dist_out, file.path(out_dir, paste0(prefix, "_distance_groups.csv")), row.names = FALSE)
  write.csv(ranked_out, file.path(out_dir, paste0(prefix, "_ranked_log2fc.csv")), row.names = FALSE)
  write.csv(curve_out, file.path(out_dir, paste0(prefix, "_enrichment_curve.csv")), row.names = FALSE)

  data.frame(
    method = "CytoSPACE",
    slide = "ST_mel1_rep2",
    cell_type = cell_type,
    pathway = fg$pathway,
    ES = fg$ES,
    NES = fg$NES,
    pval = fg$pval,
    padj = fg$padj,
    n_mapped_cells = nrow(target_meta),
    n_close = sum(distance_group == "close"),
    n_far = sum(distance_group == "far"),
    n_exhaustion_genes_used = length(exhaustion_genes),
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, lapply(c("CD4 T cells", "CD8 T cells"), compute_for_type))
write.csv(results, file.path(out_dir, "fig2c_official_enrichment_summary.csv"), row.names = FALSE)

message("[OK] wrote: ", file.path(out_dir, "fig2c_official_enrichment_summary.csv"))
print(results)
