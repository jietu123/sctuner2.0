#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(fgsea)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else "."
project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)

cytospace_dir <- if (length(args) >= 2) {
  normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
} else {
  file.path(project_root, "result", "cytospace_fig2c_official_slide1_readme_default", "cytospace_output")
}

out_dir <- if (length(args) >= 3) {
  args[[3]]
} else {
  file.path(project_root, "result", "cytospace_fig2c_official_slide1_readme_default", "fig2c_enrichment_seurat410_fgsea114")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tumor_labels <- if (length(args) >= 4 && nzchar(args[[4]])) {
  strsplit(args[[4]], "\\|")[[1]]
} else {
  c("Melanoma cells", "Melanoma")
}

slide_label <- if (length(args) >= 5 && nzchar(args[[5]])) {
  args[[5]]
} else {
  "ST_mel1_rep2"
}

expr_dir <- file.path(cytospace_dir, "assigned_expression")
loc_path <- file.path(cytospace_dir, "assigned_locations.csv")
gene_set_path <- file.path(
  project_root,
  "data", "raw", "cytospace_fig2c_melanoma", "prepared", "gene_sets",
  "tcell_exhaustion_zheng_cell2017.txt"
)

stopifnot(file.exists(file.path(expr_dir, "matrix.mtx")))
stopifnot(file.exists(file.path(expr_dir, "genes.tsv")))
stopifnot(file.exists(file.path(expr_dir, "barcodes.tsv")))
stopifnot(file.exists(loc_path))
stopifnot(file.exists(gene_set_path))

message("[Fig2c exact] Seurat ", as.character(packageVersion("Seurat")))
message("[Fig2c exact] fgsea ", as.character(packageVersion("fgsea")))
message("[Fig2c exact] reading assigned expression")

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

melanoma <- meta[meta$CellType %in% tumor_labels, , drop = FALSE]
if (nrow(melanoma) < 5) {
  stop("Need at least five mapped melanoma cells; found ", nrow(melanoma))
}
coord_cols <- if (all(c("X", "Y") %in% colnames(meta))) {
  c("X", "Y")
} else if (all(c("row", "col") %in% colnames(meta))) {
  c("row", "col")
} else {
  stop("assigned_locations.csv must contain X/Y or row/col coordinate columns")
}
tumor_xy <- as.matrix(melanoma[, coord_cols])

running_enrichment <- function(stats, geneset) {
  stats <- sort(stats[is.finite(stats)], decreasing = TRUE)
  hits <- names(stats) %in% geneset
  n_hits <- sum(hits)
  n_miss <- length(stats) - n_hits
  if (n_hits == 0 || n_miss == 0) {
    stop("Invalid ranked list for enrichment curve")
  }

  hit_norm <- sum(abs(stats[hits]))
  increments <- ifelse(hits, abs(stats) / hit_norm, -1 / n_miss)
  data.frame(
    rank = seq_along(stats),
    gene = names(stats),
    stat = as.numeric(stats),
    hit = hits,
    running_es = as.numeric(cumsum(increments)),
    stringsAsFactors = FALSE
  )
}

compute_for_type <- function(cell_type) {
  message("[Fig2c exact] processing ", cell_type)
  target_meta <- meta[meta$CellType == cell_type, , drop = FALSE]
  if (nrow(target_meta) < 20) {
    stop("Too few mapped cells for ", cell_type, ": ", nrow(target_meta))
  }

  target_xy <- as.matrix(target_meta[, coord_cols])
  dmat <- as.matrix(dist(rbind(target_xy, tumor_xy)))
  dmat <- dmat[
    seq_len(nrow(target_xy)),
    nrow(target_xy) + seq_len(nrow(tumor_xy)),
    drop = FALSE
  ]
  mean_nearest5 <- apply(dmat, 1, function(x) mean(sort(x, partial = 5)[1:5]))
  cutoff <- median(mean_nearest5)
  distance_group <- ifelse(mean_nearest5 <= cutoff, "close", "far")
  names(distance_group) <- rownames(target_meta)

  sub_counts <- counts[, rownames(target_meta), drop = FALSE]
  obj <- CreateSeuratObject(counts = sub_counts, meta.data = target_meta)
  obj$distance_group <- factor(distance_group[colnames(obj)], levels = c("close", "far"))
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  Idents(obj) <- "distance_group"

  fc <- FoldChange(
    obj,
    ident.1 = "close",
    ident.2 = "far",
    slot = "data",
    base = 2,
    pseudocount.use = 1
  )

  fc_col <- grep("avg_log", colnames(fc), value = TRUE)[1]
  if (is.na(fc_col)) {
    stop("Could not find avg_log fold-change column in FoldChange output")
  }
  stats <- fc[[fc_col]]
  names(stats) <- rownames(fc)
  stats <- sort(stats[is.finite(stats)], decreasing = TRUE)

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
    X = target_meta[[coord_cols[1]]],
    Y = target_meta[[coord_cols[2]]],
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
    slide = slide_label,
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
    seurat_version = as.character(packageVersion("Seurat")),
    fgsea_version = as.character(packageVersion("fgsea")),
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, lapply(c("CD4 T cells", "CD8 T cells"), compute_for_type))
write.csv(results, file.path(out_dir, "fig2c_official_enrichment_summary.csv"), row.names = FALSE)

message("[OK] wrote: ", file.path(out_dir, "fig2c_official_enrichment_summary.csv"))
print(results)
