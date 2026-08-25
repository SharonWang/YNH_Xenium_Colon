#!/usr/bin/env Rscript

# Build a deterministic, balanced, sparse 500-cell subset for each region.
# Raw inputs are read-only; every output and temporary file is constrained to
# the supplied colon_analysis directory.
suppressPackageStartupMessages(library(Matrix))

select_balanced_cells <- function(counts, cells, target = 500L, seed = 20260814L, eos_genes = character()) {
  stopifnot(ncol(counts) == nrow(cells), identical(colnames(counts), cells$cell_id), target > 0L)
  n_count <- as.numeric(Matrix::colSums(counts))
  n_feature <- as.numeric(Matrix::colSums(counts > 0))
  selected <- integer()
  labels <- rep(NA_character_, nrow(cells))
  add <- function(indices, label, limit) {
    indices <- setdiff(unique(as.integer(indices)), selected)
    indices <- indices[is.finite(indices) & indices >= 1L & indices <= nrow(cells)]
    if (length(indices) > limit) indices <- indices[seq_len(limit)]
    selected <<- c(selected, indices)
    labels[indices] <<- label
  }
  nearest <- function(values, boundary, number = 15L) order(abs(values - boundary), seq_along(values))[seq_len(min(number, length(values)))]
  for (boundary in c(5, 200)) add(nearest(n_feature, boundary), sprintf("nFeature_near_%s", boundary), 15L)
  for (boundary in c(10, 1000)) add(nearest(n_count, boundary), sprintf("nCount_near_%s", boundary), 15L)
  if ("nucleus_count" %in% names(cells)) {
    add(which(cells$nucleus_count == 0), "nucleus_missing", 25L)
    add(which(cells$nucleus_count > 1), "multiple_nuclei", 25L)
  }
  control_columns <- intersect(c("control_probe_counts", "genomic_control_counts", "control_codeword_counts"), names(cells))
  if (length(control_columns)) {
    control <- rowSums(cells[, control_columns, drop = FALSE], na.rm = TRUE)
    add(order(control, decreasing = TRUE), "high_control", 30L)
  }
  coordinate_columns <- intersect(c("x_centroid", "y_centroid"), names(cells))
  if (length(coordinate_columns) == 2L) {
    x_bin <- cut(cells[[coordinate_columns[[1L]]]], breaks = 5L, include.lowest = TRUE, labels = FALSE)
    y_bin <- cut(cells[[coordinate_columns[[2L]]]], breaks = 5L, include.lowest = TRUE, labels = FALSE)
    spatial_groups <- split(seq_len(nrow(cells)), interaction(x_bin, y_bin, drop = TRUE))
    spatial_pick <- unlist(lapply(spatial_groups, function(indices) indices[[ceiling(length(indices) / 2)]]), use.names = FALSE)
    add(spatial_pick, "spatial_coverage", 50L)
  }
  eos_rows <- which(rownames(counts) %in% eos_genes)
  if (length(eos_rows)) {
    eos_signal <- as.numeric(Matrix::colSums(counts[eos_rows, , drop = FALSE]))
    add(order(eos_signal, decreasing = TRUE), "eos_panel_signal", 50L)
  }
  set.seed(seed)
  remaining <- setdiff(seq_len(nrow(cells)), selected)
  needed <- max(0L, min(as.integer(target), nrow(cells)) - length(selected))
  if (needed) add(sample(remaining, needed, replace = FALSE), "random_background", needed)
  selected <- selected[seq_len(min(length(selected), as.integer(target)))]
  labels[is.na(labels)] <- "not_selected"
  data.frame(index = selected, cell_id = cells$cell_id[selected], sampling_stratum = labels[selected], stringsAsFactors = FALSE)
}

gzip_file <- function(input, output) {
  input_connection <- file(input, "rb"); on.exit(close(input_connection), add = TRUE)
  output_connection <- gzfile(output, "wb"); on.exit(close(output_connection), add = TRUE)
  repeat {
    bytes <- readBin(input_connection, what = "raw", n = 1024L * 1024L)
    if (!length(bytes)) break
    writeBin(bytes, output_connection)
  }
  invisible(output)
}

build_local_subset <- function(raw_root, analysis_root, repo_root, target = 500L, seed = 20260814L) {
  raw_root <- normalizePath(raw_root, winslash = "/", mustWork = TRUE)
  analysis_root <- normalizePath(analysis_root, winslash = "/", mustWork = TRUE)
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  output_root <- file.path(analysis_root, "subset_input", "colon_data")
  if (grepl("^[Cc]:/", output_root) || !startsWith(output_root, paste0(analysis_root, "/"))) stop("Unsafe subset output root.", call. = FALSE)
  if (grepl("adipose_analysis", output_root, fixed = TRUE)) stop("Colon subset cannot use adipose_analysis.", call. = FALSE)
  source(file.path(repo_root, "R", "source.R"))
  regions <- discover_xenium_sections(raw_root, expected_section_count = 6L)
  eos <- utils::read.delim(file.path(repo_root, "config", "eos_gene_sets.tsv"), check.names = FALSE)
  eos_column <- intersect(c("gene", "gene_symbol", "symbol"), names(eos))[[1L]]
  manifest_rows <- list()
  for (row in seq_len(nrow(regions))) {
    region_id <- regions$region_id[[row]]
    source_dir <- regions$region_dir[[row]]
    destination_dir <- file.path(output_root, basename(source_dir))
    matrix_dir <- file.path(destination_dir, "cell_feature_matrix")
    dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
    imported <- import_xenium_mex(source_dir)
    selection <- select_balanced_cells(imported$counts, imported$cells, target, seed + row, eos[[eos_column]])
    selected_matrix <- imported$counts[, selection$index, drop = FALSE]
    selected_cells <- imported$cells[selection$index, , drop = FALSE]
    stopifnot(inherits(selected_matrix, "sparseMatrix"), identical(colnames(selected_matrix), selected_cells$cell_id))
    small_files <- c("experiment.xenium", "metrics_summary.csv", "analysis_summary.html", "gene_panel.json")
    copied <- file.copy(file.path(source_dir, small_files), file.path(destination_dir, small_files), overwrite = TRUE)
    if (!all(copied)) stop(sprintf("Failed to copy required small files for %s.", region_id), call. = FALSE)
    # import_xenium_mex returns the biological gene-expression rows only, so
    # the subset feature table must use the same rows as selected_matrix.
    gene_features <- imported$features[imported$features$feature_type == "Gene Expression", , drop = FALSE]
    stopifnot(nrow(gene_features) == nrow(selected_matrix))
    feature_connection <- gzfile(file.path(matrix_dir, "features.tsv.gz"), "wt")
    utils::write.table(gene_features, feature_connection, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    close(feature_connection)
    barcode_connection <- gzfile(file.path(matrix_dir, "barcodes.tsv.gz"), "wt")
    writeLines(selected_cells$cell_id, barcode_connection, useBytes = TRUE)
    close(barcode_connection)
    cells_connection <- gzfile(file.path(destination_dir, "cells.csv.gz"), "wt")
    utils::write.csv(selected_cells, cells_connection, row.names = FALSE, quote = TRUE)
    close(cells_connection)
    temporary_matrix <- file.path(matrix_dir, "matrix.mtx")
    Matrix::writeMM(selected_matrix, temporary_matrix)
    gzip_file(temporary_matrix, file.path(matrix_dir, "matrix.mtx.gz"))
    unlink(temporary_matrix)
    selection$region_id <- region_id
    selection$source_region_dir <- source_dir
    selection$subset_region_dir <- normalizePath(destination_dir, winslash = "/", mustWork = TRUE)
    selection$seed <- seed + row
    selection$source_cells <- ncol(imported$counts)
    selection$subset_cells <- ncol(selected_matrix)
    manifest_rows[[row]] <- selection
    validate_section_integrity(destination_dir, region_id)
  }
  subset_manifest <- do.call(rbind, manifest_rows)
  manifest_path <- file.path(analysis_root, "subset_input", "subset_manifest.tsv")
  utils::write.table(subset_manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE)
  subset_manifest
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 3L) stop("Usage: create_local_subset.R RAW_ROOT ANALYSIS_ROOT REPO_ROOT", call. = FALSE)
  result <- build_local_subset(args[[1L]], args[[2L]], args[[3L]])
  cat(sprintf("Created %d subset cells across %d regions.\n", nrow(result), length(unique(result$region_id))))
}
