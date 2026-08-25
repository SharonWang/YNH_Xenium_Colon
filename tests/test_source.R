#!/usr/bin/env Rscript

# Contract tests for the colon Xenium QC source and configuration.
# All test scratch paths inherit TMPDIR from the D-drive-only launcher.
options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
test_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(test_path))

manifest_path <- file.path(repo_root, "config", "colon_sample_manifest.tsv")
threshold_path <- file.path(repo_root, "config", "fixed_cell_qc_thresholds.tsv")
extended_path <- file.path(repo_root, "config", "extended_qc_defaults.tsv")
eos_path <- file.path(repo_root, "config", "eos_gene_sets.tsv")

if (!file.exists(manifest_path)) stop("colon_sample_manifest.tsv is missing", call. = FALSE)
if (!file.exists(threshold_path)) stop("fixed_cell_qc_thresholds.tsv is missing", call. = FALSE)
if (!file.exists(extended_path)) stop("extended_qc_defaults.tsv is missing", call. = FALSE)

manifest <- utils::read.delim(manifest_path, check.names = FALSE)
expected_regions <- paste0("Region_", 1:6)
stopifnot(
  identical(as.character(manifest$region_id), expected_regions),
  identical(as.character(manifest$mouse_id), rep(c("Mouse_1", "Mouse_2"), each = 3L)),
  identical(as.character(manifest$position), rep(c("top", "middle", "bottom"), 2L)),
  identical(as.character(manifest$tissue), rep("colon", 6L)),
  identical(as.character(manifest$genotype), rep("WT", 6L)),
  all(manifest$metadata_status == "VERIFIED_USER_SUPPLIED"),
  !any(manifest$do_not_interpret)
)

thresholds <- utils::read.delim(threshold_path, check.names = FALSE)
expected_thresholds <- data.frame(
  metric = c("nFeature_Xenium", "nFeature_Xenium", "nCount_Xenium", "nCount_Xenium"),
  bound = c("lower", "upper", "lower", "upper"),
  value = c(5L, 200L, 10L, 1000L),
  inclusive = c(FALSE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)
stopifnot(identical(thresholds[, names(expected_thresholds)], expected_thresholds))

eos <- utils::read.delim(eos_path, check.names = FALSE)
stopifnot(
  nrow(eos) == 100L,
  length(unique(eos$gene)) == 100L,
  sum(eos$gene_set == "common") == 7L,
  sum(eos$gene_set == "short_lived") == 47L,
  sum(eos$gene_set == "long_lived") == 46L
)

source(file.path(repo_root, "R", "source.R"))

# Break caught: changing an exclusive threshold to inclusive, or reverting to
# distribution-derived bounds, must fail this literal boundary fixture.
make_count_column <- function(n_feature, n_count, n_genes = 250L) {
  values <- numeric(n_genes)
  values[seq_len(n_feature)] <- 1
  values[[1]] <- values[[1]] + n_count - n_feature
  values
}

boundary_cases <- data.frame(
  cell_id = paste0("cell_", 1:6),
  n_feature = c(5L, 6L, 199L, 200L, 6L, 6L),
  n_count = c(11L, 11L, 999L, 999L, 10L, 1000L),
  expected_primary = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

count_dense <- vapply(
  seq_len(nrow(boundary_cases)),
  function(i) make_count_column(boundary_cases$n_feature[[i]], boundary_cases$n_count[[i]]),
  numeric(250L)
)
rownames(count_dense) <- paste0("Gene_", seq_len(nrow(count_dense)))
colnames(count_dense) <- boundary_cases$cell_id
counts <- Matrix::Matrix(count_dense, sparse = TRUE)
cells <- data.frame(
  cell_id = boundary_cases$cell_id,
  total_counts = boundary_cases$n_count,
  control_probe_counts = 0,
  genomic_control_counts = 0,
  control_codeword_counts = 0,
  cell_area = 100,
  nucleus_count = 1L,
  x_centroid = seq_len(nrow(boundary_cases)),
  y_centroid = seq_len(nrow(boundary_cases)),
  stringsAsFactors = FALSE
)

fixed <- read_fixed_cell_qc_thresholds(threshold_path)
qc <- calculate_xenium_cell_qc(counts, cells, "Region_1", fixed)
stopifnot(identical(as.logical(qc$cell_metadata$qc_core_pass), boundary_cases$expected_primary))

# Break caught: adding segmentation/control review flags to primary_include
# violates the approved non-destructive primary mask.
mask_fixture <- qc$cell_metadata
mask_fixture$segmentation_multiplet_flag[[2]] <- TRUE
mask_fixture$high_control_flag[[3]] <- TRUE
mask_fixture$qc_review_flag <- mask_fixture$qc_review_flag |
  mask_fixture$segmentation_multiplet_flag | mask_fixture$high_control_flag
masks <- build_cell_downstream_masks(mask_fixture, provenance = "unit-test")
stopifnot(
  identical(as.logical(masks$primary_include), boundary_cases$expected_primary),
  all(!masks$strict_include | masks$primary_include),
  !masks$strict_include[[2]],
  !masks$strict_include[[3]],
  all(!masks$hotspot_sensitivity_include | masks$primary_include)
)

cat("Colon source and configuration contracts passed.\n")
