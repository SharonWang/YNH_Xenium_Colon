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

cat("Configuration contracts passed.\n")
