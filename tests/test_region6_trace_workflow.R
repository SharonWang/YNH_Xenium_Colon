#!/usr/bin/env Rscript

# Integration contract for the transparent Region 6 trace-definition and
# projection workflow. The test checks executed artifacts, not source wording.

suppressPackageStartupMessages(library(jsonlite))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
test_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(test_path))
analysis_root <- dirname(repo_root)

source(file.path(repo_root, "R", "source.R"))

config_path <- file.path(
  repo_root, "config", "unrolling", "Region_6_trace_control_points.tsv"
)
if (!file.exists(config_path)) {
  stop("Version-controlled Region 6 trace controls are missing", call. = FALSE)
}
controls <- utils::read.delim(config_path, check.names = FALSE)
control_audit <- validate_trace_control_points(controls)
stopifnot(control_audit$n_control_points == 39L)
stopifnot(identical(controls$point_order, 1:39))
stopifnot(identical(as.numeric(controls[1, c("x", "y")]), c(4500, 720)))
stopifnot(identical(as.numeric(controls[39, c("x", "y")]), c(2825, 2550)))

notebook_paths <- file.path(
  repo_root, "notebooks",
  c(
    "B3a_Region6_trace_definition.ipynb",
    "B3_Region6_Mayassi_unrolling_method_audit.ipynb"
  )
)
stopifnot(all(file.exists(notebook_paths)))
for (path in notebook_paths) {
  notebook <- jsonlite::read_json(path, simplifyVector = FALSE)
  code_cells <- notebook$cells[
    vapply(notebook$cells, `[[`, character(1), "cell_type") == "code"
  ]
  stopifnot(length(code_cells) > 0L)
  stopifnot(all(vapply(code_cells, function(cell) !is.null(cell$execution_count), logical(1))))
  output_types <- unlist(lapply(code_cells, function(cell) {
    vapply(cell$outputs, `[[`, character(1), "output_type")
  }), use.names = FALSE)
  stopifnot(!any(output_types == "error"))
}

output_dir <- file.path(analysis_root, "mayassi_unrolling_notebook", "Region_6")
trace_audit_path <- file.path(output_dir, "Region_6_trace_definition_audit.tsv")
projection_path <- file.path(output_dir, "Region_6_passQC_unrolled_subset_cells.tsv.gz")
stopifnot(file.exists(trace_audit_path), file.exists(projection_path))

trace_audit <- utils::read.delim(trace_audit_path, check.names = FALSE)
stopifnot(trace_audit$input_cells == 127935L)
stopifnot(trace_audit$reference_smooth_muscle_cells == 8059L)
stopifnot(trace_audit$control_points == 39L)

projected <- utils::read.delim(gzfile(projection_path), check.names = FALSE)
stopifnot("smooth_muscle_reference" %in% names(projected))
stopifnot(!"smooth_muscle_marker" %in% names(projected))
stopifnot(all(projected$smooth_muscle_reference %in% c(TRUE, FALSE)))

cat("Region 6 transparent trace workflow contract passed\n")
