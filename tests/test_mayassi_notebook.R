#!/usr/bin/env Rscript

# Structural and executable-code contract for the stepwise method-audit notebook.
# Full data execution is verified separately with execute_r_notebook.R.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(jsonlite))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
test_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(test_path))
notebook_path <- file.path(repo_root, "notebooks", "B3_Region6_Mayassi_unrolling_method_audit.ipynb")

if (!file.exists(notebook_path)) {
  stop("B3 Region 6 Mayassi method-audit notebook is missing", call. = FALSE)
}

notebook <- jsonlite::read_json(notebook_path, simplifyVector = FALSE)
stopifnot(identical(notebook$nbformat, 4L))
stopifnot(length(notebook$cells) >= 30L)

cell_types <- vapply(notebook$cells, `[[`, character(1), "cell_type")
stopifnot(sum(cell_types == "markdown") >= 12L)
stopifnot(sum(cell_types == "code") >= 15L)

cell_ids <- vapply(notebook$cells, `[[`, character(1), "id")
stopifnot(all(nzchar(cell_ids)), !anyDuplicated(cell_ids))

# The committed reader-facing copy is executed locally. The lightweight runner
# must emulate Jupyter auto-printing without leaking its internal withVisible()
# wrapper as duplicated `$value` / `$visible` output.
code_cells <- notebook$cells[cell_types == "code"]
stopifnot(all(vapply(code_cells, function(cell) !is.null(cell$execution_count), logical(1))))
stream_text <- unlist(lapply(code_cells, function(cell) {
  outputs <- cell$outputs
  unlist(lapply(outputs, function(output) {
    if (identical(output$output_type, "stream")) unlist(output$text) else character()
  }), use.names = FALSE)
}), use.names = FALSE)
stopifnot(!any(grepl("^\\$value|^\\$visible", stream_text)))

for (index in seq_along(notebook$cells)) {
  cell <- notebook$cells[[index]]
  source_text <- paste(unlist(cell$source), collapse = "")
  if (!nzchar(trimws(source_text))) {
    stop(sprintf("Notebook cell %d is empty", index), call. = FALSE)
  }
  if (identical(cell$cell_type, "code")) {
    parse(text = source_text, keep.source = TRUE)
  }
  if (identical(cell$cell_type, "code") && length(cell$outputs)) {
    output_types <- vapply(cell$outputs, `[[`, character(1), "output_type")
    if (any(output_types == "error")) {
      stop(sprintf("Notebook cell %d contains an execution error", index), call. = FALSE)
    }
  }
}

cat("B3 Mayassi notebook structure and R syntax contracts passed\n")
