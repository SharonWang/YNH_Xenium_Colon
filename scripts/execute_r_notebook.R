#!/usr/bin/env Rscript

# Execute R code cells sequentially without modifying the committed notebook.
# This lightweight runner avoids a Jupyter dependency while retaining the
# notebook as the reader-facing source. An executed copy is written only to
# the explicitly supplied output path.
suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: execute_r_notebook.R INPUT.ipynb OUTPUT.ipynb", call. = FALSE)
input_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output_path <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)

analysis_token <- "/colon_analysis/"
if (!grepl(analysis_token, output_path, fixed = TRUE)) stop("Executed notebook output must be below colon_analysis.", call. = FALSE)
if (grepl("^[Cc]:/", output_path)) stop("Writing executed notebooks to C: is prohibited.", call. = FALSE)
if (grepl("adipose_analysis", output_path, fixed = TRUE)) stop("Colon execution cannot use adipose_analysis.", call. = FALSE)

notebook <- jsonlite::read_json(input_path, simplifyVector = FALSE)
execution_environment <- new.env(parent = globalenv())
execution_index <- 0L

for (cell_index in seq_along(notebook$cells)) {
  cell <- notebook$cells[[cell_index]]
  if (!identical(cell$cell_type, "code")) next
  execution_index <- execution_index + 1L
  source_code <- paste(unlist(cell$source), collapse = "")
  captured <- character()
  error_message <- NULL
  captured <- tryCatch(
    capture.output({
      evaluated <- withVisible(
        eval(parse(text = source_code, keep.source = TRUE), envir = execution_environment)
      )
      # Match an interactive R/Jupyter cell: side-effect output is retained and
      # only a visible final value is auto-printed. Do not print the internal
      # list returned by withVisible(), which duplicates every table.
      if (isTRUE(evaluated$visible)) print(evaluated$value)
    }, type = "output"),
    error = function(error) {
      error_message <<- conditionMessage(error)
      character()
    }
  )
  cell$execution_count <- execution_index
  cell$outputs <- list()
  if (length(captured)) {
    cell$outputs[[1L]] <- list(name = "stdout", output_type = "stream", text = paste0(captured, "\n"))
  }
  if (!is.null(error_message)) {
    cell$outputs[[length(cell$outputs) + 1L]] <- list(
      ename = "RExecutionError", evalue = error_message,
      output_type = "error", traceback = sprintf("Cell %d: %s", cell_index, error_message)
    )
    notebook$cells[[cell_index]] <- cell
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(notebook, output_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    stop(sprintf("Notebook failed in code cell %d: %s", cell_index, error_message), call. = FALSE)
  }
  notebook$cells[[cell_index]] <- cell
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(notebook, output_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
cat(sprintf("Executed notebook written to %s\n", output_path))
