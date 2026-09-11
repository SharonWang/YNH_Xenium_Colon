options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(jsonlite))

read_notebook <- function(path) jsonlite::read_json(path, simplifyVector = FALSE)
cell_text <- function(notebook) paste(vapply(notebook$cells, function(cell) paste(unlist(cell$source), collapse = ""), character(1)), collapse = "\n")
assert_valid_execution_state <- function(notebook) {
  code <- Filter(function(cell) identical(cell$cell_type, "code"), notebook$cells)
  has_execution_count <- vapply(
    code, function(cell) !is.null(cell$execution_count), logical(1)
  )
  if (!any(has_execution_count)) {
    stopifnot(all(vapply(code, function(cell) length(cell$outputs) == 0L, logical(1))))
    return(invisible(TRUE))
  }
  # HPC-returned notebooks may be committed with outputs. Require a complete,
  # non-mixed execution state and reject any recorded cell error.
  stopifnot(all(has_execution_count))
  output_types <- unlist(lapply(code, function(cell) {
    vapply(cell$outputs, `[[`, character(1), "output_type")
  }), use.names = FALSE)
  stopifnot(!any(output_types == "error"))
  invisible(TRUE)
}

regions <- paste0("Region_", 1:6)
headings <- c("## Goal", "## Setup", "## Inputs and integrity", "## Cell QC", "## Outputs and checks", "## Next steps")
for (index in seq_along(regions)) {
  path <- file.path("notebooks", sprintf("01_QC_Region%d.ipynb", index))
  stopifnot(file.exists(path))
  raw_json <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # Jupyter nbformat requires every cell metadata field to be a JSON object.
  # jsonlite serializes an unnamed empty R list as [], which Jupyter rejects.
  stopifnot(!grepl('"metadata": []', raw_json, fixed = TRUE))
  notebook <- read_notebook(path)
  text <- cell_text(notebook)
  cell_ids <- vapply(notebook$cells, function(cell) if (is.null(cell$id)) NA_character_ else cell$id, character(1))
  stopifnot(!anyNA(cell_ids), !anyDuplicated(cell_ids), all(grepl("^[A-Za-z0-9_-]{1,64}$", cell_ids)))
  # R kernel IDs are environment-specific: the HPC-returned notebooks use
  # `ir542`, whereas a portable IRkernel installation commonly uses `ir`.
  # Require an IRkernel-style ID and R language instead of one site ID.
  stopifnot(grepl("^ir[[:alnum:]_.-]*$", notebook$metadata$kernelspec$name))
  stopifnot(identical(notebook$metadata$kernelspec$language, "R"))
  stopifnot(grepl(sprintf('REGION_ID <- "%s"', regions[[index]]), text, fixed = TRUE))
  stopifnot(all(vapply(headings, grepl, logical(1), x = text, fixed = TRUE)))
  stopifnot(grepl("FULL_HPC", text, fixed = TRUE), grepl("LOCAL_SUBSET", text, fixed = TRUE))
  stopifnot(grepl("fixed_cell_qc_thresholds.tsv", text, fixed = TRUE), grepl("primary_include", text, fixed = TRUE))
  stopifnot(grepl("Full-data HPC output - not this local subset - determines final readiness.", text, fixed = TRUE))
  stopifnot(!grepl("scWAT", text, fixed = TRUE), !grepl("adipose_analysis", text, fixed = TRUE))
  assert_valid_execution_state(notebook)
}

summary <- read_notebook(file.path("notebooks", "02_slide_QC_summary.ipynb"))
summary_raw_json <- paste(readLines(file.path("notebooks", "02_slide_QC_summary.ipynb"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
stopifnot(!grepl('"metadata": []', summary_raw_json, fixed = TRUE))
summary_text <- cell_text(summary)
summary_cell_ids <- vapply(summary$cells, function(cell) if (is.null(cell$id)) NA_character_ else cell$id, character(1))
stopifnot(!anyNA(summary_cell_ids), !anyDuplicated(summary_cell_ids), all(grepl("^[A-Za-z0-9_-]{1,64}$", summary_cell_ids)))
stopifnot(grepl("expected_colon_regions()", summary_text, fixed = TRUE))
stopifnot(grepl("summarise_colon_mouse_position", summary_text, fixed = TRUE))
stopifnot(all(vapply(c("Mouse_1", "Mouse_2", "top", "middle", "bottom"), grepl, logical(1), x = summary_text, fixed = TRUE)))
stopifnot(!grepl("scWAT", summary_text, fixed = TRUE), !grepl("adipose_analysis", summary_text, fixed = TRUE))
assert_valid_execution_state(summary)

cat("Seven colon notebook contracts passed.\n")

launcher_files <- c(
  file.path("scripts", "execute_local_subset.ps1"),
  file.path("scripts", "validate_notebooks.py"),
  file.path("shell", "run_notebook_qc_hpc.sh"),
  file.path("slurm", "colon_notebook_qc.sbatch"),
  "README.md"
)
stopifnot(all(file.exists(launcher_files)))
launcher_text <- paste(vapply(launcher_files, function(path) paste(readLines(path, warn = FALSE), collapse = "\n"), character(1)), collapse = "\n")
stopifnot(grepl("D:\\\\Xiaonan\\\\CODEX_projects\\\\Yanan_Xenium\\\\colon_analysis", launcher_text))
stopifnot(grepl("/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium/colon_analysis", launcher_text, fixed = TRUE))
stopifnot(all(vapply(regions, grepl, logical(1), x = launcher_text, fixed = TRUE)))
stopifnot(all(vapply(c("Matrix", "jsonlite", "ggplot2", "arrow", "dplyr", "RANN"), grepl, logical(1), x = launcher_text, fixed = TRUE)))
stopifnot(!grepl("adipose_analysis", launcher_text, fixed = TRUE))
cat("Local and HPC launcher contracts passed.\n")
