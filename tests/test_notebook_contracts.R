options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(jsonlite))

read_notebook <- function(path) jsonlite::read_json(path, simplifyVector = FALSE)
cell_text <- function(notebook) paste(vapply(notebook$cells, function(cell) paste(unlist(cell$source), collapse = ""), character(1)), collapse = "\n")
assert_cleared <- function(notebook) {
  code <- Filter(function(cell) identical(cell$cell_type, "code"), notebook$cells)
  stopifnot(all(vapply(code, function(cell) is.null(cell$execution_count), logical(1))))
  stopifnot(all(vapply(code, function(cell) length(cell$outputs) == 0L, logical(1))))
}

regions <- paste0("Region_", 1:6)
headings <- c("## Goal", "## Setup", "## Inputs and integrity", "## Cell QC", "## Outputs and checks", "## Next steps")
for (index in seq_along(regions)) {
  path <- file.path("notebooks", sprintf("01_QC_Region%d.ipynb", index))
  stopifnot(file.exists(path))
  notebook <- read_notebook(path)
  text <- cell_text(notebook)
  stopifnot(identical(notebook$metadata$kernelspec$name, "ir"))
  stopifnot(grepl(sprintf('REGION_ID <- "%s"', regions[[index]]), text, fixed = TRUE))
  stopifnot(all(vapply(headings, grepl, logical(1), x = text, fixed = TRUE)))
  stopifnot(grepl("FULL_HPC", text, fixed = TRUE), grepl("LOCAL_SUBSET", text, fixed = TRUE))
  stopifnot(grepl("fixed_cell_qc_thresholds.tsv", text, fixed = TRUE), grepl("primary_include", text, fixed = TRUE))
  stopifnot(!grepl("scWAT", text, fixed = TRUE), !grepl("adipose_analysis", text, fixed = TRUE))
  assert_cleared(notebook)
}

summary <- read_notebook(file.path("notebooks", "02_slide_QC_summary.ipynb"))
summary_text <- cell_text(summary)
stopifnot(grepl("expected_colon_regions()", summary_text, fixed = TRUE))
stopifnot(grepl("summarise_colon_mouse_position", summary_text, fixed = TRUE))
stopifnot(all(vapply(c("Mouse_1", "Mouse_2", "top", "middle", "bottom"), grepl, logical(1), x = summary_text, fixed = TRUE)))
stopifnot(!grepl("scWAT", summary_text, fixed = TRUE), !grepl("adipose_analysis", summary_text, fixed = TRUE))
assert_cleared(summary)

cat("Seven colon notebook contracts passed.\n")

launcher_files <- c(
  file.path("scripts", "execute_local_subset.ps1"),
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
