#!/usr/bin/env Rscript

# Generate the seven committed colon QC notebooks from canonical cell lists.
# The script writes only below the repository supplied as its first argument.
suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args)) normalizePath(args[[1L]], winslash = "/", mustWork = TRUE) else normalizePath(".", winslash = "/", mustWork = TRUE)
notebook_dir <- file.path(repo_root, "notebooks")
dir.create(notebook_dir, recursive = TRUE, showWarnings = FALSE)

markdown_cell <- function(text) list(cell_type = "markdown", metadata = list(), source = text)
code_cell <- function(text) list(cell_type = "code", execution_count = NULL, metadata = list(), outputs = list(), source = text)
notebook <- function(cells) list(
  cells = cells,
  metadata = list(
    kernelspec = list(display_name = "R", language = "R", name = "ir"),
    language_info = list(name = "R", mimetype = "text/x-r-source", file_extension = ".r")
  ),
  nbformat = 4L, nbformat_minor = 5L
)

region_cells <- function(region_id) {
  region_number <- sub("Region_", "", region_id, fixed = TRUE)
  list(
    markdown_cell(sprintf(
      "# Colon Xenium QC - %s\n\n## Goal\n\nRun reproducible, non-destructive technical QC for **%s**. This notebook is one of six structurally identical region notebooks; its region identifier is locked below. Full-data results must be produced on HPC. Local runs use a deterministic small subset only and cannot establish transcript-level or final slide readiness.",
      region_id, region_id
    )),
    markdown_cell("## Setup\n\nThe cell below defines HPC defaults and accepts explicit environment overrides for local subset validation. It refuses writable paths outside `colon_analysis`. The pipeline installs nothing; missing packages cause a clear preflight stop."),
    code_cell(sprintf(
      "REGION_ID <- \"%s\"\nEXECUTION_MODE <- toupper(Sys.getenv(\"COLON_QC_MODE\", \"FULL_HPC\"))\nstopifnot(EXECUTION_MODE %%in%% c(\"FULL_HPC\", \"LOCAL_SUBSET\"))\nPROJECT_ROOT <- Sys.getenv(\"COLON_PROJECT_ROOT\", \"/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium\")\nPIPELINE_REPO <- Sys.getenv(\"COLON_PIPELINE_REPO\", file.path(PROJECT_ROOT, \"colon_analysis\", \"YNH_Xenium_Colon\"))\nINPUT_ROOT <- Sys.getenv(\"COLON_INPUT_ROOT\", file.path(PROJECT_ROOT, \"colon_data\"))\nRUN_LABEL <- Sys.getenv(\"COLON_RUN_LABEL\", \"colon_qc_hpc\")\nRUN_ROOT <- Sys.getenv(\"COLON_RUN_ROOT\", file.path(PROJECT_ROOT, \"colon_analysis\", \"colon_qc_outputs\", RUN_LABEL))\nTEMP_ROOT <- Sys.getenv(\"COLON_TEMP_ROOT\", file.path(PROJECT_ROOT, \"colon_analysis\", \"tmp\", RUN_LABEL))\ndir.create(TEMP_ROOT, recursive = TRUE, showWarnings = FALSE)\nSys.setenv(TMPDIR = TEMP_ROOT, TMP = TEMP_ROOT, TEMP = TEMP_ROOT)\nsource(file.path(PIPELINE_REPO, \"R\", \"source.R\"))\nset.seed(20260814L)", region_id
    )),
    markdown_cell("### Key assumptions\n\nAll six regions use one installed Xenium panel. Mouse is the biological replicate; this notebook performs technical QC only. The primary cell cohort is prespecified as `5 < nFeature_Xenium < 200` and `10 < nCount_Xenium < 1000`. Segmentation, nucleus, area, control, and spatial findings remain separate review flags."),
    code_cell("runtime <- validate_runtime_paths(project_root = PROJECT_ROOT, input_root = INPUT_ROOT, output_root = RUN_ROOT, temp_root = TEMP_ROOT)\nregions <- discover_xenium_sections(INPUT_ROOT, expected_section_count = 6L)\nstopifnot(identical(regions$region_id, expected_colon_regions()))\nmanifest <- utils::read.delim(file.path(PIPELINE_REPO, \"config\", \"colon_sample_manifest.tsv\"), check.names = FALSE)\nvalidate_sample_manifest(manifest, expected_colon_regions())\nfixed_thresholds <- read_fixed_cell_qc_thresholds(file.path(PIPELINE_REPO, \"config\", \"fixed_cell_qc_thresholds.tsv\"))\nmanifest[manifest$region_id == REGION_ID, , drop = FALSE]"),
    markdown_cell("## Inputs and integrity\n\nDiscover the region directory from its identifier, inventory required files, and verify matrix/metadata alignment before loading counts. Raw data are read-only. Full transcript Parquet is not collected into memory; full-HPC transcript summaries use Arrow-backed aggregation."),
    code_cell("region_record <- discover_one_section(INPUT_ROOT, REGION_ID)\nregion_dir <- region_record$region_dir[[1L]]\ninventory <- inventory_section_files(region_dir, REGION_ID, calculate_md5 = FALSE)\nintegrity <- validate_section_integrity(region_dir, REGION_ID)\nstopifnot(all(inventory$exists), isTRUE(integrity$dimension_match[[1L]]))\nbundle <- import_xenium_mex(region_dir)\nstopifnot(inherits(bundle$counts, \"sparseMatrix\"), identical(colnames(bundle$counts), bundle$cells$cell_id))"),
    markdown_cell("## Cell QC\n\nCompute targeted-panel complexity, control burden, and segmentation review fields without deleting cells. Exact boundary values 5/200 features and 10/1000 counts fail. `primary_include` contains only the fixed rule; `strict_include` additionally excludes all review flags."),
    code_cell("cell_qc <- calculate_xenium_cell_qc(bundle$counts, bundle$cells, REGION_ID, fixed_thresholds)\nmasks <- build_cell_downstream_masks(cell_qc$cell_metadata, provenance = paste(RUN_LABEL, EXECUTION_MODE, sep = \"::\"))\nstopifnot(nrow(masks) == ncol(bundle$counts), all(!masks$strict_include | masks$primary_include))\ntable(primary_include = masks$primary_include, strict_include = masks$strict_include)"),
    markdown_cell("## Outputs and checks\n\nWrite an auditable region bundle beneath the run root: sparse raw counts, all cell metadata and masks, thresholds, inventory/integrity, run configuration, readiness evidence, and session information. The source notebook stays output-free; executed copies belong under the run directory."),
    code_cell("section_output_dir <- file.path(RUN_ROOT, \"sections\", REGION_ID)\nresult <- write_colon_region_qc_bundle(\n  project_root = PROJECT_ROOT, output_dir = section_output_dir, region_id = REGION_ID,\n  run_label = RUN_LABEL, execution_mode = EXECUTION_MODE, manifest = manifest,\n  inventory = inventory, integrity = integrity, bundle = bundle, cell_qc = cell_qc, masks = masks\n)\nstopifnot(validate_colon_region_qc_bundle(section_output_dir, REGION_ID))\nresult"),
    markdown_cell("## Next steps\n\nReview region tables and figures, especially assignment/control burden, boundary counts, segmentation flags, and spatial hotspots. A `PASS` here is technical evidence only. After all six regions finish under the same run label and mode, execute `02_slide_QC_summary.ipynb`. Full-data HPC output—not this local subset—determines final readiness.")
  )
}

summary_cells <- list(
  markdown_cell("# Colon Xenium slide QC summary\n\n## Goal\n\nCombine exactly six validated region bundles from one run label and execution mode. The summary compares technical QC across Mouse_1 and Mouse_2 at top, middle, and bottom positions without treating cells as biological replicates."),
  markdown_cell("## Setup\n\nHPC paths are defaults. Local validation injects D-drive paths and `LOCAL_SUBSET`. No package installation or writes outside `colon_analysis` are permitted."),
  code_cell("EXECUTION_MODE <- toupper(Sys.getenv(\"COLON_QC_MODE\", \"FULL_HPC\"))\nPROJECT_ROOT <- Sys.getenv(\"COLON_PROJECT_ROOT\", \"/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium\")\nPIPELINE_REPO <- Sys.getenv(\"COLON_PIPELINE_REPO\", file.path(PROJECT_ROOT, \"colon_analysis\", \"YNH_Xenium_Colon\"))\nRUN_LABEL <- Sys.getenv(\"COLON_RUN_LABEL\", \"colon_qc_hpc\")\nRUN_ROOT <- Sys.getenv(\"COLON_RUN_ROOT\", file.path(PROJECT_ROOT, \"colon_analysis\", \"colon_qc_outputs\", RUN_LABEL))\nTEMP_ROOT <- Sys.getenv(\"COLON_TEMP_ROOT\", file.path(PROJECT_ROOT, \"colon_analysis\", \"tmp\", RUN_LABEL))\ndir.create(TEMP_ROOT, recursive = TRUE, showWarnings = FALSE)\nSys.setenv(TMPDIR = TEMP_ROOT, TMP = TEMP_ROOT, TEMP = TEMP_ROOT)\nsource(file.path(PIPELINE_REPO, \"R\", \"source.R\"))\nregions <- expected_colon_regions()"),
  markdown_cell("## Inputs and integrity\n\nRequire Region_1 through Region_6, reject missing/duplicate regions and mixed run labels or modes, then reload every saved sparse object and table."),
  code_cell("manifest <- utils::read.delim(file.path(PIPELINE_REPO, \"config\", \"colon_sample_manifest.tsv\"), check.names = FALSE)\nvalidate_sample_manifest(manifest, regions)\nslide_data <- read_colon_slide_qc_outputs(RUN_ROOT, expected_regions = regions, run_label = RUN_LABEL, execution_mode = EXECUTION_MODE)\nstopifnot(identical(slide_data$coverage$region_id, regions))"),
  markdown_cell("## Cell QC\n\nSummarize fixed-bound pass fractions and review flags. `primary_include` is never redefined at slide level. Region readiness is derived from current evidence rather than inherited anchor or sensitivity labels."),
  code_cell("slide_summary <- summarise_slide_qc(slide_data)\nposition_summary <- summarise_colon_mouse_position(slide_summary$section_summary, manifest)\nposition_summary$within_mouse\nposition_summary$matched_position"),
  markdown_cell("## Outputs and checks\n\nSave combined tables, descriptive within-mouse and matched-position summaries, readiness gates, figures, session information, and a reloadable slide object."),
  code_cell("summary_output_dir <- file.path(RUN_ROOT, \"slide_summary\")\nresult <- write_colon_slide_qc_bundle(PROJECT_ROOT, summary_output_dir, RUN_LABEL, EXECUTION_MODE, slide_data, slide_summary, position_summary)\nstopifnot(validate_colon_slide_qc_bundle(summary_output_dir, expected_regions = regions))\nresult"),
  markdown_cell("## Next steps\n\nUse full-HPC evidence for the final QC decision. With only two mice, anatomical-position results are descriptive. The adipose-derived eosinophil lists remain hypotheses (`AT-short-like` and `AT-long-like`) until colon-specific coherence and identity checks are completed in a later, study-specific analysis stage.")
)

for (region_id in paste0("Region_", 1:6)) {
  number <- sub("Region_", "", region_id, fixed = TRUE)
  jsonlite::write_json(notebook(region_cells(region_id)), file.path(notebook_dir, sprintf("01_QC_Region%s.ipynb", number)), auto_unbox = TRUE, pretty = TRUE, null = "null")
}
jsonlite::write_json(notebook(summary_cells), file.path(notebook_dir, "02_slide_QC_summary.ipynb"), auto_unbox = TRUE, pretty = TRUE, null = "null")
cat("Generated six region notebooks and one slide summary.\n")
