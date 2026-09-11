#!/usr/bin/env Rscript

# Build the reader-facing notebook that documents and validates the manual
# Region 6 Swiss-roll trace before any cell projection is performed.

suppressPackageStartupMessages(library(jsonlite))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(script_path))
args <- commandArgs(trailingOnly = TRUE)
output_path <- if (length(args)) args[[1L]] else file.path(
  repo_root, "notebooks", "B3a_Region6_trace_definition.ipynb"
)

empty_object <- function() setNames(list(), character())
cells <- list()
add_markdown <- function(text, id) {
  cells[[length(cells) + 1L]] <<- list(
    cell_type = "markdown", id = id, metadata = empty_object(),
    source = list(paste0(text, "\n"))
  )
}
add_code <- function(text, id) {
  cells[[length(cells) + 1L]] <<- list(
    cell_type = "code", execution_count = NULL, id = id,
    metadata = empty_object(), outputs = list(),
    source = list(paste0(text, "\n"))
  )
}

add_markdown(
"# Region 6 manual trace definition for Mayassi-style unrolling

This notebook makes the previously implicit manual annotation step explicit. It must be reviewed before `B3_Region6_Mayassi_unrolling_method_audit.ipynb` is run.

**Input population:** all cells in `Region_6_spatial_passQC.rds`, which is the returned object restricted to `primary_include_revised == TRUE`.

**Smooth-muscle guide:** `RefAll_subtype_predicted.id == \"Smooth muscle\"`. `Xenium_cluster_subtype` is not used.",
"b3a-001-title")

add_markdown(
"## What was manual and what is reproducible

The 39 control points were manually entered during the Region 6 feasibility review by visually following the outer band of cells transferred as Mayassi-reference `Smooth muscle`, from the outer free edge through successive turns to the inner hook. Coordinates were rounded to practical plotting positions; they were not inferred by an optimization algorithm and there is no hidden click log.

Exact computational replication is achieved by version-controlling those coordinates. Independent biological replication requires reviewing the numbered overlay below and editing a copied candidate TSV when the path does not follow the intended muscle/serosal boundary.",
"b3a-002-provenance")

add_code(
"suppressPackageStartupMessages({
  library(SeuratObject)
  library(ggplot2)
})
options(stringsAsFactors = FALSE, repr.plot.width = 11, repr.plot.height = 10)
cat('R version:', R.version.string, '\n')",
"b3a-003-packages")

add_code(
"PROJECT_ROOT <- Sys.getenv(
  'COLON_PROJECT_ROOT',
  unset = if (.Platform$OS.type == 'windows') {
    'D:/Xiaonan/CODEX_projects/Yanan_Xenium'
  } else {
    '/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium'
  }
)
ANALYSIS_ROOT <- file.path(PROJECT_ROOT, 'colon_analysis')
PIPELINE_REPO <- file.path(ANALYSIS_ROOT, 'YNH_Xenium_Colon')
REGION_SPATIAL_RDS <- file.path(ANALYSIS_ROOT, 'HPC_return', '20260910_HPC_return', 'colon_downstream_outputs', 'full_notebook_qc_v2', '03_Region_6_Primary479', 'Region_6_spatial_passQC.rds')
TRACE_CONFIG <- file.path(PIPELINE_REPO, 'config', 'unrolling', 'Region_6_trace_control_points.tsv')
OUT_DIR <- file.path(ANALYSIS_ROOT, 'mayassi_unrolling_notebook', 'Region_6')
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
print(data.frame(name = c('PROJECT_ROOT','REGION_SPATIAL_RDS','TRACE_CONFIG','OUT_DIR'), path = c(PROJECT_ROOT,REGION_SPATIAL_RDS,TRACE_CONFIG,OUT_DIR)))",
"b3a-004-paths")

add_code(
"if (.Platform$OS.type == 'windows' && grepl('^[Cc]:', normalizePath(tempdir(), winslash = '/'))) stop('R temporary files are on C:. Set TMPDIR, TEMP, and TMP below colon_analysis/tmp.')
required_paths <- c(PIPELINE_REPO, REGION_SPATIAL_RDS, TRACE_CONFIG)
if (any(!file.exists(required_paths))) stop('Missing required path(s): ', paste(required_paths[!file.exists(required_paths)], collapse = '; '))
source(file.path(PIPELINE_REPO, 'R', 'source.R'))
cat('Environment and path checks: PASS\n')",
"b3a-005-preflight")

add_markdown(
"## Step 1 - Load the exact pass-QC population

No QC mask is recalculated here. The returned Seurat object is the starting population requested for unrolling. The notebook verifies that every retained cell carries `primary_include_revised == TRUE`.",
"b3a-006-input-method")

add_code(
"region_spatial <- readRDS(REGION_SPATIAL_RDS)
spatial_meta <- region_spatial[[]]
required_metadata <- c('x_centroid','y_centroid','primary_include_revised','RefAll_subtype_predicted.id')
require_columns(spatial_meta, required_metadata, 'Region 6 pass-QC metadata')
reference_smooth <- as.character(spatial_meta$RefAll_subtype_predicted.id) == 'Smooth muscle'
input_audit <- data.frame(
  input_cells = nrow(spatial_meta),
  revised_primary_true = sum(spatial_meta$primary_include_revised %in% TRUE),
  revised_primary_not_true = sum(!(spatial_meta$primary_include_revised %in% TRUE)),
  reference_smooth_muscle_cells = sum(reference_smooth, na.rm = TRUE),
  missing_reference_labels = sum(is.na(spatial_meta$RefAll_subtype_predicted.id))
)
print(input_audit)
stopifnot(input_audit$revised_primary_not_true == 0L)",
"b3a-007-input-audit")

add_markdown(
"## Step 2 - Load and validate the version-controlled manual coordinates

The TSV is the immutable input for exact replication. `point_order` defines the direction from outer free edge to inner endpoint. The comments identify the two deliberate inter-turn bridges and the inner hook. All 39 rows are printed so there is no hidden coordinate state.",
"b3a-008-controls-method")

add_code(
"control_points <- utils::read.delim(TRACE_CONFIG, check.names = FALSE)
control_audit <- validate_trace_control_points(control_points)
print(control_audit)
print(control_points, row.names = FALSE)",
"b3a-009-controls-table")

add_code(
"segment_audit <- data.frame(
  from_point = head(control_points$point_order, -1),
  to_point = tail(control_points$point_order, -1),
  segment_um = sqrt(diff(control_points$x)^2 + diff(control_points$y)^2),
  destination_comment = tail(control_points$comment, -1)
)
print(segment_audit[order(segment_audit$segment_um, decreasing = TRUE), ], row.names = FALSE)",
"b3a-010-segments")

add_markdown(
"## Step 3 - Review every numbered point against the reference-transfer guide

Grey points are all pass-QC cells. Blue cells have `RefAll_subtype_predicted.id == \"Smooth muscle\"`. The green polyline is the proposed ordered boundary and red labels are the exact control-point numbers.

Review criteria: follow one continuous outer muscle/serosal band; do not jump to an adjacent coil merely because it is spatially close; inspect points 16 and 31 as explicit inter-turn bridges; and confirm that points 37-39 follow the inner hook.",
"b3a-011-overlay-method")

add_code(
"plot_data <- data.frame(
  x = as.numeric(spatial_meta$x_centroid),
  y = as.numeric(spatial_meta$y_centroid),
  reference_smooth = reference_smooth
)
numbered_trace_plot <- ggplot(plot_data, aes(x, y)) +
  geom_point(color = 'grey82', size = 0.05, alpha = 0.18) +
  geom_point(data = plot_data[plot_data$reference_smooth %in% TRUE, ], color = '#2166AC', size = 0.12, alpha = 0.65) +
  geom_path(data = control_points, aes(x, y), inherit.aes = FALSE, color = '#00A651', linewidth = 0.65) +
  geom_point(data = control_points, aes(x, y), inherit.aes = FALSE, color = '#D73027', size = 1.2) +
  geom_text(data = control_points, aes(x, y, label = point_order), inherit.aes = FALSE, color = '#7F0000', size = 2.5, nudge_y = 55) +
  coord_equal() + theme_void() +
  labs(title = 'Region 6 manual trace review', subtitle = 'blue: Mayassi-reference Smooth muscle; green: ordered path; red numbers: editable control points')
ggsave(file.path(OUT_DIR, 'Region_6_numbered_trace_definition.png'), numbered_trace_plot, width = 10, height = 10, dpi = 220, bg = 'white')
numbered_trace_plot",
"b3a-012-overlay-plot")

add_markdown(
"![Numbered Region 6 manual trace](../../mayassi_unrolling_notebook/Region_6/Region_6_numbered_trace_definition.png)",
"b3a-013-overlay-image")

add_markdown(
"## Step 4 - Save an auditable review copy

The repository TSV remains the canonical input. This step writes a dated-run-independent review copy and an audit table beside the downstream projection outputs. Editing should be performed on a copied candidate TSV, visually reviewed, and only then promoted to the repository configuration.",
"b3a-014-save-method")

add_code(
"dense_trace <- add_trace_arc_length(interpolate_trace_control_points(control_points[, c('x','y')], spacing = 20))
trace_definition_audit <- data.frame(
  input_cells = nrow(spatial_meta),
  revised_primary_true = sum(spatial_meta$primary_include_revised %in% TRUE),
  reference_smooth_muscle_cells = sum(reference_smooth, na.rm = TRUE),
  control_points = nrow(control_points),
  dense_trace_points = nrow(dense_trace),
  dense_trace_arc_length_um = max(dense_trace$roll_arc_length),
  config_md5 = unname(tools::md5sum(TRACE_CONFIG))
)
utils::write.table(control_points, file.path(OUT_DIR, 'Region_6_trace_control_points.reviewed.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
utils::write.table(trace_definition_audit, file.path(OUT_DIR, 'Region_6_trace_definition_audit.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
print(trace_definition_audit)",
"b3a-015-save-audit")

add_markdown(
"## Decision before B3

This notebook provides exact computational provenance, but it does not convert a manual biological annotation into an objective ground truth. If the numbered line crosses an unsupported tissue gap or follows the wrong coil, stop and revise a candidate trace before running B3. The current trace remains provisional until histology or expert review confirms it.",
"b3a-016-decision")

add_code("sessionInfo()", "b3a-017-session")

notebook <- list(
  cells = cells,
  metadata = list(
    kernelspec = list(display_name = "R", language = "R", name = "ir"),
    language_info = list(name = "R")
  ),
  nbformat = 4L,
  nbformat_minor = 5L
)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(notebook, output_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
cat("Notebook written to", normalizePath(output_path, winslash = "/", mustWork = TRUE), "\n")
