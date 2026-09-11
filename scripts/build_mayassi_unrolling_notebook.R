#!/usr/bin/env Rscript

# Build the reader-facing Region 6 Mayassi method-audit notebook from one
# canonical specification. The generated notebook uses short cells, explicit
# Markdown methods, visible tables, and relative links to saved diagnostic plots.

suppressPackageStartupMessages(library(jsonlite))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(script_path))

args <- commandArgs(trailingOnly = TRUE)
output_path <- if (length(args) >= 1L) args[[1L]] else file.path(
  repo_root, "notebooks", "B3_Region6_Mayassi_unrolling_method_audit.ipynb"
)

empty_object <- function() setNames(list(), character())
cells <- list()

add_markdown <- function(text, id) {
  cells[[length(cells) + 1L]] <<- list(
    cell_type = "markdown",
    id = id,
    metadata = empty_object(),
    source = list(paste0(text, "\n"))
  )
}

add_code <- function(text, id) {
  cells[[length(cells) + 1L]] <<- list(
    cell_type = "code",
    execution_count = NULL,
    id = id,
    metadata = empty_object(),
    outputs = list(),
    source = list(paste0(text, "\n"))
  )
}

add_markdown(
"# Region 6 - Mayassi-style Swiss-roll unrolling method audit

This notebook reproduces, explains, and audits the Xenium unrolling strategy described by Mayassi et al. (Nature 2024). It is intentionally organized as short, inspectable steps in the style of `B2_Region6_primary_479.ipynb`.

**Purpose:** determine whether the published manual-trace plus global-radius projection produces trustworthy unrolled coordinates for the Region 6 pass-QC Seurat object. This is a feasibility and method-QC notebook, not a biological hypothesis test.

**Required upstream review:** run `B3a_Region6_trace_definition.ipynb` first. It documents the manual coordinates and produces the numbered trace overlay.

**Important:** `roll_arc_fraction` is anatomically unoriented. It must not be called proximal-distal until the physical rolling direction is confirmed.",
"b3-001-title")

add_markdown(
"## Result in one sentence

The smooth-muscle/serosal trace can be reconstructed, but the published global-radius constraint is considered **not analysis-ready** for this irregular roll when it causes substantial assignment changes, failed projections, or endpoint pile-up. The evidence is calculated below rather than assumed.",
"b3-002-outcome")

add_code(
"suppressPackageStartupMessages({
  library(SeuratObject)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})
options(stringsAsFactors = FALSE, repr.plot.width = 12, repr.plot.height = 8)
set.seed(20260910)
cat('R version:', R.version.string, '\n')",
"b3-003-packages")

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
HPC_RETURN <- file.path(ANALYSIS_ROOT, 'HPC_return', '20260910_HPC_return')
REGION_SPATIAL_RDS <- file.path(HPC_RETURN, 'colon_downstream_outputs', 'full_notebook_qc_v2', '03_Region_6_Primary479', 'Region_6_spatial_passQC.rds')
TRACE_TSV <- file.path(PIPELINE_REPO, 'config', 'unrolling', 'Region_6_trace_control_points.tsv')
OUT_DIR <- file.path(ANALYSIS_ROOT, 'mayassi_unrolling_notebook', 'Region_6')
MAX_CELLS <- as.integer(Sys.getenv('COLON_UNROLL_MAX_CELLS', unset = if (.Platform$OS.type == 'windows') '12000' else '0'))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
print(data.frame(name = c('PROJECT_ROOT','PIPELINE_REPO','REGION_SPATIAL_RDS','TRACE_TSV','OUT_DIR'), path = c(PROJECT_ROOT,PIPELINE_REPO,REGION_SPATIAL_RDS,TRACE_TSV,OUT_DIR)))",
"b3-004-paths")

add_code(
"if (.Platform$OS.type == 'windows' && grepl('^[Cc]:', normalizePath(tempdir(), winslash = '/'))) {
  stop('R temporary files are on C:. Set TMPDIR, TEMP, and TMP below colon_analysis/tmp.')
}
required_paths <- c(PIPELINE_REPO, REGION_SPATIAL_RDS, TRACE_TSV)
if (any(!file.exists(required_paths))) stop('Missing required path(s): ', paste(required_paths[!file.exists(required_paths)], collapse = '; '))
if (!is.finite(MAX_CELLS) || MAX_CELLS < 0L || (MAX_CELLS > 0L && MAX_CELLS < 100L)) stop('MAX_CELLS must be 0 for all cells or at least 100')
source(file.path(PIPELINE_REPO, 'R', 'source.R'))
cat('Environment and path checks: PASS\n')",
"b3-005-preflight")

add_markdown(
"## Step 1 - Load the requested Region 6 pass-QC population

This analysis starts directly from `Region_6_spatial_passQC.rds`. The returned object contains 127,935 cells retained by `primary_include_revised == TRUE`. No alternative primary mask is reconstructed and no second object is joined.",
"b3-006-populations-method")

add_code(
"region_spatial <- readRDS(REGION_SPATIAL_RDS)
spatial_meta <- region_spatial[[]]
spatial_meta$cell_id <- extract_cell_ids(spatial_meta)
required_metadata <- c('x_centroid','y_centroid','primary_include_revised','RefAll_subtype_predicted.id')
require_columns(spatial_meta, required_metadata, 'Region 6 pass-QC metadata')
print(data.frame(object = 'Region_6_spatial_passQC.rds', cells = nrow(spatial_meta), columns = ncol(spatial_meta)))",
"b3-007-load")

add_code(
"population_audit <- data.frame(
  input_cells = nrow(spatial_meta),
  primary_include_revised_true = sum(spatial_meta$primary_include_revised %in% TRUE),
  primary_include_revised_not_true = sum(!(spatial_meta$primary_include_revised %in% TRUE)),
  unique_cell_ids = length(unique(spatial_meta$cell_id))
)
print(population_audit)
stopifnot(population_audit$primary_include_revised_not_true == 0L, population_audit$unique_cell_ids == population_audit$input_cells)",
"b3-008-primary-audit")

add_code(
"reference_smooth <- as.character(spatial_meta$RefAll_subtype_predicted.id) == 'Smooth muscle'
reference_label_audit <- data.frame(
  input_cells = nrow(spatial_meta),
  reference_smooth_muscle = sum(reference_smooth, na.rm = TRUE),
  other_reference_label = sum(reference_smooth %in% FALSE),
  missing_reference_label = sum(is.na(reference_smooth))
)
print(reference_label_audit)",
"b3-009-population-table")

add_code(
"print(sort(table(spatial_meta$RefAll_subtype_predicted.id, useNA = 'ifany'), decreasing = TRUE))",
"b3-010-population-contract")

add_markdown(
"### Population contract for this notebook

For this unrolling analysis, `Region_6_spatial_passQC.rds` and `primary_include_revised` are intentionally authoritative, following the user-specified downstream population. Results from this notebook therefore describe those 127,935 retained cells, not all segmented cells and not the broader fixed-threshold population.",
"b3-011-primary-correction")

add_markdown(
"## Step 2 - Define exactly what \"smooth muscle\" means here

The sole cell-type guide is the cell-level Mayassi-reference transfer stored in `RefAll_subtype_predicted.id`. A cell is smooth muscle only when this field equals `Smooth muscle` exactly.

The corresponding continuous transfer score, `RefAll_subtype_prediction.score.Smooth.muscle`, is shown as supporting evidence. Cell type does not enter the projection equation; it only guides and validates the manual trace.",
"b3-012-celltype-method")

add_code(
"reference_definition <- data.frame(
  metadata_field = 'RefAll_subtype_predicted.id',
  required_value = 'Smooth muscle',
  n_cells = sum(reference_smooth, na.rm = TRUE),
  fraction_of_passQC = mean(reference_smooth, na.rm = TRUE)
)
print(reference_definition)",
"b3-013-marker-coverage")

add_code(
"score_columns <- c('RefAll_subtype_prediction.score.Smooth.muscle','RefAll_subtype_prediction.score.max')
require_columns(spatial_meta, score_columns, 'Region 6 reference-transfer metadata')
reference_score_audit <- data.frame(
  group = c('reference Smooth muscle','all other reference labels'),
  n = c(sum(reference_smooth %in% TRUE), sum(reference_smooth %in% FALSE)),
  median_smooth_muscle_score = c(
    median(spatial_meta$RefAll_subtype_prediction.score.Smooth.muscle[reference_smooth %in% TRUE], na.rm = TRUE),
    median(spatial_meta$RefAll_subtype_prediction.score.Smooth.muscle[reference_smooth %in% FALSE], na.rm = TRUE)
  ),
  median_max_transfer_score = c(
    median(spatial_meta$RefAll_subtype_prediction.score.max[reference_smooth %in% TRUE], na.rm = TRUE),
    median(spatial_meta$RefAll_subtype_prediction.score.max[reference_smooth %in% FALSE], na.rm = TRUE)
  )
)
print(reference_score_audit)",
"b3-014-definition-agreement")

add_code(
"reference_score_quantiles <- quantile(
  spatial_meta$RefAll_subtype_prediction.score.Smooth.muscle,
  probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.99, 1),
  na.rm = TRUE
)
print(reference_score_quantiles)",
"b3-015-cluster-definition")

add_code(
"definition_plot_data <- spatial_meta %>%
  transmute(
    x = x_centroid,
    y = y_centroid,
    reference_transfer = reference_smooth,
    smooth_transfer_score = RefAll_subtype_prediction.score.Smooth.muscle,
    max_transfer_score = RefAll_subtype_prediction.score.max
  )
p_reference <- ggplot(definition_plot_data, aes(x, y, color = reference_transfer)) + geom_point(size = 0.05, alpha = 0.35) + coord_equal() + scale_color_manual(values = c('FALSE'='grey85','TRUE'='#2166AC')) + theme_void() + ggtitle('Mayassi reference-transfer label')
score_limits <- quantile(definition_plot_data$smooth_transfer_score, c(0.50, 0.995), na.rm = TRUE)
p_score <- ggplot(definition_plot_data, aes(x, y, color = smooth_transfer_score)) + geom_point(size = 0.05, alpha = 0.45) + coord_equal() + scale_color_viridis_c(option = 'magma', limits = score_limits, oob = scales::squish) + theme_void() + ggtitle('Smooth-muscle transfer score')
p_confidence <- ggplot(definition_plot_data, aes(x, y, color = max_transfer_score)) + geom_point(size = 0.05, alpha = 0.45) + coord_equal() + scale_color_viridis_c(option = 'viridis', limits = c(0,1), oob = scales::squish) + theme_void() + ggtitle('Maximum transfer score')
celltype_plot <- p_reference + p_score + p_confidence + plot_layout(ncol = 3)
ggsave(file.path(OUT_DIR, 'Region_6_reference_smooth_muscle_definition.png'), celltype_plot, width = 18, height = 6, dpi = 200, bg = 'white')
celltype_plot",
"b3-016-celltype-plot")

add_markdown(
"![Mayassi-reference smooth-muscle definition and transfer scores](../../mayassi_unrolling_notebook/Region_6/Region_6_reference_smooth_muscle_definition.png)

The transfer label is a model-derived annotation, not histological ground truth. No cell-level significance test is valid because nearby Xenium cells are spatially correlated and are not biological replicates.",
"b3-017-celltype-image")

add_markdown(
"## Step 3 - Load the reviewed manual serosal/muscle trace

The manual construction and numbered-point review are performed in `B3a_Region6_trace_definition.ipynb`. B3 reads the version-controlled TSV from `config/unrolling` and validates its order, coordinates, and annotations before interpolation.

Consecutive control points are linearly interpolated every 20 um. Cumulative Euclidean distance along this ordered line gives `roll_arc_length`; dividing by total length gives `roll_arc_fraction` from 0 to 1.",
"b3-018-trace-method")

add_code(
"control_points <- read.delim(TRACE_TSV, check.names = FALSE)
control_audit <- validate_trace_control_points(control_points)
trace <- interpolate_trace_control_points(control_points[, c('x','y')], spacing = 20)
trace <- add_trace_arc_length(trace)
trace_center <- c(x = mean(trace$x), y = mean(trace$y))
trace$radius_from_center <- sqrt((trace$x - trace_center[['x']])^2 + (trace$y - trace_center[['y']])^2)
trace_quality <- assess_trace_quality(trace)
print(cbind(control_audit, dense_points = nrow(trace), arc_length_um = max(trace$roll_arc_length), continuous = trace_quality$continuous, large_jumps = trace_quality$n_large_jumps))",
"b3-019-build-trace")

add_code(
"trace_plot <- ggplot(definition_plot_data, aes(x, y, color = smooth_transfer_score)) +
  geom_point(size = 0.05, alpha = 0.40) +
  geom_path(data = trace, aes(x, y), inherit.aes = FALSE, color = '#00A651', linewidth = 0.55) +
  geom_point(data = control_points, aes(x, y), inherit.aes = FALSE, color = '#D73027', size = 0.65) +
  coord_equal() +
  scale_color_viridis_c(option = 'magma', limits = score_limits, oob = scales::squish) +
  theme_void() +
  ggtitle('Region 6 manual trace on Mayassi smooth-muscle transfer score', subtitle = 'green: dense trace; red: ordered control points')
ggsave(file.path(OUT_DIR, 'Region_6_trace_on_smooth_score.png'), trace_plot, width = 9, height = 9, dpi = 220, bg = 'white')
trace_plot",
"b3-020-trace-plot")

add_markdown(
"![Manual trace over the Mayassi smooth-muscle transfer score](../../mayassi_unrolling_notebook/Region_6/Region_6_trace_on_smooth_score.png)

Visual review is indispensable: a mathematically smooth trace can still be biologically wrong if it crosses unsupported tissue or jumps between coils.",
"b3-021-trace-image")

add_markdown(
"## Step 4 - Reproduce the Mayassi projection

For each cell, Euclidean distance is calculated to every dense trace point. The published constraint computes a single global center as the mean trace position and excludes trace points whose center-distance is greater than the cell's center-distance. The nearest remaining point is accepted.

The longitudinal coordinate is the accepted point's cumulative arc distance. The wall-depth coordinate is the unsigned cell-to-trace distance in micrometres. An unconstrained nearest-trace projection is also calculated strictly as a diagnostic; disagreement quantifies sensitivity to the published radial rule.",
"b3-022-projection-method")

add_code(
"passqc_coordinates <- data.frame(
  cell_id = as.character(spatial_meta$cell_id),
  x = as.numeric(spatial_meta$x_centroid),
  y = as.numeric(spatial_meta$y_centroid),
  smooth_muscle_reference = reference_smooth,
  smooth_muscle_reference_score = as.numeric(spatial_meta$RefAll_subtype_prediction.score.Smooth.muscle)
)
passqc_coordinates <- passqc_coordinates[complete.cases(passqc_coordinates[, c('cell_id','x','y','smooth_muscle_reference')]), , drop = FALSE]
set.seed(20260910)
selected <- if (MAX_CELLS > 0L && nrow(passqc_coordinates) > MAX_CELLS) sort(sample.int(nrow(passqc_coordinates), MAX_CELLS)) else seq_len(nrow(passqc_coordinates))
cells_pilot <- passqc_coordinates[selected, , drop = FALSE]
print(data.frame(passQC_input_cells = nrow(passqc_coordinates), projected_run_cells = nrow(cells_pilot), run_mode = if (MAX_CELLS == 0L) 'all cells' else 'deterministic subset'))",
"b3-023-primary-sample")

add_code(
"projected <- project_cells_to_trace(
  cells_pilot,
  trace,
  inward_constraint = TRUE,
  center = trace_center
)
cat('Published constrained projection complete for', nrow(projected), 'cells\n')",
"b3-024-constrained")

add_code(
"unconstrained <- project_cells_to_trace(
  cells_pilot,
  trace,
  inward_constraint = FALSE,
  center = trace_center
)
projected$unconstrained_trace_index <- unconstrained$trace_index
projected$unconstrained_roll_arc_fraction <- unconstrained$roll_arc_fraction
projected$unconstrained_wall_distance <- unconstrained$wall_distance
projected$projection_changed_by_constraint <- projected$trace_index != projected$unconstrained_trace_index
cat('Unconstrained diagnostic projection complete\n')",
"b3-025-unconstrained")

add_markdown(
"## Step 5 - Quantify whether the transformation is trustworthy

These are descriptive geometry/QC summaries, not inferential statistics. Pre-specified warning gates for this feasibility notebook are:

- projection success below 95%;
- more than 20% of assignments changed by the global-radius constraint;
- more than 5% of assigned cells falling in the first or last 1% of the trace.

The thresholds are engineering review gates, not biological significance cut-offs.",
"b3-026-diagnostic-method")

add_code(
"valid <- !is.na(projected$trace_index)
valid_unconstrained <- !is.na(projected$unconstrained_trace_index)
endpoint <- projected$roll_arc_fraction <= 0.01 | projected$roll_arc_fraction >= 0.99
endpoint_unconstrained <- projected$unconstrained_roll_arc_fraction <= 0.01 | projected$unconstrained_roll_arc_fraction >= 0.99
metrics <- data.frame(
  metric = c('passQC_input_cells','projected_run_cells','reference_smooth_muscle_cells_in_run','projection_success_fraction','constraint_changed_fraction','constrained_wall_distance_median_um','constrained_wall_distance_p95_um','unconstrained_wall_distance_median_um','unconstrained_wall_distance_p95_um','constrained_endpoint_fraction','unconstrained_endpoint_fraction'),
  value = c(
    nrow(passqc_coordinates), nrow(projected), sum(projected$smooth_muscle_reference %in% TRUE), mean(valid),
    mean(projected$projection_changed_by_constraint[valid], na.rm = TRUE),
    median(projected$wall_distance[valid], na.rm = TRUE),
    as.numeric(quantile(projected$wall_distance[valid], 0.95, na.rm = TRUE)),
    median(projected$unconstrained_wall_distance[valid_unconstrained], na.rm = TRUE),
    as.numeric(quantile(projected$unconstrained_wall_distance[valid_unconstrained], 0.95, na.rm = TRUE)),
    mean(endpoint[valid], na.rm = TRUE), mean(endpoint_unconstrained[valid_unconstrained], na.rm = TRUE)
  )
)
print(metrics)",
"b3-027-metrics")

add_code(
"smooth_distance <- data.frame(
  definition = c('Mayassi-reference Smooth muscle','all other reference labels'),
  n = c(sum(valid & projected$smooth_muscle_reference %in% TRUE), sum(valid & projected$smooth_muscle_reference %in% FALSE)),
  median_wall_distance_um = c(
    median(projected$wall_distance[valid & projected$smooth_muscle_reference %in% TRUE], na.rm = TRUE),
    median(projected$wall_distance[valid & projected$smooth_muscle_reference %in% FALSE], na.rm = TRUE)
  )
)
print(smooth_distance)",
"b3-028-smooth-distance")

add_code(
"coordinate_filename <- if (MAX_CELLS == 0L) 'Region_6_passQC_unrolled_cells.tsv.gz' else 'Region_6_passQC_unrolled_subset_cells.tsv.gz'
write_tsv_gz(projected, file.path(OUT_DIR, coordinate_filename))
write.table(trace, file.path(OUT_DIR, 'Region_6_ordered_dense_trace.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(metrics, file.path(OUT_DIR, 'Region_6_method_audit_metrics.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
cat('Coordinate table:', file.path(OUT_DIR, coordinate_filename), '\n')
cat('Audit tables written to:', OUT_DIR, '\n')",
"b3-029-write-tables")

add_code(
"rolled_plot <- ggplot(projected[valid, ], aes(x, y, color = roll_arc_fraction)) +
  geom_point(size = 0.16, alpha = 0.65) +
  geom_path(data = trace, aes(x, y), inherit.aes = FALSE, color = 'black', linewidth = 0.3) +
  coord_equal() + scale_color_viridis_c(option = 'turbo', name = 'roll arc\nfraction') +
  theme_void() + ggtitle('Published constrained projection on rolled coordinates')
ggsave(file.path(OUT_DIR, 'Region_6_passQC_rolled_arc.png'), rolled_plot, width = 9, height = 9, dpi = 220, bg = 'white')
rolled_plot",
"b3-030-rolled-plot")

add_markdown(
"![Rolled coordinates coloured by the constrained longitudinal coordinate](../../mayassi_unrolling_notebook/Region_6/Region_6_passQC_rolled_arc.png)",
"b3-031-rolled-image")

add_code(
"p_constrained <- ggplot(projected[valid, ], aes(roll_arc_fraction, wall_distance)) + geom_point(size = 0.18, alpha = 0.35, color = '#B2182B') + theme_bw() + labs(x = 'Constrained roll arc fraction', y = 'Wall distance (um)', title = 'Published global-radius constraint')
p_unconstrained <- ggplot(projected[valid_unconstrained, ], aes(unconstrained_roll_arc_fraction, unconstrained_wall_distance)) + geom_point(size = 0.18, alpha = 0.35, color = '#2166AC') + theme_bw() + labs(x = 'Unconstrained roll arc fraction', y = 'Nearest-trace distance (um)', title = 'Unconstrained diagnostic')
unrolled_comparison <- p_constrained + p_unconstrained
ggsave(file.path(OUT_DIR, 'Region_6_unrolled_projection_comparison.png'), unrolled_comparison, width = 14, height = 6, dpi = 220, bg = 'white')
unrolled_comparison",
"b3-032-unrolled-plot")

add_markdown(
"![Constrained and unconstrained unrolled-coordinate diagnostics](../../mayassi_unrolling_notebook/Region_6/Region_6_unrolled_projection_comparison.png)",
"b3-033-unrolled-image")

add_code(
"agreement_plot_data <- projected[valid & valid_unconstrained, ]
constraint_plot <- ggplot(agreement_plot_data, aes(unconstrained_roll_arc_fraction, roll_arc_fraction)) +
  geom_point(size = 0.18, alpha = 0.35) + geom_abline(slope = 1, intercept = 0, color = '#D73027') +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) + theme_bw() +
  labs(x = 'Unconstrained arc fraction', y = 'Constrained arc fraction', title = 'Sensitivity to the global-radius rule')
ggsave(file.path(OUT_DIR, 'Region_6_constraint_assignment_comparison.png'), constraint_plot, width = 7, height = 7, dpi = 220, bg = 'white')
constraint_plot",
"b3-034-constraint-plot")

add_markdown(
"![Cell assignments with and without the global-radius constraint](../../mayassi_unrolling_notebook/Region_6/Region_6_constraint_assignment_comparison.png)",
"b3-035-constraint-image")

add_markdown(
"## Step 6 - Statistical and scientific interpretation

- Cells are spatial observations, not independent biological replicates; therefore this notebook does not attach cell-level p-values to transfer-label or wall-distance differences.
- The manual trace is a geometric annotation and must be reviewed visually. Optimizing the trace to improve metrics while crossing unsupported tissue is invalid.
- The Mayassi reference-transfer label and its score are model-derived annotations and are not histological ground truth.
- A high constraint-change fraction or endpoint pile-up indicates model misspecification, not a biological gradient.
- The current longitudinal coordinate remains unoriented until the rolling direction is documented.",
"b3-036-statistical-audit")

add_code(
"metric_value <- setNames(metrics$value, metrics$metric)
method_gates <- data.frame(
  diagnostic = c('projection success >= 0.95','constraint-changed fraction <= 0.20','endpoint fraction <= 0.05'),
  observed = c(metric_value[['projection_success_fraction']], metric_value[['constraint_changed_fraction']], metric_value[['constrained_endpoint_fraction']]),
  pass = c(metric_value[['projection_success_fraction']] >= 0.95, metric_value[['constraint_changed_fraction']] <= 0.20, metric_value[['constrained_endpoint_fraction']] <= 0.05)
)
method_status <- if (all(method_gates$pass)) 'PASS_FOR_SCALE_UP' else 'NEEDS_REVISION_DO_NOT_SCALE'
print(method_gates)
cat('Method status:', method_status, '\n')",
"b3-037-method-gates")

add_markdown(
"## Final decision

If any gate above fails, the current Mayassi global-radius implementation must remain a feasibility result. Do not propagate it to Regions 1-5 and do not use its arc coordinate for eosinophil or other biological inference.

The appropriate next method comparison is a topology-aware local projection that retains the same manually reviewed trace while reducing cross-coil ambiguity. Morphology images and physical proximal-distal orientation remain valuable external validation evidence.",
"b3-038-decision")

add_code(
"sessionInfo()",
"b3-039-session")

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
