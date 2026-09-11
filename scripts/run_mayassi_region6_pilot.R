#!/usr/bin/env Rscript

# Region 6 feasibility pilot for the Swiss-roll coordinate transformation used
# by Mayassi et al. (Nature 2024; doi:10.1038/s41586-024-08216-z).
#
# This script deliberately analyzes a deterministic subset locally. It does
# not modify the Seurat object or any HPC-returned notebook. Paths can be
# overridden with environment variables when the same code is run on HPC.

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
repo_root_default <- dirname(dirname(script_path))
repo_root <- Sys.getenv("COLON_REPO_ROOT", unset = repo_root_default)

if (.Platform$OS.type == "windows" && grepl("^[Cc]:", normalizePath(tempdir(), winslash = "/"))) {
  stop(
    "R temporary files are on C:. Set TMPDIR, TEMP, and TMP to a D: project path before running.",
    call. = FALSE
  )
}

analysis_root_default <- if (.Platform$OS.type == "windows") {
  "D:/Xiaonan/CODEX_projects/Yanan_Xenium/colon_analysis"
} else {
  "/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium/colon_analysis"
}
analysis_root <- Sys.getenv("COLON_ANALYSIS_ROOT", unset = analysis_root_default)

input_default <- file.path(
  analysis_root,
  "HPC_return", "20260910_HPC_return", "colon_downstream_outputs",
  "full_notebook_qc_v2", "03_Region_6_Primary479",
  "Region_6_spatial_passQC.rds"
)
input_rds <- Sys.getenv("COLON_REGION6_RDS", unset = input_default)
trace_tsv <- Sys.getenv(
  "COLON_REGION6_TRACE",
  unset = file.path(repo_root, "config", "unrolling", "Region_6_trace_control_points.tsv")
)
output_dir <- Sys.getenv(
  "COLON_UNROLL_OUTPUT",
  unset = file.path(analysis_root, "mayassi_unrolling_spike", "Region_6")
)
max_cells <- as.integer(Sys.getenv(
  "COLON_UNROLL_MAX_CELLS",
  unset = if (.Platform$OS.type == "windows") "12000" else "0"
))

required_paths <- c(repo_root, input_rds, trace_tsv)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("Required path(s) missing: ", paste(missing_paths, collapse = "; "), call. = FALSE)
}
if (!is.finite(max_cells) || max_cells < 0L || (max_cells > 0L && max_cells < 100L)) {
  stop("COLON_UNROLL_MAX_CELLS must be 0 for all cells or at least 100", call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(repo_root, "R", "source.R"))
suppressPackageStartupMessages(library(ggplot2))

message("Reading Region 6 object: ", input_rds)
region6 <- readRDS(input_rds)
metadata <- region6[[]]

required_metadata <- c(
  "x_centroid", "y_centroid", "primary_include_revised",
  "RefAll_subtype_predicted.id"
)
missing_columns <- setdiff(required_metadata, names(metadata))
if (length(missing_columns) > 0L) {
  stop(
    "Region 6 metadata is missing: ", paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}
if (any(!(metadata[["primary_include_revised"]] %in% TRUE))) {
  stop("Region_6_spatial_passQC.rds contains cells not passing primary_include_revised", call. = FALSE)
}

cells <- data.frame(
  cell_id = extract_cell_ids(metadata),
  x = as.numeric(metadata[["x_centroid"]]),
  y = as.numeric(metadata[["y_centroid"]]),
  smooth_muscle_reference =
    as.character(metadata[["RefAll_subtype_predicted.id"]]) == "Smooth muscle"
)
cells <- cells[
  stats::complete.cases(cells[, c("cell_id", "x", "y", "smooth_muscle_reference")]),
  , drop = FALSE
]

# A fixed seed makes the local feasibility subset exactly reproducible.
set.seed(20260910)
if (max_cells > 0L && nrow(cells) > max_cells) {
  selected <- sort(sample.int(nrow(cells), max_cells, replace = FALSE))
  cells_pilot <- cells[selected, , drop = FALSE]
} else {
  cells_pilot <- cells
}

# Control points were manually placed, outer edge to inner endpoint, on the
# smooth-muscle guide. Linear densification is the coordinate-space equivalent
# of extracting the ordered green-pixel path in the authors' released script.
control <- utils::read.delim(trace_tsv, check.names = FALSE)
validate_trace_control_points(control)
trace <- interpolate_trace_control_points(control[, c("x", "y")], spacing = 20)
trace <- add_trace_arc_length(trace)

# Mayassi et al. used the mean trace position as the center for their inward
# projection constraint. Retaining this detail helps prevent projection onto a
# spatially close but more external coil.
trace_center <- c(x = mean(trace$x), y = mean(trace$y))
trace$radius_from_center <- sqrt(
  (trace$x - trace_center[["x"]])^2 +
    (trace$y - trace_center[["y"]])^2
)

message("Projecting ", nrow(cells_pilot), " cells onto ", nrow(trace), " trace points")
projected <- project_cells_to_trace(
  cells_pilot,
  trace,
  inward_constraint = TRUE,
  center = trace_center
)

# An unconstrained projection is computed only as a diagnostic. Large changes
# between the two versions identify locations where adjacent coils are close
# enough to create projection ambiguity.
unconstrained <- project_cells_to_trace(
  cells_pilot,
  trace,
  inward_constraint = FALSE,
  center = trace_center
)
projected$unconstrained_trace_index <- unconstrained$trace_index
projected$projection_changed_by_constraint <-
  projected$trace_index != projected$unconstrained_trace_index

valid <- !is.na(projected$trace_index)
smooth_valid <- valid & projected$smooth_muscle_reference
other_valid <- valid & !projected$smooth_muscle_reference
trace_quality <- assess_trace_quality(trace)

metrics <- data.frame(
  metric = c(
    "input_cells", "pilot_cells", "trace_control_points", "dense_trace_points",
    "trace_arc_length", "projection_success_fraction",
    "constraint_changed_fraction", "wall_distance_median",
    "wall_distance_p95", "smooth_muscle_wall_distance_median",
    "other_cell_wall_distance_median", "trace_endpoint_fraction",
    "trace_continuous", "trace_large_jumps"
  ),
  value = c(
    nrow(cells), nrow(cells_pilot), nrow(control), nrow(trace),
    max(trace$roll_arc_length), mean(valid),
    mean(projected$projection_changed_by_constraint[valid], na.rm = TRUE),
    stats::median(projected$wall_distance[valid], na.rm = TRUE),
    as.numeric(stats::quantile(projected$wall_distance[valid], 0.95, na.rm = TRUE)),
    stats::median(projected$wall_distance[smooth_valid], na.rm = TRUE),
    stats::median(projected$wall_distance[other_valid], na.rm = TRUE),
    mean(
      projected$roll_arc_fraction[valid] <= 0.01 |
        projected$roll_arc_fraction[valid] >= 0.99,
      na.rm = TRUE
    ),
    trace_quality$continuous,
    trace_quality$n_large_jumps
  )
)

write_tsv_gz(
  projected,
  file.path(
    output_dir,
    if (max_cells == 0L) {
      "Region_6_passQC_unrolled_cells.tsv.gz"
    } else {
      "Region_6_passQC_unrolled_subset_cells.tsv.gz"
    }
  )
)
utils::write.table(
  trace,
  file.path(output_dir, "Region_6_ordered_dense_trace.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
utils::write.table(
  metrics,
  file.path(output_dir, "Region_6_unrolling_pilot_metrics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

guide_plot <- ggplot(cells, aes(x = x, y = y)) +
  geom_point(color = "grey82", size = 0.06, alpha = 0.20) +
  geom_point(
    data = cells[cells$smooth_muscle_reference, , drop = FALSE],
    color = "black", size = 0.08, alpha = 0.35
  ) +
  geom_path(data = trace, color = "#00A651", linewidth = 0.55) +
  geom_point(data = control, aes(x = x, y = y), color = "#D73027", size = 0.7) +
  coord_equal() +
  theme_void() +
  ggtitle("Region 6 Mayassi-style manual trace", subtitle = "green: dense trace; red: ordered control points")

rolled_plot <- ggplot(projected[valid, , drop = FALSE], aes(x = x, y = y)) +
  geom_point(aes(color = roll_arc_fraction), size = 0.18, alpha = 0.65) +
  geom_path(data = trace, color = "black", linewidth = 0.35, inherit.aes = FALSE, aes(x = x, y = y)) +
  coord_equal() +
  scale_color_viridis_c(option = "turbo", name = "roll arc\nfraction") +
  theme_void() +
  ggtitle("Rolled Region 6 colored by projected arc position")

unrolled_plot <- ggplot(
  projected[valid, , drop = FALSE],
  aes(x = roll_arc_fraction, y = wall_distance)
) +
  geom_point(aes(color = smooth_muscle_reference), size = 0.20, alpha = 0.40) +
  scale_color_manual(values = c(`FALSE` = "grey65", `TRUE` = "#111111")) +
  labs(
    x = "Unoriented roll arc fraction (outer trace start to inner endpoint)",
    y = "Distance to accepted trace point",
    color = "Mayassi-reference\nSmooth muscle",
    title = "Region 6 unrolled coordinate map"
  ) +
  theme_bw(base_size = 10)

ggsave(file.path(output_dir, "Region_6_trace_overlay.png"), guide_plot, width = 9, height = 9, dpi = 220, bg = "white")
ggsave(file.path(output_dir, "Region_6_rolled_arc_diagnostic.png"), rolled_plot, width = 9, height = 9, dpi = 220, bg = "white")
ggsave(file.path(output_dir, "Region_6_unrolled_map.png"), unrolled_plot, width = 11, height = 6, dpi = 220, bg = "white")

writeLines(
  c(
    "Mayassi-style Region 6 feasibility pilot",
    paste("Completed:", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)),
    paste("Input:", input_rds),
    paste("Trace control points:", trace_tsv),
    paste("Pilot cells:", nrow(cells_pilot)),
    paste("Projection success fraction:", signif(mean(valid), 5)),
    paste("Median wall distance:", signif(stats::median(projected$wall_distance[valid]), 5)),
    "Coordinate orientation is intentionally anatomical-unresolved."
  ),
  file.path(output_dir, "Region_6_unrolling_pilot_run.txt")
)

print(metrics)
message("Pilot outputs written to: ", output_dir)
