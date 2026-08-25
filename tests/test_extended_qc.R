options(stringsAsFactors = FALSE)

source(file.path("R", "source.R"))

expect_error <- function(expr, pattern) {
  message <- tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
  stopifnot(!is.na(message), grepl(pattern, message, ignore.case = TRUE))
}

regions <- paste0("Region_", 1:6)
stopifnot(identical(expected_colon_regions(), regions))
stopifnot(identical(names(section_palette()), regions))
stopifnot(identical(section_downstream_status(regions), rep("QC_EVIDENCE_PENDING", 6L)))

index <- data.frame(
  region_id = regions,
  run_label = rep("unit_run", 6L),
  execution_mode = rep("LOCAL_SUBSET", 6L),
  section_output_dir = file.path("sections", regions)
)
validated <- validate_section_bundle_index(index, regions)
stopifnot(identical(validated$region_id, regions))

expect_error(validate_section_bundle_index(index[-6, ], regions), "exactly six")
duplicate <- index; duplicate$region_id[[6L]] <- "Region_5"
expect_error(validate_section_bundle_index(duplicate, regions), "duplicate|exactly six")
mixed_run <- index; mixed_run$run_label[[6L]] <- "other"
expect_error(validate_section_bundle_index(mixed_run, regions), "run label")
mixed_mode <- index; mixed_mode$execution_mode[[6L]] <- "FULL_HPC"
expect_error(validate_section_bundle_index(mixed_mode, regions), "execution mode")

edges <- data.frame(
  region_id_from = c("Region_1", "Region_2"),
  region_id_to = c("Region_1", "Region_2"),
  from_cell_id = c("a", "b"), to_cell_id = c("c", "d")
)
stopifnot(isTRUE(validate_within_region_spatial_edges(edges)))
cross_region <- edges; cross_region$region_id_to[[2L]] <- "Region_3"
expect_error(validate_within_region_spatial_edges(cross_region), "cross-region")

manifest <- utils::read.delim(file.path("config", "colon_sample_manifest.tsv"), check.names = FALSE)
section_summary <- data.frame(region_id = regions, input_cells = 101:106, core_qc_pass = 91:96)
position_summary <- summarise_colon_mouse_position(section_summary, manifest)
stopifnot(identical(position_summary$region_summary$region_id, regions))
stopifnot(nrow(position_summary$within_mouse) == 6L)
stopifnot(nrow(position_summary$matched_position) == 3L)
stopifnot(all(position_summary$matched_position$n_mice == 2L))

gates <- expand.grid(region_id = regions, gate = "overall", stringsAsFactors = FALSE)
gates$status <- c("PASS", "WARN", "HOLD", "PASS", "PENDING", "PASS")
readiness <- derive_colon_region_readiness(gates, regions)
stopifnot(identical(readiness$region_id, regions))
stopifnot(identical(readiness$status, gates$status))
stopifnot(!any(readiness$status %in% c("PRIMARY", "PRIMARY_CONDITIONAL", "SENSITIVITY_ONLY")))

cat("Six-region colon aggregation contracts passed.\n")
