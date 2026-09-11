#!/usr/bin/env Rscript

# Unit tests for the Mayassi-style Swiss-roll coordinate transformation.
# These tests use small hand-checked geometries so expected values do not
# depend on the implementation under test.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
test_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)
repo_root <- dirname(dirname(test_path))

source(file.path(repo_root, "R", "source.R"))

testthat::test_that("fixed primary audit detects exact exclusive-boundary agreement", {
  cells <- data.frame(
    nFeature_Xenium = c(5, 6, 199, 200),
    nCount_Xenium = c(11, 11, 999, 999),
    primary_include = c(FALSE, TRUE, TRUE, FALSE)
  )

  observed <- audit_fixed_primary_include(cells)

  testthat::expect_equal(observed$n_fixed_primary, 2L)
  testthat::expect_equal(observed$n_disagree, 0L)
  testthat::expect_true(observed$stored_matches_fixed_rule)
})

testthat::test_that("marker coverage retains absent markers in the denominator", {
  markers <- data.frame(
    Gene_Symbol = c("g1", "g2", "g3", "g4"),
    CellType_main = c("A", "A", "B", "B"),
    CellType_subtype = c("A1", "A1", "B1", "B1")
  )

  observed <- calculate_marker_definition_coverage(markers, c("g1", "g3", "g4"))

  testthat::expect_equal(observed$n_defined, c(2L, 2L))
  testthat::expect_equal(observed$n_present, c(1L, 2L))
  testthat::expect_equal(observed$coverage, c(0.5, 1))
})

testthat::test_that("binary cell-type agreement reports overlap without pseudoreplication", {
  observed <- compare_binary_cell_definitions(
    c(TRUE, TRUE, FALSE, FALSE),
    c(TRUE, FALSE, TRUE, FALSE),
    first_name = "marker",
    second_name = "reference"
  )

  testthat::expect_equal(observed$n_intersection, 1L)
  testthat::expect_equal(observed$n_union, 3L)
  testthat::expect_equal(observed$jaccard, 1 / 3)
})

testthat::test_that("gzip TSV writer creates a compressed round-trip artifact", {
  path <- tempfile(fileext = ".tsv.gz")
  on.exit(unlink(path), add = TRUE)
  expected <- data.frame(cell_id = c("a", "b"), value = c(1, 2))

  write_tsv_gz(expected, path)
  observed <- utils::read.delim(gzfile(path), check.names = FALSE)
  magic <- readBin(path, what = "raw", n = 2)

  testthat::expect_identical(as.integer(magic), c(31L, 139L))
  testthat::expect_equal(observed, expected)
})

testthat::test_that("cell IDs fall back to unique non-empty metadata row names", {
  metadata <- data.frame(x_centroid = c(1, 2), row.names = c("cell-a", "cell-b"))

  observed <- extract_cell_ids(metadata)

  testthat::expect_identical(observed, c("cell-a", "cell-b"))
})

testthat::test_that("ordered control points are interpolated at the requested spacing", {
  control <- data.frame(x = c(0, 3), y = c(0, 0))
  observed <- interpolate_trace_control_points(control, spacing = 1)

  testthat::expect_equal(observed$x, c(0, 1, 2, 3))
  testthat::expect_equal(observed$y, c(0, 0, 0, 0))
})

testthat::test_that("manual trace controls require finite sequential reviewed points", {
  controls <- data.frame(
    point_order = 1:3,
    x = c(10, 20, 30),
    y = c(40, 50, 60),
    comment = c("outer edge", "middle turn", "inner endpoint")
  )

  observed <- validate_trace_control_points(controls)

  testthat::expect_equal(observed$n_control_points, 3L)
  testthat::expect_equal(observed$start_x, 10)
  testthat::expect_equal(observed$end_y, 60)
  broken <- controls[c(1, 3), ]
  testthat::expect_error(
    validate_trace_control_points(broken),
    "point_order must be consecutive"
  )
})

testthat::test_that("cumulative arc length starts at zero and follows the ordered trace", {
  trace <- data.frame(x = c(0, 3, 3), y = c(0, 0, 4))
  observed <- add_trace_arc_length(trace)

  testthat::expect_equal(observed$roll_arc_length, c(0, 3, 7))
  testthat::expect_equal(observed$roll_arc_fraction, c(0, 3 / 7, 1))
})

testthat::test_that("projection returns longitudinal position and distance to trace", {
  trace <- data.frame(
    x = c(0, 10, 20),
    y = c(0, 0, 0),
    roll_arc_length = c(0, 10, 20),
    roll_arc_fraction = c(0, 0.5, 1),
    radius_from_center = c(20, 10, 0)
  )
  cells <- data.frame(cell_id = c("a", "b"), x = c(1, 19), y = c(2, 3))

  observed <- project_cells_to_trace(cells, trace, inward_constraint = FALSE)

  testthat::expect_equal(observed$trace_index, c(1L, 3L))
  testthat::expect_equal(observed$roll_arc_length, c(0, 20))
  testthat::expect_equal(observed$wall_distance, c(sqrt(5), sqrt(10)))
})

testthat::test_that("inward constraint prevents projection to a more external coil", {
  trace <- data.frame(
    x = c(10, 2), y = c(0, 0),
    roll_arc_length = c(0, 1), roll_arc_fraction = c(0, 1),
    radius_from_center = c(10, 2)
  )
  cells <- data.frame(cell_id = "inside", x = 3, y = 0)

  observed <- project_cells_to_trace(
    cells,
    trace,
    inward_constraint = TRUE,
    center = c(x = 0, y = 0)
  )

  testthat::expect_equal(observed$trace_index, 2L)
})

testthat::test_that("trace quality flags discontinuous jumps", {
  trace <- data.frame(x = c(0, 1, 2, 20), y = c(0, 0, 0, 0))
  quality <- assess_trace_quality(trace, jump_multiplier = 5)

  testthat::expect_false(quality$continuous)
  testthat::expect_equal(quality$n_large_jumps, 1L)
})
