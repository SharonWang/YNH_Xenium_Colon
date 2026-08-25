options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(Matrix))
source(file.path("scripts", "create_local_subset.R"))

counts <- Matrix::sparseMatrix(
  i = rep(1:20, 30), j = rep(1:30, each = 20), x = 1,
  dims = c(20, 30), dimnames = list(paste0("Gene", 1:20), paste0("cell", 1:30))
)
cells <- data.frame(
  cell_id = colnames(counts), nucleus_count = rep(c(0, 1, 2), 10),
  control_probe_counts = 1:30, genomic_control_counts = 0,
  control_codeword_counts = 0, x_centroid = 1:30, y_centroid = 30:1
)
first <- select_balanced_cells(counts, cells, target = 20L, seed = 42L, eos_genes = "Gene1")
second <- select_balanced_cells(counts, cells, target = 20L, seed = 42L, eos_genes = "Gene1")
stopifnot(identical(first, second), nrow(first) == 20L, !anyDuplicated(first$cell_id))
stopifnot(all(first$index >= 1L & first$index <= ncol(counts)))
cat("Deterministic balanced subset selection contract passed.\n")
