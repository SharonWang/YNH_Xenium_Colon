# YNH Xenium Colon

Reproducible R and Jupyter workflow for quality control and spatial analysis of a 10x Genomics Xenium experiment containing six wild-type mouse colon Swiss-roll regions.

## Project abstract

This project investigates spatial cell organization in colon tissue from two wild-type mice. Each mouse contributes three positions—top, middle, and bottom—placed as six regions on one Xenium slide. The repository provides a reproducible pipeline for per-region quality control (QC), slide-level QC synthesis, Region 6 cell-type analysis, and an auditable pilot for converting rolled Swiss-roll coordinates into an approximate one-dimensional unrolled coordinate.

The panel design combines the Xenium mouse base panel (380 genes) with a custom 100-gene eosinophil panel. The installed assay contains 479 biological targets after panel reconciliation. The custom panel includes seven shared eosinophil genes, 47 genes associated with a short-lived eosinophil program, and 46 genes associated with a long-lived eosinophil program. Because these programs were derived from adipose tissue, their use in colon is explicitly hypothesis-generating rather than a validated definition of eosinophil lifespan state.

The primary cell-level QC rule is fixed across all six regions:

```r
nFeature_Xenium > 5 & nFeature_Xenium < 200 &
  nCount_Xenium > 10 & nCount_Xenium < 1000
```

All cells and raw sparse counts are retained. Segmentation, nucleus, cell-area, control-burden, and spatial flags are recorded as review evidence and do not alter `primary_include`; `strict_include` is retained as a sensitivity-analysis mask.

## Xenium sample metadata

| Xenium region | Mouse | Colon position | Genotype | Biological replicate | Technical section |
|---|---|---|---|---|---|
| `Region_1` | Mouse 1 | Top | WT | Mouse 1 | Region 1 |
| `Region_2` | Mouse 1 | Middle | WT | Mouse 1 | Region 2 |
| `Region_3` | Mouse 1 | Bottom | WT | Mouse 1 | Region 3 |
| `Region_4` | Mouse 2 | Top | WT | Mouse 2 | Region 4 |
| `Region_5` | Mouse 2 | Middle | WT | Mouse 2 | Region 5 |
| `Region_6` | Mouse 2 | Bottom | WT | Mouse 2 | Region 6 |

Treatment, experimental condition, and age have not yet been provided. The machine-readable record is [`config/colon_sample_manifest.tsv`](config/colon_sample_manifest.tsv). With only two biological replicates, mouse and position comparisons are descriptive and are not treated as adequately powered population-level inference.

## Notebook guide

The notebooks are written as step-by-step analytical records: methods, code, statistics, diagnostic plots, and outputs are presented together so that each decision can be reviewed before the next stage is run.

### Initial QC

| Notebook | Scope and purpose |
|---|---|
| [`01_QC_Region1.ipynb`](notebooks/01_QC_Region1.ipynb) | Region 1 / Mouse 1 top: input audit, panel reconciliation, fixed cell QC, transcript and spatial QC, plots, and export of the region QC bundle. |
| [`01_QC_Region2.ipynb`](notebooks/01_QC_Region2.ipynb) | Region 2 / Mouse 1 middle: the same reproducible region-level QC workflow. |
| [`01_QC_Region3.ipynb`](notebooks/01_QC_Region3.ipynb) | Region 3 / Mouse 1 bottom: the same reproducible region-level QC workflow. |
| [`01_QC_Region4.ipynb`](notebooks/01_QC_Region4.ipynb) | Region 4 / Mouse 2 top: the same reproducible region-level QC workflow. |
| [`01_QC_Region5.ipynb`](notebooks/01_QC_Region5.ipynb) | Region 5 / Mouse 2 middle: the same reproducible region-level QC workflow. |
| [`01_QC_Region6.ipynb`](notebooks/01_QC_Region6.ipynb) | Region 6 / Mouse 2 bottom: the same reproducible region-level QC workflow. |
| [`02_slide_QC_summary.ipynb`](notebooks/02_slide_QC_summary.ipynb) | Combines validated outputs from all six regions; compares technical QC by region, mouse, and colon position; and creates the slide-level readiness summary. |

### Region 6 downstream and Swiss-roll analysis

| Notebook | Scope and purpose |
|---|---|
| [`B2_Region6_primary_479.ipynb`](notebooks/B2_Region6_primary_479.ipynb) | Builds and audits the Region 6 spatial object using the reconciled 479-gene panel, QC-passed cells, transferred reference labels, marker evidence, and refined cell-type annotations. Its pass-QC object is the input to the unrolling workflow. |
| [`B3a_Region6_trace_definition.ipynb`](notebooks/B3a_Region6_trace_definition.ipynb) | Defines the Swiss-roll centreline transparently from version-controlled control points. It starts from `Region_6_spatial_passQC.rds`, uses `RefAll_subtype_predicted.id == "Smooth muscle"` as the anatomical guide, numbers every control point, and produces diagnostics for manual review. |
| [`B3_Region6_Mayassi_unrolling_method_audit.ipynb`](notebooks/B3_Region6_Mayassi_unrolling_method_audit.ipynb) | Audits a Mayassi-style centreline projection and generates candidate unrolled coordinates, projection diagnostics, quality gates, and export tables. The current result is methodological development and must not be treated as a validated anatomical coordinate until the trace and orientation checks pass. |

## Source code and configuration

### Core R functions

[`R/source.R`](R/source.R) is the shared, documented function library used by the notebooks. Keeping reusable logic here ensures that every region uses the same implementation. It contains functions for:

- safe path and runtime validation;
- Xenium file discovery, integrity checks, MEX import, and panel reconciliation;
- fixed and extended cell-, transcript-, gene-, and spatial-QC calculations;
- region-level and slide-level summaries, figures, output bundles, and validation;
- Seurat object construction, reference-label transfer, marker scoring, and eosinophil annotation audits;
- Swiss-roll trace validation, interpolation, cell-to-trace projection, and Mayassi-style coordinate transfer.

### Configuration files

| File | Purpose |
|---|---|
| [`config/colon_sample_manifest.tsv`](config/colon_sample_manifest.tsv) | Authoritative mapping from region to mouse, colon position, replicate, genotype, and metadata status. |
| [`config/eos_gene_sets.tsv`](config/eos_gene_sets.tsv) | Custom eosinophil gene sets: shared, short-lived-associated, and long-lived-associated genes. |
| [`config/fixed_cell_qc_thresholds.tsv`](config/fixed_cell_qc_thresholds.tsv) | Fixed primary thresholds for Xenium feature and transcript counts. |
| [`config/extended_qc_defaults.tsv`](config/extended_qc_defaults.tsv) | Parameters for evidence-only transcript, spatial, control, and sensitivity QC. |
| [`config/subset_qc_reference.tsv`](config/subset_qc_reference.tsv) | Expected results used to validate the small local test subset. |
| [`config/unrolling/Region_6_trace_control_points.tsv`](config/unrolling/Region_6_trace_control_points.tsv) | Version-controlled Region 6 centreline control points used by B3a and B3. |
| [`config/unrolling/README.md`](config/unrolling/README.md) | Provenance, coordinate convention, and review instructions for the trace file. |

### Execution and notebook-building scripts

| File | Purpose |
|---|---|
| [`scripts/build_notebooks.R`](scripts/build_notebooks.R) | Generates the six region QC notebooks and the slide-level summary notebook from consistent templates. |
| [`scripts/build_region6_trace_definition_notebook.R`](scripts/build_region6_trace_definition_notebook.R) | Rebuilds B3a, including its methods text, trace checks, figures, and outputs. |
| [`scripts/build_mayassi_unrolling_notebook.R`](scripts/build_mayassi_unrolling_notebook.R) | Rebuilds the auditable B3 unrolling notebook. |
| [`scripts/create_local_subset.R`](scripts/create_local_subset.R) | Creates a small, balanced input subset for code-path validation without attempting full local analysis. |
| [`scripts/execute_local_subset.ps1`](scripts/execute_local_subset.ps1) | Runs the Windows subset validation with the project-specific R installation and D-drive-only paths. |
| [`scripts/execute_r_notebook.R`](scripts/execute_r_notebook.R) | Executes an R Jupyter notebook cell by cell and writes a valid executed notebook. |
| [`scripts/run_mayassi_region6_pilot.R`](scripts/run_mayassi_region6_pilot.R) | Runs a bounded Region 6 unrolling pilot outside Jupyter for rapid validation. |
| [`scripts/validate_notebooks.py`](scripts/validate_notebooks.py) | Validates notebook JSON and, when available, the Jupyter `nbformat` schema. |
| [`shell/run_notebook_qc_hpc.sh`](shell/run_notebook_qc_hpc.sh) | HPC entry point for one region or the complete six-region QC workflow. |
| [`slurm/colon_notebook_qc.sbatch`](slurm/colon_notebook_qc.sbatch) | SLURM submission wrapper for the HPC QC workflow. |

### Automated checks

The [`tests`](tests) directory contains R and Python contract tests for source functions, subset generation, QC outputs, notebook structure, trace provenance, and unrolling calculations. The local subset validates code paths only; scientific readiness is assessed from the complete HPC outputs.

## Output locations

The repository contains code and configuration, while generated data products remain outside Git:

- region QC: `colon_analysis/colon_qc_outputs/<RUN_LABEL>/sections/Region_N/`;
- slide summary: `colon_analysis/colon_qc_outputs/<RUN_LABEL>/slide_summary/`;
- executed notebook copies: `colon_analysis/colon_qc_outputs/<RUN_LABEL>/executed_notebooks/`;
- Region 6 unrolling: `colon_analysis/mayassi_unrolling_notebook/Region_6/`.

Source notebooks are not overwritten during routine QC execution.

## Reproducibility notes

- The production workflow is designed for the HPC clone because the complete Xenium files are too large for local analysis.
- Local runs use only a deliberately constructed subset and must not be interpreted as estimates of full-region QC rates.
- Region 6 unrolling begins with cells for which `primary_include_revised == TRUE` in `Region_6_spatial_passQC.rds`.
- The unrolling guide is the transferred reference label `RefAll_subtype_predicted.id == "Smooth muscle"`; `Xenium_cluster_subtype` is not used to define the trace.
- The final unrolled axis remains anatomically unoriented until the physical proximal–distal direction is confirmed.

## Related work

The README organization was informed by the clear project summary and linked notebook catalog in the [Patel Study repository](https://github.com/SharonWang/Patel_Study). The Region 6 pilot adapts the centreline-based conceptual strategy reported by Mayassi *et al.* (Nature, 2024), with project-specific diagnostics and explicit quality gates.
