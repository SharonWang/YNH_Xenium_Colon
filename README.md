# YNH Xenium Colon QC

Independent, reproducible R pipeline for initial QC of six WT mouse colon Xenium regions. Region 1-3 are Mouse 1 top/middle/bottom; Region 4-6 are Mouse 2 top/middle/bottom. The pipeline preserves raw sparse counts and all cells. Primary inclusion is fixed and exclusive:

```r
nFeature_Xenium > 5 & nFeature_Xenium < 200 &
  nCount_Xenium > 10 & nCount_Xenium < 1000
```

Segmentation, nucleus, cell-area, control-burden, and spatial flags are review evidence. They do not change `primary_include`; `strict_include` is the documented sensitivity mask.

## Required R packages

Core/subset: `Matrix`, `jsonlite`, `ggplot2`. Full HPC transcript/spatial QC: `arrow`, `dplyr`, `RANN`. The launchers install nothing.

## Local subset execution

The local launcher is locked to `D:/Xiaonan/CODEX_projects/Yanan_Xenium/colon_analysis`, R 4.6.1 in `D:/Programs/R-4.6.1`, and `D:/Programs/R_library`.

```powershell
Set-Location 'D:\Xiaonan\CODEX_projects\Yanan_Xenium\colon_analysis\YNH_Xenium_Colon'
& '.\scripts\execute_local_subset.ps1' -Mode Region_1 -RunLabel local_subset_validation
```

After the Region 1 pilot succeeds:

```powershell
& '.\scripts\execute_local_subset.ps1' -Mode ALL -RunLabel local_subset_validation
```

## HPC preflight and pilot

Chunk 1 — enter the cloned repository and create colon-only log/temp directories:

```bash
PROJECT_ROOT=/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium
PIPELINE_REPO=${PROJECT_ROOT}/colon_analysis/YNH_Xenium_Colon
mkdir -p ${PROJECT_ROOT}/colon_analysis/logs ${PROJECT_ROOT}/colon_analysis/tmp
cd ${PIPELINE_REPO}
```

Chunk 2 — confirm the site-provided runtime and packages (no installation):

```bash
command -v Rscript
command -v python3
python3 scripts/validate_notebooks.py --require-nbformat notebooks/*.ipynb
Rscript -e 'required <- c("Matrix","jsonlite","ggplot2","arrow","dplyr","RANN"); print(data.frame(package=required, available=vapply(required, requireNamespace, logical(1), quietly=TRUE))); stopifnot(all(vapply(required, requireNamespace, logical(1), quietly=TRUE)))'
```

Chunk 3 — run Region 1 interactively as the full-data pilot:

```bash
RUN_LABEL=colon_qc_hpc_pilot_$(date +%Y%m%d)
bash shell/run_notebook_qc_hpc.sh Region_1 ${RUN_LABEL}
```

Chunk 4 — inspect the pilot outputs and executed notebook before submitting all regions:

```bash
find ${PROJECT_ROOT}/colon_analysis/colon_qc_outputs/${RUN_LABEL}/sections/Region_1 -maxdepth 1 -type f -printf '%f\n' | sort
test -s ${PROJECT_ROOT}/colon_analysis/colon_qc_outputs/${RUN_LABEL}/executed_notebooks/01_QC_Region1.ipynb
```

Chunk 5 — submit all regions plus the summary after pilot review. Add site-specific account/partition flags if required:

```bash
mkdir -p ${PROJECT_ROOT}/colon_analysis/logs
sbatch --export=ALL,MODE=ALL,RUN_LABEL=colon_qc_hpc_full_$(date +%Y%m%d) slurm/colon_notebook_qc.sbatch
```

Alternatively, run sequentially in an allocated interactive job:

```bash
bash shell/run_notebook_qc_hpc.sh ALL colon_qc_hpc_full_$(date +%Y%m%d)
```

## Outputs

Region artifacts are written to `colon_analysis/colon_qc_outputs/<RUN_LABEL>/sections/Region_N/`; the summary is written to `.../slide_summary/`; executed notebook copies are written to `.../executed_notebooks/`. Source notebooks are never overwritten.

Local subset results validate code paths only. Final technical readiness must be based on full HPC output. With two mice, top/middle/bottom comparisons are descriptive, and adipose-derived eosinophil programs are treated as colon hypotheses rather than validated lifespan identities.

The expected local code-path reference is versioned in `config/subset_qc_reference.tsv`. Because the subset deliberately enriches boundary cells and QC anomalies, its pass fractions must not be interpreted as full-region cell-quality estimates.

## Region 6 Swiss-roll trace and unrolling on HPC

This workflow starts from `Region_6_spatial_passQC.rds` and therefore analyzes the 127,935 cells retained by `primary_include_revised == TRUE`. The trace guide is defined only by `RefAll_subtype_predicted.id == "Smooth muscle"`. The cluster-derived `Xenium_cluster_subtype` label is not used.

Chunk 1 - set the exact paths, force temporary files below the colon analysis directory, and request all pass-QC cells:

```bash
export COLON_PROJECT_ROOT=/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium
export COLON_ANALYSIS_ROOT=${COLON_PROJECT_ROOT}/colon_analysis
export COLON_REPO_ROOT=${COLON_ANALYSIS_ROOT}/YNH_Xenium_Colon
export TMPDIR=${COLON_ANALYSIS_ROOT}/tmp
export TEMP=${TMPDIR}
export TMP=${TMPDIR}
export COLON_UNROLL_MAX_CELLS=0
mkdir -p ${TMPDIR} ${COLON_ANALYSIS_ROOT}/mayassi_unrolling_notebook/Region_6
cd ${COLON_REPO_ROOT}
```

`COLON_UNROLL_MAX_CELLS=0` means all 127,935 pass-QC cells. Use a positive integer such as `12000` only for a validation subset.

Chunk 2 - check the R kernel packages and required inputs. The commands install nothing:

```bash
Rscript -e 'required <- c("SeuratObject","ggplot2","dplyr","patchwork","jsonlite","testthat"); print(data.frame(package=required,available=vapply(required,requireNamespace,logical(1),quietly=TRUE))); stopifnot(all(vapply(required,requireNamespace,logical(1),quietly=TRUE)))'
test -s ${COLON_ANALYSIS_ROOT}/HPC_return/20260910_HPC_return/colon_downstream_outputs/full_notebook_qc_v2/03_Region_6_Primary479/Region_6_spatial_passQC.rds
test -s ${COLON_REPO_ROOT}/config/unrolling/Region_6_trace_control_points.tsv
```

Chunk 3 - regenerate and execute the trace-definition notebook:

```bash
Rscript scripts/build_region6_trace_definition_notebook.R notebooks/B3a_Region6_trace_definition.ipynb
Rscript scripts/execute_r_notebook.R notebooks/B3a_Region6_trace_definition.ipynb notebooks/B3a_Region6_trace_definition.ipynb
```

Open `B3a_Region6_trace_definition.ipynb` and inspect `Region_6_numbered_trace_definition.png`. In particular, review the two labelled inter-turn bridges and points 37-39 at the inner hook. Stop here if the path follows the wrong coil or crosses unsupported tissue.

Chunk 4 - after trace approval, regenerate and execute the full projection notebook:

```bash
Rscript scripts/build_mayassi_unrolling_notebook.R notebooks/B3_Region6_Mayassi_unrolling_method_audit.ipynb
Rscript scripts/execute_r_notebook.R notebooks/B3_Region6_Mayassi_unrolling_method_audit.ipynb notebooks/B3_Region6_Mayassi_unrolling_method_audit.ipynb
```

Chunk 5 - verify the executed artifacts:

```bash
Rscript tests/test_unrolling.R
Rscript tests/test_mayassi_notebook.R
Rscript tests/test_region6_trace_workflow.R
test -s ${COLON_ANALYSIS_ROOT}/mayassi_unrolling_notebook/Region_6/Region_6_passQC_unrolled_cells.tsv.gz
```

The HPC run with `COLON_UNROLL_MAX_CELLS=0` writes `Region_6_passQC_unrolled_cells.tsv.gz`; a bounded validation run writes `Region_6_passQC_unrolled_subset_cells.tsv.gz`. The unrolled coordinate remains anatomically unoriented until the physical proximal-distal direction is confirmed.
