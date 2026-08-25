# Colon Xenium QC Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and locally validate a colon-only, six-region Xenium QC notebook pipeline with fixed non-destructive cell-inclusion masks and an HPC-ready execution path.

**Architecture:** One documented R function library provides path safety, six-region discovery, metadata validation, sparse import, fixed-threshold cell QC, transcript/spatial diagnostics, artifact writing, and slide aggregation. Six parameter-locked region notebooks call the same API; one summary notebook consumes the six validated region bundles. Local tests run on a deterministic 500-cell-per-region subset, while HPC mode reads full `transcripts.parquet` with Arrow.

**Tech Stack:** R 4.3+, Matrix, jsonlite, ggplot2, arrow, dplyr, RANN, IRkernel/Jupyter, Python 3 notebook tooling, PowerShell local launcher, Bash/SLURM HPC launcher.

**Spec:** `D:/Xiaonan/CODEX_projects/Yanan_Xenium/colon_analysis/COLON_XENIUM_QC_RESEARCH_PLAN_AND_PROGRESS.md`

## Global Constraints

- Work only below `D:/Xiaonan/CODEX_projects/Yanan_Xenium/colon_analysis` locally.
- Write no project data, temporary files, caches, logs, or libraries to `C:`.
- Do not read from or write to `adipose_analysis` at runtime.
- Treat `D:/Xiaonan/CODEX_projects/Yanan_Xenium/colon_data` as read-only input.
- Use exactly six regions and the approved mouse/position mapping.
- Preserve raw sparse counts and never destructively delete cells during QC.
- Define `primary_include` using exclusive bounds: `nFeature_Xenium > 5`, `< 200`, `nCount_Xenium > 10`, `< 1000`.
- Keep transcript tables on disk/chunked in full HPC mode.
- Update the living plan after every completed or changed stage.

---

### Task 1: Preserve baseline and add colon configuration contracts

**Files:**
- Create: `config/colon_sample_manifest.tsv`
- Create: `config/fixed_cell_qc_thresholds.tsv`
- Create: `config/extended_qc_defaults.tsv`
- Create: `.gitignore`
- Create: `tests/test_source.R`

**Interfaces:**
- Consumes: approved metadata and path contracts from the spec.
- Produces: six-row manifest; versioned threshold/config tables; failing R tests that define the colon API.

- [ ] Record a Git baseline of the current untracked starting files before modifying them.
- [ ] Write `tests/test_source.R` with assertions for six-region discovery, exact metadata order, unique eos genes, D/HPC path safety, and exclusive threshold boundaries.
- [ ] Run the test with `TMPDIR`, `TMP`, and `TEMP` below `colon_analysis/tmp`; confirm it fails because the colon API/configuration is absent.
- [ ] Add the three configuration files with exact values from the approved design.
- [ ] Add ignore rules for subset data, outputs, temporary files, R libraries, logs, and executed notebooks.
- [ ] Re-run the configuration-only assertions and record the result in the living plan.

### Task 2: Implement the documented fixed-threshold colon QC API

**Files:**
- Modify: `R/source.R`
- Modify: `tests/test_source.R`

**Interfaces:**
- Consumes: sparse counts, Xenium cell metadata, `region_id`, and versioned thresholds.
- Produces: `calculate_xenium_cell_qc()` and `build_cell_downstream_masks()` outputs with aligned metadata, exclusive primary bounds, non-destructive review fields, and provenance.

- [ ] Add failing tests where values exactly 5/200 features and 10/1000 counts fail, while 6/199 and 11/999 pass.
- [ ] Add failing tests showing `primary_include` is unaffected by segmentation/control review flags and `strict_include` is a subset of `primary_include`.
- [ ] Replace robust nFeature/nCount filtering in the active QC path with the four exact exclusive bounds; retain descriptive distribution summaries only.
- [ ] Remove scWAT synthetic metadata and four-region assumptions from active discovery/manifest/palette functions.
- [ ] Add roxygen-style documentation to every public function used by notebooks and focused comments for assay-specific logic.
- [ ] Run `tests/test_source.R`; require zero failures and no writable `C:` path in outputs.

### Task 3: Generalize six-region artifacts, spatial QC, and slide aggregation

**Files:**
- Modify: `R/source.R`
- Create: `tests/test_extended_qc.R`

**Interfaces:**
- Consumes: six validated region bundles with one run label and one execution mode.
- Produces: general region artifact validators/readers, six-region summary tables, within-mouse position comparisons, matched-position summaries, and readiness gates.

- [ ] Write failing tests for missing/duplicate/mixed-mode region outputs, no cross-region spatial edges, and exactly six ordered regions.
- [ ] Add tests that reject inherited Region 3 anchor and Region 4 sensitivity-only policies in colon QC outputs.
- [ ] Generalize artifact validators/readers to accept `expected_regions = paste0("Region_", 1:6)` without functions named `validate_four_*` in the notebook execution path.
- [ ] Derive region readiness from current files/integrity/panel/alerts rather than hard-coded region status.
- [ ] Implement top/middle/bottom summaries within mouse and same-position descriptive summaries across mice.
- [ ] Run both R test files and record pass/fail counts.

### Task 4: Generate six clean region notebooks and one clean summary notebook

**Files:**
- Create: `scripts/build_notebooks.py`
- Create: `scripts/render_notebooks.py`
- Create: `scripts/execute_r_notebook.R`
- Modify: `notebooks/01_QC_Region1.ipynb`
- Modify: `notebooks/01_QC_Region2.ipynb`
- Modify: `notebooks/01_QC_Region3.ipynb`
- Modify: `notebooks/01_QC_Region4.ipynb`
- Create: `notebooks/01_QC_Region5.ipynb`
- Create: `notebooks/01_QC_Region6.ipynb`
- Modify: `notebooks/02_slide_QC_summary.ipynb`
- Create: `tests/test_notebook_contracts.py`

**Interfaces:**
- Consumes: documented R API and colon configuration.
- Produces: seven source notebooks with cleared outputs, locked region IDs, HPC defaults, detailed Markdown, and identical region-notebook structure.

- [ ] Write failing structural tests for seven filenames, six locked region IDs, six-section expectation, colon paths, required headings, cleared outputs, and prohibited `scWAT`/`adipose_analysis` strings.
- [ ] Implement `build_notebooks.py` to render one canonical region structure into six committed notebooks and one six-region summary.
- [ ] Include Markdown explaining inputs, path safety, scientific intent, fixed thresholds, flags, interpretation, outputs, and local-versus-HPC limitations.
- [ ] Clear all execution counts and outputs in source notebooks.
- [ ] Run notebook contract tests and JSON parsing checks; require zero failures.

### Task 5: Add D-only local and HPC launchers

**Files:**
- Create: `scripts/execute_local_subset.ps1`
- Create: `shell/run_notebook_qc_hpc.sh`
- Create: `slurm/colon_notebook_qc.sbatch`
- Create: `README.md`
- Modify: `tests/test_notebook_contracts.py`

**Interfaces:**
- Consumes: source notebooks, renderer/executor, configuration, subset/full input roots.
- Produces: local `ALL|SUMMARY|Region_1..Region_6` execution and equivalent HPC execution with preflight and reload validation.

- [ ] Add failing launcher-contract tests for exact local/HPC roots, six regions, required packages, D-only temp paths, and prohibited adipose paths.
- [ ] Implement local launcher with explicit D-drive runtime paths and no dependency installation.
- [ ] Implement Bash and SLURM launchers with the approved HPC paths, package/kernel/input/capacity preflight, Region 1 pilot, ALL, and SUMMARY modes.
- [ ] Document copy/paste-ready HPC chunks and output contracts in `README.md`.
- [ ] Run Python contract tests and PowerShell parser validation; run Bash syntax validation when Bash is available and otherwise record the local environment blocker.

### Task 6: Create the deterministic balanced local subset

**Files:**
- Create: `scripts/create_local_subset.R`
- Create: `tests/test_subset_builder.R`
- Generate: `colon_analysis/subset_input/colon_data/**`
- Generate: `colon_analysis/subset_input/subset_manifest.tsv`
- Generate: `colon_analysis/subset_input/README_SUBSET_ONLY.txt`

**Interfaces:**
- Consumes: six read-only raw region bundles.
- Produces: 500 aligned cells per region with sparse MEX data, cell metadata, panel/metrics inputs, and provenance.

- [ ] Write failing tests for deterministic selection, 500 cells/region, alignment, boundary strata, spatial coverage, eos-marker enrichment, and D-only outputs.
- [ ] Implement a two-pass sampler that guarantees fixed-boundary/review strata first and fills remaining slots reproducibly.
- [ ] Write sparse subset matrices without densification and copy only small required run/panel metadata files.
- [ ] Run the builder once and verify 3,000 total unique region-prefixed cells, aligned MEX dimensions, hashes, and zero raw-input modifications.
- [ ] Record subset evidence and reference burden in the living plan.

### Task 7: Execute local end-to-end notebook validation

**Files:**
- Generate: `colon_analysis/colon_qc_outputs/local_subset_validation/**`
- Modify: `config/subset_qc_reference.tsv`
- Modify: `COLON_XENIUM_QC_RESEARCH_PLAN_AND_PROGRESS.md`

**Interfaces:**
- Consumes: validated subset, six source notebooks, summary notebook, R API, and launchers.
- Produces: six executed notebook copies, six region artifact bundles, slide summary, validation logs, and subset reference table.

- [ ] Run all R and Python tests in clean processes with D-only temp variables.
- [ ] Execute Region 1 locally and validate/reload every expected artifact.
- [ ] Execute Regions 2-6 under the same run label.
- [ ] Execute the summary notebook and validate six-region counts, masks, sparse objects, metadata, and figures.
- [ ] Scan generated text/JSON artifacts for `C:` writable paths, `adipose_analysis`, `scWAT`, stale four-region assumptions, errors, and warnings.
- [ ] Record exact commands, package versions, counts, failures, and limitations in the living plan.

### Task 8: Prepare the verified HPC handoff

**Files:**
- Modify: `README.md`
- Modify: `COLON_XENIUM_QC_RESEARCH_PLAN_AND_PROGRESS.md`

**Interfaces:**
- Consumes: passing local evidence and fixed HPC path contracts.
- Produces: Region 1 pilot, Regions 2-6, ALL, SUMMARY, and post-run validation command chunks.

- [ ] Run final source/notebook/launcher/path validation from the repository root.
- [ ] Confirm no HPC script installs software or writes outside `colon_analysis`.
- [ ] Document required site-supplied module/account/partition details without inventing them.
- [ ] Update progress, decisions, and change log with the exact local validation evidence and remaining full-HPC responsibility.
- [ ] Commit the validated implementation in focused commits after reviewing the complete diff.

