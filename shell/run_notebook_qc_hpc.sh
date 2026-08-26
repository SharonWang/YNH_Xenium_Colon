#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-Region_1}"
RUN_LABEL="${2:-colon_qc_hpc_$(date +%Y%m%d)}"
PROJECT_ROOT="/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium"
PIPELINE_REPO="${PROJECT_ROOT}/colon_analysis/YNH_Xenium_Colon"
INPUT_ROOT="${PROJECT_ROOT}/colon_data"
RUN_ROOT="${PROJECT_ROOT}/colon_analysis/colon_qc_outputs/${RUN_LABEL}"
TEMP_ROOT="${PROJECT_ROOT}/colon_analysis/tmp/${RUN_LABEL}"
R_LIBS_USER="${R_LIBS_USER:-${PROJECT_ROOT}/colon_analysis/R_libs}"
export PROJECT_ROOT PIPELINE_REPO INPUT_ROOT RUN_ROOT R_LIBS_USER
export COLON_PROJECT_ROOT="$PROJECT_ROOT" COLON_PIPELINE_REPO="$PIPELINE_REPO"
export COLON_INPUT_ROOT="$INPUT_ROOT" COLON_RUN_ROOT="$RUN_ROOT"
export COLON_RUN_LABEL="$RUN_LABEL" COLON_TEMP_ROOT="$TEMP_ROOT" COLON_QC_MODE="FULL_HPC"
export TMPDIR="$TEMP_ROOT" TMP="$TEMP_ROOT" TEMP="$TEMP_ROOT" LC_ALL=C

case "$MODE" in ALL|SUMMARY|Region_1|Region_2|Region_3|Region_4|Region_5|Region_6) ;; *) echo "Invalid mode: $MODE" >&2; exit 2;; esac
[[ "$RUN_ROOT" == "$PROJECT_ROOT/colon_analysis/"* ]] || { echo "Unsafe run root" >&2; exit 2; }
[[ "$TEMP_ROOT" == "$PROJECT_ROOT/colon_analysis/"* ]] || { echo "Unsafe temp root" >&2; exit 2; }
command -v Rscript >/dev/null
command -v python3 >/dev/null
mkdir -p "$RUN_ROOT/executed_notebooks" "$TEMP_ROOT"

python3 "$PIPELINE_REPO/scripts/validate_notebooks.py" --require-nbformat "$PIPELINE_REPO"/notebooks/*.ipynb
Rscript -e 'required <- c("Matrix","jsonlite","ggplot2","arrow","dplyr","RANN"); missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]; if(length(missing)) stop(paste("Missing R packages:", paste(missing, collapse=", ")))' 
for region in {1..6}; do
  compgen -G "$INPUT_ROOT/output-*__Region_${region}__*" >/dev/null || { echo "Missing Region_${region} input" >&2; exit 2; }
done

run_one() {
  local notebook="$1"
  Rscript "$PIPELINE_REPO/scripts/execute_r_notebook.R" \
    "$PIPELINE_REPO/notebooks/$notebook" \
    "$RUN_ROOT/executed_notebooks/$notebook"
}

if [[ "$MODE" == "ALL" ]]; then
  for region in {1..6}; do run_one "01_QC_Region${region}.ipynb"; done
  run_one "02_slide_QC_summary.ipynb"
elif [[ "$MODE" == "SUMMARY" ]]; then
  run_one "02_slide_QC_summary.ipynb"
else
  run_one "01_QC_Region${MODE#Region_}.ipynb"
fi
