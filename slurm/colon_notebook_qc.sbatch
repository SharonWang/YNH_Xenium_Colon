#!/usr/bin/env bash
#SBATCH --job-name=colon_xenium_qc
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium/colon_analysis/logs/%x_%j.out
#SBATCH --error=/dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium/colon_analysis/logs/%x_%j.err

set -euo pipefail
MODE="${MODE:-Region_1}"
RUN_LABEL="${RUN_LABEL:-colon_qc_hpc_${SLURM_JOB_ID}}"
bash /dssg/home/acct-svetoslav_chakarov/svetoslav_chakarov/Lab_members/Yanan_Hu/YNH_Xenium/colon_analysis/YNH_Xenium_Colon/shell/run_notebook_qc_hpc.sh "$MODE" "$RUN_LABEL"
