#!/bin/bash
# ============================================================
# Step 19_cnv - CopyKAT SLURM array job
# HGSC Single-Cell Atlas  ·  Cook Lab  ·  2026
# ============================================================
#
# One array task per sample. Each task runs CopyKAT + subclone
# clustering on ONE sample. Resumable: already-done samples
# (DONE.txt present) exit immediately. Tune --array range and
# %N (max concurrent) + --cpus-per-task to your cluster limits.
#
# Before submitting:
#   1. Upload bundle (inputs + script + manifest) to cluster
#   2. Edit the --account / --partition lines for your allocation
#   3. Build pending_samples.txt:
#        ls -d per_sample/*/ | xargs -n1 basename | sort > pending_samples.txt
#      (or copy the file produced by 19_cnv_bundle.sh)
#   4. Make sure R has copykat + Matrix + cluster packages:
#        Rscript -e 'if(!requireNamespace("copykat")) remotes::install_github("navinlabcode/copykat")'
#   5. Submit:
#        N=$(wc -l < pending_samples.txt)
#        sbatch --array=1-${N}%50 scripts/19_cnv_slurm.sh
#      (%50 = max 50 concurrent array tasks; raise/lower as allowed)
#
# After the array finishes:
#   rsync -av <cluster>:/path/per_sample/ output/19_cnv/per_sample/
#   python scripts/19_cnv_status.py
#   python scripts/19_cnv_independence.py
#   python scripts/19_cnv_coexistence.py

#SBATCH --job-name=19_cnv
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=logs/copykat_%A_%a.out
#SBATCH --error=logs/copykat_%A_%a.err
# # Uncomment / edit for your cluster:
# #SBATCH --account=YOUR_ACCOUNT
# #SBATCH --partition=YOUR_PARTITION

set -euo pipefail

# Resolve project root: this script lives at <root>/scripts/, so root=../
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

MANIFEST="${ROOT}/pending_samples.txt"   # one sample_id per line
PER_SAMPLE="${ROOT}/per_sample"
R_SCRIPT="${SCRIPT_DIR}/19_cnv_copykat.R"
LOG_DIR="${ROOT}/logs"
mkdir -p "${LOG_DIR}"

export BATCH_LOG="${LOG_DIR}/batch_progress.log"
if [ ! -s "${BATCH_LOG}" ]; then
    printf "timestamp\tsample_id\tstatus\telapsed\textras\n" > "${BATCH_LOG}"
fi

# Pull the sample_id for this array index (1-based)
if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID not set. Submit with sbatch --array=..."
    exit 1
fi
SAMPLE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")
if [ -z "${SAMPLE_ID}" ]; then
    echo "ERROR: no sample_id at line ${SLURM_ARRAY_TASK_ID} of ${MANIFEST}"
    exit 1
fi

SAMPLE_DIR="${PER_SAMPLE}/${SAMPLE_ID}"
if [ ! -d "${SAMPLE_DIR}" ]; then
    echo "ERROR: sample dir not found: ${SAMPLE_DIR}"
    exit 1
fi

echo "[$(date)] array=${SLURM_ARRAY_TASK_ID} sample=${SAMPLE_ID} cpus=${SLURM_CPUS_PER_TASK:-8}"

Rscript "${R_SCRIPT}" "${SAMPLE_DIR}" "${SLURM_CPUS_PER_TASK:-8}"
rc=$?

echo "[$(date)] array=${SLURM_ARRAY_TASK_ID} sample=${SAMPLE_ID} rc=${rc}"
exit ${rc}
