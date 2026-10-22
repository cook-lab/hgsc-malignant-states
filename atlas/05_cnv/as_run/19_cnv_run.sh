#!/bin/bash
# Step 19_cnv — CopyKAT batch runner
# ===================================
# HGSC Single-Cell Atlas  ·  Cook Lab  ·  2026
#
# Drives 19_cnv_copykat.R across extracted samples with parallelism
# control. Skips already-done samples (resumable after crash / Ctrl-C).
#
# Usage:
#   bash scripts/19_cnv_run.sh [n_parallel] [n_cores] [limit]
#
# Defaults:
#   n_parallel = 2   (2 samples at a time)
#   n_cores    = 4   (4 R cores per CopyKAT run)
#   limit      = 0   (0 = all pending; N = pilot: first N pending)
#
# Total CPU: n_parallel * n_cores (default 2 * 4 = 8 of 14 on M3 Max,
# leaves 6 free for other work).

set -euo pipefail

N_PARALLEL=${1:-2}
N_CORES=${2:-4}
LIMIT=${3:-0}

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT_DIR="/Volumes/CookLab/Sarah/2026_final_atlas/output/19_cnv"
SAMPLE_DIR="${OUT_DIR}/per_sample"
R_SCRIPT="${SCRIPT_DIR}/19_cnv_copykat.R"
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "${LOG_DIR}"
export BATCH_LOG="${LOG_DIR}/batch_progress.log"

echo "========================================================"
echo "  19_cnv CopyKAT batch runner"
echo "  Parallel: ${N_PARALLEL} samples x ${N_CORES} cores"
echo "  Limit:    ${LIMIT}  (0 = all)"
echo "  Log:      ${BATCH_LOG}"
echo "========================================================"

TOTAL=$(ls -d "${SAMPLE_DIR}"/*/ 2>/dev/null | wc -l | tr -d ' ')
DONE=$(find "${SAMPLE_DIR}" -name "DONE.txt" 2>/dev/null | wc -l | tr -d ' ')
FAILED=$(find "${SAMPLE_DIR}" -name "FAILED.txt" 2>/dev/null | wc -l | tr -d ' ')
PENDING=$((TOTAL - DONE - FAILED))

echo "  Total extracted:   ${TOTAL}"
echo "  Already done:      ${DONE}"
echo "  Previously failed: ${FAILED}"
echo "  Pending:           ${PENDING}"
echo ""

if [ "${PENDING}" -eq 0 ]; then
    echo "  All samples processed. Nothing to do."
    exit 0
fi

# Deterministic pending list, sorted largest-first by total cell count
# (so big samples run early when RAM headroom is highest)
PENDING_FILE=$(mktemp)
MANIFEST="${OUT_DIR}/tables/sample_manifest.csv"
for dir in $(ls -d "${SAMPLE_DIR}"/*/ | sort); do
    if [ ! -f "${dir}/DONE.txt" ] && [ ! -f "${dir}/FAILED.txt" ]; then
        sid=$(basename "${dir}")
        # Look up total cells (n_epithelial + n_reference) from manifest
        total=$(awk -F',' -v id="${sid}" 'NR>1 && $1==id {print $4+$5}' "${MANIFEST}")
        echo "${total:-0} ${dir}" >> "${PENDING_FILE}"
    fi
done
# Sort descending by cell count, then strip the count column
sort -t' ' -k1 -nr "${PENDING_FILE}" | awk '{print $2}' > "${PENDING_FILE}.sorted"
mv "${PENDING_FILE}.sorted" "${PENDING_FILE}"

if [ "${LIMIT}" -gt 0 ]; then
    head -n "${LIMIT}" "${PENDING_FILE}" > "${PENDING_FILE}.limit"
    mv "${PENDING_FILE}.limit" "${PENDING_FILE}"
    echo "  Running first ${LIMIT} pending samples only."
fi

N_TO_RUN=$(wc -l < "${PENDING_FILE}" | tr -d ' ')
echo "  Will run: ${N_TO_RUN} samples"
echo ""

# Header line if log is new
if [ ! -s "${BATCH_LOG}" ]; then
    echo -e "timestamp\tsample_id\tstatus\telapsed\textras" > "${BATCH_LOG}"
fi

cat "${PENDING_FILE}" | xargs -P "${N_PARALLEL}" -I {} \
    Rscript "${R_SCRIPT}" {} "${N_CORES}"

rm -f "${PENDING_FILE}"

DONE_AFTER=$(find "${SAMPLE_DIR}" -name "DONE.txt" 2>/dev/null | wc -l | tr -d ' ')
FAILED_AFTER=$(find "${SAMPLE_DIR}" -name "FAILED.txt" 2>/dev/null | wc -l | tr -d ' ')
NEW_DONE=$((DONE_AFTER - DONE))
NEW_FAILED=$((FAILED_AFTER - FAILED))

echo ""
echo "========================================================"
echo "  Batch complete"
echo "  Newly completed: ${NEW_DONE}"
echo "  Newly failed:    ${NEW_FAILED}"
echo "  Total done:      ${DONE_AFTER} / ${TOTAL}"
echo "  Progress log:    ${BATCH_LOG}"
echo "========================================================"
