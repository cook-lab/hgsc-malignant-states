#!/bin/bash
# ============================================================
# Step 19_cnv - Build an upload bundle for the HPC cluster
# HGSC Single-Cell Atlas  ·  Cook Lab  ·  2026
# ============================================================
#
# Creates a tarball containing ONLY what the cluster needs to
# run CopyKAT + subclone clustering. Strips CopyKAT's large
# intermediate outputs from any already-done samples so the
# upload stays small (inputs only: counts.mtx.gz + barcodes.csv
# + genes.txt + ref_barcodes.txt).
#
# Usage:
#   bash scripts/19_cnv_bundle.sh [out_tar]
#
# Default output: 19_cnv_bundle.tar.gz at the project root.
# Expected size: ~11 GB for 251 samples (mostly counts matrices).

set -euo pipefail

OUT_TAR=${1:-"19_cnv_bundle.tar.gz"}

ROOT="/Volumes/CookLab/Sarah/2026_final_atlas"
STAGE="$(mktemp -d)"
trap 'rm -rf "${STAGE}"' EXIT

echo "Staging at ${STAGE}..."
mkdir -p "${STAGE}/per_sample" "${STAGE}/scripts" "${STAGE}/logs"

# Copy only the 4 input files per sample (skip CopyKAT outputs, DONE/FAILED)
N=0
for d in "${ROOT}"/output/19_cnv/per_sample/*/; do
    sid=$(basename "${d}")
    dst="${STAGE}/per_sample/${sid}"
    mkdir -p "${dst}"
    for f in counts.mtx.gz barcodes.csv genes.txt ref_barcodes.txt; do
        if [ -f "${d}${f}" ]; then
            cp "${d}${f}" "${dst}/${f}"
        fi
    done
    N=$((N + 1))
done
echo "Staged ${N} sample input dirs"

# Scripts
cp "${ROOT}/scripts/19_cnv_copykat.R"       "${STAGE}/scripts/"
cp "${ROOT}/scripts/19_cnv_slurm.sh"        "${STAGE}/scripts/"
cp "${ROOT}/scripts/19_cnv_status.py"       "${STAGE}/scripts/"

# Sample manifest (for reference) + pending list (array index source)
cp "${ROOT}/output/19_cnv/tables/sample_manifest.csv" \
   "${STAGE}/sample_manifest.csv"
ls -d "${STAGE}"/per_sample/*/ | xargs -n1 basename | sort \
   > "${STAGE}/pending_samples.txt"
echo "pending_samples.txt: $(wc -l < "${STAGE}/pending_samples.txt") samples"

# Inline README
cat > "${STAGE}/README.md" <<'EOF'
# Step 19 CNV - Cluster Run Bundle

## Layout after unpacking
```
<root>/
  per_sample/<sample_id>/counts.mtx.gz, barcodes.csv, genes.txt, ref_barcodes.txt
  scripts/19_cnv_copykat.R       # R worker (1 sample per invocation)
  scripts/19_cnv_slurm.sh        # SLURM array job template
  scripts/19_cnv_status.py       # optional: live progress checker
  pending_samples.txt            # array index -> sample_id
  sample_manifest.csv            # original manifest for reference
  logs/                          # per-job stdout/err will land here
```

## One-time setup on the cluster
```bash
# Load R (your module name may differ)
module load R/4.3.0

# Install copykat (once, into your user lib)
Rscript -e 'install.packages(c("remotes","Matrix","cluster"), repos="https://cloud.r-project.org")'
Rscript -e 'remotes::install_github("navinlabcode/copykat")'
```

## Submit the full batch
```bash
cd <root>

# Edit scripts/19_cnv_slurm.sh first:
#   #SBATCH --account=YOUR_ACCOUNT
#   #SBATCH --partition=YOUR_PARTITION
#   #SBATCH --cpus-per-task=8   (raise if your nodes have more)
#   #SBATCH --mem=24G
#   #SBATCH --time=04:00:00     (raise for ref-heavy samples)

N=$(wc -l < pending_samples.txt)
sbatch --array=1-${N}%50 scripts/19_cnv_slurm.sh
# %50 = max 50 concurrent array tasks. Match your QoS limits.
```

## Monitor
```bash
squeue -u $USER
python scripts/19_cnv_status.py      # reads DONE.txt / FAILED.txt / QC CSVs
tail -f logs/batch_progress.log
```

## Download results
Only need the per_sample/ outputs back. Rsync over:
```bash
# from your local machine
rsync -av <cluster_host>:/path/to/<root>/per_sample/ \
  /Volumes/CookLab/Sarah/2026_final_atlas/output/19_cnv/per_sample/
rsync -av <cluster_host>:/path/to/<root>/logs/ \
  /Volumes/CookLab/Sarah/2026_final_atlas/output/19_cnv/logs/
```

Then locally:
```bash
python scripts/19_cnv_status.py
python scripts/19_cnv_independence.py
python scripts/19_cnv_coexistence.py
```

## Troubleshooting
- **FAILED samples**: open `per_sample/<sid>/FAILED.txt` for the R error; array tasks auto-remove FAILED.txt on next retry so just resubmit the same sbatch.
- **Partial re-run**: `sbatch --array=<comma-list-of-ids>%N scripts/19_cnv_slurm.sh` to re-run specific indices.
- **Timeouts**: bump `--time=` in slurm.sh for ref-heavy samples (aneu > 3000).
EOF

# Build tar
echo "Creating ${OUT_TAR}..."
(cd "${STAGE}" && tar -czf "${ROOT}/${OUT_TAR}" .)

SIZE=$(du -sh "${ROOT}/${OUT_TAR}" | awk '{print $1}')
echo ""
echo "========================================================"
echo "  Bundle:  ${ROOT}/${OUT_TAR}"
echo "  Size:    ${SIZE}"
echo "  Samples: ${N}"
echo "========================================================"
echo "  Upload with:"
echo "    rsync -avP ${OUT_TAR} <cluster_host>:<workdir>/"
echo "  Unpack on cluster:"
echo "    cd <workdir> && tar -xzf ${OUT_TAR}"
echo "========================================================"
