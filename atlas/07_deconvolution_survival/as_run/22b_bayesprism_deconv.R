# Trimmed copy (2026-10-21) of 202601_final_atlas/scripts/22b_bayesprism_deconv.R, the script that ran.
# Kept: the BayesPrism deconvolution of TCGA-OV with the atlas reference, written as CIBERSORTx_Results.txt (epithelial fraction
# read by 22d). Removed: the outputs nothing reads (posterior RDS, fractions CSV, run summary) and a closing pointer to 22c.
# The untrimmed script is in the lab project folder.
# 22b_bayesprism_deconv.R
# BayesPrism deconvolution of TCGA-OV bulk RNA-seq
# ====================================================================
# HGSC Single-Cell Atlas  ·  Cook Lab  ·  2026
#
# Local alternative to CIBERSORTx. BayesPrism is a Bayesian approach
# specifically designed for tumor deconvolution from bulk RNA-seq using
# scRNA-seq references. It properly handles:
#   - Tumor cell heterogeneity (SecA/Trans/SecB as separate types)
#   - TME contamination (the exact problem we're solving)
#   - Cross-platform batch effects (scRNA-seq → bulk RNA-seq)
#
# Input:
#   - output/22_tcga_deconvolution/cibersortx_sc_reference_v2.txt
#   - output/22_tcga_deconvolution/cibersortx_phenotypes_v2.txt
#   - output/22_tcga_deconvolution/cibersortx_mixture.txt
#
# Output (output/22_tcga_deconvolution/):
#   - bayesprism_fractions.csv            (estimated fractions, samples x cell types)
#   - bayesprism_fractions_cibersort_fmt.txt  (tab-delimited, CIBERSORTx-compatible)
#   - bayesprism_theta.rds                (full posterior, for advanced analyses)
#   - bayesprism_summary.txt              (run summary)
#
# Usage:
#   Rscript scripts/22b_bayesprism_deconv.R
# ====================================================================

cat("=================================================================\n")
cat("  Step 22b — BayesPrism deconvolution of TCGA-OV\n")
cat("=================================================================\n\n")

# ── Paths ──────────────────────────────────────────────────────
project_dir <- "/Volumes/CookLab/Sarah/2026_final_atlas"
out_dir     <- file.path(project_dir, "output", "22_tcga_deconvolution")

ref_file    <- file.path(out_dir, "cibersortx_sc_reference_v2.txt")
pheno_file  <- file.path(out_dir, "cibersortx_phenotypes_v2.txt")
mix_file    <- file.path(out_dir, "cibersortx_mixture.txt")

# ── Install BayesPrism if needed ───────────────────────────────
cat("[0] Checking / installing BayesPrism...\n")

if (!requireNamespace("snowfall", quietly = TRUE)) {
  install.packages("snowfall", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BayesPrism", quietly = TRUE)) {
  cat("    Installing BayesPrism from GitHub...\n")
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cloud.r-project.org")
  }
  devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")
}

library(BayesPrism)
library(data.table)
cat("    BayesPrism loaded successfully\n")

# ============================================================================
# 1. LOAD REFERENCE (scRNA-seq)
# ============================================================================
cat("\n[1] Loading scRNA-seq reference...\n")

ref_dt <- fread(ref_file, header = TRUE, sep = "\t", check.names = FALSE)
ref_genes <- ref_dt$Gene
ref_mat   <- as.matrix(ref_dt[, -1, with = FALSE])
rownames(ref_mat) <- ref_genes
cat("    Reference matrix:", nrow(ref_mat), "genes x", ncol(ref_mat), "cells\n")

# Load phenotypes
pheno <- fread(pheno_file, header = FALSE, sep = "\t")
cell_types <- pheno$V2
names(cell_types) <- pheno$V1
cat("    Cell types:", length(unique(cell_types)), "groups\n")
cat("    Groups:", paste(sort(unique(cell_types)), collapse = ", "), "\n")

# BayesPrism expects: cells x genes (transposed from our genes x cells)
sc_ref <- t(ref_mat)
cat("    Transposed to:", nrow(sc_ref), "cells x", ncol(sc_ref), "genes\n")

# Align cell type labels with reference matrix rows
cell_type_labels <- cell_types[rownames(sc_ref)]
stopifnot(all(!is.na(cell_type_labels)))
cat("    All cell labels matched\n")

rm(ref_dt, ref_mat); gc()

# ============================================================================
# 2. LOAD MIXTURE (TCGA bulk)
# ============================================================================
cat("\n[2] Loading TCGA mixture...\n")

mix_dt   <- fread(mix_file, header = TRUE, sep = "\t", check.names = FALSE)
mix_genes <- mix_dt$Gene
mix_mat   <- as.matrix(mix_dt[, -1, with = FALSE])
rownames(mix_mat) <- mix_genes
cat("    Mixture matrix:", nrow(mix_mat), "genes x", ncol(mix_mat), "samples\n")

# BayesPrism expects: samples x genes (transposed)
bulk_mix <- t(mix_mat)
cat("    Transposed to:", nrow(bulk_mix), "samples x", ncol(bulk_mix), "genes\n")

rm(mix_dt, mix_mat); gc()

# ============================================================================
# 3. FIND COMMON GENES
# ============================================================================
cat("\n[3] Intersecting gene spaces...\n")

common_genes <- intersect(colnames(sc_ref), colnames(bulk_mix))
cat("    Common genes:", length(common_genes), "\n")
cat("    Reference genes:", ncol(sc_ref), "\n")
cat("    Mixture genes:", ncol(bulk_mix), "\n")

sc_ref   <- sc_ref[, common_genes]
bulk_mix <- bulk_mix[, common_genes]

cat("    After subsetting:\n")
cat("      Reference:", nrow(sc_ref), "x", ncol(sc_ref), "\n")
cat("      Mixture:  ", nrow(bulk_mix), "x", ncol(bulk_mix), "\n")

# ============================================================================
# 4. DEFINE TUMOR vs ENVIRONMENT CELL TYPES
# ============================================================================
cat("\n[4] Defining cell type categories...\n")

# BayesPrism distinguishes tumor from environment for better estimation
tumor_types <- c("SecA_epithelium", "Transitioning_epithelium",
                 "SecB_epithelium", "Ciliated_epithelium")
env_types   <- setdiff(unique(cell_type_labels), tumor_types)

cat("    Tumor types (", length(tumor_types), "):",
    paste(tumor_types, collapse = ", "), "\n")
cat("    Environment types (", length(env_types), "):",
    paste(env_types, collapse = ", "), "\n")

# ============================================================================
# 5. CREATE BayesPrism OBJECTS
# ============================================================================
cat("\n[5] Creating BayesPrism reference object...\n")

# BayesPrism new.prism() creates the reference
# cell.type.labels: character vector of cell type for each cell
# cell.state.labels: can be same as cell.type.labels for our purpose
# (cell states would be finer resolution within each type)

bp_ref <- new.prism(
  reference = sc_ref,
  mixture   = bulk_mix,
  input.type = "count.matrix",
  cell.type.labels  = as.character(cell_type_labels),
  cell.state.labels = as.character(cell_type_labels),
  key = NULL,
  outlier.cut = 0.01,
  outlier.fraction = 0.1
)

cat("    BayesPrism reference object created\n")

# ============================================================================
# 6. RUN DECONVOLUTION
# ============================================================================
cat("\n[6] Running BayesPrism deconvolution...\n")
cat("    This may take 10-30 minutes for 429 samples...\n")

t_start <- Sys.time()

bp_result <- run.prism(
  prism   = bp_ref,
  n.cores = 4,
  update.gibbs = TRUE
)

t_end <- Sys.time()
runtime <- difftime(t_end, t_start, units = "mins")
cat("    Deconvolution complete in", round(as.numeric(runtime), 1), "minutes\n")

# ============================================================================
# 7. EXTRACT FRACTIONS
# ============================================================================
cat("\n[7] Extracting estimated fractions...\n")

# get.fraction() returns samples x cell_types matrix of estimated fractions
theta <- get.fraction(bp = bp_result,
                      which.theta = "final",
                      state.or.type = "type")

cat("    Fractions matrix:", nrow(theta), "samples x", ncol(theta), "cell types\n")
cat("    Column sums (should be ~1):", round(range(rowSums(theta)), 4), "\n")

# Summary per cell type
cat("\n    Cell type fraction summary:\n")
for (ct in colnames(theta)) {
  vals <- theta[, ct]
  cat(sprintf("      %-30s  mean=%.4f  median=%.4f  range=[%.4f, %.4f]\n",
              ct, mean(vals), median(vals), min(vals), max(vals)))
}

# ============================================================================
# 8. SAVE RESULTS
# ============================================================================
cat("\n[8] Saving results...\n")

# CIBERSORTx-compatible format (so 22c_survival_analysis.py can read it)
# Format: tab-delimited, first column = "Mixture" (sample IDs)
cibersort_fmt <- data.frame(Mixture = rownames(theta), theta, check.names = FALSE)
write.table(cibersort_fmt,
            file = file.path(out_dir, "CIBERSORTx_Results.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("    Saved: CIBERSORTx_Results.txt (compatible format for 22c)\n")

cat("\n", rep("=", 65), "\n", sep = "")
cat("  Step 22b complete!\n")
cat(rep("=", 65), "\n", sep = "")
cat("  Output:", out_dir, "\n")
cat("  Runtime:", round(as.numeric(runtime), 1), "minutes\n")
cat(rep("=", 65), "\n", sep = "")
