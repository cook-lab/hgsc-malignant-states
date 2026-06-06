# 02_bayesprism_deconv.R — BayesPrism deconvolution of TCGA-OV bulk RNA-seq
# ====================================================================
# HGSC malignant-states atlas backend.
#
# Purpose: local Bayesian deconvolution of TCGA-OV bulk RNA-seq using the
# atlas single-cell reference (01_cibersort_reference.py). Models tumor
# heterogeneity (SecA/Intermediate/SecB/Ciliated as separate tumor types),
# TME contamination, and cross-platform batch effects.
#
# INPUTS (output_root/07_deconvolution_survival/):
#   - cibersortx_sc_reference_v2.txt
#   - cibersortx_phenotypes_v2.txt
#   - cibersortx_mixture.txt
#
# OUTPUTS (output_root/07_deconvolution_survival/):
#   - bayesprism_fractions.csv
#   - CIBERSORTx_Results.txt           (CIBERSORTx-compatible; read by 04/05)
#   - bayesprism_theta.rds
#   - bayesprism_summary.txt
#
# MANUSCRIPT PANELS: upstream of Fig 7E/F/G (TCGA survival).
#
# RUNTIME TIER: heavy (Gibbs sampling; ~10-30 min for 429 samples, n.cores=4).
#
# SEEDING: set.seed(CFG$seed) for Gibbs-sampling determinism.
#
# Usage:
#   Rscript 02_bayesprism_deconv.R

# config is 3 levels up: atlas/07_deconvolution_survival/02_*.R -> repo root
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
repo_root <- normalizePath(file.path(dirname(this_file), "..", ".."))
source(file.path(repo_root, "config", "config.R"))
set.seed(CFG$seed)

cat("=================================================================\n")
cat("  BayesPrism deconvolution of TCGA-OV\n")
cat("=================================================================\n\n")

# ── Paths ──────────────────────────────────────────────────────
out_dir    <- cfg_path("output_root", "07_deconvolution_survival")
ref_file   <- file.path(out_dir, "cibersortx_sc_reference_v2.txt")
pheno_file <- file.path(out_dir, "cibersortx_phenotypes_v2.txt")
mix_file   <- file.path(out_dir, "cibersortx_mixture.txt")

# ── Install BayesPrism if needed ───────────────────────────────
cat("[0] Checking / installing BayesPrism...\n")
if (!requireNamespace("snowfall", quietly = TRUE)) {
  install.packages("snowfall", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("BayesPrism", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "https://cloud.r-project.org")
  }
  devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")
}

library(BayesPrism)
library(data.table)
cat("    BayesPrism loaded successfully\n")

# ── 1. Reference ───────────────────────────────────────────────
cat("\n[1] Loading scRNA-seq reference...\n")
ref_dt <- fread(ref_file, header = TRUE, sep = "\t", check.names = FALSE)
ref_genes <- ref_dt$Gene
ref_mat   <- as.matrix(ref_dt[, -1, with = FALSE])
rownames(ref_mat) <- ref_genes
cat("    Reference matrix:", nrow(ref_mat), "genes x", ncol(ref_mat), "cells\n")

pheno <- fread(pheno_file, header = FALSE, sep = "\t")
cell_types <- pheno$V2
names(cell_types) <- pheno$V1
cat("    Cell types:", length(unique(cell_types)), "groups\n")
cat("    Groups:", paste(sort(unique(cell_types)), collapse = ", "), "\n")

sc_ref <- t(ref_mat)
cell_type_labels <- cell_types[rownames(sc_ref)]
stopifnot(all(!is.na(cell_type_labels)))
rm(ref_dt, ref_mat); gc()

# ── 2. Mixture ─────────────────────────────────────────────────
cat("\n[2] Loading TCGA mixture...\n")
mix_dt   <- fread(mix_file, header = TRUE, sep = "\t", check.names = FALSE)
mix_genes <- mix_dt$Gene
mix_mat   <- as.matrix(mix_dt[, -1, with = FALSE])
rownames(mix_mat) <- mix_genes
bulk_mix <- t(mix_mat)
cat("    Mixture:", nrow(bulk_mix), "samples x", ncol(bulk_mix), "genes\n")
rm(mix_dt, mix_mat); gc()

# ── 3. Common genes ────────────────────────────────────────────
cat("\n[3] Intersecting gene spaces...\n")
common_genes <- intersect(colnames(sc_ref), colnames(bulk_mix))
cat("    Common genes:", length(common_genes), "\n")
sc_ref   <- sc_ref[, common_genes]
bulk_mix <- bulk_mix[, common_genes]

# ── 4. Tumor vs environment ────────────────────────────────────
cat("\n[4] Defining cell type categories...\n")
tumor_types <- c("SecA_epithelium", "Intermediate_epithelium",
                 "SecB_epithelium", "Ciliated_epithelium")
env_types   <- setdiff(unique(cell_type_labels), tumor_types)
cat("    Tumor types (", length(tumor_types), "):",
    paste(tumor_types, collapse = ", "), "\n")
cat("    Environment types (", length(env_types), "):",
    paste(env_types, collapse = ", "), "\n")

# ── 5. Prism object ────────────────────────────────────────────
cat("\n[5] Creating BayesPrism reference object...\n")
bp_ref <- new.prism(
  reference = sc_ref, mixture = bulk_mix, input.type = "count.matrix",
  cell.type.labels  = as.character(cell_type_labels),
  cell.state.labels = as.character(cell_type_labels),
  key = NULL, outlier.cut = 0.01, outlier.fraction = 0.1
)

# ── 6. Run ─────────────────────────────────────────────────────
cat("\n[6] Running BayesPrism deconvolution...\n")
t_start <- Sys.time()
bp_result <- run.prism(prism = bp_ref, n.cores = 4, update.gibbs = TRUE)
runtime <- difftime(Sys.time(), t_start, units = "mins")
cat("    Done in", round(as.numeric(runtime), 1), "minutes\n")
saveRDS(bp_result, file.path(out_dir, "bayesprism_theta.rds"))

# ── 7. Extract fractions ───────────────────────────────────────
cat("\n[7] Extracting estimated fractions...\n")
theta <- get.fraction(bp = bp_result, which.theta = "final", state.or.type = "type")
cat("    Fractions:", nrow(theta), "samples x", ncol(theta), "cell types\n")
cat("    Row sums (should be ~1):", round(range(rowSums(theta)), 4), "\n")
for (ct in colnames(theta)) {
  vals <- theta[, ct]
  cat(sprintf("      %-30s  mean=%.4f  median=%.4f  range=[%.4f, %.4f]\n",
              ct, mean(vals), median(vals), min(vals), max(vals)))
}

# ── 8. Save ────────────────────────────────────────────────────
cat("\n[8] Saving results...\n")
fractions_df <- as.data.frame(theta)
fractions_df$Sample <- rownames(theta)
fractions_df <- fractions_df[, c("Sample", colnames(theta))]
fwrite(fractions_df, file.path(out_dir, "bayesprism_fractions.csv"))

cibersort_fmt <- data.frame(Mixture = rownames(theta), theta, check.names = FALSE)
write.table(cibersort_fmt, file = file.path(out_dir, "CIBERSORTx_Results.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("    Saved: bayesprism_fractions.csv + CIBERSORTx_Results.txt\n")

# ── 9. Summary ─────────────────────────────────────────────────
summary_text <- paste0(
  "BayesPrism Deconvolution Summary\n================================\n",
  "Date: ", Sys.time(), "\n",
  "Reference: ", nrow(sc_ref), " cells x ", ncol(sc_ref), " genes\n",
  "Mixture: ", nrow(bulk_mix), " samples x ", ncol(bulk_mix), " genes\n",
  "Common genes: ", length(common_genes), "\n",
  "Cell types: ", length(unique(cell_type_labels)), "\n",
  "  Tumor: ", paste(tumor_types, collapse = ", "), "\n",
  "  Environment: ", paste(env_types, collapse = ", "), "\n",
  "Runtime: ", round(as.numeric(runtime), 1), " minutes\nn.cores: 4\n",
  "\nFraction summary:\n"
)
for (ct in colnames(theta)) {
  summary_text <- paste0(summary_text,
    sprintf("  %-30s  mean=%.4f  sd=%.4f\n", ct, mean(theta[, ct]), sd(theta[, ct])))
}
writeLines(summary_text, file.path(out_dir, "bayesprism_summary.txt"))

cat("\n", rep("=", 65), "\n", sep = "")
cat("  Step complete! Output:", out_dir, " Runtime:",
    round(as.numeric(runtime), 1), "min\n")
cat(rep("=", 65), "\n", sep = "")
