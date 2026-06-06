# 03_bayesprism_merged.R — BayesPrism with merged Secretory epithelium (validation)
# ====================================================================
# HGSC malignant-states atlas backend.
#
# Purpose: validation run that merges SecA + Intermediate + SecB into a single
# Secretory_epithelium type (10 cell types instead of 12). Cross-checks the
# epithelial fraction estimates against the full 02_bayesprism_deconv.R run
# (compared in 06_validate_survival.py).
#
# INPUTS (output_root/07_deconvolution_survival/):
#   - cibersortx_sc_reference_v2.txt
#   - cibersortx_phenotypes_merged.txt   (merged-label phenotypes)
#   - cibersortx_mixture.txt
#
# OUTPUTS (output_root/07_deconvolution_survival/):
#   - bayesprism_merged_theta.rds
#   - bayesprism_merged_results.txt
#   - bayesprism_merged_fractions.csv
#
# MANUSCRIPT PANELS: robustness check for Fig 7E/F/G.
#
# RUNTIME TIER: heavy (Gibbs sampling).
#
# SEEDING: set.seed(CFG$seed) for Gibbs-sampling determinism.
#
# Usage:
#   Rscript 03_bayesprism_merged.R

this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
repo_root <- normalizePath(file.path(dirname(this_file), "..", ".."))
source(file.path(repo_root, "config", "config.R"))
set.seed(CFG$seed)

cat("=================================================================\n")
cat("  BayesPrism (merged Secretory) validation\n")
cat("=================================================================\n\n")

out_dir    <- cfg_path("output_root", "07_deconvolution_survival")
ref_file   <- file.path(out_dir, "cibersortx_sc_reference_v2.txt")
pheno_file <- file.path(out_dir, "cibersortx_phenotypes_merged.txt")
mix_file   <- file.path(out_dir, "cibersortx_mixture.txt")

library(BayesPrism)
library(data.table)
cat("    BayesPrism loaded\n")

# 1. REFERENCE
cat("\n[1] Loading scRNA-seq reference...\n")
ref_dt <- fread(ref_file, header = TRUE, sep = "\t", check.names = FALSE)
ref_genes <- ref_dt$Gene
ref_mat   <- as.matrix(ref_dt[, -1, with = FALSE])
rownames(ref_mat) <- ref_genes
cat("    Reference:", nrow(ref_mat), "genes x", ncol(ref_mat), "cells\n")

pheno <- fread(pheno_file, header = FALSE, sep = "\t")
cell_types <- pheno$V2
names(cell_types) <- pheno$V1
cat("    Cell types:", length(unique(cell_types)), "groups\n")
cat("    Groups:", paste(sort(unique(cell_types)), collapse = ", "), "\n")

sc_ref <- t(ref_mat)
cell_type_labels <- cell_types[rownames(sc_ref)]
stopifnot(all(!is.na(cell_type_labels)))
rm(ref_dt, ref_mat); gc()

# 2. MIXTURE
cat("\n[2] Loading TCGA mixture...\n")
mix_dt   <- fread(mix_file, header = TRUE, sep = "\t", check.names = FALSE)
mix_genes <- mix_dt$Gene
mix_mat   <- as.matrix(mix_dt[, -1, with = FALSE])
rownames(mix_mat) <- mix_genes
bulk_mix <- t(mix_mat)
cat("    Mixture:", nrow(bulk_mix), "samples x", ncol(bulk_mix), "genes\n")
rm(mix_dt, mix_mat); gc()

# 3. INTERSECT GENES
cat("\n[3] Gene intersection...\n")
common_genes <- intersect(colnames(sc_ref), colnames(bulk_mix))
cat("    Common genes:", length(common_genes), "\n")
sc_ref   <- sc_ref[, common_genes]
bulk_mix <- bulk_mix[, common_genes]

# 4. CELL TYPE CATEGORIES
cat("\n[4] Cell type categories...\n")
tumor_types <- c("Secretory_epithelium", "Ciliated_epithelium")
env_types   <- setdiff(unique(cell_type_labels), tumor_types)
cat("    Tumor:", paste(tumor_types, collapse = ", "), "\n")
cat("    Environment:", paste(env_types, collapse = ", "), "\n")

# 5. PRISM OBJECT
cat("\n[5] Creating BayesPrism object...\n")
bp_ref <- new.prism(
  reference = sc_ref, mixture = bulk_mix, input.type = "count.matrix",
  cell.type.labels  = as.character(cell_type_labels),
  cell.state.labels = as.character(cell_type_labels),
  key = NULL, outlier.cut = 0.01, outlier.fraction = 0.1
)

# 6. RUN
cat("\n[6] Running BayesPrism...\n")
t_start <- Sys.time()
bp_result <- run.prism(prism = bp_ref, n.cores = 4, update.gibbs = TRUE)
runtime <- difftime(Sys.time(), t_start, units = "mins")
cat("    Done in", round(as.numeric(runtime), 1), "minutes\n")
saveRDS(bp_result, file.path(out_dir, "bayesprism_merged_theta.rds"))

# 7. EXTRACT FRACTIONS
cat("\n[7] Extracting fractions...\n")
theta <- get.fraction(bp = bp_result, which.theta = "final", state.or.type = "type")
cat("    Fractions:", nrow(theta), "x", ncol(theta), "\n")
cat("    Row sums:", round(range(rowSums(theta)), 4), "\n")
for (ct in colnames(theta)) {
  vals <- theta[, ct]
  cat(sprintf("    %-25s mean=%.4f  median=%.4f  range=[%.4f, %.4f]\n",
              ct, mean(vals), median(vals), min(vals), max(vals)))
}

# 8. SAVE
cibersort_fmt <- data.frame(Mixture = rownames(theta), theta, check.names = FALSE)
write.table(cibersort_fmt, file = file.path(out_dir, "bayesprism_merged_results.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
fractions_df <- as.data.frame(theta)
fractions_df$Sample <- rownames(theta)
fractions_df <- fractions_df[, c("Sample", colnames(theta))]
fwrite(fractions_df, file.path(out_dir, "bayesprism_merged_fractions.csv"))

cat("\n=================================================================\n")
cat("  Step complete! Runtime:", round(as.numeric(runtime), 1), "min\n")
cat("=================================================================\n")
