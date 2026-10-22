# ============================================================================
# 03_score_signatures.R
# Score SecA, SecB, and polarization per sample using:
#   1. UCell::ScoreSignatures_UCell  (signature-rank, robust)
#   2. z-scored mean of log-expression (sanity check)
# Polarization = z(SecB) - z(SecA) computed across full cohort.
#
# Outputs:
#   data/processed/sample_metadata_scored.tsv  (n=380 with scores)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(UCell)
})

ROOT <- "/Users/dpcook/My Drive/Lab/Manuscripts/hgsc_epitypes/supplemental_analysis/ICON7"
EXPR_FILE <- file.path(ROOT, "data/processed/expr_gene_symbol.rds")
META_FILE <- file.path(ROOT, "data/processed/sample_metadata.tsv")
OUT_FILE  <- file.path(ROOT, "data/processed/sample_metadata_scored.tsv")

gene_mat <- readRDS(EXPR_FILE)
meta <- fread(META_FILE)

# ---- Signatures -----------------------------------------------------------
sigs <- list(
  SecA = c("MECOM","FBXO21","LGR5","LPAR3","PBX1","SOX17","RCN2"),
  SecB = c("KRT17","KRT19","KRT7","LCN2","PRSS22","SLPI","TACSTD2")
)
# Restrict each signature to genes actually present on the array
sigs_present <- lapply(sigs, function(g) intersect(g, rownames(gene_mat)))
for (nm in names(sigs)) {
  miss <- setdiff(sigs[[nm]], sigs_present[[nm]])
  message(sprintf("%s: %d/%d genes on array; missing: %s",
                  nm, length(sigs_present[[nm]]), length(sigs[[nm]]),
                  if (length(miss)) paste(miss, collapse = ",") else "none"))
}

# ---- UCell scores ---------------------------------------------------------
# UCell expects genes x cells/samples matrix (which is what we have)
message("Running UCell ...")
uc_scores <- UCell::ScoreSignatures_UCell(gene_mat, features = sigs_present, ncores = 1)
# Result columns are <Signature>_UCell
uc_df <- as.data.frame(uc_scores)
uc_df$sample_title <- rownames(uc_df)
setnames(uc_df,
         old = c("SecA_UCell", "SecB_UCell"),
         new = c("SecA_ucell", "SecB_ucell"),
         skip_absent = TRUE)

# ---- z-scored mean of log-expression --------------------------------------
zscore_score <- function(mat, genes) {
  g_used <- intersect(genes, rownames(mat))
  if (length(g_used) == 0) return(rep(NA_real_, ncol(mat)))
  # z-score each gene across samples, then average across genes
  z <- t(scale(t(mat[g_used, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}

zsec_a <- zscore_score(gene_mat, sigs_present$SecA)
zsec_b <- zscore_score(gene_mat, sigs_present$SecB)

z_df <- data.frame(
  sample_title = colnames(gene_mat),
  SecA_zmean = zsec_a,
  SecB_zmean = zsec_b
)

# ---- Polarization scores --------------------------------------------------
# z-score the UCell outputs across the cohort to get matched scales
uc_df$SecA_ucell_z <- as.numeric(scale(uc_df$SecA_ucell))
uc_df$SecB_ucell_z <- as.numeric(scale(uc_df$SecB_ucell))
uc_df$polarization_ucell <- uc_df$SecB_ucell_z - uc_df$SecA_ucell_z

# For zmean: SecA_zmean / SecB_zmean are already z-scored per-gene then
# averaged. To match the manuscript Polarization = z(SecB) - z(SecA), we
# additionally z-score the per-sample averages across the cohort.
z_df$SecA_zmean_z <- as.numeric(scale(z_df$SecA_zmean))
z_df$SecB_zmean_z <- as.numeric(scale(z_df$SecB_zmean))
z_df$polarization_zmean <- z_df$SecB_zmean_z - z_df$SecA_zmean_z

# ---- Merge into metadata --------------------------------------------------
scored <- merge(meta, uc_df, by = "sample_title", all.x = TRUE)
scored <- merge(scored, z_df, by = "sample_title", all.x = TRUE)

# Score-method concordance check
if (all(c("polarization_ucell","polarization_zmean") %in% names(scored))) {
  rho <- cor(scored$polarization_ucell, scored$polarization_zmean,
             use = "complete.obs", method = "spearman")
  message(sprintf("Polarization UCell vs zmean Spearman rho = %.3f", rho))
}

# ---- Save -----------------------------------------------------------------
fwrite(scored, OUT_FILE, sep = "\t")
message(sprintf("Wrote: %s  (n=%d rows)", OUT_FILE, nrow(scored)))
