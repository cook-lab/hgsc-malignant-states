#!/usr/bin/env Rscript
# ============================================================================
# CCLE/DepMap — score HGSC cell lines for SecA / SecB / polarization   [PLACEHOLDER]
# ----------------------------------------------------------------------------
# *** PLACEHOLDER PANEL — final figure/panel letter NOT yet assigned. ***
# Slated as a supplementary / mechanistic side-quest (the SMAC-mimetic
# pharmacological vulnerability of SecB-polarized HGSC lines, sibling to the
# in vivo niche analyses). Rename this directory + scripts and set MANUSCRIPT
# PANEL(S) once the panel is decided. See ./README.md and the source module
# ../../2026_final_xenium_analysis/davids side quests/ccle_depmap/ (+ its
# report.html / README.md).
#
# PURPOSE
#   PREREQUISITE stage. Score DepMap HGSC cell lines for SecA / SecB / SecB
#   polarization using UCell on log2(TPM+1) expression, plus a per-gene
#   z-scored-mean variant (primary for bulk RNA-seq, because UCell SecA rank
#   scores saturate at ~0 for ~half of CCLE lines where LGR5/LPAR3/SOX17 are at
#   log2 TPM+1 ~ 0). Defines the per-line polarization axis (z(SecB) - z(SecA))
#   that scripts 02 (PRISM) and 03 (CRISPR) correlate against.
#
# INPUTS
#   cfg_obj("ccle_model_meta")  -> Model.csv (DepMap 24Q4 sample metadata)
#   cfg_obj("ccle_expression")  -> OmicsExpressionProteinCodingGenesTPMLogp1.csv
#   Signatures: shared/signatures.yml (SecA / SecB 7-gene noBCAM sets)
#   Shared helpers: config/config.R
#
#   NOTE: the CANONICAL per-line scores consumed by scripts 02 and 03 are the
#   DEPOSITED cache cfg_obj("ccle_line_scores"). This script REGENERATES that
#   exact file (hgsc_line_scores.tsv) from the raw matrices for
#   provenance/verification; the scoring is deterministic so it reproduces the
#   deposited cache identically.
#
# OUTPUTS (all via cfg_path("figures_dir","figure_ccle_smac_mimetics", ...))
#   hgsc_line_scores.tsv              (one row per HGSC line; SecA/SecB/polarization)
#   hgsc_model_meta.tsv               (metadata for the lines kept)
#   01_polarization_ranking.pdf       (z-mean ranking; primary)
#   01_polarization_ranking_UCell.pdf (UCell ranking; secondary)
#
# MANUSCRIPT PANEL(S): TBD (placeholder — supplementary / mechanistic; panel not yet assigned)
# RUNTIME TIER: moderate (loads the full TPM matrix; ~30 s to read; fast scoring)
# ============================================================================

.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))

set.seed(CFG$seed)

suppressPackageStartupMessages({
  library(data.table)
  library(UCell)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

FIG_DIR <- cfg_path("figures_dir", "figure_ccle_smac_mimetics")

# ---- signatures (single source of truth: shared/signatures.yml) -------------
.sig <- yaml::read_yaml(file.path(.here, "..", "..", "shared", "signatures.yml"))
sigs <- list(SecA = .sig$SecA, SecB = .sig$SecB)

# ---- Model.csv: filter to HGSC (permissive + strict) ------------------------
cat("Loading Model.csv...\n")
model <- fread(cfg_obj("ccle_model_meta"))

# Permissive: OncotreePrimaryDisease == Ovarian Epithelial Tumor AND HGSOC code
perm <- model[OncotreePrimaryDisease == "Ovarian Epithelial Tumor" &
              OncotreeCode == "HGSOC"]
cat(sprintf("Permissive HGSOC set: %d lines\n", nrow(perm)))

# Strict Domcke-2013 set (9 lines)
domcke_names <- c("KURAMOCHI","OVSAHO","COV362","COV318","OVCAR4",
                  "OVCAR3","OVCAR8","CAOV3","SNU119")
# NIH:OVCAR-3 stripped name is "NIHOVCAR3" — handle that
model[, StrippedNorm := toupper(gsub("[-_ :]","", StrippedCellLineName))]
model[StrippedNorm == "NIHOVCAR3", StrippedNorm := "OVCAR3"]
strict <- model[StrippedNorm %in% domcke_names]
cat(sprintf("Strict Domcke-2013 set: %d/9 lines\n", nrow(strict)))
print(strict[, .(StrippedCellLineName, ModelID, OncotreeSubtype, OncotreeCode)])

# Union of model IDs we need to score
all_hgsc_ids <- union(perm$ModelID, strict$ModelID)
cat(sprintf("Total unique HGSC lines to score: %d\n", length(all_hgsc_ids)))

# ---- expression matrix -------------------------------------------------------
cat("\nLoading expression matrix (this takes ~30 s)...\n")
# columns: 'GeneSymbol (entrezID)', rows: ModelID
expr <- fread(cfg_obj("ccle_expression"), header = TRUE)
# first column is ModelID (no header label). data.table sets V1.
setnames(expr, 1, "ModelID")
cat(sprintf("Expression matrix: %d models x %d genes\n",
            nrow(expr), ncol(expr) - 1L))

# Subset to HGSC lines
expr_hgsc <- expr[ModelID %in% all_hgsc_ids]
cat(sprintf("HGSC lines with expression: %d / %d\n",
            nrow(expr_hgsc), length(all_hgsc_ids)))

missing_expr <- setdiff(all_hgsc_ids, expr_hgsc$ModelID)
if (length(missing_expr)) {
  cat("Lines missing expression:\n")
  print(model[ModelID %in% missing_expr, .(StrippedCellLineName, ModelID)])
}

# Clean column names: extract gene symbol (before space-paren)
old_names  <- colnames(expr_hgsc)[-1]
gene_syms  <- sub("\\s*\\(.+\\)\\s*$", "", old_names)
# Detect duplicates after parsing
dup_syms <- gene_syms[duplicated(gene_syms)]
if (length(dup_syms)) {
  cat("Duplicated gene symbols after parsing:", length(unique(dup_syms)),
      "(will keep first)\n")
}

# Build gene x sample numeric matrix
mat <- as.matrix(expr_hgsc[, -1, with = FALSE])
rownames(mat) <- expr_hgsc$ModelID
colnames(mat) <- gene_syms
# Transpose: genes (rows) x samples (cols), as UCell expects
mat <- t(mat)
# De-duplicate genes by keeping the first occurrence
mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
cat(sprintf("Final matrix: %d genes x %d cell lines\n", nrow(mat), ncol(mat)))

# Confirm signature gene coverage
for (s in names(sigs)) {
  present <- intersect(sigs[[s]], rownames(mat))
  cat(sprintf("  %s coverage: %d/%d  (%s)\n",
              s, length(present), length(sigs[[s]]),
              paste(setdiff(sigs[[s]], rownames(mat)), collapse = ",")))
}

# ---- UCell scoring -----------------------------------------------------------
cat("\nRunning UCell scoring (signature rank approach)...\n")
ucell <- ScoreSignatures_UCell(mat, features = sigs, name = "")
ucell <- as.data.frame(ucell)
ucell$ModelID <- rownames(ucell)
setDT(ucell)

# Also compute simple z-scored mean of log-expression as sanity check
z_score_mean <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (!length(g)) return(rep(NA_real_, ncol(mat)))
  zmat <- t(scale(t(mat[g, , drop = FALSE])))  # z across samples per gene
  colMeans(zmat, na.rm = TRUE)
}
scores <- data.table(
  ModelID  = colnames(mat),
  SecA_UCell = ucell$SecA,
  SecB_UCell = ucell$SecB,
  SecA_zmean = z_score_mean(mat, sigs$SecA),
  SecB_zmean = z_score_mean(mat, sigs$SecB)
)

# Z-score across the HGSC line set (per the manuscript convention)
# Note: UCell SecA scores are near-zero for many CCLE lines because of
# rank saturation when several signature genes have low/zero TPM. Report
# both UCell and z-scored mean; treat z-mean as primary for CCLE bulk
# (the manuscript uses UCell on scRNA-seq, where rank scoring is well-defined).
scores[, SecA_UCell_z := scale(SecA_UCell)[,1]]
scores[, SecB_UCell_z := scale(SecB_UCell)[,1]]
scores[, polarization_UCell := SecB_UCell_z - SecA_UCell_z]
scores[, polarization_zmean := SecB_zmean - SecA_zmean]
# Primary polarization = z-mean (more robust for bulk RNA-seq with sparse
# expression of LGR5 / LPAR3 / SOX17 in some lines)
scores[, polarization := polarization_zmean]

# Merge metadata
meta <- model[, .(ModelID, StrippedCellLineName, OncotreeSubtype,
                  OncotreeCode, Sex, PrimaryOrMetastasis,
                  SampleCollectionSite, Age, RRID)]
scores <- merge(scores, meta, by = "ModelID", all.x = TRUE)

# Annotate which sets each line belongs to
scores[, in_permissive := ModelID %in% perm$ModelID]
scores[, in_strict     := ModelID %in% strict$ModelID]

scores <- scores[order(-polarization)]
fwrite(scores, file.path(FIG_DIR, "hgsc_line_scores.tsv"), sep = "\t")
fwrite(meta[ModelID %in% all_hgsc_ids],
       file.path(FIG_DIR, "hgsc_model_meta.tsv"), sep = "\t")

cat("\n--- ranking (top 10 SecB-polarized; primary = z-mean) ---\n")
print(scores[, .(StrippedCellLineName, polarization,
                 polarization_UCell, in_strict)][1:10])
cat("\n--- ranking (top 10 SecA-polarized) ---\n")
print(scores[, .(StrippedCellLineName, polarization,
                 polarization_UCell, in_strict)][(.N-9):.N])

# Sanity check on known reference lines
ref_lines <- c("KURAMOCHI","SKOV3","OVCAR4","OVSAHO","NIHOVCAR3","CAOV3")
cat("\n--- reference line ranking ---\n")
print(scores[toupper(gsub("[-:]","",StrippedCellLineName)) %in% ref_lines,
             .(StrippedCellLineName, polarization, polarization_UCell)])

# Concordance UCell vs z-mean
cor_method <- cor(scores$polarization_UCell, scores$polarization_zmean,
                  method = "spearman", use = "pairwise.complete.obs")
cat(sprintf("\nUCell vs z-mean polarization Spearman rho: %.3f\n", cor_method))

# ---- ranking figure ---------------------------------------------------------
p <- ggplot(scores,
            aes(x = reorder(StrippedCellLineName, polarization),
                y = polarization,
                fill = polarization)) +
  geom_col(alpha = 0.9, colour = "grey20", linewidth = 0.2) +
  geom_point(data = scores[in_strict == TRUE],
             aes(x = StrippedCellLineName, y = polarization),
             colour = "black", shape = 18, size = 2.5) +
  coord_flip() +
  scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c",
                       midpoint = 0, name = "Polarization\n(z SecB - z SecA)") +
  labs(x = NULL, y = "SecB polarization (z-mean of SecB minus z-mean of SecA)",
       title = "HGSC cell lines ranked by SecB polarization",
       subtitle = sprintf("DepMap 24Q4; n=%d HGSC lines; diamonds = Domcke-9", nrow(scores))) +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(size = 7))

ggsave(file.path(FIG_DIR, "01_polarization_ranking.pdf"), p,
       width = 5.5, height = 7)

# Also show UCell ranking
p_u <- ggplot(scores,
            aes(x = reorder(StrippedCellLineName, polarization_UCell),
                y = polarization_UCell,
                fill = polarization_UCell)) +
  geom_col(alpha = 0.9, colour = "grey20", linewidth = 0.2) +
  geom_point(data = scores[in_strict == TRUE],
             aes(x = StrippedCellLineName, y = polarization_UCell),
             colour = "black", shape = 18, size = 2.5) +
  coord_flip() +
  scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c",
                       midpoint = 0, name = "UCell\npolarization (z)") +
  labs(x = NULL, y = "SecB polarization (UCell)",
       title = "HGSC cell lines ranked by SecB polarization (UCell)",
       subtitle = sprintf("DepMap 24Q4; n=%d HGSC lines; diamonds = Domcke-9", nrow(scores))) +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(size = 7))

ggsave(file.path(FIG_DIR, "01_polarization_ranking_UCell.pdf"), p_u,
       width = 5.5, height = 7)

cat("\nDone. Wrote:\n")
cat("  ", file.path(FIG_DIR, "hgsc_line_scores.tsv"), "\n")
cat("  ", file.path(FIG_DIR, "01_polarization_ranking.pdf"), "\n")
cat("[PLACEHOLDER — assign panel letter, then rename per ./README.md]\n")
