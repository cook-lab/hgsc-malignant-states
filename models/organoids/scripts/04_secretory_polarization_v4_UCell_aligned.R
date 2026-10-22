# Trimmed copy (2026-10-21) of 202603_organoids/scripts/04_secretory_polarization_v4_UCell_aligned.R, the script
# that ran. Kept: organoid_ucell_scores_v4_aligned.csv (read by 04v4 and 04v5) and
# seurat_secretory_scored_v4_aligned.rds. Removed: the atlas-score load, the per-model summary table and the
# figures. The untrimmed script, which the header below describes, is in the lab project folder.

# 04_secretory_polarization_v4_UCell_aligned.R
# v4: Re-score organoid cells with UCell using the ALIGNED gene space
#     (intersection of atlas 14,470 genes and organoid 34,526 genes = 14,271 shared)
#     so that scores are directly comparable to atlas UCell scores from step 18b.
#
# Changes from v3:
#   - Subsets the count matrix to 14,271 shared genes BEFORE UCell scoring
#   - Adds atlas-derived summary statistics for contextual interpretation
#   - Generates an HTML report comparing organoid scores against atlas distributions
#
# Input:
#   - output/01_Data_Processing_and_QC/seurat_untreated_baseline.rds
#   - /Volumes/CookLab/Sarah/2026_final_atlas/output/18_ucell_atlas/shared_gene_list.txt
#   - /Volumes/CookLab/Sarah/2026_final_atlas/output/18_ucell_atlas/atlas_ucell_scores.csv
#   - /Volumes/CookLab/Sarah/2026_final_atlas/output/18_ucell_atlas/atlas_secretory_metadata.csv
#
# Output: output/02_Secretory_Polarization_v4_UCell_aligned/
#         figs/02_Secretory_Polarization_v4_UCell_aligned/

source("scripts/00_setup.R")
library(UCell)

# --- Override SecA genes: remove BCAM (same as v3) ---
secA_genes <- setdiff(secA_genes, "BCAM")
cat("v4 SecA genes (BCAM removed):", paste(secA_genes, collapse = ", "), "\n")
cat("v4 SecB genes:", paste(secB_genes, collapse = ", "), "\n")
cat("Scoring method: UCell on ALIGNED gene space (14,271 shared genes)\n")

# --- Paths ---
input_rds    <- file.path(output_dir, "01_Data_Processing_and_QC/seurat_untreated_baseline.rds")
out_dir      <- file.path(output_dir, "02_Secretory_Polarization_v4_UCell_aligned")
shared_genes_file <- "/Volumes/CookLab/Sarah/2026_final_atlas/output/18_ucell_atlas/shared_gene_list.txt"

if (!dir.exists(out_dir))  dir.create(out_dir, recursive = TRUE)

# ===========================================================================
# 1. Load data
# ===========================================================================
cat("\n=== [1] Loading Seurat object ===\n")
seurat_obj <- readRDS(input_rds)
cat("  ", ncol(seurat_obj), "cells x", nrow(seurat_obj), "genes (full matrix)\n")

# ===========================================================================
# 2. Load shared gene list & subset
# ===========================================================================
cat("\n=== [2] Loading shared gene list and subsetting ===\n")
shared_genes <- readLines(shared_genes_file)
cat("  Shared genes from atlas intersection:", length(shared_genes), "\n")

# Check overlap
available <- shared_genes[shared_genes %in% rownames(seurat_obj)]
missing   <- shared_genes[!shared_genes %in% rownames(seurat_obj)]
cat("  Available in organoid object:", length(available), "\n")
if (length(missing) > 0) {
  cat("  Missing:", length(missing), "\n")
  cat("    ", paste(head(missing, 10), collapse = ", "), "\n")
}

# Subset to shared genes
seurat_aligned <- subset(seurat_obj, features = available)
cat("  Aligned matrix:", ncol(seurat_aligned), "cells x", nrow(seurat_aligned), "genes\n")

# ===========================================================================
# 3. Check signature gene availability
# ===========================================================================
cat("\n=== [3] Checking SecA/SecB markers in aligned matrix ===\n")
aligned_genes <- rownames(seurat_aligned)

secA_present <- secA_genes[secA_genes %in% aligned_genes]
secA_missing <- secA_genes[!secA_genes %in% aligned_genes]
secB_present <- secB_genes[secB_genes %in% aligned_genes]
secB_missing <- secB_genes[!secB_genes %in% aligned_genes]

cat("  SecA genes present:", length(secA_present), "/", length(secA_genes), "\n")
if (length(secA_missing) > 0) cat("    Missing:", paste(secA_missing, collapse = ", "), "\n")
cat("  SecB genes present:", length(secB_present), "/", length(secB_genes), "\n")
if (length(secB_missing) > 0) cat("    Missing:", paste(secB_missing, collapse = ", "), "\n")

# ===========================================================================
# 4. UCell scoring on aligned matrix
# ===========================================================================
cat("\n=== [4] Computing UCell scores on aligned gene space ===\n")
cat("  Parameters: maxRank=1500, slot='counts' (identical to atlas 18b)\n")

signatures <- list(
  SecA = secA_present,
  SecB = secB_present
)

seurat_aligned <- AddModuleScore_UCell(seurat_aligned, features = signatures,
                                        maxRank = 1500)

# Copy scores back to original object metadata for convenience
seurat_obj$SecA_score_aligned <- seurat_aligned$SecA_UCell
seurat_obj$SecB_score_aligned <- seurat_aligned$SecB_UCell
seurat_obj$sec_polarization_aligned <- seurat_aligned$SecB_UCell - seurat_aligned$SecA_UCell

cat("  SecA UCell range:", round(range(seurat_obj$SecA_score_aligned), 4), "\n")
cat("  SecB UCell range:", round(range(seurat_obj$SecB_score_aligned), 4), "\n")
cat("  Polarization range:", round(range(seurat_obj$sec_polarization_aligned), 4), "\n")

# ===========================================================================
# 8. Save updated object
# ===========================================================================
cat("\n=== [8] Saving updated Seurat object ===\n")
saveRDS(seurat_obj, file.path(out_dir, "seurat_secretory_scored_v4_aligned.rds"))

# ===========================================================================
# 9. Export scores CSV for Python report
# ===========================================================================
cat("\n=== [9] Exporting scores for report ===\n")
scores_export <- data.frame(
  barcode = colnames(seurat_obj),
  SecA_UCell = seurat_obj$SecA_score_aligned,
  SecB_UCell = seurat_obj$SecB_score_aligned,
  sec_polarization = seurat_obj$sec_polarization_aligned,
  model = seurat_obj$pdo_model,
  stringsAsFactors = FALSE
)
write.csv(scores_export, file.path(out_dir, "organoid_ucell_scores_v4_aligned.csv"),
          row.names = FALSE)

cat("\n=== DONE (v4: aligned gene space UCell) ===\n")
cat("Scores:   ", file.path(out_dir, "organoid_ucell_scores_v4_aligned.csv"), "\n")
cat("Object:   ", file.path(out_dir, "seurat_secretory_scored_v4_aligned.rds"), "\n")
