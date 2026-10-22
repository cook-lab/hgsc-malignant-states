# 18b_ucell_atlas_scoring.R
# UCell scoring on atlas secretory epithelial cells
# ====================================================================
# HGSC Single-Cell Atlas  ·  Cook Lab  ·  2026
#
# Loads atlas secretory counts (from 18a), intersects the gene space with
# the organoid Seurat object to ensure identical ranking universes, then
# runs UCell with the same parameters as the organoid pipeline.
#
# Input:
#   - output/18_ucell_atlas/atlas_secretory_counts.mtx.gz
#   - output/18_ucell_atlas/atlas_secretory_barcodes.tsv
#   - output/18_ucell_atlas/atlas_secretory_genes.tsv
#   - /Volumes/CookLab/Sarah/organoids/seurat_untreated_baseline.rds
#
# Output (output/18_ucell_atlas/):
#   - atlas_ucell_scores.csv
#   - shared_gene_list.txt
#
# Usage:
#   Rscript scripts/18b_ucell_atlas_scoring.R

library(Matrix)
library(Seurat)
library(UCell)

# ── Paths ──────────────────────────────────────────────────────
project_dir  <- "/Volumes/CookLab/Sarah/2026_final_atlas"
out_dir      <- file.path(project_dir, "output", "18_ucell_atlas")
organoid_rds <- "/Volumes/CookLab/Sarah/organoids/seurat_untreated_baseline.rds"

# ── Gene signatures (matching organoid v3: BCAM removed from SecA) ──
secA_genes <- c("MECOM", "FBXO21", "LGR5", "LPAR3", "PBX1", "SOX17", "RCN2")
secB_genes <- c("KRT17", "KRT19", "KRT7", "LCN2", "PRSS22", "SLPI", "TACSTD2")

# ============================================================================
# 1. LOAD ORGANOID GENE LIST
# ============================================================================
cat("=" , rep("=", 59), "\n", sep = "")
cat("  Step 18b — UCell scoring on atlas secretory cells\n")
cat("=" , rep("=", 59), "\n", sep = "")

cat("\n[1] Loading organoid Seurat object for gene list...\n")
organoid_obj <- readRDS(organoid_rds)
organoid_genes <- rownames(organoid_obj)
cat("    Organoid genes:", length(organoid_genes), "\n")
rm(organoid_obj)
gc()

# ============================================================================
# 2. LOAD ATLAS EXPORTED DATA
# ============================================================================
cat("\n[2] Loading atlas secretory counts...\n")

barcodes <- read.table(file.path(out_dir, "atlas_secretory_barcodes.tsv"),
                        stringsAsFactors = FALSE)$V1
genes    <- read.table(file.path(out_dir, "atlas_secretory_genes.tsv"),
                        stringsAsFactors = FALSE)$V1
counts   <- readMM(file.path(out_dir, "atlas_secretory_counts.mtx.gz"))

# ReadMM gives genes × cells (as exported by 18a)
rownames(counts) <- genes
colnames(counts) <- barcodes
cat("    Atlas matrix:", nrow(counts), "genes ×", ncol(counts), "cells\n")

# ============================================================================
# 3. INTERSECT GENE SPACES
# ============================================================================
cat("\n[3] Computing gene intersection...\n")

atlas_genes  <- genes
shared_genes <- intersect(atlas_genes, organoid_genes)
cat("    Atlas genes:    ", length(atlas_genes), "\n")
cat("    Organoid genes: ", length(organoid_genes), "\n")
cat("    Shared genes:   ", length(shared_genes), "\n")

# Check signature genes
secA_in_shared <- secA_genes[secA_genes %in% shared_genes]
secA_missing   <- secA_genes[!secA_genes %in% shared_genes]
secB_in_shared <- secB_genes[secB_genes %in% shared_genes]
secB_missing   <- secB_genes[!secB_genes %in% shared_genes]

cat("    SecA genes in shared set:", length(secA_in_shared), "/", length(secA_genes), "\n")
if (length(secA_missing) > 0) cat("      Missing:", paste(secA_missing, collapse = ", "), "\n")
cat("    SecB genes in shared set:", length(secB_in_shared), "/", length(secB_genes), "\n")
if (length(secB_missing) > 0) cat("      Missing:", paste(secB_missing, collapse = ", "), "\n")

if (length(secA_in_shared) < length(secA_genes) || length(secB_in_shared) < length(secB_genes)) {
  warning("Some signature genes are missing from the shared gene space!")
}

# Subset counts to shared genes
counts_shared <- counts[shared_genes, ]
cat("    Subsetted matrix:", nrow(counts_shared), "genes ×", ncol(counts_shared), "cells\n")

# Save shared gene list for organoid re-scoring
writeLines(shared_genes, file.path(out_dir, "shared_gene_list.txt"))
cat("    Saved: shared_gene_list.txt (", length(shared_genes), " genes)\n")

# Clean up
rm(counts)
gc()

# ============================================================================
# 4. CREATE SEURAT OBJECT
# ============================================================================
cat("\n[4] Creating Seurat object...\n")

seurat_obj <- CreateSeuratObject(counts = counts_shared, project = "atlas_secretory")
cat("    Seurat object:", ncol(seurat_obj), "cells ×", nrow(seurat_obj), "genes\n")

rm(counts_shared)
gc()

# ============================================================================
# 5. RUN UCELL
# ============================================================================
cat("\n[5] Running UCell scoring...\n")
cat("    Parameters: maxRank=1500, slot='counts' (raw), defaults for w_neg & ties\n")

signatures <- list(
  SecA = secA_in_shared,
  SecB = secB_in_shared
)

cat("    SecA signature:", paste(secA_in_shared, collapse = ", "), "\n")
cat("    SecB signature:", paste(secB_in_shared, collapse = ", "), "\n")

seurat_obj <- AddModuleScore_UCell(seurat_obj, features = signatures,
                                    maxRank = 1500)

# Extract scores
seurat_obj$SecA_score <- seurat_obj$SecA_UCell
seurat_obj$SecB_score <- seurat_obj$SecB_UCell
seurat_obj$sec_polarization <- seurat_obj$SecB_score - seurat_obj$SecA_score

cat("    SecA UCell score range:", round(range(seurat_obj$SecA_score), 4), "\n")
cat("    SecB UCell score range:", round(range(seurat_obj$SecB_score), 4), "\n")
cat("    Polarization range:    ", round(range(seurat_obj$sec_polarization), 4), "\n")

# ============================================================================
# 6. EXPORT SCORES
# ============================================================================
cat("\n[6] Exporting scores...\n")

scores_df <- data.frame(
  barcode          = colnames(seurat_obj),
  SecA_UCell       = seurat_obj$SecA_score,
  SecB_UCell       = seurat_obj$SecB_score,
  sec_polarization = seurat_obj$sec_polarization,
  stringsAsFactors = FALSE
)

scores_path <- file.path(out_dir, "atlas_ucell_scores.csv")
write.csv(scores_df, scores_path, row.names = FALSE)
cat("    Saved:", scores_path, "\n")
cat("    Rows:", nrow(scores_df), "\n")

# Quick summary
cat("\n    === Quick Summary ===\n")
cat("    SecA score: mean=", round(mean(scores_df$SecA_UCell), 4),
    " median=", round(median(scores_df$SecA_UCell), 4), "\n")
cat("    SecB score: mean=", round(mean(scores_df$SecB_UCell), 4),
    " median=", round(median(scores_df$SecB_UCell), 4), "\n")
cat("    Polarization: mean=", round(mean(scores_df$sec_polarization), 4),
    " median=", round(median(scores_df$sec_polarization), 4), "\n")

cat("\n", rep("=", 60), "\n", sep = "")
cat("  Step 18b complete!\n")
cat("  Output: ", out_dir, "\n")
cat(rep("=", 60), "\n", sep = "")
