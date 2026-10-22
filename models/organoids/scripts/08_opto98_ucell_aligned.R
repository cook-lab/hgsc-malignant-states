# Trimmed copy (2026-10-21) of 202603_organoids/scripts/08_opto98_ucell_aligned.R, the script that ran. Kept:
# opto98_ucell_scores_aligned.csv (Fig 3c-e; the 11,543 cells of the Methods) and seurat_opto98_scored.rds. Removed:
# the creation of an unused figure folder. The untrimmed script is in the lab project folder.

# 08_opto98_ucell_aligned.R
# UCell scoring of OPTO-098 across all conditions (Growth, Media, Treatment)
# using the aligned gene space — identical method to v4 (step 04).
#
# Input:
#   - data/seurat_opto98.rds
#   - /Volumes/CookLab/Sarah/2026_final_atlas/output/18_ucell_atlas/shared_gene_list.txt
#
# Output: output/08_OPTO98_SecB_Conditions/
#         figs/08_OPTO98_SecB_Conditions/

source("scripts/00_setup.R")
library(UCell)

# --- Fix paths (setup.R project_dir is stale) ---
project_dir <- "/Volumes/CookLab/Sarah/2026_organoids"
output_dir  <- file.path(project_dir, "output")

# --- Override SecA genes: remove BCAM (same as v3/v4) ---
secA_genes <- setdiff(secA_genes, "BCAM")
cat("SecA genes (BCAM removed):", paste(secA_genes, collapse = ", "), "\n")
cat("SecB genes:", paste(secB_genes, collapse = ", "), "\n")
cat("Scoring method: UCell on ALIGNED gene space (identical to v4)\n")

# --- Paths ---
input_rds       <- file.path(project_dir, "data/seurat_opto98.rds")
out_dir         <- file.path(output_dir, "08_OPTO98_SecB_Conditions")
shared_genes_file <- "/Volumes/CookLab/Sarah/2026_final_atlas/output/18_ucell_atlas/shared_gene_list.txt"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ===========================================================================
# 1. Load data
# ===========================================================================
cat("\n=== [1] Loading OPTO-098 Seurat object ===\n")
seurat_obj <- readRDS(input_rds)
cat("  Full object:", ncol(seurat_obj), "cells x", nrow(seurat_obj), "genes\n")

# ===========================================================================
# 2. Filter multiplets / negatives
# ===========================================================================
cat("\n=== [2] Filtering multiplets and negatives ===\n")
keep <- !seurat_obj$barcode %in% c("multiplet", "negative") & !is.na(seurat_obj$experiment)
seurat_obj <- seurat_obj[, keep]
cat("  After filtering:", ncol(seurat_obj), "cells\n")
cat("  Experiments:", paste(sort(unique(seurat_obj$experiment)), collapse = ", "), "\n")
cat("  Samples:", paste(sort(unique(seurat_obj$sample_id)), collapse = ", "), "\n")

# ===========================================================================
# 3. Load shared gene list & build aligned count matrix
# ===========================================================================
cat("\n=== [3] Building aligned count matrix ===\n")
shared_genes <- readLines(shared_genes_file)
cat("  Shared genes from atlas intersection:", length(shared_genes), "\n")

# Work with RNA assay counts directly (avoids Seurat v5 Assay5/SCT subsetting issues)
DefaultAssay(seurat_obj) <- "RNA"
if (inherits(seurat_obj[["RNA"]], "Assay5")) {
  seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
}
rna_counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")

available <- shared_genes[shared_genes %in% rownames(rna_counts)]
missing   <- shared_genes[!shared_genes %in% rownames(rna_counts)]
cat("  Available in OPTO-098:", length(available), "\n")
if (length(missing) > 0) {
  cat("  Missing:", length(missing), "\n")
}

# Create a minimal Seurat object with only the aligned genes for UCell
aligned_counts <- rna_counts[available, ]
seurat_aligned <- CreateSeuratObject(counts = aligned_counts, meta.data = seurat_obj@meta.data)
cat("  Aligned matrix:", ncol(seurat_aligned), "cells x", nrow(seurat_aligned), "genes\n")

# ===========================================================================
# 4. Check signature gene availability
# ===========================================================================
cat("\n=== [4] Checking SecA/SecB markers in aligned matrix ===\n")
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
# 5. UCell scoring — identical parameters to v4/atlas
# ===========================================================================
cat("\n=== [5] Computing UCell scores on aligned gene space ===\n")
cat("  Parameters: maxRank=1500, slot='counts' (identical to atlas 18b / organoid v4)\n")

signatures <- list(
  SecA = secA_present,
  SecB = secB_present
)

seurat_aligned <- AddModuleScore_UCell(seurat_aligned, features = signatures,
                                        maxRank = 1500)

# Copy scores back to original object metadata
seurat_obj$SecA_UCell <- seurat_aligned$SecA_UCell
seurat_obj$SecB_UCell <- seurat_aligned$SecB_UCell
seurat_obj$sec_polarization <- seurat_aligned$SecB_UCell - seurat_aligned$SecA_UCell

cat("  SecA UCell range:", round(range(seurat_obj$SecA_UCell), 4), "\n")
cat("  SecB UCell range:", round(range(seurat_obj$SecB_UCell), 4), "\n")
cat("  Polarization range:", round(range(seurat_obj$sec_polarization), 4), "\n")

# ===========================================================================
# 6. Export scores CSV
# ===========================================================================
cat("\n=== [6] Exporting scores ===\n")
scores_export <- data.frame(
  barcode      = colnames(seurat_obj),
  experiment   = seurat_obj$experiment,
  sample_id    = seurat_obj$sample_id,
  SecA_UCell   = seurat_obj$SecA_UCell,
  SecB_UCell   = seurat_obj$SecB_UCell,
  sec_polarization = seurat_obj$sec_polarization,
  phase        = seurat_obj$Phase,
  seurat_clusters = seurat_obj$seurat_clusters,
  stringsAsFactors = FALSE
)
write.csv(scores_export, file.path(out_dir, "opto98_ucell_scores_aligned.csv"),
          row.names = FALSE)

# ===========================================================================
# 7. Save scored object
# ===========================================================================
cat("\n=== [7] Saving scored Seurat object ===\n")
saveRDS(seurat_obj, file.path(out_dir, "seurat_opto98_scored.rds"))

cat("\n=== DONE (OPTO-098 UCell scoring — aligned gene space) ===\n")
cat("Scores CSV:", file.path(out_dir, "opto98_ucell_scores_aligned.csv"), "\n")
cat("Object:    ", file.path(out_dir, "seurat_opto98_scored.rds"), "\n")
