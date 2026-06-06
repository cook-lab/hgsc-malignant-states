# 04_ucell_atlas_scoring.R
# UCell scoring on atlas secretory epithelial cells (cross-platform, noBCAM)
# ====================================================================
# HGSC malignant-states atlas backend.
#
# Purpose: load the atlas secretory counts exported by 03_ucell_atlas_export.py,
# intersect the gene space with the organoid Seurat object so the UCell ranking
# universe is identical, then score SecA/SecB with the canonical noBCAM 7-gene
# signatures (shared/signatures.yml). This is the OLD 18b scoring that matches
# the xenium noBCAM annotation chain (NOT the divergent 18b_v2).
#
# INPUTS:
#   - output_root/03_epithelial_nmf/ucell_atlas/atlas_secretory_{counts.mtx.gz,
#     barcodes.tsv,genes.tsv}
#   - <data_root>/2026_organoids/output/01_Data_Processing_and_QC/
#       seurat_untreated_baseline.rds   (gene-space reference)
#   - shared/signatures.yml             (SecA/SecB noBCAM 7-gene sets)
#
# OUTPUTS (output_root/03_epithelial_nmf/ucell_atlas/):
#   - atlas_ucell_scores.csv     (barcode, SecA_UCell, SecB_UCell, sec_polarization)
#   - shared_gene_list.txt
#
# MANUSCRIPT PANELS: atlas-side cross-platform UCell scoring underpinning the
#   organoid SecB comparison (Fig 3B) and the WT/TMA polarization strip (SF11).
#
# RUNTIME TIER: moderate (Seurat object load + UCell on secretory subset).
#
# SEEDING: UCell ranking is deterministic; set.seed(CFG$seed) for safety.
#
# Usage:
#   Rscript 04_ucell_atlas_scoring.R

# config is 3 levels up: atlas/03_epithelial_nmf/04_*.R -> repo root
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
repo_root <- normalizePath(file.path(dirname(this_file), "..", ".."))
source(file.path(repo_root, "config", "config.R"))
set.seed(CFG$seed)

library(Matrix)
library(Seurat)
library(UCell)
library(yaml)

# ── Paths ──────────────────────────────────────────────────────
out_dir      <- cfg_path("output_root", "03_epithelial_nmf", "ucell_atlas")
organoid_rds <- file.path(path.expand(CFG$paths$data_root), "2026_organoids",
                          "output", "01_Data_Processing_and_QC",
                          "seurat_untreated_baseline.rds")

# ── Gene signatures (shared/signatures.yml — noBCAM 7-gene set) ──
sigs <- yaml::read_yaml(file.path(repo_root, "shared", "signatures.yml"))
secA_genes <- as.character(sigs$SecA)
secB_genes <- as.character(sigs$SecB)

cat(strrep("=", 60), "\n", sep = "")
cat("  UCell scoring on atlas secretory cells (noBCAM)\n")
cat(strrep("=", 60), "\n", sep = "")

# ── 1. Organoid gene list ──────────────────────────────────────
cat("\n[1] Loading organoid Seurat object for gene list...\n")
organoid_obj <- readRDS(organoid_rds)
organoid_genes <- rownames(organoid_obj)
cat("    Organoid genes:", length(organoid_genes), "\n")
rm(organoid_obj); gc()

# ── 2. Atlas exported data ─────────────────────────────────────
cat("\n[2] Loading atlas secretory counts...\n")
barcodes <- read.table(file.path(out_dir, "atlas_secretory_barcodes.tsv"),
                       stringsAsFactors = FALSE)$V1
genes    <- read.table(file.path(out_dir, "atlas_secretory_genes.tsv"),
                       stringsAsFactors = FALSE)$V1
counts   <- readMM(file.path(out_dir, "atlas_secretory_counts.mtx.gz"))
rownames(counts) <- genes      # ReadMM gives genes x cells (as exported by 03)
colnames(counts) <- barcodes
cat("    Atlas matrix:", nrow(counts), "genes x", ncol(counts), "cells\n")

# ── 3. Intersect gene spaces ───────────────────────────────────
cat("\n[3] Computing gene intersection...\n")
shared_genes <- intersect(genes, organoid_genes)
cat("    Atlas genes:    ", length(genes), "\n")
cat("    Organoid genes: ", length(organoid_genes), "\n")
cat("    Shared genes:   ", length(shared_genes), "\n")

secA_in_shared <- secA_genes[secA_genes %in% shared_genes]
secA_missing   <- secA_genes[!secA_genes %in% shared_genes]
secB_in_shared <- secB_genes[secB_genes %in% shared_genes]
secB_missing   <- secB_genes[!secB_genes %in% shared_genes]
cat("    SecA genes in shared set:", length(secA_in_shared), "/", length(secA_genes), "\n")
if (length(secA_missing) > 0) cat("      Missing:", paste(secA_missing, collapse = ", "), "\n")
cat("    SecB genes in shared set:", length(secB_in_shared), "/", length(secB_genes), "\n")
if (length(secB_missing) > 0) cat("      Missing:", paste(secB_missing, collapse = ", "), "\n")
if (length(secA_in_shared) < length(secA_genes) ||
    length(secB_in_shared) < length(secB_genes)) {
  warning("Some signature genes are missing from the shared gene space!")
}

counts_shared <- counts[shared_genes, ]
cat("    Subsetted matrix:", nrow(counts_shared), "genes x", ncol(counts_shared), "cells\n")
writeLines(shared_genes, file.path(out_dir, "shared_gene_list.txt"))
rm(counts); gc()

# ── 4. Seurat object ───────────────────────────────────────────
cat("\n[4] Creating Seurat object...\n")
seurat_obj <- CreateSeuratObject(counts = counts_shared, project = "atlas_secretory")
cat("    Seurat object:", ncol(seurat_obj), "cells x", nrow(seurat_obj), "genes\n")
rm(counts_shared); gc()

# ── 5. UCell ───────────────────────────────────────────────────
cat("\n[5] Running UCell scoring (maxRank=1500, raw counts)...\n")
signatures <- list(SecA = secA_in_shared, SecB = secB_in_shared)
cat("    SecA signature:", paste(secA_in_shared, collapse = ", "), "\n")
cat("    SecB signature:", paste(secB_in_shared, collapse = ", "), "\n")

seurat_obj <- AddModuleScore_UCell(seurat_obj, features = signatures, maxRank = 1500)
seurat_obj$SecA_score <- seurat_obj$SecA_UCell
seurat_obj$SecB_score <- seurat_obj$SecB_UCell
seurat_obj$sec_polarization <- seurat_obj$SecB_score - seurat_obj$SecA_score

cat("    SecA UCell range:", round(range(seurat_obj$SecA_score), 4), "\n")
cat("    SecB UCell range:", round(range(seurat_obj$SecB_score), 4), "\n")
cat("    Polarization range:", round(range(seurat_obj$sec_polarization), 4), "\n")

# ── 6. Export scores ───────────────────────────────────────────
cat("\n[6] Exporting scores...\n")
scores_df <- data.frame(
  barcode          = colnames(seurat_obj),
  SecA_UCell       = seurat_obj$SecA_score,
  SecB_UCell       = seurat_obj$SecB_score,
  sec_polarization = seurat_obj$sec_polarization,
  stringsAsFactors = FALSE
)
write.csv(scores_df, file.path(out_dir, "atlas_ucell_scores.csv"), row.names = FALSE)
cat("    Saved: atlas_ucell_scores.csv (", nrow(scores_df), "rows )\n")

cat("\n", rep("=", 60), "\n", sep = "")
cat("  Step complete! Output:", out_dir, "\n")
cat(rep("=", 60), "\n", sep = "")
