# Trimmed copy (2026-10-21) of 202603_organoids/scripts/10a_geo_cluster_annotate.R, the script that ran. Kept:
# 10a_seurat_epithelial.rds (read by 12), its cached clustering and the per-cluster module scores that the Section 8
# labels were checked against. Removed: the QC table, dot plot, FindAllMarkers and first-pass label tables, the
# UMAPs, and the full annotated object and its label table. The untrimmed script is in the lab project folder.

# 10a_geo_cluster_annotate.R
# =============================================================================
# Step 10a — GSE299200 (public HGSC tumor-immune scRNA-seq):
#   load -> QC -> integrate -> cluster -> annotate to the authors' scheme ->
#   isolate the epithelial/HGSC compartment for SecA/SecB scoring (10b).
#
# Dataset: primary HGSC tumor cells + autologous PBMCs, 3 samples:
#   PBMCs_control (normoxia, immune-only)  <- internal control for immune calls
#   Tumor_PBMC_hypoxia  (tECMHyp)
#   Tumor_PBMC_physoxia (nECMPhys)
#
# QC matches the organoid baseline (scripts/03): drop lncRNA/miRNA/snoRNA/
# ribosomal (keep mito) -> scDblFinder per sample ->
# nCount>500 & nFeature>300 & percent.mt<20 & singlet.
# Verified pass rates: PBMC 95.9%, hypoxia 87.5%, physoxia 95.7% (no depth bias).
#
# Annotation: cell types matched to Plesselova et al. (bioRxiv 2025.06.23.661098) —
# Epithelial(HGSC)/CAF/CD4 T/CD8 T/NK/B/Myeloid(+Cycling). Markers in
# scripts/10_markers_config.R. Epithelial anchored on EPCAM/PAX8/WT1 (keratins
# reserved for SecB scoring).
#
# Output: output/10_HGSC_public_SecAB_scoring/ , figs/10_HGSC_public_SecAB_scoring/
#
# NOTE: cluster -> cell-type assignment REQUIRES manual review (Section 8).
#       Run through Section 7, inspect the diagnostics, edit `cluster_labels`,
#       then re-run Sections 8-10.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(dplyr)
})

set.seed(42)

# --- Paths (self-contained) ---
# Source 00_setup FIRST for theme_lab + palettes, THEN set our paths — 00_setup
# resets project_dir to a stale value, so it must not run after our definitions.
source("/Volumes/CookLab/Sarah/202603_organoids/scripts/00_setup.R")

project_dir <- "/Volumes/CookLab/Sarah/202603_organoids"
geo_dir     <- file.path(project_dir, "data/GSE299200/filtered")
out_dir     <- file.path(project_dir, "output/10_HGSC_public_SecAB_scoring")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

source(file.path(project_dir, "scripts/10_markers_config.R"))

# --- Sample map: folder -> (sample id, condition) ---
sample_map <- data.frame(
  folder    = c("PBMCs_control", "Tumor_PBMC_hypoxia", "Tumor_PBMC_physoxia"),
  sample    = c("PBMCs_control", "Tumor_PBMC_hypoxia", "Tumor_PBMC_physoxia"),
  condition = c("normoxia_ctrl", "hypoxia",            "physoxia"),
  gsm       = c("GSM9035123",    "GSM9035124",         "GSM9035125"),
  stringsAsFactors = FALSE
)

RESOLUTION     <- 1.0
unlabeled_path <- file.path(out_dir, "10a_seurat_clustered_UNLABELED.rds")

# =============================================================================
# 1-6. Build the integrated object ONCE and cache it. Reusing the cached graph
#      guarantees reproducible cluster numbering (so the reviewed label map below
#      stays valid) and makes re-labeling fast.
# =============================================================================
if (!file.exists(unlabeled_path)) {

# =============================================================================
# 1. Load 10x matrices -> per-sample Seurat objects -> merge
# =============================================================================
cat("=== [1] Loading GSE299200 10x matrices ===\n")

obj_list <- lapply(seq_len(nrow(sample_map)), function(i) {
  s <- sample_map[i, ]
  mtx <- Read10X(data.dir = file.path(geo_dir, s$folder))
  o <- CreateSeuratObject(counts = mtx, project = s$sample, min.cells = 0, min.features = 0)
  o$sample    <- s$sample
  o$condition <- s$condition
  o$gsm       <- s$gsm
  cat(sprintf("  %-22s %6d cells x %d genes\n", s$sample, ncol(o), nrow(o)))
  o
})
names(obj_list) <- sample_map$sample

seurat_obj <- merge(obj_list[[1]], y = obj_list[-1],
                    add.cell.ids = sample_map$sample, project = "GSE299200")
cat("  Merged:", ncol(seurat_obj), "cells x", nrow(seurat_obj), "genes\n")

# =============================================================================
# 2. Gene filtering (same logic as baseline script 03)
# =============================================================================
cat("\n=== [2] Gene filtering ===\n")
g <- rownames(seurat_obj)
has_symbol    <- !grepl("^ENSG[0-9]+$", g)
is_lncrna     <- grepl("^LINC|^LOC[0-9]", g)
is_mirna      <- grepl("^MIR[0-9]", g)
is_snorna     <- grepl("^SNOR|^RNU[0-9]|^SCARNA", g)
is_ribosomal  <- grepl("^RP[SL][0-9]", g)
keep_genes    <- has_symbol & !is_lncrna & !is_mirna & !is_snorna & !is_ribosomal
cat("  Genes before:", length(g), " -> kept:", sum(keep_genes),
    " (mito retained:", sum(grepl("^MT-", g)), ")\n")
seurat_obj <- subset(seurat_obj, features = g[keep_genes])

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# =============================================================================
# 3. Doublet detection (scDblFinder, per sample)
# =============================================================================
cat("\n=== [3] Doublet detection (scDblFinder per sample) ===\n")
sce <- as.SingleCellExperiment(JoinLayers(seurat_obj))
sce <- scDblFinder(sce, samples = "sample")
# scDblFinder outputs are UNNAMED vectors in cell order -> name them before
# matching by barcode (indexing an unnamed vector by name yields all NA).
dbl_class <- setNames(as.character(sce$scDblFinder.class), colnames(sce))
dbl_score <- setNames(sce$scDblFinder.score,             colnames(sce))
seurat_obj$scDblFinder_class <- dbl_class[colnames(seurat_obj)]
seurat_obj$scDblFinder_score <- dbl_score[colnames(seurat_obj)]
stopifnot(!anyNA(seurat_obj$scDblFinder_class))
print(seurat_obj@meta.data %>% group_by(sample) %>%
        summarize(n = n(), doublets = sum(scDblFinder_class == "doublet"),
                  pct = round(100 * doublets / n, 1), .groups = "drop") %>% as.data.frame())

# =============================================================================
# 4. Cell-level QC filtering (baseline thresholds)
# =============================================================================
cat("\n=== [4] Cell QC filtering ===\n")
NCOUNT_MIN <- 500; NFEATURE_MIN <- 300; PCT_MT_MAX <- 20
n_before <- ncol(seurat_obj)
keep <- seurat_obj$nCount_RNA   >  NCOUNT_MIN &
        seurat_obj$nFeature_RNA >  NFEATURE_MIN &
        seurat_obj$percent.mt   <  PCT_MT_MAX &
        seurat_obj$scDblFinder_class == "singlet"
seurat_obj <- seurat_obj[, keep]
cat(sprintf("  Retained %d / %d cells (%.1f%%)\n", ncol(seurat_obj), n_before,
            100 * ncol(seurat_obj) / n_before))
cat("  Per sample after QC:\n"); print(table(seurat_obj$sample))

# =============================================================================
# 5. Normalize, HVGs, PCA, integrate across samples (Seurat v5 layers)
# =============================================================================
cat("\n=== [5] Normalize / HVG / PCA / integrate ===\n")
seurat_obj <- JoinLayers(seurat_obj)   # collapse merge-split layers first
seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$sample)
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst",
                                   nfeatures = 2000, verbose = FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)

# Integrate to remove condition-driven splitting of shared cell types.
# CCA is built into Seurat (no extra deps). Swap to HarmonyIntegration if preferred.
seurat_obj <- IntegrateLayers(seurat_obj, method = CCAIntegration,
                              orig.reduction = "pca",
                              new.reduction = "integrated.cca", verbose = FALSE)
RED <- "integrated.cca"
seurat_obj <- JoinLayers(seurat_obj)   # rejoin for DE / plotting

# =============================================================================
# 6. Neighbors + UMAP (part of the reusable integrated object)
# =============================================================================
cat("\n=== [6] Neighbors + UMAP ===\n")
seurat_obj <- FindNeighbors(seurat_obj, reduction = "integrated.cca", dims = 1:30, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, reduction = "integrated.cca", dims = 1:30, verbose = FALSE)
saveRDS(seurat_obj, unlabeled_path)
cat("  Cached integrated object:", unlabeled_path, "\n")

} else {
  cat("=== [1-6] Loading cached integrated object:", unlabeled_path, "===\n")
  seurat_obj <- readRDS(unlabeled_path)
}

# --- Cluster at the chosen resolution (deterministic on the fixed graph) ---
cat("\n=== [6b] Clustering at resolution", RESOLUTION, "===\n")
seurat_obj <- FindClusters(seurat_obj, resolution = RESOLUTION, verbose = FALSE)
cat("  Clusters found:", length(levels(seurat_obj$seurat_clusters)), "\n")

# =============================================================================
# 7. Annotation diagnostics  (INSPECT THESE before assigning labels)
# =============================================================================
cat("\n=== [7] Annotation diagnostics ===\n")

# 7a. Per-cell compartment module scores (mean of markers present)
for (ct in names(annotation_markers)) {
  present <- intersect(annotation_markers[[ct]], rownames(seurat_obj))
  if (length(present) == 0) next
  seurat_obj <- AddModuleScore(seurat_obj, features = list(present),
                               name = paste0("mod_", ct), seed = 42)
  # AddModuleScore appends "1" -> rename to mod_<ct>
  nm <- paste0("mod_", ct, "1")
  seurat_obj[[paste0("modscore_", ct)]] <- seurat_obj[[nm]]; seurat_obj[[nm]] <- NULL
}

modcols <- paste0("modscore_", names(annotation_markers))
modcols <- modcols[modcols %in% colnames(seurat_obj@meta.data)]

# 7b. Per-cluster mean module score -> heatmap + table (evidence for labels)
clust_mod <- seurat_obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarize(across(all_of(modcols), ~ round(mean(.x), 3)), n = n(), .groups = "drop")
write.csv(clust_mod, file.path(out_dir, "10a_cluster_module_scores.csv"), row.names = FALSE)
cat("  Wrote per-cluster module scores.\n"); print(as.data.frame(clust_mod))

# =============================================================================
# 8. Cluster -> cell-type assignment (reviewed from Section 7, resolution 1.0)
# =============================================================================
# Reviewed against per-cluster marker means (2026-07-14). 15 clusters at res 1.0:
#   CD4_T: CD3D+/IL7R+/CD4+    CD8_T: CD3D+/CD8A/B+/GZMB   NK: GNLY/NKG7 hi, CD3 lo
#   B: MS4A1/CD79A/CD19        Myeloid: LYZ/CD14/CD68/FCGR3A
#   CAF: COL1A1/DCN/ACTA2, PTPRC~0    Epithelial: EPCAM/KRT8/PAX8/WT1/MUC16, PTPRC~0
cluster_labels <- c(
  "0"="CD4_T", "1"="Myeloid",    "2"="CD8_T",   "3"="B",   "4"="CD4_T",
  "5"="B",     "6"="Myeloid",    "7"="CAF",     "8"="CD4_T","9"="NK",
  "10"="CD4_T","11"="Epithelial","12"="Myeloid","13"="CAF","14"="B"
)
stopifnot(all(levels(seurat_obj$seurat_clusters) %in% names(cluster_labels)))
# unname(): a named vector assigned to $celltype is matched to cell barcodes by
# Seurat and errors ("No cell overlap"); we want positional assignment.
seurat_obj$celltype <- factor(unname(cluster_labels[as.character(seurat_obj$seurat_clusters)]))

# Guard against cluster-number drift if the object is ever rebuilt: the cluster
# we call Epithelial must actually be EPCAM-high / PTPRC-low.
epi_chk <- AverageExpression(seurat_obj, features = c("EPCAM","PTPRC"),
             group.by = "celltype", assays = "RNA", layer = "data")$RNA
cat("\n  Epithelial sanity (EPCAM high, PTPRC low expected):\n")
print(round(epi_chk[, "Epithelial", drop = FALSE], 2))
if (epi_chk["EPCAM","Epithelial"] < 0.5 || epi_chk["PTPRC","Epithelial"] > 3)
  warning("Epithelial cluster failed EPCAM-high/PTPRC-low check — RE-REVIEW labels!")
cat("\n=== [8] Cell types assigned ===\n"); print(table(seurat_obj$celltype, seurat_obj$condition))

# 8b. Validation: PBMCs_control should be immune-only (no Epithelial/CAF)
pbmc_ct <- table(seurat_obj$celltype[seurat_obj$sample == "PBMCs_control"])
cat("\n  PBMCs_control composition (expect ~0 Epithelial/CAF):\n"); print(pbmc_ct)

# =============================================================================
# 9. Isolate epithelial / HGSC compartment
# =============================================================================
cat("\n=== [9] Subsetting epithelial compartment ===\n")
epi <- subset(seurat_obj, subset = celltype == "Epithelial")
cat("  Epithelial cluster cells:", ncol(epi), "\n")
cat("  By condition (raw):\n"); print(table(epi$condition))
# Purity: the PBMC-only control cannot contain tumor epithelium, so any
# control-derived cells here are ambient/doublets -> drop them, leaving only
# genuine tumor cells (hypoxia + physoxia) to feed 10b.
epi <- subset(epi, subset = condition != "normoxia_ctrl")
epi$condition <- droplevels(factor(epi$condition))
cat("  After dropping control-derived (ambient) cells:", ncol(epi), "\n")
cat("  Epithelial for scoring (feeds 10b), by condition:\n"); print(table(epi$condition))

# =============================================================================
# 10. Save
# =============================================================================
cat("\n=== [10] Saving ===\n")
saveRDS(epi,        file.path(out_dir, "10a_seurat_epithelial.rds"))

cat("\n=== 10a DONE ===\n")
cat("  Epithelial subset:     10a_seurat_epithelial.rds  ->  feeds 10b\n")
