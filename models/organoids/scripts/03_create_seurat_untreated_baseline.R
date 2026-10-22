# 03_create_seurat_untreated_baseline.R
# Build a new Seurat object from raw DGE_filtered data, including only
# untreated/baseline samples:
#   OPTO98_C1_D0, PDO66_C1_D0, OCAD93, OCAD96, OCAD97, OCAD106, OPTO129, OPTO112
#
# Pipeline: read raw matrix → filter to 8 samples → gene filtering (same as
# 01_create_seurat.R) → QC metrics → doublet detection (scDblFinder) →
# cell-level QC filtering → normalization → HVGs → PCA → UMAP → save
#
# Output: output/01_Data_Processing_and_QC/seurat_untreated_baseline.rds

library(Seurat)
library(Matrix)
library(scDblFinder)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)

# --- Paths ---
data_dir   <- "data/single_cell/all-sample/DGE_filtered"
output_dir <- "output/01_Data_Processing_and_QC"
output_file <- file.path(output_dir, "seurat_untreated_baseline.rds")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- Samples to include ---
samples_keep <- c(
  "OPTO98_C1_D0",
  "PDO66_C1_D0",
  "OCAD93",
  "OCAD96",
  "OCAD97",
  "OCAD106",
  "OPTO129",
  "OPTO112"
)

# ===========================================================================
# 1. Read raw data
# ===========================================================================
cat("=== Step 1: Reading raw data ===\n")

cat("Reading count matrix...\n")
counts <- readMM(file.path(data_dir, "count_matrix.mtx"))

cat("Reading gene information...\n")
genes <- read.csv(file.path(data_dir, "all_genes.csv"))

cat("Reading cell metadata...\n")
cell_meta <- read.csv(file.path(data_dir, "cell_metadata_updated.csv"))
cell_meta$sample <- trimws(cell_meta$sample)

cat("  Total cells in raw data:", nrow(cell_meta), "\n")
cat("  Total genes in raw data:", nrow(genes), "\n")

# ===========================================================================
# 2. Filter to target samples (before transposing, for efficiency)
# ===========================================================================
cat("\n=== Step 2: Filtering to target samples ===\n")

# Verify all requested samples exist
available <- unique(cell_meta$sample)
missing <- setdiff(samples_keep, available)
if (length(missing) > 0) {
  stop("Samples not found in metadata: ", paste(missing, collapse = ", "))
}

# Cell indices to keep (1-based row indices in the cells x genes matrix)
cell_idx <- which(cell_meta$sample %in% samples_keep)
cat("  Cells matching target samples:", length(cell_idx), "\n")

# Subset matrix rows (cells) before transpose
counts <- counts[cell_idx, ]
cell_meta <- cell_meta[cell_idx, ]
rownames(cell_meta) <- cell_meta$bc_wells

cat("  Cells per sample:\n")
print(table(cell_meta$sample))

# ===========================================================================
# 3. Transpose and set dimension names
# ===========================================================================
cat("\n=== Step 3: Transposing matrix (cells x genes -> genes x cells) ===\n")
counts <- t(counts)

gene_names <- make.unique(genes$gene_name)
rownames(counts) <- gene_names
colnames(counts) <- cell_meta$bc_wells

# ===========================================================================
# 4. Gene filtering (same logic as 01_create_seurat.R)
# ===========================================================================
cat("\n=== Step 4: Gene filtering ===\n")
cat("  Total genes before filtering:", nrow(counts), "\n")

has_symbol   <- !grepl("^ENSG[0-9]+$", gene_names)
is_lncrna    <- grepl("^LINC|^LOC[0-9]", gene_names)
is_mirna     <- grepl("^MIR[0-9]", gene_names)
is_snorna    <- grepl("^SNOR|^RNU[0-9]|^SCARNA", gene_names)
is_pseudogene <- grepl("-AS[0-9]|P[0-9]+$", gene_names) & !grepl("^ATP|^SLC|^MAP", gene_names)
is_ribosomal <- grepl("^RP[SL][0-9]", gene_names)
is_mito      <- grepl("^MT-", gene_names)

# Keep genes with symbols, exclude lncRNA/miRNA/snoRNA/ribosomal (retain mito)
keep_genes <- has_symbol & !is_lncrna & !is_mirna & !is_snorna & !is_ribosomal

cat("  Genes without symbols (ENSG only):", sum(!has_symbol), "\n")
cat("  lncRNAs (LINC/LOC):", sum(is_lncrna & has_symbol), "\n")
cat("  miRNAs:", sum(is_mirna & has_symbol), "\n")
cat("  snoRNAs/snRNAs:", sum(is_snorna & has_symbol), "\n")
cat("  Ribosomal proteins:", sum(is_ribosomal & has_symbol), "\n")
cat("  Mitochondrial (retained):", sum(is_mito & has_symbol), "\n")
cat("  Genes kept:", sum(keep_genes), "\n")

counts <- counts[keep_genes, ]
genes  <- genes[keep_genes, ]

# Convert to CsparseMatrix
counts <- as(counts, "CsparseMatrix")

# ===========================================================================
# 5. Create Seurat object
# ===========================================================================
cat("\n=== Step 5: Creating Seurat object ===\n")
seurat_obj <- CreateSeuratObject(
  counts   = counts,
  meta.data = cell_meta,
  project  = "cyclic_chemo_untreated_baseline"
)

# Add model and QC metadata
seurat_obj$pdo_model  <- gsub("_.*", "", seurat_obj$sample)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

cat("  Object created:", ncol(seurat_obj), "cells x", nrow(seurat_obj), "genes\n")

# ===========================================================================
# 6. Doublet detection (scDblFinder, per sample)
# ===========================================================================
cat("\n=== Step 6: Doublet detection (scDblFinder) ===\n")

sce <- as.SingleCellExperiment(seurat_obj)
set.seed(42)
sce <- scDblFinder(sce, samples = "sample")

seurat_obj$scDblFinder_class <- sce$scDblFinder.class
seurat_obj$scDblFinder_score <- sce$scDblFinder.score

n_doublet <- sum(seurat_obj$scDblFinder_class == "doublet")
pct_doublet <- round(100 * n_doublet / ncol(seurat_obj), 1)
cat("  Doublets detected:", n_doublet, "(", pct_doublet, "%)\n")

cat("  Per-sample doublet counts:\n")
dbl_stats <- seurat_obj@meta.data %>%
  group_by(sample) %>%
  summarize(
    total = n(),
    doublets = sum(scDblFinder_class == "doublet"),
    pct_doublet = round(100 * doublets / total, 1),
    .groups = "drop"
  )
print(as.data.frame(dbl_stats))

# ===========================================================================
# 7. Cell-level QC filtering
# ===========================================================================
cat("\n=== Step 7: Cell-level QC filtering ===\n")

NCOUNT_MIN   <- 500
NFEATURE_MIN <- 300
PCT_MT_MAX   <- 20

n_before <- ncol(seurat_obj)

fail_ncount  <- sum(seurat_obj$nCount_RNA   <= NCOUNT_MIN)
fail_nfeat   <- sum(seurat_obj$nFeature_RNA <= NFEATURE_MIN)
fail_mito    <- sum(seurat_obj$percent.mt   >= PCT_MT_MAX)
fail_doublet <- sum(seurat_obj$scDblFinder_class == "doublet")

cat(sprintf("  Cells failing nCount_RNA <= %d:    %d\n",  NCOUNT_MIN,   fail_ncount))
cat(sprintf("  Cells failing nFeature_RNA <= %d:  %d\n",  NFEATURE_MIN, fail_nfeat))
cat(sprintf("  Cells failing percent.mt >= %.0f%%: %d\n", PCT_MT_MAX,   fail_mito))
cat(sprintf("  Doublets (scDblFinder):              %d\n", fail_doublet))

keep <- seurat_obj$nCount_RNA   > NCOUNT_MIN  &
        seurat_obj$nFeature_RNA > NFEATURE_MIN &
        seurat_obj$percent.mt   < PCT_MT_MAX   &
        seurat_obj$scDblFinder_class == "singlet"

seurat_obj <- seurat_obj[, keep]
n_after  <- ncol(seurat_obj)
n_removed <- n_before - n_after

cat(sprintf("  Total cells removed: %d (%.1f%%)\n", n_removed, 100 * n_removed / n_before))
cat(sprintf("  Cells retained:      %d (%.1f%%)\n", n_after,   100 * n_after  / n_before))

cat("\n  Cells per sample after filtering:\n")
print(table(seurat_obj$sample))

cat("\n  Cells per PDO model after filtering:\n")
print(table(seurat_obj$pdo_model))

# ===========================================================================
# 8. Standard Seurat workflow: normalize, HVGs, PCA, UMAP
# ===========================================================================
cat("\n=== Step 8: Normalization, HVGs, PCA, UMAP ===\n")

seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst",
                                   nfeatures = 2000, verbose = FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)

cat("  Normalization, 2000 HVGs (VST), 30 PCA dims, UMAP complete\n")

# ===========================================================================
# 9. Save
# ===========================================================================
cat("\n=== Step 9: Saving ===\n")
cat("Saving to:", output_file, "\n")
saveRDS(seurat_obj, output_file)

cat("\n=== DONE ===\n")
cat("Final object:", nrow(seurat_obj), "genes x", ncol(seurat_obj), "cells\n")
cat("Samples:", length(unique(seurat_obj$sample)), "\n")
cat("Output:", output_file, "\n")
