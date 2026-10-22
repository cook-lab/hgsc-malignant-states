# ============================================================================
# 03_ucell_scoring_noBCAM.R — UCell SecA / SecB polarization (Xenium ∩ atlas genes)
# ============================================================================
# PURPOSE: Score every SingleR "Secretory epithelium" cell with UCell for the
#   SecA and SecB signatures, in the gene space shared by the Xenium panel and
#   the atlas (457 genes), so that Xenium and atlas scores are comparable:
#   polarization_UCell = SecB_UCell - SecA_UCell. The atlas secretory cells are
#   scored on the same genes; 04_reclassification_polarization.R derives the
#   SecA / Intermediate / SecB thresholds from those atlas scores.
#
# SAMPLES: every SFE in sfe_names_all.
#
# SIGNATURES: the 7 + 7 SecA / SecB genes of shared/signatures.yml via 00_setup.R.
#   UCell runs on raw counts with maxRank = the number of shared genes (457): a
#   full ranking, below UCell's default of 1500, so nothing is truncated.
#
# NAMING: atlas NMF labels "Transitioning epithelium" are read as
#   "Intermediate epithelium".
#
# INPUTS:
#   - <sfe_dir>/sfe_*  (SFEs with singler_label)
#   - atlas UCell export: atlas_gene_list.txt, atlas_secretory_counts.mtx.gz,
#       atlas_secretory_barcodes.tsv, atlas_secretory_genes.tsv,
#       atlas_secretory_metadata.csv. Folder: ATLAS_UCELL_DIR if set, otherwise the
#       deposited data_root/2026_final_atlas/output/18_ucell_atlas/ (written by
#       atlas/03_epithelial_nmf/as_run/18a_ucell_atlas_export.py + 18b_ucell_atlas_scoring.R),
#       which the final run used.
#
# OUTPUTS:
#   - SFEs updated with SecA_UCell, SecB_UCell, polarization_UCell (NA for cells
#     that are not SingleR "Secretory epithelium")
#   - <output_root>/06d_annotation_noBCAM/
#       xenium_ucell_scores.csv, atlas_ucell_scores.csv  (atlas scores read by
#         04_reclassification_polarization.R)
#       shared_gene_list_xenium_atlas.txt  (the 457 genes; read by
#         tables/06_tableS6_secA_secB_signatures.py and by the count export,
#         spatial/05_count_export/as_run/split_recompute_stage1_polarization.R)
#       checkpoint_{xenium,atlas}_scores.rds  (resume points: when they, or both
#         score CSVs, exist, UCell is not re-run)
#
# MANUSCRIPT PANEL(S): polarization_UCell, the score behind every Xenium SecA /
#   Intermediate / SecB call (Figs 4–7, Figures S11–S18, Tables S7–S11, through
#   04_reclassification_polarization.R and 06_extract_singler.R); the 457-gene
#   Xenium space of the Methods (Polarization scoring) and Table S6.
#
# RUNTIME TIER: heavy (UCell over all secretory cells + atlas secretory pool)
# ============================================================================

source("spatial/00_setup/00_setup.R")
library(UCell)
library(Seurat)
library(Matrix)

message("\n=== UCell Secretory Scoring (noBCAM) — Shared Gene Space ===")

# --- Gene signatures (noBCAM, shared) ----------------------------------------

secA_genes <- SECA_GENES
secB_genes <- SECB_GENES

message("SecA genes (", length(secA_genes), ", noBCAM): ", paste(secA_genes, collapse = ", "))
message("SecB genes (", length(secB_genes), "): ", paste(secB_genes, collapse = ", "))

# Standardize atlas NMF label naming.
remap_nmf_label <- function(x) {
  x <- as.character(x)
  x[x == "Transitioning epithelium"] <- "Intermediate epithelium"
  x
}

# --- Paths --------------------------------------------------------------------

out_path <- file.path(out_dir, "06d_annotation_noBCAM")
fig_path <- file.path(fig_dir, "06d_annotation_noBCAM")
for (d in c(out_path, fig_path)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Atlas UCell export folder: ATLAS_UCELL_DIR if set, otherwise the deposited
# export (data_root/2026_final_atlas/output/18_ucell_atlas/), which the final run read.
atlas_dir <- Sys.getenv("ATLAS_UCELL_DIR", unset = "")
if (!nzchar(atlas_dir)) {
  atlas_dir <- file.path(path.expand(CFG$paths$data_root),
                         "2026_final_atlas", "output", "18_ucell_atlas")
}
message("Atlas UCell export dir: ", atlas_dir)

sfe_names <- sfe_names_all   # every SFE: TMA A, the cohort.tma_new arrays, WT01–WT08

# ===========================================================================
# PART 1: Compute shared gene space & extract secretory cells
# ===========================================================================

message("\n=== PART 1: Shared gene space + secretory cell extraction ===")

message("\n--- [1a] Computing gene intersection ---")
atlas_genes <- readLines(file.path(atlas_dir, "atlas_gene_list.txt"))

message("\n--- [1b] Extracting secretory cell counts ---")

count_list   <- list()
barcode_meta <- list()

for (sname in sfe_names) {
  message("  Loading ", sname, " ...")
  sfe <- load_sfe(sname)

  is_sec <- sfe$singler_label == "Secretory epithelium"
  is_sec[is.na(is_sec)] <- FALSE
  n_sec <- sum(is_sec)

  if (n_sec > 0) {
    cts <- counts(sfe[, is_sec])
    cts <- as(cts, "dgCMatrix")

    barcodes <- paste0(sname, "_", colnames(cts))
    colnames(cts) <- barcodes

    count_list[[sname]] <- cts
    barcode_meta[[sname]] <- data.frame(
      barcode_unique = barcodes,
      barcode_orig   = colnames(sfe)[is_sec],
      sample         = sname,
      stringsAsFactors = FALSE
    )
  }

  message("    ", format(n_sec, big.mark = ","), " secretory cells")
  rm(sfe); gc(verbose = FALSE)
}

message("\n--- Merging secretory counts ---")
merged_counts <- do.call(cbind, count_list)
meta_df <- do.call(rbind, barcode_meta)
rownames(meta_df) <- meta_df$barcode_unique

panel_genes <- rownames(merged_counts)
message("  Merged: ", ncol(merged_counts), " secretory cells x ",
        nrow(merged_counts), " genes")

shared_genes <- sort(intersect(panel_genes, atlas_genes))
message("  Atlas genes: ", length(atlas_genes))
message("  Xenium genes: ", length(panel_genes))
message("  Shared genes (Xenium ∩ Atlas): ", length(shared_genes))

secA_present <- secA_genes[secA_genes %in% shared_genes]
secB_present <- secB_genes[secB_genes %in% shared_genes]
secA_missing <- secA_genes[!secA_genes %in% shared_genes]
secB_missing <- secB_genes[!secB_genes %in% shared_genes]

message("  SecA in shared: ", length(secA_present), "/", length(secA_genes))
if (length(secA_missing) > 0) message("    Missing: ", paste(secA_missing, collapse = ", "))
message("  SecB in shared: ", length(secB_present), "/", length(secB_genes))
if (length(secB_missing) > 0) message("    Missing: ", paste(secB_missing, collapse = ", "))

writeLines(shared_genes, file.path(out_path, "shared_gene_list_xenium_atlas.txt"))

message("\n--- [1c] Subsetting to shared gene space ---")
merged_counts <- merged_counts[shared_genes, ]
message("  Xenium matrix (shared space): ", ncol(merged_counts), " cells x ",
        nrow(merged_counts), " genes")

# ===========================================================================
# PART 2: UCell scoring on shared gene space
# ===========================================================================

message("\n=== PART 2: UCell scoring (shared gene space) ===")

# Resume from checkpoints/CSV exports if present (UCell is expensive).
checkpoint_xenium <- file.path(out_path, "checkpoint_xenium_scores.rds")
checkpoint_atlas  <- file.path(out_path, "checkpoint_atlas_scores.rds")
csv_xenium <- file.path(out_path, "xenium_ucell_scores.csv")
csv_atlas  <- file.path(out_path, "atlas_ucell_scores.csv")

if (file.exists(checkpoint_xenium) && file.exists(checkpoint_atlas)) {
  message("  Found RDS checkpoints — loading saved UCell scores (skipping UCell)")
  scores_df    <- readRDS(checkpoint_xenium)
  atlas_scores <- readRDS(checkpoint_atlas)
  rm(merged_counts, count_list); gc(verbose = FALSE)
  message("  Xenium: ", format(nrow(scores_df), big.mark = ","), " cells")
  message("  Atlas:  ", format(nrow(atlas_scores), big.mark = ","), " cells")
} else if (file.exists(csv_xenium) && file.exists(csv_atlas)) {
  message("  Found CSV exports — loading saved UCell scores (skipping UCell)")
  scores_df    <- read.csv(csv_xenium, stringsAsFactors = FALSE)
  atlas_scores <- read.csv(csv_atlas, stringsAsFactors = FALSE)
  rm(merged_counts, count_list); gc(verbose = FALSE)
  message("  Xenium: ", format(nrow(scores_df), big.mark = ","), " cells")
  message("  Atlas:  ", format(nrow(atlas_scores), big.mark = ","), " cells")
} else {

message("\n--- [2a] UCell on Xenium secretory cells ---")

seurat_sec <- CreateSeuratObject(counts = merged_counts, meta.data = meta_df)
rm(merged_counts, count_list); gc(verbose = FALSE)

signatures <- list(
  SecA = secA_present,
  SecB = secB_present
)

# NB: UCell scoring here is single-threaded on purpose. Neither
# MulticoreParam/ncores (forking fails on this R build's S4 method tables — same
# issue as SingleR — and silently falls back to serial) nor SnowParam (socket
# workers: no speedup) parallelizes usefully: the cost is UCell's per-cell
# ranking at maxRank=nrow, which is serial regardless of BPPARAM. To speed this
# stage up, split the cells into chunks across independent Rscript processes
# and combine the scores.
n_genes <- nrow(seurat_sec)
message("  Running UCell (maxRank = ", n_genes, " [all shared genes], counts slot) ...")
seurat_sec <- AddModuleScore_UCell(seurat_sec, features = signatures,
                                    maxRank = n_genes)

seurat_sec$polarization_UCell <- seurat_sec$SecB_UCell - seurat_sec$SecA_UCell

message("  SecA UCell range: ", paste(round(range(seurat_sec$SecA_UCell), 4), collapse = " — "))
message("  SecB UCell range: ", paste(round(range(seurat_sec$SecB_UCell), 4), collapse = " — "))
message("  Polarization range: ", paste(round(range(seurat_sec$polarization_UCell), 4), collapse = " — "))

scores_df <- data.frame(
  barcode_unique     = colnames(seurat_sec),
  barcode_orig       = seurat_sec$barcode_orig,
  sample             = seurat_sec$sample,
  SecA_UCell         = seurat_sec$SecA_UCell,
  SecB_UCell         = seurat_sec$SecB_UCell,
  polarization_UCell = seurat_sec$polarization_UCell,
  stringsAsFactors   = FALSE
)

rm(seurat_sec); gc(verbose = FALSE)

saveRDS(scores_df, file.path(out_path, "checkpoint_xenium_scores.rds"))
message("  Checkpoint saved: checkpoint_xenium_scores.rds")

message("\n--- [2b] UCell on atlas secretory cells (same shared space) ---")

atlas_mtx <- readMM(file.path(atlas_dir, "atlas_secretory_counts.mtx.gz"))
atlas_barcodes <- readLines(file.path(atlas_dir, "atlas_secretory_barcodes.tsv"))
atlas_gene_names <- readLines(file.path(atlas_dir, "atlas_secretory_genes.tsv"))

rownames(atlas_mtx) <- atlas_gene_names
colnames(atlas_mtx) <- atlas_barcodes
atlas_mtx <- as(atlas_mtx, "dgCMatrix")

message("  Atlas secretory matrix: ", ncol(atlas_mtx), " cells x ", nrow(atlas_mtx), " genes")

atlas_shared <- shared_genes[shared_genes %in% rownames(atlas_mtx)]
atlas_mtx <- atlas_mtx[atlas_shared, ]
message("  Atlas (shared space): ", ncol(atlas_mtx), " cells x ",
        nrow(atlas_mtx), " genes")

seurat_atlas <- CreateSeuratObject(counts = atlas_mtx)
rm(atlas_mtx); gc(verbose = FALSE)

signatures_atlas <- list(
  SecA = secA_present[secA_present %in% rownames(seurat_atlas)],
  SecB = secB_present[secB_present %in% rownames(seurat_atlas)]
)

n_genes_atlas <- nrow(seurat_atlas)
message("  Running UCell on atlas (maxRank = ", n_genes_atlas, " [all shared genes]) ...")
seurat_atlas <- AddModuleScore_UCell(seurat_atlas, features = signatures_atlas,
                                      maxRank = n_genes_atlas)

atlas_meta <- read.csv(file.path(atlas_dir, "atlas_secretory_metadata.csv"),
                       row.names = 1)

atlas_scores <- data.frame(
  barcode      = colnames(seurat_atlas),
  SecA_UCell   = seurat_atlas$SecA_UCell,
  SecB_UCell   = seurat_atlas$SecB_UCell,
  polarization = seurat_atlas$SecB_UCell - seurat_atlas$SecA_UCell,
  stringsAsFactors = FALSE
)

atlas_scores <- merge(atlas_scores, atlas_meta,
                       by.x = "barcode", by.y = "row.names", all.x = TRUE)

message("  Atlas SecA range: ",
        paste(round(range(atlas_scores$SecA_UCell), 4), collapse = " — "))
message("  Atlas SecB range: ",
        paste(round(range(atlas_scores$SecB_UCell), 4), collapse = " — "))

rm(seurat_atlas); gc(verbose = FALSE)

saveRDS(atlas_scores, file.path(out_path, "checkpoint_atlas_scores.rds"))
message("  Checkpoint saved: checkpoint_atlas_scores.rds")

} # end else (no checkpoints)

# Standardize atlas NMF label naming (Transitioning -> Intermediate).
if ("celltype_nmf" %in% colnames(atlas_scores)) {
  atlas_scores$celltype_nmf <- remap_nmf_label(atlas_scores$celltype_nmf)
}

# --- [2c] Export scores -------------------------------------------------------

message("\n--- [2c] Exporting scores ---")

write.csv(scores_df, file.path(out_path, "xenium_ucell_scores.csv"),
          row.names = FALSE)
write.csv(atlas_scores, file.path(out_path, "atlas_ucell_scores.csv"),
          row.names = FALSE)

# ===========================================================================
# PART 3: Write UCell scores back to SFEs
# ===========================================================================

message("\n=== PART 3: Writing UCell scores back to SFEs ===")

for (sname in sfe_names) {
  message("  Updating ", sname, " ...")
  sfe <- load_sfe(sname)

  is_sec <- sfe$singler_label == "Secretory epithelium"
  is_sec[is.na(is_sec)] <- FALSE

  sfe$SecA_UCell         <- NA_real_
  sfe$SecB_UCell         <- NA_real_
  sfe$polarization_UCell <- NA_real_

  if (sum(is_sec) > 0) {
    sample_scores <- scores_df[scores_df$sample == sname, ]
    m <- match(colnames(sfe)[is_sec], sample_scores$barcode_orig)

    sfe$SecA_UCell[is_sec]         <- sample_scores$SecA_UCell[m]
    sfe$SecB_UCell[is_sec]         <- sample_scores$SecB_UCell[m]
    sfe$polarization_UCell[is_sec] <- sample_scores$polarization_UCell[m]

    n_matched <- sum(!is.na(m))
    message("    Matched ", n_matched, " / ", sum(is_sec), " secretory cells")
  }

  for (a in assayNames(sfe)) {
    assay(sfe, a) <- as(assay(sfe, a), "dgCMatrix")
  }

  save_sfe(sfe, sname)
  rm(sfe); gc(verbose = FALSE)
}

# ===========================================================================
# Summary
# ===========================================================================

message("\n=== UCell scoring complete ===")
message("Shared genes (Xenium ∩ Atlas): ", length(shared_genes))
message("Xenium secretory cells scored: ", format(nrow(scores_df), big.mark = ","))
message("Atlas secretory cells scored: ", format(nrow(atlas_scores), big.mark = ","))
message("\nOutput:  ", out_path)
message("SFEs updated with: SecA_UCell, SecB_UCell, polarization_UCell")

log_session()
