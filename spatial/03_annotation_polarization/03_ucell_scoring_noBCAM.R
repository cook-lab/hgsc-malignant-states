# ============================================================================
# 03_ucell_scoring_noBCAM.R — UCell secretory scoring (noBCAM) + atlas context
# ============================================================================
# PURPOSE: Score all SingleR "Secretory epithelium" cells with UCell using the
#   shared noBCAM signatures, on the Xenium-atlas gene intersection so Xenium
#   and atlas scores are directly comparable. Computes the intersection (~457
#   genes), scores Xenium secretory cells, scores atlas secretory cells on the
#   same gene space, writes scores back into the SFEs, and renders comparison
#   figures.
#
# COHORT PIN: operates over sfe_tma + the 8 published whole tissues.
#
# SIGNATURES: SecA/SecB from shared/signatures.yml (noBCAM) via 00_setup.R.
#   UCell maxRank = nrow() (full ranking on the shared ~457-gene space; both
#   panels are well below UCell's default 1500, so no truncation occurs).
#
# NAMING: atlas NMF labels remapped on read "Transitioning epithelium" ->
#   "Intermediate epithelium" for consistency with the manuscript schema.
#
# INPUTS:
#   - <sfe_dir>/sfe_tma, <sfe_dir>/sfe_<wt>  (annotated SFEs)
#   - atlas UCell exports (from atlas/03_epithelial_nmf 03_ucell_atlas_export.py
#       + 04_ucell_atlas_scoring.R), read from
#       output_root/03_epithelial_nmf/ucell_atlas/ (preferred) or the deposited
#       data_root/2026_final_atlas/output/18_ucell_atlas/ (fallback);
#       override with ATLAS_UCELL_DIR:
#       atlas_gene_list.txt, atlas_secretory_counts.mtx.gz,
#       atlas_secretory_barcodes.tsv, atlas_secretory_genes.tsv,
#       atlas_secretory_metadata.csv
#
# OUTPUTS:
#   - SFEs updated with SecA_UCell, SecB_UCell, polarization_UCell
#   - <output_root>/06d_annotation_noBCAM/{xenium,atlas}_ucell_scores.csv (+summaries)
#   - <output_root>/figures/06d_annotation_noBCAM/*.pdf
#
# MANUSCRIPT PANEL(S): polarization-chain backend; produces the UCell scores
#   that 04_reclassification_polarization.R / 05_clean_split_rctd.R consume
#   (Fig 4–6).
#
# RUNTIME TIER: heavy (UCell over all secretory cells + atlas secretory pool)
# ============================================================================

source("spatial/00_setup/00_setup.R")
library(UCell)
library(Seurat)
library(Matrix)
library(dplyr)

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

# Atlas UCell export dir. The atlas backend (atlas/03_epithelial_nmf/
# 03_ucell_atlas_export.py + 04_ucell_atlas_scoring.R) writes these exports to
# output_root/03_epithelial_nmf/ucell_atlas/. Prefer that regenerated copy on a
# clean run; fall back to the deposited input under data_root
# (2026_final_atlas/output/18_ucell_atlas/) — same precedence as load_sfe().
# ATLAS_UCELL_DIR overrides both if set.
.atlas_ucell_candidates <- c(
  file.path(out_dir, "03_epithelial_nmf", "ucell_atlas"),
  file.path(path.expand(CFG$paths$data_root),
            "2026_final_atlas", "output", "18_ucell_atlas")
)
atlas_dir <- Sys.getenv("ATLAS_UCELL_DIR", unset = "")
if (!nzchar(atlas_dir)) {
  hit <- .atlas_ucell_candidates[file.exists(
    file.path(.atlas_ucell_candidates, "atlas_gene_list.txt"))]
  atlas_dir <- if (length(hit) > 0) hit[1] else .atlas_ucell_candidates[1]
}
message("Atlas UCell export dir: ", atlas_dir)

sfe_names <- sfe_names_all   # sfe_tma + published 8 whole tissues (cohort PIN)

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

sample_summary <- scores_df |>
  dplyr::group_by(sample) |>
  dplyr::summarize(
    n_cells = dplyr::n(),
    mean_SecA = round(mean(SecA_UCell), 4),
    median_SecA = round(median(SecA_UCell), 4),
    mean_SecB = round(mean(SecB_UCell), 4),
    median_SecB = round(median(SecB_UCell), 4),
    mean_polarization = round(mean(polarization_UCell), 4),
    pct_SecA_dominant = round(100 * mean(SecA_UCell > SecB_UCell), 1),
    pct_SecB_dominant = round(100 * mean(SecB_UCell > SecA_UCell), 1),
    .groups = "drop"
  )
write.csv(sample_summary, file.path(out_path, "per_sample_summary.csv"),
          row.names = FALSE)
message("  Per-sample summary:")
print(as.data.frame(sample_summary))

atlas_summary <- atlas_scores |>
  dplyr::group_by(celltype_nmf) |>
  dplyr::summarize(
    n_cells = dplyr::n(),
    mean_SecA = round(mean(SecA_UCell), 4),
    median_SecA = round(median(SecA_UCell), 4),
    mean_SecB = round(mean(SecB_UCell), 4),
    median_SecB = round(median(SecB_UCell), 4),
    mean_polarization = round(mean(polarization), 4),
    .groups = "drop"
  )
write.csv(atlas_summary, file.path(out_path, "atlas_summary.csv"),
          row.names = FALSE)
message("\n  Atlas summary (shared gene space):")
print(as.data.frame(atlas_summary))

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
# PART 4: Figures
# ===========================================================================

message("\n=== PART 4: Figures ===")

message("\n--- [4a] Score histograms ---")

p_hist_secA <- ggplot(scores_df, aes(x = SecA_UCell)) +
  geom_histogram(bins = 50, fill = ref_palette["SecA epithelium"],
                 color = "white", linewidth = 0.2) +
  facet_wrap(~sample, scales = "free_y") +
  labs(title = "SecA UCell score per Xenium sample [noBCAM, shared gene space]",
       x = "SecA UCell score", y = "Cell count") +
  theme_lab()
ggsave(file.path(fig_path, "hist_SecA_score_per_sample.pdf"), p_hist_secA,
       width = 12, height = 10)

p_hist_secB <- ggplot(scores_df, aes(x = SecB_UCell)) +
  geom_histogram(bins = 50, fill = ref_palette["SecB epithelium"],
                 color = "white", linewidth = 0.2) +
  facet_wrap(~sample, scales = "free_y") +
  labs(title = "SecB UCell score per Xenium sample [shared gene space]",
       x = "SecB UCell score", y = "Cell count") +
  theme_lab()
ggsave(file.path(fig_path, "hist_SecB_score_per_sample.pdf"), p_hist_secB,
       width = 12, height = 10)

message("--- [4b] Overlaid histograms ---")

score_long <- tidyr::pivot_longer(
  scores_df,
  cols = c("SecA_UCell", "SecB_UCell"),
  names_to = "signature", values_to = "score"
)
score_long$signature <- gsub("_UCell", "", score_long$signature)

p_hist_overlay <- ggplot(score_long, aes(x = score, fill = signature)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity",
                 color = "white", linewidth = 0.2) +
  facet_wrap(~sample, scales = "free_y") +
  scale_fill_manual(values = c("SecA" = ref_palette["SecA epithelium"],
                                "SecB" = ref_palette["SecB epithelium"])) +
  labs(title = "SecA vs SecB UCell distributions per sample [noBCAM, shared gene space]",
       x = "UCell score", y = "Cell count", fill = "Signature") +
  theme_lab()
ggsave(file.path(fig_path, "hist_SecA_vs_SecB_overlay_per_sample.pdf"), p_hist_overlay,
       width = 12, height = 10)

message("--- [4c] Scatter plots ---")

p_scatter <- ggplot(scores_df, aes(x = SecA_UCell, y = SecB_UCell, color = sample)) +
  geom_point(size = 0.1, alpha = 0.3) +
  scale_color_manual(values = setNames(
    colorRampPalette(okabe_ito)(length(unique(scores_df$sample))),
    sort(unique(scores_df$sample)))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "SecA vs SecB UCell scores [noBCAM, shared gene space]",
       subtitle = "Above diagonal = SecB-biased; below = SecA-biased",
       x = "SecA UCell score", y = "SecB UCell score", color = "Sample") +
  theme_lab()
ggsave(file.path(fig_path, "scatter_SecA_vs_SecB.pdf"), p_scatter,
       width = 8, height = 7)

p_scatter_facet <- ggplot(scores_df, aes(x = SecA_UCell, y = SecB_UCell)) +
  geom_point(size = 0.1, alpha = 0.3, color = "grey30") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~sample) +
  labs(title = "SecA vs SecB per sample [noBCAM, UCell, shared gene space]",
       x = "SecA UCell score", y = "SecB UCell score") +
  theme_lab()
ggsave(file.path(fig_path, "scatter_SecA_vs_SecB_per_sample.pdf"), p_scatter_facet,
       width = 12, height = 10)

message("--- [4d] Polarization violins ---")

p_violin <- ggplot(scores_df, aes(x = sample, y = polarization_UCell, fill = sample)) +
  geom_violin(scale = "width", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = setNames(
    colorRampPalette(okabe_ito)(length(unique(scores_df$sample))),
    sort(unique(scores_df$sample)))) +
  labs(title = "Secretory polarization (SecB − SecA) per sample [noBCAM, UCell]",
       x = NULL, y = "Polarization (SecB − SecA)") +
  theme_lab() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave(file.path(fig_path, "violin_polarization_per_sample.pdf"), p_violin,
       width = 10, height = 6)

message("--- [4e] Atlas context overlay ---")

xenium_dens <- data.frame(
  SecA = scores_df$SecA_UCell,
  SecB = scores_df$SecB_UCell,
  Polarization = scores_df$polarization_UCell,
  source = "Xenium"
)
atlas_dens <- data.frame(
  SecA = atlas_scores$SecA_UCell,
  SecB = atlas_scores$SecB_UCell,
  Polarization = atlas_scores$polarization,
  source = "Atlas"
)
combined_dens <- rbind(xenium_dens, atlas_dens)
combined_long <- tidyr::pivot_longer(combined_dens, cols = c("SecA", "SecB", "Polarization"),
                                      names_to = "score_type", values_to = "value")
combined_long$score_type <- factor(combined_long$score_type,
                                    levels = c("SecA", "SecB", "Polarization"))

p_atlas_density <- ggplot(combined_long, aes(x = value, fill = source)) +
  geom_density(alpha = 0.5, linewidth = 0.3) +
  facet_wrap(~score_type, scales = "free") +
  scale_fill_manual(values = c("Atlas" = "grey60", "Xenium" = "#E6A141")) +
  labs(title = "UCell score distributions: Xenium vs Atlas [shared gene space]",
       x = "UCell score", y = "Density", fill = "Source") +
  theme_lab()
ggsave(file.path(fig_path, "atlas_context_density_overlay.pdf"), p_atlas_density,
       width = 12, height = 5)

message("--- [4f] Atlas context scatter ---")

p_atlas_scatter <- ggplot() +
  geom_point(data = atlas_scores,
             aes(x = SecA_UCell, y = SecB_UCell),
             color = "grey80", size = 0.05, alpha = 0.2) +
  geom_point(data = scores_df,
             aes(x = SecA_UCell, y = SecB_UCell),
             color = "#E6A141", size = 0.05, alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "SecA vs SecB: Xenium (gold) over Atlas (grey) [shared gene space]",
       subtitle = "Both scored on identical gene universe",
       x = "SecA UCell score", y = "SecB UCell score") +
  theme_lab()
ggsave(file.path(fig_path, "atlas_context_scatter.pdf"), p_atlas_scatter,
       width = 8, height = 7)

p_atlas_nmf <- ggplot(atlas_scores,
                       aes(x = SecA_UCell, y = SecB_UCell, color = celltype_nmf)) +
  geom_point(size = 0.05, alpha = 0.3) +
  scale_color_manual(values = ref_palette) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "Atlas secretory cells: SecA vs SecB UCell [shared gene space]",
       x = "SecA UCell score", y = "SecB UCell score", color = "NMF label") +
  theme_lab()
ggsave(file.path(fig_path, "atlas_context_scatter_nmf.pdf"), p_atlas_nmf,
       width = 9, height = 7)

message("--- [4g] Spatial maps ---")
options(ggrastr.default.dpi = 300, ggrastr.default.dev = "ragg")

for (sname in sfe_names) {
  message("  Spatial plots for ", sname, " ...")
  sfe <- load_sfe(sname)

  is_sec <- sfe$singler_label == "Secretory epithelium"
  is_sec[is.na(is_sec)] <- FALSE

  if (sum(is_sec) == 0) {
    rm(sfe); gc(verbose = FALSE)
    next
  }

  sp_df <- data.frame(
    x = spatialCoords(sfe)[is_sec, 1],
    y = spatialCoords(sfe)[is_sec, 2],
    SecA_UCell = sfe$SecA_UCell[is_sec],
    SecB_UCell = sfe$SecB_UCell[is_sec],
    polarization = sfe$polarization_UCell[is_sec]
  )

  p_sp_secA <- ggplot(sp_df, aes(x = x, y = y, color = SecA_UCell)) +
    spatial_point_layer() +
    scale_color_gradientn(colors = expr_spatial) +
    coord_fixed() +
    labs(title = paste0(sname, ": SecA UCell [noBCAM]"), color = "SecA") +
    theme_lab() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), axis.line = element_blank())
  ggsave(file.path(fig_path, paste0("spatial_SecA_score_", sname, ".pdf")),
         p_sp_secA, width = 8, height = 7)

  p_sp_secB <- ggplot(sp_df, aes(x = x, y = y, color = SecB_UCell)) +
    spatial_point_layer() +
    scale_color_gradientn(colors = expr_spatial) +
    coord_fixed() +
    labs(title = paste0(sname, ": SecB UCell"), color = "SecB") +
    theme_lab() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), axis.line = element_blank())
  ggsave(file.path(fig_path, paste0("spatial_SecB_score_", sname, ".pdf")),
         p_sp_secB, width = 8, height = 7)

  p_sp_pol <- ggplot(sp_df, aes(x = x, y = y, color = polarization)) +
    spatial_point_layer() +
    scale_color_gradient2(low = ref_palette["SecA epithelium"],
                          mid = "grey90",
                          high = ref_palette["SecB epithelium"],
                          midpoint = 0) +
    coord_fixed() +
    labs(title = paste0(sname, ": Polarization (SecB − SecA)"),
         color = "Polarization") +
    theme_lab() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), axis.line = element_blank())
  ggsave(file.path(fig_path, paste0("spatial_polarization_", sname, ".pdf")),
         p_sp_pol, width = 8, height = 7)

  rm(sfe, sp_df); gc(verbose = FALSE)
}

# ===========================================================================
# Summary
# ===========================================================================

message("\n=== UCell scoring complete ===")
message("Shared genes (Xenium ∩ Atlas): ", length(shared_genes))
message("Xenium secretory cells scored: ", format(nrow(scores_df), big.mark = ","))
message("Atlas secretory cells scored: ", format(nrow(atlas_scores), big.mark = ","))
message("\nOutput:  ", out_path)
message("Figures: ", fig_path)
message("SFEs updated with: SecA_UCell, SecB_UCell, polarization_UCell")

log_session()
