# ============================================================================
# 01_colocalization.R
# ----------------------------------------------------------------------------
# PURPOSE: Pairwise spatial colocalization / avoidance of cell types per sample & TMA core.
#
# INPUTS:
#   - SFEs (load_sfe): sfe_tma_filtered + 8 whole-tissue, with cell_label (06g/06f)
#
# OUTPUTS:
#   - output/08_colocalization/  colocalization tables + figures
#
# MANUSCRIPT PANEL(S): Backend for Fig 4 composition/niche context (no single published panel).
# RUNTIME TIER: moderate
#
# Migrated from 2026_final_xenium_analysis/scripts/. Analytical logic preserved;
# paths routed through central config, seed from CFG$seed, epithelial label
# "Transitioning" -> "Intermediate", SecA/SecB from shared/signatures.yml.
# ============================================================================

# --- Config + shared setup (replaces hardcoded /Volumes/CookLab/Sarah paths) ---
here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
source(file.path(here, "..", "..", "config", "config.R"))   # CFG, cfg_obj, cfg_path
source(file.path(here, "..", "00_setup", "00_setup.R"))      # load_sfe, save_sfe, theme_lab, nb_names, palettes
set.seed(CFG$seed)

library(FNN)
library(dbscan)
library(spatstat.geom)
library(spatstat.explore)
library(spatstat.random)
library(circlize)
library(RColorBrewer)
library(ggrepel)

message("\n", strrep("=", 70))
message("08_colocalization.R — Spatial Co-localization Analysis")
message(strrep("=", 70))

# --- Output directories -----------------------------------------------------

colocal_dir <- file.path(out_dir, "08_colocalization")
if (!dir.exists(colocal_dir)) dir.create(colocal_dir, recursive = TRUE)

# --- Constants & palettes ----------------------------------------------------

# 18 cell types in canonical order (by lineage)
celltype_order <- c(
  "Ciliated epithelium",
  "SecA epithelium",

  "Intermediate epithelium",
  "SecB epithelium",
  "Mesothelial",
  "Fibroblast", "Smooth muscle", "Pericyte", "Endothelial",
  "T cell", "NK cell", "B cell", "Plasma cell",
  "Macrophage", "Conventional dendritic cell", "Plasmacytoid dendritic cell",
  "Neutrophil", "Mast cell"
)

# Lineage grouping
celltype_lineage <- c(
  rep("Epithelial", 5), rep("Stromal", 4), rep("Lymphoid", 4),
  rep("Myeloid", 5)
)
names(celltype_lineage) <- celltype_order

lineage_colors <- c(
  "Epithelial" = "#E6A141",
  "Stromal"    = "#D14E6C",
  "Lymphoid"   = "#87CEFA",
  "Myeloid"    = "#8FBC8F"
)

# Secretory subtypes
secretory_subtypes <- c("SecA epithelium", "Intermediate epithelium",
                         "SecB epithelium")

# 12 focal pairs for Cross-K analysis
focal_pairs <- data.frame(
  typeA = c("Macrophage", "Macrophage", "Fibroblast", "T cell",
            "B cell", "Endothelial", "NK cell", "Fibroblast",
            "SecB epithelium", "SecB epithelium",
            "SecA epithelium", "Intermediate epithelium"),
  typeB = c("T cell", "SecA epithelium", "SecA epithelium", "SecA epithelium",
            "Plasma cell", "Pericyte", "SecB epithelium", "Macrophage",
            "Macrophage", "T cell",
            "SecB epithelium", "SecB epithelium"),
  stringsAsFactors = FALSE
)
focal_pairs$pair_label <- paste0(focal_pairs$typeA, " \u2194 ", focal_pairs$typeB)

# Subtype palette for focused analysis
subtype_pal <- c(
  "SecA epithelium"          = "#E6A141",
  "Intermediate epithelium" = "#C08E48",
  "SecB epithelium"          = "#9A7D55"
)

# Whole tissue samples for Cross-K
wt_crossk_samples <- c("sfe_OTB_2461", "sfe_OTB_2454", "sfe_SP24_24824")

# --- Helper functions --------------------------------------------------------

#' Extract coordinates and cell_label from an SFE object
extract_coords_labels <- function(sfe, label_col = "cell_label") {
  xy <- spatialCoords(sfe)
  labs <- as.character(colData(sfe)[[label_col]])
  # Rename-mismatch fix: deposited SFEs still carry the legacy
  # "Transitioning epithelium"; standardize to the canonical "Intermediate
  # epithelium" right at the read point, before any downstream match / color /
  # filter on it (idempotent — a no-op if the legacy value is absent).
  labs[labs == "Transitioning epithelium"] <- "Intermediate epithelium"
  data.frame(x = xy[, 1], y = xy[, 2], label = labs, stringsAsFactors = FALSE)
}

#' Compute kNN enrichment for a single tissue unit
#'
#' @param coords N x 2 matrix of cell coordinates
#' @param labels Character vector of cell type labels (length N)
#' @param k Number of neighbors
#' @param min_cells Minimum cells of a type to include
#' @return Named list: enrichment, observed, expected, n_cells, type_counts
compute_knn_enrichment <- function(coords, labels, k = 20, min_cells = 20) {
  types_present <- sort(unique(labels))
  n <- length(labels)

  # Global type proportions (expected under CSR)
  type_freq <- table(labels)
  type_prop <- type_freq / n

  # kNN search
  nn <- FNN::get.knnx(data = coords, query = coords, k = k + 1)
  nn_idx <- nn$nn.index[, -1, drop = FALSE]  # remove self

  # Build neighbor type labels per cell

  nn_labels <- matrix(labels[nn_idx], nrow = n, ncol = k)

  # Observed: for each type A, fraction of A's neighbors that are type B
  obs_mat <- matrix(0, nrow = length(types_present), ncol = length(types_present),
                    dimnames = list(types_present, types_present))

  for (i in seq_along(types_present)) {
    type_a <- types_present[i]
    mask <- labels == type_a
    n_a <- sum(mask)
    if (n_a < min_cells) next
    nn_of_a <- nn_labels[mask, , drop = FALSE]
    for (j in seq_along(types_present)) {
      obs_mat[i, j] <- sum(nn_of_a == types_present[j]) / (n_a * k)
    }
  }

  # Expected: global proportions
  exp_mat <- matrix(rep(as.numeric(type_prop[types_present]),
                        each = length(types_present)),
                    nrow = length(types_present), ncol = length(types_present),
                    dimnames = list(types_present, types_present))

  # Log2 enrichment with pseudocount
  pseudo <- 1 / (n * k)
  enrich_mat <- log2((obs_mat + pseudo) / (exp_mat + pseudo))

  list(
    enrichment  = enrich_mat,
    observed    = obs_mat,
    expected    = exp_mat,
    n_cells     = n,
    type_counts = type_freq[types_present]
  )
}

#' Compute Cross-K for a pair of cell types in a tissue unit
#'
#' @param coords N x 2 matrix
#' @param labels Character vector of cell type labels
#' @param typeA, typeB The two focal types
#' @param rmax Maximum radius in microns
#' @param nsim Number of simulations for envelope
#' @return data.frame with r, L_obs, L_hi, L_lo, L_theo; or NULL
compute_crossk <- function(coords, labels, typeA, typeB,
                            rmax = 200, nsim = 199) {
  types <- factor(labels)
  if (!(typeA %in% levels(types)) || !(typeB %in% levels(types))) return(NULL)

  n_a <- sum(labels == typeA)
  n_b <- sum(labels == typeB)
  if (n_a < 10 || n_b < 10) return(NULL)

  win <- owin(range(coords[, 1]), range(coords[, 2]))
  pp <- ppp(coords[, 1], coords[, 2], window = win, marks = types)

  env <- tryCatch(
    envelope(pp, Kcross, i = typeA, j = typeB,
             r = seq(0, rmax, length.out = 100),
             nsim = nsim, correction = "border", verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(env)) return(NULL)

  env_df <- as.data.frame(env)
  env_df$L_obs  <- sqrt(env_df$obs / pi) - env_df$r
  env_df$L_hi   <- sqrt(env_df$hi / pi)  - env_df$r
  env_df$L_lo   <- sqrt(env_df$lo / pi)  - env_df$r
  env_df$L_theo <- 0
  env_df
}

#' Build consensus heatmap from a list of enrichment matrices
#' @param enrich_list Named list of enrichment matrices
#' @param types Character vector of cell types (row/col order)
#' @return Median consensus matrix aligned to types
build_consensus <- function(enrich_list, types) {
  arr <- array(NA, dim = c(length(types), length(types), length(enrich_list)),
               dimnames = list(types, types, names(enrich_list)))
  for (i in seq_along(enrich_list)) {
    mat <- enrich_list[[i]]
    shared <- intersect(types, rownames(mat))
    arr[shared, shared, i] <- mat[shared, shared]
  }
  consensus <- apply(arr, c(1, 2), median, na.rm = TRUE)
  consensus[is.nan(consensus)] <- 0
  consensus
}

#' Draw and save a ComplexHeatmap consensus heatmap to PDF
draw_consensus_heatmap <- function(mat, title, pdf_path, types, width = 10,
                                    height = 9) {
  col_fun <- colorRamp2(
    breaks = seq(-3, 3, length.out = 11),
    colors = rev(brewer.pal(11, "RdBu"))
  )
  row_ha <- HeatmapAnnotation(
    Lineage = celltype_lineage[types],
    col = list(Lineage = lineage_colors),
    which = "row",
    show_legend = TRUE,
    show_annotation_name = FALSE
  )
  ht <- Heatmap(
    mat,
    name = "log2 enrichment",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 45,
    left_annotation = row_ha,
    cell_fun = function(j, i, x, y, width, height, fill) {
      val <- mat[i, j]
      if (!is.na(val) && abs(val) > 0.5) {
        grid.text(sprintf("%.1f", val), x, y,
                  gp = gpar(fontsize = 6,
                            col = ifelse(abs(val) > 2, "white", "black")))
      }
    },
    column_title = title,
    column_title_gp = gpar(fontsize = 11),
    heatmap_legend_param = list(direction = "horizontal")
  )
  pdf(pdf_path, width = width, height = height)
  draw(ht, heatmap_legend_side = "bottom")
  dev.off()
  message("  Saved: ", pdf_path)
}


# =============================================================================
# PHASE 1: Extract spatial data
# =============================================================================

message("\n--- Phase 1: Extracting spatial data ---")

# --- 1a. TMA ----------------------------------------------------------------

message("Loading sfe_tma_filtered ...")
sfe_tma <- load_sfe("sfe_tma_filtered")

tma_dat <- extract_coords_labels(sfe_tma)
tma_dat$core_id     <- sfe_tma$core_id
tma_dat$sample_type <- sfe_tma$sample_type
tma_dat$patient_id  <- sfe_tma$patient_id

# Remove off-core cells
tma_dat <- tma_dat[!is.na(tma_dat$core_id) & tma_dat$core_id != "Off core", ]

# Core-level metadata
core_meta <- unique(tma_dat[, c("core_id", "sample_type", "patient_id")])
n_cores <- length(unique(tma_dat$core_id))
message("TMA: ", format(nrow(tma_dat), big.mark = ","), " cells in ",
        n_cores, " cores (",
        sum(core_meta$sample_type == "tumour"), " tumour, ",
        sum(core_meta$sample_type == "fallopian"), " fallopian)")

rm(sfe_tma); gc(verbose = FALSE)

# --- 1b. Whole tissue -------------------------------------------------------

wt_names <- sfe_names_wt

wt_dat <- list()
for (sname in wt_names) {
  message("Extracting ", sname, " ...")
  sfe <- load_sfe(sname)
  df <- extract_coords_labels(sfe)
  df$sample_id <- sname
  wt_dat[[sname]] <- df
  message("  ", format(nrow(df), big.mark = ","), " cells, ",
          length(unique(df$label)), " types")
  rm(sfe); gc(verbose = FALSE)
}

wt_cells <- sum(sapply(wt_dat, nrow))
message("Whole tissue total: ", format(wt_cells, big.mark = ","), " cells")


# =============================================================================
# PHASE 2: kNN Enrichment — TMA Cores
# =============================================================================

message("\n--- Phase 2: kNN enrichment — TMA cores (k=20) ---")

core_ids <- sort(unique(tma_dat$core_id))

tma_enrichments <- list()
tma_ncells      <- integer(0)
skipped_cores   <- character(0)

for (cid in core_ids) {
  sub <- tma_dat[tma_dat$core_id == cid, ]
  if (nrow(sub) < 100) {
    skipped_cores <- c(skipped_cores, cid)
    next
  }
  res <- compute_knn_enrichment(
    coords = as.matrix(sub[, c("x", "y")]),
    labels = sub$label,
    k = 20
  )
  tma_enrichments[[cid]] <- res$enrichment
  tma_ncells[cid] <- res$n_cells
}

n_computed <- length(tma_enrichments)
message("Computed enrichment for ", n_computed, " cores; skipped ",
        length(skipped_cores), " (< 100 cells)")

# Standard cell type set present in data
all_types <- celltype_order[celltype_order %in% unique(tma_dat$label)]

# --- 2a. TMA consensus heatmap ----------------------------------------------

consensus_tma <- build_consensus(tma_enrichments, all_types)

write.csv(consensus_tma,
          file.path(colocal_dir, "tma_consensus_enrichment_k20.csv"))

draw_consensus_heatmap(
  mat      = consensus_tma,
  title    = paste0("TMA Consensus kNN Enrichment (k=20, median of ",
                     n_computed, " cores)"),
  pdf_path = file.path(colocal_dir, "tma_consensus_heatmap_k20.pdf"),
  types    = all_types
)

# --- 2b. Tumour vs FT separate consensus ------------------------------------

tumour_cores <- core_meta$core_id[core_meta$sample_type == "tumour"]
ft_cores     <- core_meta$core_id[core_meta$sample_type == "fallopian"]

tumour_enrich <- tma_enrichments[intersect(names(tma_enrichments), tumour_cores)]
ft_enrich     <- tma_enrichments[intersect(names(tma_enrichments), ft_cores)]

consensus_tumour <- build_consensus(tumour_enrich, all_types)
consensus_ft     <- build_consensus(ft_enrich, all_types)

write.csv(consensus_tumour,
          file.path(colocal_dir, "tma_tumour_consensus_enrichment_k20.csv"))
write.csv(consensus_ft,
          file.path(colocal_dir, "tma_ft_consensus_enrichment_k20.csv"))

# Side-by-side heatmaps: Tumour, FT, Difference
col_fun <- colorRamp2(
  breaks = seq(-3, 3, length.out = 11),
  colors = rev(brewer.pal(11, "RdBu"))
)
diff_col <- colorRamp2(
  breaks = seq(-2, 2, length.out = 11),
  colors = rev(brewer.pal(11, "PiYG"))
)

row_ha <- HeatmapAnnotation(
  Lineage = celltype_lineage[all_types],
  col = list(Lineage = lineage_colors),
  which = "row",
  show_legend = TRUE,
  show_annotation_name = FALSE
)

diff_mat <- consensus_tumour - consensus_ft

ht_tumour <- Heatmap(
  consensus_tumour, name = "Tumour", col = col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
  column_names_rot = 45, left_annotation = row_ha,
  column_title = paste0("Tumour cores (n=", length(tumour_enrich), ")"),
  column_title_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(direction = "horizontal")
)
ht_ft <- Heatmap(
  consensus_ft, name = "FT", col = col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
  column_names_rot = 45,
  column_title = paste0("Fallopian tube cores (n=", length(ft_enrich), ")"),
  column_title_gp = gpar(fontsize = 10),
  show_heatmap_legend = FALSE
)
ht_diff <- Heatmap(
  diff_mat, name = "Difference", col = diff_col,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
  column_names_rot = 45,
  column_title = "Tumour - FT",
  column_title_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(direction = "horizontal")
)

pdf(file.path(colocal_dir, "tma_tumour_vs_ft_heatmaps.pdf"),
    width = 18, height = 9)
draw(ht_tumour + ht_ft + ht_diff,
     heatmap_legend_side = "bottom",
     column_title = "kNN Enrichment: Tumour vs. Fallopian Tube",
     column_title_gp = gpar(fontsize = 12))
dev.off()
message("  Saved: tma_tumour_vs_ft_heatmaps.pdf")


# =============================================================================
# PHASE 3: kNN Enrichment — Whole Tissue
# =============================================================================

message("\n--- Phase 3: kNN enrichment — whole tissue (k=20) ---")

wt_enrichments <- list()

for (sname in names(wt_dat)) {
  message("  ", sname, " ...")
  df <- wt_dat[[sname]]
  res <- compute_knn_enrichment(
    coords = as.matrix(df[, c("x", "y")]),
    labels = df$label,
    k = 20
  )
  wt_enrichments[[sname]] <- res
  message("    ", format(res$n_cells, big.mark = ","), " cells, ",
          length(res$type_counts), " types present")
}

# --- 3a. WT consensus -------------------------------------------------------

wt_enrich_mats <- lapply(wt_enrichments, `[[`, "enrichment")
consensus_wt <- build_consensus(wt_enrich_mats, all_types)

write.csv(consensus_wt,
          file.path(colocal_dir, "wt_consensus_enrichment_k20.csv"))

draw_consensus_heatmap(
  mat      = consensus_wt,
  title    = paste0("Whole Tissue Consensus kNN Enrichment (k=20, median of ",
                     length(wt_enrichments), " samples)"),
  pdf_path = file.path(colocal_dir, "wt_consensus_heatmap_k20.pdf"),
  types    = all_types
)

# --- 3b. Individual WT heatmaps (multi-page PDF) ----------------------------

pdf(file.path(colocal_dir, "wt_individual_heatmaps_k20.pdf"),
    width = 10, height = 9)
for (sname in names(wt_enrichments)) {
  mat <- wt_enrichments[[sname]]$enrichment
  shared <- intersect(all_types, rownames(mat))
  plot_mat <- matrix(0, nrow = length(all_types), ncol = length(all_types),
                     dimnames = list(all_types, all_types))
  plot_mat[shared, shared] <- mat[shared, shared]

  ht <- Heatmap(
    plot_mat, name = "log2 enrichment", col = col_fun,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
    column_names_rot = 45,
    column_title = paste0(gsub("sfe_", "", sname), " (n=",
                           format(wt_enrichments[[sname]]$n_cells,
                                  big.mark = ","), " cells)"),
    column_title_gp = gpar(fontsize = 11),
    heatmap_legend_param = list(direction = "horizontal")
  )
  draw(ht, heatmap_legend_side = "bottom")
}
dev.off()
message("  Saved: wt_individual_heatmaps_k20.pdf")


# =============================================================================
# PHASE 4: Tumour vs FT Wilcoxon Comparison
# =============================================================================

message("\n--- Phase 4: Tumour vs FT Wilcoxon rank-sum tests ---")

test_results <- list()
pair_idx <- 0

for (i in seq_along(all_types)) {
  for (j in seq_along(all_types)) {
    typeA <- all_types[i]
    typeB <- all_types[j]

    tumour_vals <- sapply(
      tma_enrichments[intersect(names(tma_enrichments), tumour_cores)],
      function(mat) {
        if (typeA %in% rownames(mat) && typeB %in% colnames(mat))
          mat[typeA, typeB] else NA
      })
    ft_vals <- sapply(
      tma_enrichments[intersect(names(tma_enrichments), ft_cores)],
      function(mat) {
        if (typeA %in% rownames(mat) && typeB %in% colnames(mat))
          mat[typeA, typeB] else NA
      })

    tumour_vals <- tumour_vals[!is.na(tumour_vals)]
    ft_vals     <- ft_vals[!is.na(ft_vals)]

    if (length(tumour_vals) < 5 || length(ft_vals) < 5) next

    wt <- wilcox.test(tumour_vals, ft_vals, conf.int = TRUE)
    pair_idx <- pair_idx + 1
    test_results[[pair_idx]] <- data.frame(
      typeA         = typeA,
      typeB         = typeB,
      pair          = paste0(typeA, " -> ", typeB),
      median_tumour = median(tumour_vals),
      median_ft     = median(ft_vals),
      diff          = median(tumour_vals) - median(ft_vals),
      hl_estimate   = wt$estimate,
      p_value       = wt$p.value,
      n_tumour      = length(tumour_vals),
      n_ft          = length(ft_vals),
      stringsAsFactors = FALSE
    )
  }
}

test_df <- do.call(rbind, test_results)
test_df$p_adj <- p.adjust(test_df$p_value, method = "BH")
test_df$sig   <- test_df$p_adj < 0.05
test_df <- test_df[order(test_df$p_adj), ]

message("Tested ", nrow(test_df), " cell type pairs; ",
        sum(test_df$sig), " significant at FDR < 0.05")

write.csv(test_df,
          file.path(colocal_dir, "tma_tumour_vs_ft_wilcoxon.csv"),
          row.names = FALSE)

# --- 4a. Volcano plot -------------------------------------------------------

test_df$neg_log10_p <- -log10(test_df$p_adj)
test_df$neg_log10_p[is.infinite(test_df$neg_log10_p)] <-
  max(test_df$neg_log10_p[is.finite(test_df$neg_log10_p)]) + 1

top_pairs <- head(test_df, 20)

p_volcano <- ggplot(test_df, aes(x = diff, y = neg_log10_p)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(aes(color = sig), size = 1.5, alpha = 0.6) +
  scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "grey70"),
                     name = "FDR < 0.05") +
  geom_text_repel(
    data = top_pairs,
    aes(label = pair),
    size = 2.5, max.overlaps = 20, segment.size = 0.3
  ) +
  labs(x = "Enrichment difference (Tumour - FT)",
       y = "-log10(adjusted p-value)",
       title = "Tumour vs. Fallopian Tube \u2014 kNN Enrichment Differences",
       subtitle = paste0(sum(test_df$sig), " / ", nrow(test_df),
                         " pairs significant (BH-FDR < 0.05)")) +
  theme_lab()

pdf(file.path(colocal_dir, "tumour_vs_ft_volcano.pdf"), width = 10, height = 8)
print(p_volcano)
dev.off()
message("  Saved: tumour_vs_ft_volcano.pdf")

# --- 4b. Focal pair violin plots --------------------------------------------

pair_values <- list()
for (p in seq_len(nrow(focal_pairs))) {
  a <- focal_pairs$typeA[p]
  b <- focal_pairs$typeB[p]
  vals <- sapply(tma_enrichments, function(mat) {
    if (a %in% rownames(mat) && b %in% colnames(mat)) mat[a, b] else NA
  })
  st <- core_meta$sample_type[match(names(vals), core_meta$core_id)]
  pair_values[[p]] <- data.frame(
    pair        = focal_pairs$pair_label[p],
    core_id     = names(vals),
    enrichment  = as.numeric(vals),
    sample_type = st,
    stringsAsFactors = FALSE
  )
}
pair_df <- do.call(rbind, pair_values)
pair_df <- pair_df[!is.na(pair_df$enrichment), ]

p_violin <- ggplot(pair_df, aes(x = pair, y = enrichment, fill = sample_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_violin(alpha = 0.6, scale = "width", position = position_dodge(0.7)) +
  geom_boxplot(width = 0.15, position = position_dodge(0.7),
               outlier.size = 0.5, alpha = 0.8) +
  scale_fill_manual(values = c("tumour" = "#D55E00", "fallopian" = "#56B4E9"),
                    name = "Tissue type") +
  labs(x = NULL, y = "log2 kNN enrichment (k=20)",
       title = "Per-core kNN enrichment for focal cell type pairs",
       subtitle = "TMA cores, tumour vs. fallopian tube") +
  theme_lab() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

pdf(file.path(colocal_dir, "focal_pairs_violin.pdf"), width = 14, height = 8)
print(p_violin)
dev.off()
message("  Saved: focal_pairs_violin.pdf")

write.csv(pair_df,
          file.path(colocal_dir, "focal_pairs_per_core_enrichment.csv"),
          row.names = FALSE)


# =============================================================================
# PHASE 5: Cross-K Function Analysis
# =============================================================================

message("\n--- Phase 5: Cross-K function analysis ---")

# --- 5a. Select representative TMA cores ------------------------------------

core_sizes <- tma_ncells[names(tma_enrichments)]
core_st <- setNames(core_meta$sample_type, core_meta$core_id)

tumour_cores_avail <- names(core_sizes)[core_st[names(core_sizes)] == "tumour"]
ft_cores_avail     <- names(core_sizes)[core_st[names(core_sizes)] == "fallopian"]

# Top 3 largest tumour cores
tumour_large <- head(sort(core_sizes[tumour_cores_avail], decreasing = TRUE), 3)
# 3 random diverse tumour cores
set.seed(CFG$seed)
tumour_diverse <- sample(setdiff(tumour_cores_avail, names(tumour_large)), 3)
# Top 2 largest FT cores
ft_select <- head(sort(core_sizes[ft_cores_avail], decreasing = TRUE), 2)

crossk_cores <- c(names(tumour_large), tumour_diverse, names(ft_select))
message("Selected ", length(crossk_cores), " cores for Cross-K: ",
        paste(crossk_cores, collapse = ", "))

# --- 5b. TMA Cross-K --------------------------------------------------------

message("Computing Cross-K for TMA cores (nsim=199, rmax=200) ...")
rmax <- 200

crossk_results <- list()
for (cid in crossk_cores) {
  sub <- tma_dat[tma_dat$core_id == cid, ]
  for (p in seq_len(nrow(focal_pairs))) {
    key <- paste0(cid, "|", focal_pairs$pair_label[p])
    message("  Cross-K: core ", cid, ", ", focal_pairs$pair_label[p])
    res <- compute_crossk(
      coords = as.matrix(sub[, c("x", "y")]),
      labels = sub$label,
      typeA  = focal_pairs$typeA[p],
      typeB  = focal_pairs$typeB[p],
      rmax   = rmax,
      nsim   = 199
    )
    if (!is.null(res)) {
      res$core_id     <- cid
      res$pair        <- focal_pairs$pair_label[p]
      res$sample_type <- core_st[cid]
      crossk_results[[key]] <- res
    }
  }
}

crossk_df <- do.call(rbind, crossk_results)
message("Cross-K computed for ", length(crossk_results), " core-pair combinations")

write.csv(crossk_df,
          file.path(colocal_dir, "crossk_tma_results.csv"),
          row.names = FALSE)

# Cross-K L(r)-r plots — TMA
if (!is.null(crossk_df) && nrow(crossk_df) > 0) {
  p_crossk_tma <- ggplot(crossk_df, aes(x = r)) +
    geom_ribbon(aes(ymin = L_lo, ymax = L_hi), fill = "grey85", alpha = 0.5) +
    geom_line(aes(y = L_obs, color = sample_type), linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    facet_grid(pair ~ core_id, scales = "free_y") +
    scale_color_manual(values = c("tumour" = "#D55E00", "fallopian" = "#56B4E9")) +
    labs(x = expression(r ~ "(µm)"), y = "L(r) - r",
         title = "Cross-K L(r)-r curves for focal cell type pairs",
         subtitle = "Grey envelope = 95% CI under CSR (199 simulations)") +
    theme_lab(base_size = 7) +
    theme(strip.text.y = element_text(size = 6, angle = 0))

  pdf(file.path(colocal_dir, "crossk_tma_curves.pdf"), width = 16, height = 22)
  print(p_crossk_tma)
  dev.off()
  message("  Saved: crossk_tma_curves.pdf")
}

# --- 5c. Whole tissue Cross-K (inhomogeneous, subsampled) --------------------

message("Computing inhomogeneous Cross-K for whole tissue (nsim=19) ...")

wt_crossk <- list()

for (sname in wt_crossk_samples) {
  df <- wt_dat[[sname]]
  message("  ", sname, " ...")

  for (p in seq_len(nrow(focal_pairs))) {
    typeA <- focal_pairs$typeA[p]
    typeB <- focal_pairs$typeB[p]

    n_a <- sum(df$label == typeA)
    n_b <- sum(df$label == typeB)
    if (n_a < 50 || n_b < 50) {
      message("    Skipping ", focal_pairs$pair_label[p],
              " (n_A=", n_a, ", n_B=", n_b, ")")
      next
    }

    # Subsample: 250 per focal type + 1000 background
    set.seed(CFG$seed)
    idx_a <- which(df$label == typeA)
    idx_b <- which(df$label == typeB)
    idx_other <- which(!df$label %in% c(typeA, typeB))
    keep <- c(
      if (length(idx_a) > 250) sample(idx_a, 250) else idx_a,
      if (length(idx_b) > 250) sample(idx_b, 250) else idx_b,
      sample(idx_other, min(1000, length(idx_other)))
    )
    sub <- df[keep, ]

    types <- factor(sub$label)
    win <- owin(range(sub$x), range(sub$y))
    pp <- ppp(sub$x, sub$y, window = win, marks = types)

    env <- tryCatch(
      envelope(pp, Kcross.inhom, i = typeA, j = typeB,
               r = seq(0, rmax, length.out = 100),
               nsim = 19, correction = "border", verbose = FALSE),
      error = function(e) {
        message("    Error: ", e$message)
        NULL
      }
    )

    if (!is.null(env)) {
      env_df <- as.data.frame(env)
      env_df$L_obs  <- sqrt(env_df$obs / pi) - env_df$r
      env_df$L_hi   <- sqrt(env_df$hi / pi)  - env_df$r
      env_df$L_lo   <- sqrt(env_df$lo / pi)  - env_df$r
      env_df$sample <- sname
      env_df$pair   <- focal_pairs$pair_label[p]
      wt_crossk[[paste0(sname, "|", focal_pairs$pair_label[p])]] <- env_df
    }
  }
}

wt_crossk_df <- do.call(rbind, wt_crossk)
message("Whole tissue Cross-K computed for ", length(wt_crossk),
        " sample-pair combinations")

write.csv(wt_crossk_df,
          file.path(colocal_dir, "crossk_wt_results.csv"),
          row.names = FALSE)

# Cross-K L(r)-r plots — WT
if (!is.null(wt_crossk_df) && nrow(wt_crossk_df) > 0) {
  p_crossk_wt <- ggplot(wt_crossk_df, aes(x = r)) +
    geom_ribbon(aes(ymin = L_lo, ymax = L_hi), fill = "grey85", alpha = 0.5) +
    geom_line(aes(y = L_obs), color = "#D55E00", linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    facet_grid(pair ~ sample, scales = "free_y") +
    labs(x = expression(r ~ "(µm)"), y = "L(r) - r (inhomogeneous)",
         title = "Inhomogeneous Cross-K for focal cell type pairs \u2014 whole tissue",
         subtitle = "Grey envelope = 95% CI under CSR (19 simulations)") +
    theme_lab(base_size = 7) +
    theme(strip.text.y = element_text(size = 6, angle = 0))

  pdf(file.path(colocal_dir, "crossk_wt_curves.pdf"), width = 14, height = 20)
  print(p_crossk_wt)
  dev.off()
  message("  Saved: crossk_wt_curves.pdf")
}


# =============================================================================
# PHASE 6: Cross-Context Consistency (TMA vs Whole Tissue)
# =============================================================================

message("\n--- Phase 6: Cross-context consistency ---")

shared_types <- intersect(rownames(consensus_tma), rownames(consensus_wt))
tma_vec <- as.vector(consensus_tma[shared_types, shared_types])
wt_vec  <- as.vector(consensus_wt[shared_types, shared_types])

pair_labels <- expand.grid(row = shared_types, col = shared_types,
                           stringsAsFactors = FALSE)
pair_labels$label <- paste0(pair_labels$row, " -> ", pair_labels$col)

context_df <- data.frame(
  tma_enrichment = tma_vec,
  wt_enrichment  = wt_vec,
  pair_label     = pair_labels$label,
  is_diagonal    = pair_labels$row == pair_labels$col,
  stringsAsFactors = FALSE
)
context_df <- context_df[!is.nan(context_df$tma_enrichment) &
                           !is.nan(context_df$wt_enrichment), ]

cor_val <- cor(context_df$tma_enrichment, context_df$wt_enrichment,
               method = "spearman", use = "complete.obs")
message("TMA vs WT Spearman rho = ", round(cor_val, 3))

# Concordance classification
context_df$tma_sign    <- sign(context_df$tma_enrichment)
context_df$wt_sign     <- sign(context_df$wt_enrichment)
context_df$concordance <- ifelse(context_df$tma_sign == context_df$wt_sign,
                                  "Concordant", "Discordant")
context_df$magnitude   <- abs(context_df$tma_enrichment) +
                            abs(context_df$wt_enrichment)

write.csv(context_df,
          file.path(colocal_dir, "cross_context_consistency.csv"),
          row.names = FALSE)

# Scatter plot
p_context <- ggplot(context_df, aes(x = tma_enrichment, y = wt_enrichment)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_vline(xintercept = 0, color = "grey80") +
  geom_point(aes(color = is_diagonal), size = 1.5, alpha = 0.6) +
  scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2"),
                     labels = c("TRUE" = "Self-enrichment", "FALSE" = "Cross-type"),
                     name = "Pair type") +
  annotate("text",
           x = min(context_df$tma_enrichment, na.rm = TRUE) + 0.5,
           y = max(context_df$wt_enrichment, na.rm = TRUE) - 0.5,
           label = paste0("Spearman rho = ", round(cor_val, 3)),
           size = 4, hjust = 0) +
  labs(x = "TMA consensus enrichment",
       y = "Whole tissue consensus enrichment",
       title = "Cross-context consistency: TMA vs. whole tissue",
       subtitle = "Each point = one cell type pair") +
  theme_lab()

pdf(file.path(colocal_dir, "cross_context_scatter.pdf"), width = 8, height = 8)
print(p_context)
dev.off()
message("  Saved: cross_context_scatter.pdf")


# =============================================================================
# PHASE 7: Secretory Subtype Focused Neighbor Analysis
# =============================================================================

message("\n--- Phase 7: Secretory subtype neighbor analysis ---")

# --- 7a. TMA consensus neighbor composition comparison -----------------------

neighbor_rows <- list()
for (st in secretory_subtypes) {
  if (st %in% rownames(consensus_tma)) {
    vals <- consensus_tma[st, ]
    neighbor_rows[[st]] <- data.frame(
      cell_type  = names(vals),
      enrichment = as.numeric(vals),
      subtype    = st,
      stringsAsFactors = FALSE
    )
  }
}
neighbor_comp <- do.call(rbind, neighbor_rows)
# Remove self-comparisons within secretory subtypes
neighbor_comp <- neighbor_comp[!neighbor_comp$cell_type %in% secretory_subtypes, ]

p_neighbor <- ggplot(neighbor_comp,
                      aes(x = reorder(cell_type, enrichment),
                          y = enrichment, fill = subtype)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = subtype_pal) +
  coord_flip() +
  labs(x = NULL, y = "log2 kNN enrichment (k=20)",
       title = "Neighbor profiles by secretory subtype",
       subtitle = "TMA consensus (median across cores)",
       fill = "Secretory subtype") +
  theme_lab()

pdf(file.path(colocal_dir, "secretory_neighbor_profiles.pdf"),
    width = 10, height = 7)
print(p_neighbor)
dev.off()
message("  Saved: secretory_neighbor_profiles.pdf")

write.csv(neighbor_comp,
          file.path(colocal_dir, "secretory_neighbor_enrichment.csv"),
          row.names = FALSE)

# --- 7b. Secretory self- and cross-enrichment heatmap -----------------------

sec_present <- intersect(secretory_subtypes, rownames(consensus_tma))
if (length(sec_present) >= 2) {
  sec_mat <- consensus_tma[sec_present, sec_present]

  sec_col <- colorRamp2(
    breaks = seq(-3, 3, length.out = 11),
    colors = rev(brewer.pal(11, "RdBu"))
  )

  ht_sec <- Heatmap(
    sec_mat,
    name = "log2 enrichment",
    col = sec_col,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 9),
    column_names_rot = 30,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", sec_mat[i, j]), x, y,
                gp = gpar(fontsize = 8))
    },
    column_title = "Secretory subtype self- and cross-enrichment (TMA consensus)",
    column_title_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(direction = "horizontal")
  )

  pdf(file.path(colocal_dir, "secretory_cross_enrichment.pdf"),
      width = 6, height = 5)
  draw(ht_sec, heatmap_legend_side = "bottom")
  dev.off()
  message("  Saved: secretory_cross_enrichment.pdf")
}

# --- 7c. Per-sample SecB vs SecA differential neighbor enrichment ------------

subtype_diffs <- list()
for (sname in names(wt_enrichments)) {
  mat <- wt_enrichments[[sname]]$enrichment
  present <- intersect(secretory_subtypes, rownames(mat))
  if (!all(c("SecB epithelium", "SecA epithelium") %in% present)) next

  secb_row <- mat["SecB epithelium", ]
  seca_row <- mat["SecA epithelium", ]
  shared <- intersect(names(secb_row), names(seca_row))
  shared <- setdiff(shared, secretory_subtypes)

  subtype_diffs[[sname]] <- data.frame(
    sample    = sname,
    cell_type = shared,
    secb      = secb_row[shared],
    seca      = seca_row[shared],
    diff      = secb_row[shared] - seca_row[shared],
    stringsAsFactors = FALSE
  )
}

if (length(subtype_diffs) > 0) {
  sd_df <- do.call(rbind, subtype_diffs)
  sd_df$sample_short <- gsub("sfe_", "", sd_df$sample)

  p_subdiff <- ggplot(sd_df,
                       aes(x = cell_type, y = diff, fill = sample_short)) +
    geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = NULL, y = "Enrichment diff (SecB - SecA)",
         title = "SecB vs. SecA: differential neighbor enrichment",
         subtitle = "Positive = SecB more enriched near this type; whole tissue samples",
         fill = "Sample") +
    theme_lab() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

  pdf(file.path(colocal_dir, "secretory_secb_vs_seca_diff.pdf"),
      width = 14, height = 8)
  print(p_subdiff)
  dev.off()
  message("  Saved: secretory_secb_vs_seca_diff.pdf")

  write.csv(sd_df,
            file.path(colocal_dir, "secretory_secb_vs_seca_diff.csv"),
            row.names = FALSE)
}


# =============================================================================
# Summary
# =============================================================================

message("\n", strrep("=", 70))
message("08_colocalization.R complete.")
message(strrep("=", 70))
message("Output directory: ", colocal_dir)
message("Key outputs:")
message("  CSVs:")
message("    - tma_consensus_enrichment_k20.csv")
message("    - wt_consensus_enrichment_k20.csv")
message("    - tma_tumour_consensus_enrichment_k20.csv")
message("    - tma_ft_consensus_enrichment_k20.csv")
message("    - tma_tumour_vs_ft_wilcoxon.csv")
message("    - focal_pairs_per_core_enrichment.csv")
message("    - crossk_tma_results.csv")
message("    - crossk_wt_results.csv")
message("    - cross_context_consistency.csv")
message("    - secretory_neighbor_enrichment.csv")
message("    - secretory_secb_vs_seca_diff.csv")
message("  PDFs:")
message("    - tma_consensus_heatmap_k20.pdf")
message("    - wt_consensus_heatmap_k20.pdf")
message("    - wt_individual_heatmaps_k20.pdf")
message("    - tma_tumour_vs_ft_heatmaps.pdf")
message("    - tumour_vs_ft_volcano.pdf")
message("    - focal_pairs_violin.pdf")
message("    - crossk_tma_curves.pdf")
message("    - crossk_wt_curves.pdf")
message("    - cross_context_scatter.pdf")
message("    - secretory_neighbor_profiles.pdf")
message("    - secretory_cross_enrichment.pdf")
message("    - secretory_secb_vs_seca_diff.pdf")

log_session()
