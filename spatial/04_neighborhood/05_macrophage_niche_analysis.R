# ============================================================================
# 05_macrophage_niche_analysis.R
# ----------------------------------------------------------------------------
# PURPOSE: Context-dependent Macrophage characterization: cells in the immune-rich niche vs cells infiltrating SecB-enriched epithelium (niche-conditioned DEG + correlation).
#
# INPUTS:
#   - SFEs (load_sfe) with cell_label + neighborhood
#
# OUTPUTS:
#   - output/13_macrophage_niche/ macrophage niche tables + figures
#
# MANUSCRIPT PANEL(S): Backend for Fig 6 lymphocyte-excluded niche narrative.
# RUNTIME TIER: heavy
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

# --- Output directory --------------------------------------------------------

out_13 <- file.path(out_dir, "13_macrophage_niche")
if (!dir.exists(out_13)) dir.create(out_13, recursive = TRUE)

fig_13 <- file.path(fig_dir, "13_macrophage_niche")
if (!dir.exists(fig_13)) dir.create(fig_13, recursive = TRUE)

message("\n=== Script 13: Macrophage Niche Analysis ===")

# --- Load neighborhood assignments (stored externally) -----------------------

nb_assign <- read.csv(file.path(out_dir, "09_neighborhood", "neighborhood_assignments.csv"),
                      stringsAsFactors = FALSE)
nb_assign$niche_name <- nb_names[nb_assign$neighborhood]
message(sprintf("Loaded %s neighborhood assignments",
                format(nrow(nb_assign), big.mark = ",")))

# --- SFE names ---------------------------------------------------------------

sfe_names <- c("sfe_tma_filtered", "sfe_OTB_2384", "sfe_OTB_2417", "sfe_OTB_2432",
               "sfe_OTB_2454", "sfe_OTB_2457", "sfe_OTB_2461",
               "sfe_SP24_24824", "sfe_SP24_25573")

# --- Niche definitions -------------------------------------------------------

niche_a <- "Immune niche"                # nb_2 — dedicated immune neighborhood
niche_b <- "SecB-enriched epithelium"    # nb_6 — SecB end of trajectory, immune cells in tumor context

# --- Macrophage polarization gene sets (curated for Xenium 541-gene panel) ---

# M1-like (pro-inflammatory / anti-tumour)
m1_genes <- c("CD86", "IRF1", "STAT1", "CXCL9", "CXCL10", "CXCL11",
              "IDO1", "NFKB1", "NFKB2", "TNF", "IL15", "ICAM1",
              "COTL1", "TAP1", "CIITA")

# M2-like (tissue remodelling / pro-tumour)
m2_genes <- c("MRC1", "C1QA", "C1QB", "C1QC", "TGFB1", "TGFBI",
              "CD14", "TREM2", "INHBA", "VEGFA", "MMP11")

# TAM-associated (tumour-associated macrophage markers)
tam_genes <- c("HAVCR2", "CD274", "LGALS9", "ADAM17", "ADAM10",
               "CTSS", "ADAMDEC1", "FCGR3A", "SLAMF7", "CCR1")

polarization_sets <- list(m1_like = m1_genes, m2_like = m2_genes, tam = tam_genes)

# ============================================================================
# STEP 1: Collect macrophage data across all SFEs
# ============================================================================

message("\n--- Step 1: Collecting macrophage data from all SFEs ---")

mac_meta_list  <- list()
mac_expr_list  <- list()
gradient_list  <- list()  # for spatial gradient analysis

# Epithelial neighborhoods for gradient analysis
epi_nbs <- c("nb_4", "nb_7", "nb_9", "nb_3")  # nb_10 is epi-stroma interface
epi_nb_labels <- c(
  nb_4  = "Ciliated-mesenchymal",
  nb_7  = "Intermediate epithelium",
  nb_9  = "Vascular-stromal",
  nb_3  = "Fibroblast-rich stroma"
)

for (sname in sfe_names) {

  message("  Loading ", sname, " ...")
  sfe <- load_sfe(sname)

  cd <- as.data.frame(colData(sfe))
  cd$cell_id <- colnames(sfe)

  # Join neighborhood from external assignments
  nb_match <- nb_assign[match(cd$cell_id, nb_assign$cell_id), ]
  cd$neighborhood <- nb_match$neighborhood
  cd$niche_name   <- nb_match$niche_name
  cd$nb_id        <- nb_match$neighborhood

  # --- Macrophage subset for DEG / scoring ---
  keep <- !is.na(cd$niche_name) &
          cd$cell_label == "Macrophage" &
          cd$niche_name %in% c(niche_a, niche_b)

  if (sum(keep) < 10) {
    message("    Skipping: only ", sum(keep), " macrophages in target niches")
    rm(sfe); gc(verbose = FALSE)
    next
  }

  sfe_mac <- sfe[, keep]
  cd_mac  <- cd[keep, ]

  # Expression matrix (logcounts)
  lc <- as.matrix(logcounts(sfe_mac))

  # Build metadata
  meta <- data.frame(
    cell_id    = colnames(sfe_mac),
    sample     = sname,
    niche      = cd_mac$niche_name,
    stringsAsFactors = FALSE
  )

  # Grab pathway scores if present
  pw_cols <- grep("^pathway_", colnames(cd_mac), value = TRUE)
  if (length(pw_cols) > 0) {
    meta <- cbind(meta, cd_mac[, pw_cols, drop = FALSE])
  }

  # For TMA, include core_id for pseudobulk grouping
  if ("core_id" %in% colnames(cd_mac) && !all(is.na(cd_mac$core_id))) {
    meta$group_id <- paste0(sname, "_core", cd_mac$core_id, "_", cd_mac$niche_name)
  } else {
    meta$group_id <- paste0(sname, "_", cd_mac$niche_name)
  }

  mac_meta_list[[sname]] <- meta
  mac_expr_list[[sname]] <- lc

  n_a <- sum(cd_mac$niche_name == niche_a)
  n_b <- sum(cd_mac$niche_name == niche_b)
  message(sprintf("    %s: %d macrophages (%d in nb_2, %d in nb_6)",
                  sname, ncol(sfe_mac), n_a, n_b))

  # --- Spatial gradient: macrophage coords + distance to epithelial niches ---
  coords <- spatialCoords(sfe)
  mac_idx <- which(keep)
  mac_xy  <- coords[mac_idx, , drop = FALSE]

  # For each epithelial neighborhood, compute distance from each macrophage
  # to the nearest cell in that neighborhood (subsample epi cells for speed)
  max_epi_sample <- 20000

  for (enb in epi_nbs) {
    epi_idx <- which(!is.na(cd$nb_id) & cd$nb_id == enb)
    if (length(epi_idx) < 10) next

    # Subsample epithelial cells for speed
    if (length(epi_idx) > max_epi_sample) {
      epi_idx <- sample(epi_idx, max_epi_sample)
    }
    epi_xy <- coords[epi_idx, , drop = FALSE]

    # Compute nearest-neighbor distance from each macrophage to epi cells
    # Process in chunks to avoid memory issues
    chunk_size <- 5000
    nn_dists <- numeric(nrow(mac_xy))
    for (ci in seq(1, nrow(mac_xy), by = chunk_size)) {
      ci_end <- min(ci + chunk_size - 1, nrow(mac_xy))
      mac_chunk <- mac_xy[ci:ci_end, , drop = FALSE]
      # Euclidean distance matrix (chunk x epi), take row min
      dx <- outer(mac_chunk[, 1], epi_xy[, 1], "-")
      dy <- outer(mac_chunk[, 2], epi_xy[, 2], "-")
      d_mat <- sqrt(dx^2 + dy^2)
      nn_dists[ci:ci_end] <- apply(d_mat, 1, min)
    }

    grad_df <- data.frame(
      cell_id     = colnames(sfe_mac),
      sample      = sname,
      niche       = cd_mac$niche_name,
      epi_nb      = epi_nb_labels[enb],
      dist_to_epi = nn_dists,
      stringsAsFactors = FALSE
    )
    gradient_list[[paste0(sname, "_", enb)]] <- grad_df
    message(sprintf("      Distance to %s: median = %.0f um", epi_nb_labels[enb],
                    median(nn_dists)))
  }

  rm(sfe, sfe_mac, lc, cd, cd_mac, coords, mac_xy); gc(verbose = FALSE)
}

# Combine across samples
mac_meta <- do.call(rbind, mac_meta_list)
mac_expr <- do.call(cbind, mac_expr_list)

message(sprintf("\nTotal macrophages collected: %s",
                format(nrow(mac_meta), big.mark = ",")))
message(sprintf("  %s (nb_2): %s",
                niche_a, format(sum(mac_meta$niche == niche_a), big.mark = ",")))
message(sprintf("  %s (nb_6): %s",
                niche_b, format(sum(mac_meta$niche == niche_b), big.mark = ",")))

# ============================================================================
# STEP 2: Pseudobulk DEG analysis (macrophages only, nb_2 vs nb_6)
# ============================================================================

message("\n--- Step 2: Pseudobulk DEG analysis ---")

# Aggregate to pseudobulk by group_id
groups <- unique(mac_meta$group_id)
genes  <- rownames(mac_expr)

# Count cells per group and filter groups with < 20 macrophages
group_sizes <- table(mac_meta$group_id)
valid_groups <- names(group_sizes[group_sizes >= 20])
message(sprintf("  Pseudobulk groups: %d total, %d with >= 20 cells",
                length(groups), length(valid_groups)))

# Aggregate counts (sum) per group
pb_counts <- matrix(0, nrow = length(genes), ncol = length(valid_groups),
                    dimnames = list(genes, valid_groups))

for (grp in valid_groups) {
  idx <- which(mac_meta$group_id == grp)
  pb_counts[, grp] <- rowSums(mac_expr[, idx, drop = FALSE])
}

# Pseudobulk group metadata
pb_meta <- data.frame(
  group_id = valid_groups,
  niche    = sapply(valid_groups, function(g) {
    mac_meta$niche[mac_meta$group_id == g][1]
  }),
  n_cells  = as.integer(group_sizes[valid_groups]),
  stringsAsFactors = FALSE
)

# Log2-CPM normalization
lib_sizes <- colSums(pb_counts)
pb_expr   <- log2(t(t(pb_counts) / lib_sizes) * 1e6 + 1)

# Wilcoxon rank-sum test per gene
idx_a <- which(pb_meta$niche == niche_a)
idx_b <- which(pb_meta$niche == niche_b)

message(sprintf("  Replicates: %d (%s) vs %d (%s)",
                length(idx_a), niche_a, length(idx_b), niche_b))

deg_results <- data.frame(
  gene    = genes,
  log2FC  = NA_real_,
  mean_nb2 = NA_real_,
  mean_nb6 = NA_real_,
  p_value = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(genes)) {
  vals_a <- pb_expr[i, idx_a]
  vals_b <- pb_expr[i, idx_b]

  deg_results$mean_nb2[i] <- mean(vals_a)
  deg_results$mean_nb6[i] <- mean(vals_b)
  deg_results$log2FC[i]   <- mean(vals_a) - mean(vals_b)

  if (length(vals_a) >= 3 && length(vals_b) >= 3) {
    wt <- tryCatch(
      wilcox.test(vals_a, vals_b, exact = FALSE),
      error = function(e) list(p.value = NA_real_)
    )
    deg_results$p_value[i] <- wt$p.value
  }
}

deg_results$p_adj <- p.adjust(deg_results$p_value, method = "BH")
deg_results$sig   <- !is.na(deg_results$p_adj) &
                     deg_results$p_adj < 0.05 &
                     abs(deg_results$log2FC) > 0.25

# Sort by absolute log2FC
deg_results <- deg_results[order(-abs(deg_results$log2FC)), ]

write.csv(deg_results,
          file.path(out_13, "deg_macrophage_nb2_vs_nb6.csv"),
          row.names = FALSE)

n_sig <- sum(deg_results$sig, na.rm = TRUE)
n_up  <- sum(deg_results$sig & deg_results$log2FC > 0, na.rm = TRUE)
n_dn  <- sum(deg_results$sig & deg_results$log2FC < 0, na.rm = TRUE)

message(sprintf("  DEGs (padj < 0.05, |log2FC| > 0.25): %d total (%d up in nb_2, %d up in nb_6)",
                n_sig, n_up, n_dn))

# ============================================================================
# STEP 3: M1/M2/TAM Polarization Scoring
# ============================================================================

message("\n--- Step 3: Polarization scoring ---")

# Score each macrophage for M1, M2, TAM gene sets
panel_genes <- rownames(mac_expr)
polarization_scores <- data.frame(niche = mac_meta$niche)

for (set_name in names(polarization_sets)) {
  gs <- intersect(polarization_sets[[set_name]], panel_genes)
  if (length(gs) >= 2) {
    polarization_scores[[set_name]] <- colMeans(mac_expr[gs, , drop = FALSE])
    message(sprintf("  %s: %d/%d genes on panel", set_name, length(gs),
                    length(polarization_sets[[set_name]])))
  } else {
    polarization_scores[[set_name]] <- NA_real_
    warning(sprintf("  %s: only %d genes on panel, skipping", set_name, length(gs)))
  }
}

# M1/M2 ratio
if (!is.na(polarization_scores$m1_like[1]) && !is.na(polarization_scores$m2_like[1])) {
  polarization_scores$m1_m2_ratio <- (polarization_scores$m1_like + 0.01) /
                                     (polarization_scores$m2_like + 0.01)
}

# Statistical tests (Wilcoxon) per score
polar_stats <- list()
for (score_name in setdiff(colnames(polarization_scores), "niche")) {
  vals_a <- polarization_scores[[score_name]][polarization_scores$niche == niche_a]
  vals_b <- polarization_scores[[score_name]][polarization_scores$niche == niche_b]

  vals_a <- vals_a[!is.na(vals_a)]
  vals_b <- vals_b[!is.na(vals_b)]
  wt <- if (length(vals_a) >= 3 && length(vals_b) >= 3) {
    tryCatch(wilcox.test(vals_a, vals_b, exact = FALSE),
             error = function(e) list(p.value = NA_real_))
  } else {
    list(p.value = NA_real_)
  }
  polar_stats[[score_name]] <- data.frame(
    score       = score_name,
    mean_nb2    = mean(vals_a, na.rm = TRUE),
    mean_nb6    = mean(vals_b, na.rm = TRUE),
    median_nb2  = median(vals_a, na.rm = TRUE),
    median_nb6  = median(vals_b, na.rm = TRUE),
    log2FC      = log2((mean(vals_a, na.rm = TRUE) + 0.01) /
                       (mean(vals_b, na.rm = TRUE) + 0.01)),
    p_value     = wt$p.value,
    stringsAsFactors = FALSE
  )
}

polar_stats_df <- do.call(rbind, polar_stats)
polar_stats_df$p_adj <- p.adjust(polar_stats_df$p_value, method = "BH")
write.csv(polar_stats_df,
          file.path(out_13, "polarization_scores_nb2_vs_nb6.csv"),
          row.names = FALSE)

message("  Polarization score comparison:")
for (i in seq_len(nrow(polar_stats_df))) {
  message(sprintf("    %s: nb_2 = %.4f, nb_6 = %.4f, p = %.2e",
                  polar_stats_df$score[i],
                  polar_stats_df$mean_nb2[i],
                  polar_stats_df$mean_nb6[i],
                  polar_stats_df$p_value[i]))
}

# ============================================================================
# STEP 4: Pathway Activity Comparison (macrophages only)
# ============================================================================

message("\n--- Step 4: Pathway activity comparison ---")

pw_cols <- grep("^pathway_", colnames(mac_meta), value = TRUE)

if (length(pw_cols) > 0) {

  pw_stats <- list()

  for (pw in pw_cols) {
    vals_a <- mac_meta[[pw]][mac_meta$niche == niche_a]
    vals_b <- mac_meta[[pw]][mac_meta$niche == niche_b]

    vals_a <- vals_a[!is.na(vals_a)]
    vals_b <- vals_b[!is.na(vals_b)]

    if (length(vals_a) >= 10 && length(vals_b) >= 10) {
      wt <- tryCatch(
        wilcox.test(vals_a, vals_b, exact = FALSE),
        error = function(e) list(p.value = NA_real_)
      )
      eff <- (mean(vals_a) - mean(vals_b)) / sqrt((sd(vals_a)^2 + sd(vals_b)^2) / 2)

      pw_stats[[pw]] <- data.frame(
        pathway     = sub("^pathway_", "", pw),
        mean_nb2    = mean(vals_a),
        mean_nb6    = mean(vals_b),
        diff        = mean(vals_a) - mean(vals_b),
        cohens_d    = eff,
        p_value     = wt$p.value,
        stringsAsFactors = FALSE
      )
    }
  }

  pw_stats_df <- do.call(rbind, pw_stats)
  if (is.null(pw_stats_df) || nrow(pw_stats_df) == 0) {
    message("  No pathways had sufficient observations for comparison — skipping")
    pw_stats_df <- data.frame()
  } else {
    pw_stats_df$p_adj <- p.adjust(pw_stats_df$p_value, method = "BH")
    pw_stats_df <- pw_stats_df[order(pw_stats_df$p_adj), ]

    write.csv(pw_stats_df,
              file.path(out_13, "pathway_comparison_macrophage_nb2_vs_nb6.csv"),
              row.names = FALSE)
  }

  n_sig_pw <- if (nrow(pw_stats_df) > 0) sum(pw_stats_df$p_adj < 0.05, na.rm = TRUE) else 0L
  message(sprintf("  %d / %d pathways significantly different (padj < 0.05)",
                  n_sig_pw, nrow(pw_stats_df)))

  # Top pathways by effect size
  top_pw <- head(pw_stats_df[order(-abs(pw_stats_df$cohens_d)), ], 10)
  for (i in seq_len(nrow(top_pw))) {
    direction <- ifelse(top_pw$diff[i] > 0, "higher in nb_2", "higher in nb_6")
    message(sprintf("    %s: Cohen's d = %.3f (%s), padj = %.2e",
                    top_pw$pathway[i], top_pw$cohens_d[i], direction,
                    top_pw$p_adj[i]))
  }

} else {
  message("  No pathway scores found in colData. Run 9b_pathway_scoring.R first.")
}

# ============================================================================
# STEP 5: Gene-Gene Correlation Structure Comparison
# ============================================================================

message("\n--- Step 5: Gene-gene correlation structure ---")

# Use macrophage-relevant genes for correlation analysis
mac_relevant_genes <- unique(c(
  # Canonical macrophage markers
  "CD14", "CD68", "C1QA", "C1QB", "C1QC", "MRC1", "ITGAM",
  # M1-associated
  "CD86", "IRF1", "STAT1", "CXCL9", "CXCL10", "CXCL11", "IDO1",
  "TNF", "NFKB1", "NFKB2", "IL15", "ICAM1",
  # M2-associated
  "TGFB1", "TGFBI", "VEGFA", "MMP11", "TREM2", "INHBA",
  # TAM / immune checkpoint
  "HAVCR2", "CD274", "LGALS9", "CTSS", "ADAMDEC1",
  # Complement
  "C1S", "C3", "C7", "CLU",
  # Chemokines
  "CCR1", "CCR5", "CX3CR1", "CXCL3",
  # Antigen presentation
  "CIITA", "TAP1", "TAP2", "TAPBP",
  # Signalling
  "JAK1", "JAK2", "STAT3", "STAT4", "STAT6"
))

# Filter to genes on panel
cor_genes <- intersect(mac_relevant_genes, panel_genes)
message(sprintf("  Correlation analysis: %d macrophage-relevant genes on panel", length(cor_genes)))

# Compute Spearman correlations separately for each niche
idx_a <- which(mac_meta$niche == niche_a)
idx_b <- which(mac_meta$niche == niche_b)

if (length(cor_genes) >= 5 && length(idx_a) >= 50 && length(idx_b) >= 50) {

  # Subsample for efficiency if very large
  max_cells <- 50000
  if (length(idx_a) > max_cells) idx_a <- sample(idx_a, max_cells)
  if (length(idx_b) > max_cells) idx_b <- sample(idx_b, max_cells)

  cor_nb2 <- cor(t(mac_expr[cor_genes, idx_a]), method = "spearman")
  cor_nb6 <- cor(t(mac_expr[cor_genes, idx_b]), method = "spearman")

  # Correlation difference matrix
  cor_diff <- cor_nb2 - cor_nb6

  # Save correlation matrices
  write.csv(as.data.frame(cor_nb2),
            file.path(out_13, "correlation_macrophage_nb2.csv"))
  write.csv(as.data.frame(cor_nb6),
            file.path(out_13, "correlation_macrophage_nb6.csv"))
  write.csv(as.data.frame(cor_diff),
            file.path(out_13, "correlation_difference_nb2_minus_nb6.csv"))

  # Top differential correlations
  upper_idx <- which(upper.tri(cor_diff), arr.ind = TRUE)
  diff_df <- data.frame(
    gene_a   = cor_genes[upper_idx[, 1]],
    gene_b   = cor_genes[upper_idx[, 2]],
    rho_nb2  = cor_nb2[upper_idx],
    rho_nb6  = cor_nb6[upper_idx],
    rho_diff = cor_diff[upper_idx],
    stringsAsFactors = FALSE
  )
  diff_df <- diff_df[order(-abs(diff_df$rho_diff)), ]

  write.csv(diff_df,
            file.path(out_13, "differential_correlations_nb2_vs_nb6.csv"),
            row.names = FALSE)

  message("  Top differential correlations (nb_2 - nb_6):")
  top_dc <- head(diff_df, 10)
  for (i in seq_len(nrow(top_dc))) {
    message(sprintf("    %s ~ %s: rho_nb2 = %.3f, rho_nb6 = %.3f, diff = %.3f",
                    top_dc$gene_a[i], top_dc$gene_b[i],
                    top_dc$rho_nb2[i], top_dc$rho_nb6[i], top_dc$rho_diff[i]))
  }
}

# ============================================================================
# STEP 6: Spatial Gradient Analysis — Polarization vs Distance to Epithelium
# ============================================================================

message("\n--- Step 6: Spatial gradient analysis ---")

gradient_df <- do.call(rbind, gradient_list)
rownames(gradient_df) <- NULL

# Merge polarization scores into gradient data
polar_for_merge <- polarization_scores[, setdiff(colnames(polarization_scores), "niche"),
                                       drop = FALSE]
polar_for_merge$cell_id <- mac_meta$cell_id
gradient_df <- merge(gradient_df, polar_for_merge, by = "cell_id", all.x = TRUE)

message(sprintf("  Gradient data: %s rows across %d epithelial neighborhoods",
                format(nrow(gradient_df), big.mark = ","),
                length(unique(gradient_df$epi_nb))))

# Bin distances (microns) and compute mean scores per bin x niche x epi_nb
gradient_df$dist_bin <- cut(gradient_df$dist_to_epi,
                            breaks = c(0, 25, 50, 100, 150, 200, 300, 500, 1000, Inf),
                            labels = c("0-25", "25-50", "50-100", "100-150",
                                       "150-200", "200-300", "300-500",
                                       "500-1000", ">1000"),
                            include.lowest = TRUE)

# Compute bin midpoints for continuous plotting
bin_midpoints <- c("0-25" = 12.5, "25-50" = 37.5, "50-100" = 75,
                   "100-150" = 125, "150-200" = 175, "200-300" = 250,
                   "300-500" = 400, "500-1000" = 750, ">1000" = 1250)
gradient_df$dist_mid <- bin_midpoints[as.character(gradient_df$dist_bin)]

# Score columns to analyze
score_cols <- c("m1_like", "m2_like", "tam", "m1_m2_ratio")
score_cols <- intersect(score_cols, colnames(gradient_df))

# Aggregate: mean score per distance bin x niche x epi_nb
grad_agg_list <- list()
for (sc in score_cols) {
  dt <- data.table(
    score     = gradient_df[[sc]],
    dist_bin  = gradient_df$dist_bin,
    dist_mid  = gradient_df$dist_mid,
    mac_niche = gradient_df$niche,
    epi_nb    = gradient_df$epi_nb
  )
  agg <- dt[!is.na(score),
            .(mean_score = mean(score),
              se         = sd(score) / sqrt(.N),
              n_cells    = .N),
            by = .(dist_bin, dist_mid, mac_niche, epi_nb)]
  setnames(agg, "mac_niche", "niche")
  agg$score_type <- sc
  grad_agg_list[[sc]] <- as.data.frame(agg)
}

grad_agg <- do.call(rbind, grad_agg_list)
rownames(grad_agg) <- NULL

write.csv(grad_agg,
          file.path(out_13, "spatial_gradient_polarization_by_epi_distance.csv"),
          row.names = FALSE)

# Spearman correlation: score vs distance, per niche x epi_nb x score
grad_cor_list <- list()
for (sc in score_cols) {
  for (nch in c(niche_a, niche_b)) {
    for (enb in unique(gradient_df$epi_nb)) {
      sub <- gradient_df[gradient_df$niche == nch &
                         gradient_df$epi_nb == enb &
                         !is.na(gradient_df[[sc]]), ]
      if (nrow(sub) < 30) next
      ct <- cor.test(sub$dist_to_epi, sub[[sc]], method = "spearman", exact = FALSE)
      grad_cor_list[[length(grad_cor_list) + 1]] <- data.frame(
        score_type = sc,
        niche      = nch,
        epi_nb     = enb,
        rho        = ct$estimate,
        p_value    = ct$p.value,
        n_cells    = nrow(sub),
        stringsAsFactors = FALSE
      )
    }
  }
}

grad_cor_df <- do.call(rbind, grad_cor_list)
rownames(grad_cor_df) <- NULL
grad_cor_df$p_adj <- p.adjust(grad_cor_df$p_value, method = "BH")

write.csv(grad_cor_df,
          file.path(out_13, "spatial_gradient_correlations.csv"),
          row.names = FALSE)

message("  Gradient correlations (score ~ distance to epithelial niche):")
for (i in seq_len(nrow(grad_cor_df))) {
  sig_flag <- ifelse(grad_cor_df$p_adj[i] < 0.05, " *", "")
  message(sprintf("    %s | %s | %s: rho = %.3f, padj = %.2e, n = %d%s",
                  grad_cor_df$score_type[i],
                  ifelse(grad_cor_df$niche[i] == niche_a, "nb_2", "nb_6"),
                  grad_cor_df$epi_nb[i],
                  grad_cor_df$rho[i],
                  grad_cor_df$p_adj[i],
                  grad_cor_df$n_cells[i],
                  sig_flag))
}

# Save raw gradient data
write.csv(gradient_df[, c("cell_id", "sample", "niche", "epi_nb", "dist_to_epi",
                           score_cols)],
          file.path(out_13, "spatial_gradient_raw.csv"),
          row.names = FALSE)

# ============================================================================
# STEP 7: Visualizations
# ============================================================================

# Override ggrastr device to ragg (Cairo unavailable on this system)
options(ggrastr.default.dpi = 300)

message("\n--- Step 7: Generating figures ---")

# --- 6a. Volcano plot: DEG macrophages nb_2 vs nb_6 -------------------------

deg_plot <- deg_results[!is.na(deg_results$p_adj), ]
deg_plot$neg_log10_p <- -log10(deg_plot$p_adj)
deg_plot$neg_log10_p[deg_plot$neg_log10_p > 50] <- 50  # cap for display
deg_plot$category <- ifelse(!deg_plot$sig, "NS",
                     ifelse(deg_plot$log2FC > 0, "Up in nb_2", "Up in nb_6"))

# Label top genes
top_label <- head(deg_results[deg_results$sig == TRUE, ], 30)

p_volcano <- ggplot(deg_plot, aes(x = log2FC, y = neg_log10_p, color = category)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_manual(values = c("NS" = "grey70",
                                "Up in nb_2" = "#8FBC8F",
                                "Up in nb_6" = "#87CEFA"),
                     name = "") +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  ggrepel::geom_text_repel(
    data = top_label,
    aes(label = gene, y = pmin(-log10(p_adj), 50)),
    size = 2.2, max.overlaps = 25, color = "black", segment.color = "grey50"
  ) +
  labs(title = "Macrophage DEGs: Immune-rich vs SecB-enriched epithelium",
       subtitle = sprintf("%d significant (padj < 0.05, |log2FC| > 0.25): %d up nb_2, %d up nb_6",
                          n_sig, n_up, n_dn),
       x = expression(log[2]~fold~change~"(nb_2 / nb_6)"),
       y = expression(-log[10]~adjusted~p)) +
  theme_lab(base_size = 8) +
  theme(legend.position = "top")

ggsave(file.path(fig_13, "volcano_macrophage_nb2_vs_nb6.pdf"),
       p_volcano, width = 7, height = 6)

# --- 6b. Polarization score violin plots ------------------------------------

polar_long <- tidyr::pivot_longer(
  polarization_scores,
  cols = -niche,
  names_to = "score_type",
  values_to = "score"
)

p_polar <- ggplot(polar_long, aes(x = score_type, y = score, fill = niche)) +
  geom_violin(scale = "width", alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.9, linewidth = 0.3,
               position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Immune niche" = "#8FBC8F",
                                "SecB-enriched epithelium" = "#87CEFA"),
                    name = "Niche") +
  labs(title = "Macrophage polarization scores by niche",
       x = NULL, y = "Mean logcounts score") +
  theme_lab(base_size = 8) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "top")

ggsave(file.path(fig_13, "polarization_violin_nb2_vs_nb6.pdf"),
       p_polar, width = 6, height = 5)

# --- 6c. Pathway activity dot plot ------------------------------------------

if (exists("pw_stats_df")) {

  pw_stats_df$sig_label <- ifelse(pw_stats_df$p_adj < 0.001, "***",
                           ifelse(pw_stats_df$p_adj < 0.01, "**",
                           ifelse(pw_stats_df$p_adj < 0.05, "*", "ns")))

  p_pathway <- ggplot(pw_stats_df,
                      aes(x = cohens_d,
                          y = reorder(pathway, cohens_d),
                          size = -log10(p_adj),
                          color = cohens_d)) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    scale_color_gradient2(low = "#87CEFA", mid = "grey80", high = "#8FBC8F",
                          midpoint = 0, name = "Cohen's d") +
    scale_size_continuous(name = expression(-log[10]~p[adj]), range = c(1, 6)) +
    labs(title = "Pathway activity: macrophages in nb_2 vs nb_6",
         subtitle = "Positive = higher in Immune-rich niche",
         x = "Cohen's d (effect size)", y = NULL) +
    theme_lab(base_size = 8)

  ggsave(file.path(fig_13, "pathway_dotplot_macrophage_nb2_vs_nb6.pdf"),
         p_pathway, width = 7, height = 6)
}

# --- 6d. Correlation heatmap (difference) -----------------------------------

if (exists("cor_diff") && length(cor_genes) >= 5) {

  pdf(file.path(fig_13, "correlation_heatmap_diff_nb2_minus_nb6.pdf"),
      width = 10, height = 9)
  Heatmap(cor_diff,
          name = "rho diff\n(nb8 - nb9)",
          col = circlize::colorRamp2(c(-0.5, 0, 0.5),
                                     c("#87CEFA", "white", "#8FBC8F")),
          row_names_gp = gpar(fontsize = 5),
          column_names_gp = gpar(fontsize = 5),
          column_title = "Differential correlation: Immune-rich minus SecB-enriched epithelium",
          show_row_dend = TRUE,
          show_column_dend = TRUE)
  dev.off()
}

# --- 7e. Spatial gradient: polarization vs distance to epithelium ------------

# Faceted line plot: score vs distance, colored by niche, faceted by epi_nb
for (sc in score_cols) {
  sub_agg <- grad_agg[grad_agg$score_type == sc & grad_agg$n_cells >= 20, ]
  if (nrow(sub_agg) < 5) next

  sc_label <- gsub("_", " ", sc)
  sc_label <- paste0(toupper(substring(sc_label, 1, 1)), substring(sc_label, 2))

  p_grad <- ggplot(sub_agg,
                   aes(x = dist_mid, y = mean_score,
                       color = niche, group = niche)) +
    geom_line(linewidth = 0.8) +
    geom_point(aes(size = n_cells), alpha = 0.8) +
    geom_ribbon(aes(ymin = mean_score - se, ymax = mean_score + se, fill = niche),
                alpha = 0.15, color = NA) +
    facet_wrap(~ epi_nb, scales = "free_y", ncol = 2) +
    scale_color_manual(values = c("Macrophage-dominant niche" = "#8FBC8F",
                                  "Lymphoid-rich niche" = "#87CEFA"),
                       name = "Niche") +
    scale_fill_manual(values = c("Immune niche" = "#8FBC8F",
                                 "SecB-enriched epithelium" = "#87CEFA"),
                      guide = "none") +
    scale_size_continuous(name = "Cells", range = c(0.8, 3.5)) +
    scale_x_continuous(trans = "log10",
                       breaks = c(25, 50, 100, 200, 500, 1000),
                       labels = c("25", "50", "100", "200", "500", "1000")) +
    labs(title = sprintf("Spatial gradient: %s score vs distance to epithelial niche", sc_label),
         subtitle = "Macrophages in nb_2 and nb_6, binned by distance to nearest epithelial cell",
         x = expression("Distance to epithelial niche ("*mu*"m, log scale)"),
         y = paste0("Mean ", sc_label, " score")) +
    theme_lab(base_size = 8) +
    theme(legend.position = "top",
          strip.text = element_text(size = rel(0.85)))

  ggsave(file.path(fig_13, paste0("spatial_gradient_", sc, "_vs_epi_distance.pdf")),
         p_grad, width = 8, height = 7)
}

# Combined heatmap: rho values for all gradient correlations
if (nrow(grad_cor_df) > 0) {
  grad_cor_df$niche_short <- ifelse(grad_cor_df$niche == niche_a, "nb_2", "nb_6")
  grad_cor_df$label <- paste0(grad_cor_df$score_type, " | ", grad_cor_df$niche_short)
  grad_cor_df$sig_star <- ifelse(grad_cor_df$p_adj < 0.001, "***",
                          ifelse(grad_cor_df$p_adj < 0.01, "**",
                          ifelse(grad_cor_df$p_adj < 0.05, "*", "")))

  p_grad_heat <- ggplot(grad_cor_df,
                        aes(x = epi_nb, y = label, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f%s", rho, sig_star)), size = 2.5) +
    scale_fill_gradient2(low = "#87CEFA", mid = "white", high = "#8FBC8F",
                         midpoint = 0, name = "Spearman rho",
                         limits = c(-0.3, 0.3)) +
    labs(title = "Polarization-distance correlations by niche and epithelial neighborhood",
         subtitle = "Spearman rho (score ~ distance); * padj<0.05, ** <0.01, *** <0.001",
         x = NULL, y = NULL) +
    theme_lab(base_size = 8) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  ggsave(file.path(fig_13, "spatial_gradient_correlation_heatmap.pdf"),
         p_grad_heat, width = 8, height = 5)
}

# ============================================================================
# STEP 8: Summary table
# ============================================================================

message("\n--- Step 8: Summary ---")

summary_df <- data.frame(
  metric = c(
    "Total macrophages analyzed",
    "Macrophages in nb_2 (Immune niche)",
    "Macrophages in nb_6 (SecB-enriched epithelium)",
    "Pseudobulk groups (>= 20 cells)",
    "DEGs (padj < 0.05, |log2FC| > 0.25)",
    "DEGs upregulated in nb_2",
    "DEGs upregulated in nb_6",
    "Significant pathways (padj < 0.05)",
    "Genes in correlation analysis",
    "Significant gradient correlations (padj < 0.05)"
  ),
  value = c(
    nrow(mac_meta),
    sum(mac_meta$niche == niche_a),
    sum(mac_meta$niche == niche_b),
    length(valid_groups),
    n_sig,
    n_up,
    n_dn,
    ifelse(exists("pw_stats_df"), sum(pw_stats_df$p_adj < 0.05, na.rm = TRUE), NA),
    length(cor_genes),
    ifelse(exists("grad_cor_df"), sum(grad_cor_df$p_adj < 0.05, na.rm = TRUE), NA)
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_df,
          file.path(out_13, "analysis_summary.csv"),
          row.names = FALSE)

print(summary_df)

# ============================================================================
# Done
# ============================================================================

message("\n=== Script 13 Complete ===")
message("Outputs saved to: ", out_13)
message("Figures saved to: ", fig_13)
message("\nFiles generated:")
message("  - deg_macrophage_nb2_vs_nb6.csv")
message("  - polarization_scores_nb2_vs_nb6.csv")
message("  - pathway_comparison_macrophage_nb2_vs_nb6.csv")
message("  - correlation_macrophage_nb2.csv / nb9.csv / difference.csv")
message("  - differential_correlations_nb2_vs_nb6.csv")
message("  - spatial_gradient_polarization_by_epi_distance.csv")
message("  - spatial_gradient_correlations.csv")
message("  - spatial_gradient_raw.csv")
message("  - analysis_summary.csv")
message("  - volcano_macrophage_nb2_vs_nb6.pdf")
message("  - polarization_violin_nb2_vs_nb6.pdf")
message("  - pathway_dotplot_macrophage_nb2_vs_nb6.pdf")
message("  - correlation_heatmap_diff_nb2_minus_nb6.pdf")
message("  - spatial_gradient_*_vs_epi_distance.pdf (per score)")
message("  - spatial_gradient_correlation_heatmap.pdf")

log_session()
