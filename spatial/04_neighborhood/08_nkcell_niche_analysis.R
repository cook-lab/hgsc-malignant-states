# ============================================================================
# 08_nkcell_niche_analysis.R
# ----------------------------------------------------------------------------
# PURPOSE: Context-dependent NK cell characterization: cells in the immune-rich niche vs cells infiltrating SecB-enriched epithelium (niche-conditioned DEG + correlation).
#
# INPUTS:
#   - SFEs (load_sfe) with cell_label + neighborhood
#
# OUTPUTS:
#   - output/13_macrophage_niche/ nkcell niche tables + figures
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

out_13d <- file.path(out_dir, "13d_nkcell_niche")
if (!dir.exists(out_13d)) dir.create(out_13d, recursive = TRUE)

fig_13d <- file.path(fig_dir, "13d_nkcell_niche")
if (!dir.exists(fig_13d)) dir.create(fig_13d, recursive = TRUE)

message("\n=== Script 13d: NK Cell Niche Analysis ===")

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

# --- NK cell functional state gene sets (curated for Xenium 541-gene panel) --

# CD56dim cytotoxic NK cells
cd56dim_genes <- c("FCGR3A", "FGFBP2", "CST7", "KLRF1", "PRF1",
                   "GZMA", "GZMB", "GZMH", "CX3CR1")

# CD56bright regulatory / cytokine-producing NK cells
cd56bright_genes <- c("IGFBP2", "KRT81", "SELL", "NCR1", "CD7",
                      "KLRC1", "KLRC2")

# NK activation / effector
activation_genes <- c("NKG7", "GNLY", "HCST", "KLRD1", "CD69",
                      "IFNG", "TNF", "FASLG", "SH2D1A")

# IFN-responding NK cells
ifn_responding_genes <- c("IFI44L", "IFI6", "ISG20", "MX1", "IFIT1",
                          "IFIT3", "RSAD2", "OAS1")

# NK exhaustion / inhibitory receptors
exhaustion_genes <- c("TIGIT", "HAVCR2", "LAG3", "PDCD1",
                      "LILRB1", "LILRB2")

functional_sets <- list(
  cd56dim        = cd56dim_genes,
  cd56bright     = cd56bright_genes,
  activation     = activation_genes,
  ifn_responding = ifn_responding_genes,
  exhaustion     = exhaustion_genes
)

# ============================================================================
# STEP 1: Collect NK cell data across all SFEs
# ============================================================================

message("\n--- Step 1: Collecting NK cell data from all SFEs ---")

nk_meta_list  <- list()
nk_expr_list  <- list()
gradient_list <- list()

# Epithelial neighborhoods for gradient analysis
epi_nbs <- c("nb_4", "nb_7", "nb_9", "nb_3")
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

  # --- NK cell subset ---
  keep <- !is.na(cd$niche_name) &
          cd$cell_label == "NK cell" &
          cd$niche_name %in% c(niche_a, niche_b)

  if (sum(keep) < 10) {
    message("    Skipping: only ", sum(keep), " NK cells in target niches")
    rm(sfe); gc(verbose = FALSE)
    next
  }

  sfe_nk <- sfe[, keep]
  cd_nk  <- cd[keep, ]

  # Expression matrix (logcounts)
  lc <- as.matrix(logcounts(sfe_nk))

  # Build metadata
  meta <- data.frame(
    cell_id    = colnames(sfe_nk),
    sample     = sname,
    niche      = cd_nk$niche_name,
    stringsAsFactors = FALSE
  )

  # Grab pathway scores if present
  pw_cols <- grep("^pathway_", colnames(cd_nk), value = TRUE)
  if (length(pw_cols) > 0) {
    meta <- cbind(meta, cd_nk[, pw_cols, drop = FALSE])
  }

  # Pseudobulk grouping
  if ("core_id" %in% colnames(cd_nk) && !all(is.na(cd_nk$core_id))) {
    meta$group_id <- paste0(sname, "_core", cd_nk$core_id, "_", cd_nk$niche_name)
  } else {
    meta$group_id <- paste0(sname, "_", cd_nk$niche_name)
  }

  nk_meta_list[[sname]] <- meta
  nk_expr_list[[sname]] <- lc

  n_a <- sum(cd_nk$niche_name == niche_a)
  n_b <- sum(cd_nk$niche_name == niche_b)
  message(sprintf("    %s: %d NK cells (%d in nb_2, %d in nb_6)",
                  sname, ncol(sfe_nk), n_a, n_b))

  # --- Spatial gradient: NK cell coords + distance to epithelial niches ---
  coords <- spatialCoords(sfe)
  nk_idx <- which(keep)
  nk_xy  <- coords[nk_idx, , drop = FALSE]

  max_epi_sample <- 20000

  for (enb in epi_nbs) {
    epi_idx <- which(!is.na(cd$nb_id) & cd$nb_id == enb)
    if (length(epi_idx) < 10) next

    if (length(epi_idx) > max_epi_sample) {
      epi_idx <- sample(epi_idx, max_epi_sample)
    }
    epi_xy <- coords[epi_idx, , drop = FALSE]

    chunk_size <- 5000
    nn_dists <- numeric(nrow(nk_xy))
    for (ci in seq(1, nrow(nk_xy), by = chunk_size)) {
      ci_end <- min(ci + chunk_size - 1, nrow(nk_xy))
      nk_chunk <- nk_xy[ci:ci_end, , drop = FALSE]
      dx <- outer(nk_chunk[, 1], epi_xy[, 1], "-")
      dy <- outer(nk_chunk[, 2], epi_xy[, 2], "-")
      d_mat <- sqrt(dx^2 + dy^2)
      nn_dists[ci:ci_end] <- apply(d_mat, 1, min)
    }

    grad_df <- data.frame(
      cell_id     = colnames(sfe_nk),
      sample      = sname,
      niche       = cd_nk$niche_name,
      epi_nb      = epi_nb_labels[enb],
      dist_to_epi = nn_dists,
      stringsAsFactors = FALSE
    )
    gradient_list[[paste0(sname, "_", enb)]] <- grad_df
    message(sprintf("      Distance to %s: median = %.0f um", epi_nb_labels[enb],
                    median(nn_dists)))
  }

  rm(sfe, sfe_nk, lc, cd, cd_nk, coords, nk_xy); gc(verbose = FALSE)
}

# Combine across samples
nk_meta <- do.call(rbind, nk_meta_list)
nk_expr <- do.call(cbind, nk_expr_list)

message(sprintf("\nTotal NK cells collected: %s",
                format(nrow(nk_meta), big.mark = ",")))
message(sprintf("  %s (nb_2): %s",
                niche_a, format(sum(nk_meta$niche == niche_a), big.mark = ",")))
message(sprintf("  %s (nb_6): %s",
                niche_b, format(sum(nk_meta$niche == niche_b), big.mark = ",")))

# ============================================================================
# STEP 2: Pseudobulk DEG analysis (NK cells only, nb_2 vs nb_6)
# ============================================================================

message("\n--- Step 2: Pseudobulk DEG analysis ---")

groups <- unique(nk_meta$group_id)
genes  <- rownames(nk_expr)

group_sizes <- table(nk_meta$group_id)
valid_groups <- names(group_sizes[group_sizes >= 20])
message(sprintf("  Pseudobulk groups: %d total, %d with >= 20 cells",
                length(groups), length(valid_groups)))

if (length(valid_groups) >= 2) {

  pb_counts <- matrix(0, nrow = length(genes), ncol = length(valid_groups),
                      dimnames = list(genes, valid_groups))

  for (grp in valid_groups) {
    idx <- which(nk_meta$group_id == grp)
    pb_counts[, grp] <- rowSums(nk_expr[, idx, drop = FALSE])
  }

  pb_meta <- data.frame(
    group_id = valid_groups,
    niche    = sapply(valid_groups, function(g) {
      nk_meta$niche[nk_meta$group_id == g][1]
    }),
    n_cells  = as.integer(group_sizes[valid_groups]),
    stringsAsFactors = FALSE
  )

  lib_sizes <- colSums(pb_counts)
  pb_expr   <- log2(t(t(pb_counts) / lib_sizes) * 1e6 + 1)

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

  deg_results <- deg_results[order(-abs(deg_results$log2FC)), ]

} else {
  message("  WARNING: Fewer than 2 valid pseudobulk groups. Skipping DEG analysis.")
  deg_results <- data.frame(
    gene = genes, log2FC = NA_real_, mean_nb2 = NA_real_,
    mean_nb6 = NA_real_, p_value = NA_real_, p_adj = NA_real_,
    sig = FALSE, stringsAsFactors = FALSE
  )
}

write.csv(deg_results,
          file.path(out_13d, "deg_nkcell_nb2_vs_nb6.csv"),
          row.names = FALSE)

n_sig <- sum(deg_results$sig, na.rm = TRUE)
n_up  <- sum(deg_results$sig & deg_results$log2FC > 0, na.rm = TRUE)
n_dn  <- sum(deg_results$sig & deg_results$log2FC < 0, na.rm = TRUE)

message(sprintf("  DEGs (padj < 0.05, |log2FC| > 0.25): %d total (%d up in nb_2, %d up in nb_6)",
                n_sig, n_up, n_dn))

# ============================================================================
# STEP 3: NK Cell Functional State Scoring
# ============================================================================

message("\n--- Step 3: Functional state scoring ---")

panel_genes <- rownames(nk_expr)
functional_scores <- data.frame(niche = nk_meta$niche)

for (set_name in names(functional_sets)) {
  gs <- intersect(functional_sets[[set_name]], panel_genes)
  if (length(gs) >= 2) {
    functional_scores[[set_name]] <- colMeans(nk_expr[gs, , drop = FALSE])
    message(sprintf("  %s: %d/%d genes on panel (%s)", set_name, length(gs),
                    length(functional_sets[[set_name]]),
                    paste(gs, collapse = ", ")))
  } else {
    functional_scores[[set_name]] <- NA_real_
    warning(sprintf("  %s: only %d genes on panel, skipping", set_name, length(gs)))
  }
}

# CD56dim / CD56bright ratio (cytotoxic vs regulatory)
if (!is.na(functional_scores$cd56dim[1]) && !is.na(functional_scores$cd56bright[1])) {
  functional_scores$dim_bright_ratio <-
    (functional_scores$cd56dim + 0.01) / (functional_scores$cd56bright + 0.01)
}

# Statistical tests (Wilcoxon) per score
func_stats <- list()
for (score_name in setdiff(colnames(functional_scores), "niche")) {
  vals_a <- functional_scores[[score_name]][functional_scores$niche == niche_a]
  vals_b <- functional_scores[[score_name]][functional_scores$niche == niche_b]

  vals_a <- vals_a[!is.na(vals_a)]
  vals_b <- vals_b[!is.na(vals_b)]

  if (length(vals_a) >= 10 && length(vals_b) >= 10) {
    wt <- tryCatch(
      wilcox.test(vals_a, vals_b, exact = FALSE),
      error = function(e) list(p.value = NA_real_)
    )
    func_stats[[score_name]] <- data.frame(
      score       = score_name,
      mean_nb2    = mean(vals_a),
      mean_nb6    = mean(vals_b),
      median_nb2  = median(vals_a),
      median_nb6  = median(vals_b),
      log2FC      = log2((mean(vals_a) + 0.01) / (mean(vals_b) + 0.01)),
      p_value     = wt$p.value,
      stringsAsFactors = FALSE
    )
  }
}

func_stats_df <- do.call(rbind, func_stats)
func_stats_df$p_adj <- p.adjust(func_stats_df$p_value, method = "BH")
write.csv(func_stats_df,
          file.path(out_13d, "functional_scores_nb2_vs_nb6.csv"),
          row.names = FALSE)

message("  Functional score comparison:")
for (i in seq_len(nrow(func_stats_df))) {
  message(sprintf("    %s: nb_2 = %.4f, nb_6 = %.4f, p = %.2e",
                  func_stats_df$score[i],
                  func_stats_df$mean_nb2[i],
                  func_stats_df$mean_nb6[i],
                  func_stats_df$p_value[i]))
}

# ============================================================================
# STEP 4: Pathway Activity Comparison (NK cells only)
# ============================================================================

message("\n--- Step 4: Pathway activity comparison ---")

pw_cols <- grep("^pathway_", colnames(nk_meta), value = TRUE)

if (length(pw_cols) > 0) {

  pw_stats <- list()

  for (pw in pw_cols) {
    vals_a <- nk_meta[[pw]][nk_meta$niche == niche_a]
    vals_b <- nk_meta[[pw]][nk_meta$niche == niche_b]

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
  pw_stats_df$p_adj <- p.adjust(pw_stats_df$p_value, method = "BH")
  pw_stats_df <- pw_stats_df[order(pw_stats_df$p_adj), ]

  write.csv(pw_stats_df,
            file.path(out_13d, "pathway_comparison_nkcell_nb2_vs_nb6.csv"),
            row.names = FALSE)

  n_sig_pw <- sum(pw_stats_df$p_adj < 0.05, na.rm = TRUE)
  message(sprintf("  %d / %d pathways significantly different (padj < 0.05)",
                  n_sig_pw, nrow(pw_stats_df)))

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

# NK cell-relevant genes for correlation analysis
nk_relevant_genes <- unique(c(
  # NK cell markers
  "CD7", "NCR1", "KLRD1", "KLRF1", "KLRC1", "KLRC2", "HCST", "NCAM1",
  # CD56dim / cytotoxic
  "FCGR3A", "FGFBP2", "CST7", "PRF1", "GZMA", "GZMB", "GZMH", "CX3CR1",
  # CD56bright
  "IGFBP2", "KRT81", "SELL",
  # Activation / effector
  "NKG7", "GNLY", "IFNG", "TNF", "FASLG", "SH2D1A", "CD69",
  # IFN responding
  "IFI44L", "IFI6", "ISG20", "MX1", "IFIT1", "IFIT3", "RSAD2", "OAS1",
  # Exhaustion / checkpoint
  "TIGIT", "HAVCR2", "LAG3", "PDCD1", "LILRB1", "LILRB2",
  # Signalling
  "STAT1", "STAT3", "STAT4", "JAK1", "JAK3", "IRF1", "FYN", "LCK",
  # Chemokines
  "CXCR4", "CXCR6", "CCR5", "CXCL9", "CXCL10"
))

cor_genes <- intersect(nk_relevant_genes, panel_genes)
message(sprintf("  Correlation analysis: %d NK cell-relevant genes on panel", length(cor_genes)))

idx_a <- which(nk_meta$niche == niche_a)
idx_b <- which(nk_meta$niche == niche_b)

if (length(cor_genes) >= 5 && length(idx_a) >= 50 && length(idx_b) >= 50) {

  max_cells <- 50000
  if (length(idx_a) > max_cells) idx_a <- sample(idx_a, max_cells)
  if (length(idx_b) > max_cells) idx_b <- sample(idx_b, max_cells)

  cor_nb2 <- cor(t(nk_expr[cor_genes, idx_a]), method = "spearman")
  cor_nb6 <- cor(t(nk_expr[cor_genes, idx_b]), method = "spearman")

  cor_diff <- cor_nb2 - cor_nb6

  write.csv(as.data.frame(cor_nb2),
            file.path(out_13d, "correlation_nkcell_nb2.csv"))
  write.csv(as.data.frame(cor_nb6),
            file.path(out_13d, "correlation_nkcell_nb6.csv"))
  write.csv(as.data.frame(cor_diff),
            file.path(out_13d, "correlation_difference_nb2_minus_nb6.csv"))

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
            file.path(out_13d, "differential_correlations_nb2_vs_nb6.csv"),
            row.names = FALSE)

  message("  Top differential correlations (nb_2 - nb_6):")
  top_dc <- head(diff_df, 10)
  for (i in seq_len(nrow(top_dc))) {
    message(sprintf("    %s ~ %s: rho_nb2 = %.3f, rho_nb6 = %.3f, diff = %.3f",
                    top_dc$gene_a[i], top_dc$gene_b[i],
                    top_dc$rho_nb2[i], top_dc$rho_nb6[i], top_dc$rho_diff[i]))
  }
} else {
  message("  Insufficient cells for correlation analysis in one or both niches")
  message(sprintf("  nb_2: %d cells, nb_6: %d cells, genes: %d",
                  length(idx_a), length(idx_b), length(cor_genes)))
}

# ============================================================================
# STEP 6: Spatial Gradient Analysis — Functional Scores vs Distance to Epithelium
# ============================================================================

message("\n--- Step 6: Spatial gradient analysis ---")

gradient_df <- do.call(rbind, gradient_list)
rownames(gradient_df) <- NULL

# Merge functional scores into gradient data
func_for_merge <- functional_scores[, setdiff(colnames(functional_scores), "niche"),
                                     drop = FALSE]
func_for_merge$cell_id <- nk_meta$cell_id
gradient_df <- merge(gradient_df, func_for_merge, by = "cell_id", all.x = TRUE)

message(sprintf("  Gradient data: %s rows across %d epithelial neighborhoods",
                format(nrow(gradient_df), big.mark = ","),
                length(unique(gradient_df$epi_nb))))

# Bin distances
gradient_df$dist_bin <- cut(gradient_df$dist_to_epi,
                            breaks = c(0, 25, 50, 100, 150, 200, 300, 500, 1000, Inf),
                            labels = c("0-25", "25-50", "50-100", "100-150",
                                       "150-200", "200-300", "300-500",
                                       "500-1000", ">1000"),
                            include.lowest = TRUE)

bin_midpoints <- c("0-25" = 12.5, "25-50" = 37.5, "50-100" = 75,
                   "100-150" = 125, "150-200" = 175, "200-300" = 250,
                   "300-500" = 400, "500-1000" = 750, ">1000" = 1250)
gradient_df$dist_mid <- bin_midpoints[as.character(gradient_df$dist_bin)]

# Score columns to analyze
score_cols <- c("cd56dim", "cd56bright", "activation", "ifn_responding",
                "exhaustion", "dim_bright_ratio")
score_cols <- intersect(score_cols, colnames(gradient_df))

# Aggregate
grad_agg_list <- list()
for (sc in score_cols) {
  dt <- data.table(
    score     = gradient_df[[sc]],
    dist_bin  = gradient_df$dist_bin,
    dist_mid  = gradient_df$dist_mid,
    nk_niche  = gradient_df$niche,
    epi_nb    = gradient_df$epi_nb
  )
  agg <- dt[!is.na(score),
            .(mean_score = mean(score),
              se         = sd(score) / sqrt(.N),
              n_cells    = .N),
            by = .(dist_bin, dist_mid, nk_niche, epi_nb)]
  setnames(agg, "nk_niche", "niche")
  agg$score_type <- sc
  grad_agg_list[[sc]] <- as.data.frame(agg)
}

grad_agg <- do.call(rbind, grad_agg_list)
rownames(grad_agg) <- NULL

write.csv(grad_agg,
          file.path(out_13d, "spatial_gradient_functional_by_epi_distance.csv"),
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
          file.path(out_13d, "spatial_gradient_correlations.csv"),
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
          file.path(out_13d, "spatial_gradient_raw.csv"),
          row.names = FALSE)

# ============================================================================
# STEP 7: Visualizations
# ============================================================================

message("\n--- Step 7: Generating figures ---")

# --- 7a. Volcano plot -------------------------------------------------------

deg_plot <- deg_results[!is.na(deg_results$p_adj), ]
deg_plot$neg_log10_p <- -log10(deg_plot$p_adj)
deg_plot$neg_log10_p[deg_plot$neg_log10_p > 50] <- 50
deg_plot$category <- ifelse(!deg_plot$sig, "NS",
                     ifelse(deg_plot$log2FC > 0, "Up in nb_2", "Up in nb_6"))

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
  labs(title = "NK Cell DEGs: Immune-rich vs SecB-enriched epithelium",
       subtitle = sprintf("%d significant (padj < 0.05, |log2FC| > 0.25): %d up nb_2, %d up nb_6",
                          n_sig, n_up, n_dn),
       x = expression(log[2]~fold~change~"(nb_2 / nb_6)"),
       y = expression(-log[10]~adjusted~p)) +
  theme_lab(base_size = 8) +
  theme(legend.position = "top")

ggsave(file.path(fig_13d, "volcano_nkcell_nb2_vs_nb6.pdf"),
       p_volcano, width = 7, height = 6)

# --- 7b. Functional score violin plots --------------------------------------

func_long <- tidyr::pivot_longer(
  functional_scores,
  cols = -niche,
  names_to = "score_type",
  values_to = "score"
)

p_func <- ggplot(func_long, aes(x = score_type, y = score, fill = niche)) +
  geom_violin(scale = "width", alpha = 0.7, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.9, linewidth = 0.3,
               position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("Immune niche" = "#8FBC8F",
                                "SecB-enriched epithelium" = "#87CEFA"),
                    name = "Niche") +
  labs(title = "NK cell functional state scores by niche",
       x = NULL, y = "Mean logcounts score") +
  theme_lab(base_size = 8) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "top")

ggsave(file.path(fig_13d, "functional_violin_nb2_vs_nb6.pdf"),
       p_func, width = 7, height = 5)

# --- 7c. Pathway activity dot plot ------------------------------------------

if (exists("pw_stats_df")) {

  pw_stats_df$sig_label <- ifelse(pw_stats_df$p_adj < 0.001, "***",
                           ifelse(pw_stats_df$p_adj < 0.01, "**",
                           ifelse(pw_stats_df$p_adj < 0.05, "*", "ns")))

  p_pathway <- ggplot(pw_stats_df,
                      aes(x = cohens_d,
                          y = reorder(pathway, cohens_d),
                          size = -log10(pmax(p_adj, 1e-300)),
                          color = cohens_d)) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    scale_color_gradient2(low = "#87CEFA", mid = "grey80", high = "#8FBC8F",
                          midpoint = 0, name = "Cohen's d") +
    scale_size_continuous(name = expression(-log[10]~p[adj]), range = c(1, 6)) +
    labs(title = "Pathway activity: NK cells in nb_2 vs nb_6",
         subtitle = "Positive = higher in Immune-rich niche",
         x = "Cohen's d (effect size)", y = NULL) +
    theme_lab(base_size = 8)

  ggsave(file.path(fig_13d, "pathway_dotplot_nkcell_nb2_vs_nb6.pdf"),
         p_pathway, width = 7, height = 6)
}

# --- 7d. Correlation heatmap (difference) -----------------------------------

if (exists("cor_diff") && length(cor_genes) >= 5) {

  pdf(file.path(fig_13d, "correlation_heatmap_diff_nb2_minus_nb6.pdf"),
      width = 10, height = 9)
  Heatmap(cor_diff,
          name = "rho diff\n(nb7 - nb1)",
          col = circlize::colorRamp2(c(-0.5, 0, 0.5),
                                     c("#87CEFA", "white", "#8FBC8F")),
          row_names_gp = gpar(fontsize = 5),
          column_names_gp = gpar(fontsize = 5),
          column_title = "NK cell differential correlation: nb_2 minus nb_6",
          show_row_dend = TRUE,
          show_column_dend = TRUE)
  dev.off()
}

# --- 7e. Spatial gradient line plots ----------------------------------------

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
         subtitle = "NK cells in nb_2 and nb_6, binned by distance to nearest epithelial cell",
         x = expression("Distance to epithelial niche ("*mu*"m, log scale)"),
         y = paste0("Mean ", sc_label, " score")) +
    theme_lab(base_size = 8) +
    theme(legend.position = "top",
          strip.text = element_text(size = rel(0.85)))

  ggsave(file.path(fig_13d, paste0("spatial_gradient_", sc, "_vs_epi_distance.pdf")),
         p_grad, width = 8, height = 7)
}

# Combined gradient correlation heatmap
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
    labs(title = "NK cell functional score-distance correlations by niche",
         subtitle = "Spearman rho (score ~ distance); * padj<0.05, ** <0.01, *** <0.001",
         x = NULL, y = NULL) +
    theme_lab(base_size = 8) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  ggsave(file.path(fig_13d, "spatial_gradient_correlation_heatmap.pdf"),
         p_grad_heat, width = 8, height = 6)
}

# ============================================================================
# STEP 8: Summary table
# ============================================================================

message("\n--- Step 8: Summary ---")

summary_df <- data.frame(
  metric = c(
    "Total NK cells analyzed",
    "NK cells in nb_2 (Immune niche)",
    "NK cells in nb_6 (SecB-enriched epithelium)",
    "Pseudobulk groups (>= 20 cells)",
    "DEGs (padj < 0.05, |log2FC| > 0.25)",
    "DEGs upregulated in nb_2",
    "DEGs upregulated in nb_6",
    "Significant pathways (padj < 0.05)",
    "Genes in correlation analysis",
    "Significant gradient correlations (padj < 0.05)"
  ),
  value = c(
    nrow(nk_meta),
    sum(nk_meta$niche == niche_a),
    sum(nk_meta$niche == niche_b),
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
          file.path(out_13d, "analysis_summary.csv"),
          row.names = FALSE)

print(summary_df)

# ============================================================================
# Done
# ============================================================================

message("\n=== Script 13d Complete ===")
message("Outputs saved to: ", out_13d)
message("Figures saved to: ", fig_13d)
message("\nFiles generated:")
message("  - deg_nkcell_nb2_vs_nb6.csv")
message("  - functional_scores_nb2_vs_nb6.csv")
message("  - pathway_comparison_nkcell_nb2_vs_nb6.csv")
message("  - correlation_nkcell_nb2.csv / nb9.csv / difference.csv")
message("  - differential_correlations_nb2_vs_nb6.csv")
message("  - spatial_gradient_functional_by_epi_distance.csv")
message("  - spatial_gradient_correlations.csv")
message("  - spatial_gradient_raw.csv")
message("  - analysis_summary.csv")
message("  - volcano_nkcell_nb2_vs_nb6.pdf")
message("  - functional_violin_nb2_vs_nb6.pdf")
message("  - pathway_dotplot_nkcell_nb2_vs_nb6.pdf")
message("  - correlation_heatmap_diff_nb2_minus_nb6.pdf")
message("  - spatial_gradient_*_vs_epi_distance.pdf (per score)")
message("  - spatial_gradient_correlation_heatmap.pdf")

log_session()
