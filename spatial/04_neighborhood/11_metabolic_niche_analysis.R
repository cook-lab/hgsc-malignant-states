# ============================================================================
# 11_metabolic_niche_analysis.R
# ----------------------------------------------------------------------------
# PURPOSE: Metabolic niche characterization (hypoxia / glycolysis pathway context) across SecA / Intermediate / SecB and surrounding cell types.
#
# INPUTS:
#   - SFEs (load_sfe) with cell_label + pathway_* (9b)
#
# OUTPUTS:
#   - output/13_macrophage_niche/ metabolic niche tables + figures
#
# MANUSCRIPT PANEL(S): Backend for Fig 6 metabolic-niche narrative.
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

out_13m <- file.path(out_dir, "13m_metabolic_niche")
if (!dir.exists(out_13m)) dir.create(out_13m, recursive = TRUE)

fig_13m <- file.path(fig_dir, "13m_metabolic_niche")
if (!dir.exists(fig_13m)) dir.create(fig_13m, recursive = TRUE)

message("\n=== Script 13m: Metabolic Niche Analysis (All Neighborhoods) ===")

# --- Load neighborhood assignments (stored externally) -----------------------

nb_assign <- read.csv(file.path(out_dir, "09_neighborhood", "neighborhood_assignments.csv"),
                      stringsAsFactors = FALSE)
nb_assign$niche_name <- nb_names[nb_assign$neighborhood]

# --- SFE names ---------------------------------------------------------------

sfe_names <- c("sfe_tma_filtered", sfe_names_wt)

# --- All 10 neighborhoods ---------------------------------------------------

all_niches <- nb_names[paste0("nb_", 1:10)]
message("Neighborhoods:")
for (i in seq_along(all_niches)) {
  message(sprintf("  nb_%d: %s", i, all_niches[i]))
}

# --- Metabolic gene sets (curated for Xenium 477-gene panel) -----------------

glycolysis_genes <- c("LDHA", "ALDOA", "ENO1", "PGK1", "SLC2A1",
                      "SLC16A3", "PDK1", "HIF1A", "EPAS1")

oxphos_genes <- c("CYC1", "NDUFV1", "IDH1", "IDH2", "CPT1A")

hypoxia_genes <- c("HIF1A", "EPAS1", "VEGFA", "SLC2A1", "LDHA", "PDK1",
                   "ENO1", "ALDOA", "PGK1", "SLC16A3", "SERPINE1")

antioxidant_genes <- c("NFE2L2", "CAT", "PRDX1", "GPX3", "TXN")

mtor_genes <- c("MTOR", "TSC2", "RHEB", "AKT1", "AKT2", "AKT3",
                "PIK3CA", "PIK3CB", "PIK3R1", "PTEN", "STK11", "RICTOR")

aa_catabolism_genes <- c("IDO1", "SLC7A2", "IDH1", "IDH2")

metabolic_sets <- list(
  glycolysis    = glycolysis_genes,
  oxphos        = oxphos_genes,
  hypoxia       = hypoxia_genes,
  antioxidant   = antioxidant_genes,
  mtor_signal   = mtor_genes,
  aa_catabolism = aa_catabolism_genes
)

# Cell types for matched analysis
target_celltypes <- c("Macrophage", "T cell", "B cell", "NK cell",
                      "Fibroblast", "Endothelial", "Plasma cell",
                      "Conventional dendritic cell",
                      "SecA epithelium", "SecB epithelium",
                      "Intermediate epithelium", "Mesothelial", "Pericyte")

# ============================================================================
# STEP 1: Collect cell data across all SFEs
# ============================================================================

message("\n--- Step 1: Collecting cell data from all SFEs ---")

all_meta_list <- list()
all_expr_list <- list()

for (sname in sfe_names) {

  message("  Loading ", sname, " ...")
  sfe <- load_sfe(sname)

  cd <- as.data.frame(colData(sfe))
  cd$cell_id <- colnames(sfe)

  # Join neighborhood from external assignments
  nb_match <- nb_assign[match(cd$cell_id, nb_assign$cell_id), ]
  cd$neighborhood <- nb_match$neighborhood
  cd$niche_name   <- nb_match$niche_name

  # Keep all assigned cells (all 10 neighborhoods)
  keep <- !is.na(cd$niche_name)

  if (sum(keep) < 10) {
    rm(sfe); gc(verbose = FALSE)
    next
  }

  sfe_sub <- sfe[, keep]
  cd_sub  <- cd[keep, ]

  # All metabolic genes on panel
  all_met_genes <- unique(unlist(metabolic_sets))
  panel_genes <- rownames(sfe_sub)
  met_on_panel <- intersect(all_met_genes, panel_genes)

  # Expression for metabolic genes
  lc <- as.matrix(logcounts(sfe_sub[met_on_panel, ]))

  # Build metadata
  meta <- data.frame(
    cell_id    = colnames(sfe_sub),
    sample     = sname,
    niche      = cd_sub$niche_name,
    nb_id      = cd_sub$neighborhood,
    cell_type  = cd_sub$cell_label,
    stringsAsFactors = FALSE
  )

  # Rename-mismatch fix: deposited SFEs still carry the legacy epithelial
  # label; downstream target_celltypes filters on "Intermediate epithelium".
  # Idempotent (harmless if the legacy value is absent).
  meta$cell_type[meta$cell_type == "Transitioning epithelium"] <- "Intermediate epithelium"

  # Grab pathway scores if present
  pw_cols <- grep("^pathway_", colnames(cd_sub), value = TRUE)
  if (length(pw_cols) > 0) {
    meta <- cbind(meta, cd_sub[, pw_cols, drop = FALSE])
  }

  all_meta_list[[sname]] <- meta
  all_expr_list[[sname]] <- lc

  message(sprintf("    %s: %d cells across %d neighborhoods",
                  sname, ncol(sfe_sub),
                  length(unique(cd_sub$niche_name[!is.na(cd_sub$niche_name)]))))

  rm(sfe, sfe_sub, lc, cd, cd_sub); gc(verbose = FALSE)
}

# Combine
all_meta <- do.call(rbind, all_meta_list)
all_expr <- do.call(cbind, all_expr_list)
rownames(all_meta) <- NULL

message(sprintf("\nTotal cells collected: %s across %d neighborhoods",
                format(nrow(all_meta), big.mark = ","),
                length(unique(all_meta$niche))))

for (nch in all_niches) {
  message(sprintf("  %s: %s cells",
                  nch, format(sum(all_meta$niche == nch), big.mark = ",")))
}

# ============================================================================
# STEP 2: Per-cell metabolic scoring
# ============================================================================

message("\n--- Step 2: Per-cell metabolic scoring ---")

panel_genes <- rownames(all_expr)
met_scores <- data.frame(row.names = seq_len(nrow(all_meta)))

for (set_name in names(metabolic_sets)) {
  gs <- intersect(metabolic_sets[[set_name]], panel_genes)
  if (length(gs) >= 2) {
    met_scores[[set_name]] <- colMeans(all_expr[gs, , drop = FALSE])
    message(sprintf("  %s: %d genes (%s)", set_name, length(gs), paste(gs, collapse = ", ")))
  } else {
    met_scores[[set_name]] <- NA_real_
    message(sprintf("  %s: only %d genes, skipping", set_name, length(gs)))
  }
}

# Warburg index
glyc_genes <- intersect(metabolic_sets$glycolysis, panel_genes)
oxph_genes <- intersect(metabolic_sets$oxphos, panel_genes)
met_scores$warburg_index <- (colMeans(all_expr[glyc_genes, , drop = FALSE]) + 0.01) /
                            (colMeans(all_expr[oxph_genes, , drop = FALSE]) + 0.01)

all_meta <- cbind(all_meta, met_scores)

score_cols <- c(names(metabolic_sets), "warburg_index")
score_cols <- intersect(score_cols, colnames(all_meta))

# ============================================================================
# STEP 3: Neighborhood-level comparison (all cells, Kruskal-Wallis)
# ============================================================================

message("\n--- Step 3: Neighborhood-level metabolic comparison ---")

# Mean scores per neighborhood
nb_means_list <- list()
for (sc in score_cols) {
  for (nch in all_niches) {
    vals <- all_meta[[sc]][all_meta$niche == nch]
    vals <- vals[!is.na(vals)]
    nb_means_list[[paste0(nch, "_", sc)]] <- data.frame(
      score       = sc,
      neighborhood = nch,
      mean        = mean(vals),
      median      = median(vals),
      sd          = sd(vals),
      n_cells     = length(vals),
      stringsAsFactors = FALSE
    )
  }
}
nb_means <- do.call(rbind, nb_means_list)
rownames(nb_means) <- NULL

write.csv(nb_means,
          file.path(out_13m, "neighborhood_metabolic_means.csv"),
          row.names = FALSE)

# Kruskal-Wallis test per score
kw_results <- list()
for (sc in score_cols) {
  vals <- all_meta[[sc]][!is.na(all_meta[[sc]])]
  grps <- all_meta$niche[!is.na(all_meta[[sc]])]
  kw <- kruskal.test(vals ~ grps)
  kw_results[[sc]] <- data.frame(
    score = sc,
    chi_sq = kw$statistic,
    df = kw$parameter,
    p_value = kw$p.value,
    stringsAsFactors = FALSE
  )
  message(sprintf("  %s: Kruskal-Wallis chi2 = %.0f, p = %.2e",
                  sc, kw$statistic, kw$p.value))
}
kw_df <- do.call(rbind, kw_results)
rownames(kw_df) <- NULL

write.csv(kw_df,
          file.path(out_13m, "kruskal_wallis_all_neighborhoods.csv"),
          row.names = FALSE)

# ============================================================================
# STEP 4: Cell-type-matched comparison across ALL neighborhoods
# ============================================================================

message("\n--- Step 4: Cell-type-matched metabolic comparison ---")

# For each cell type x score: compute mean per neighborhood, run Kruskal-Wallis,
# and pairwise Dunn tests

ct_nb_means_list <- list()
ct_kw_list <- list()
ct_pairwise_list <- list()

for (ct in target_celltypes) {

  ct_data <- all_meta[all_meta$cell_type == ct, ]

  # Which neighborhoods have enough cells?
  nb_counts <- table(ct_data$niche)
  valid_nbs <- names(nb_counts[nb_counts >= 30])

  if (length(valid_nbs) < 2) {
    message(sprintf("  %s: only %d neighborhoods with >=30 cells, skipping",
                    ct, length(valid_nbs)))
    next
  }

  ct_valid <- ct_data[ct_data$niche %in% valid_nbs, ]

  for (sc in score_cols) {
    vals <- ct_valid[[sc]]
    grps <- ct_valid$niche

    if (all(is.na(vals))) next

    # Means per neighborhood
    for (nch in valid_nbs) {
      v <- vals[grps == nch]
      v <- v[!is.na(v)]
      ct_nb_means_list[[paste0(ct, "_", nch, "_", sc)]] <- data.frame(
        cell_type    = ct,
        score        = sc,
        neighborhood = nch,
        mean         = mean(v),
        median       = median(v),
        sd           = sd(v),
        n_cells      = length(v),
        stringsAsFactors = FALSE
      )
    }

    # Kruskal-Wallis
    kw <- kruskal.test(vals[!is.na(vals)] ~ grps[!is.na(vals)])
    ct_kw_list[[paste0(ct, "_", sc)]] <- data.frame(
      cell_type = ct,
      score     = sc,
      n_neighborhoods = length(valid_nbs),
      chi_sq    = kw$statistic,
      df        = kw$parameter,
      p_value   = kw$p.value,
      stringsAsFactors = FALSE
    )

    # Pairwise Wilcoxon with BH correction
    pw_test <- tryCatch({
      pairwise.wilcox.test(vals[!is.na(vals)], grps[!is.na(vals)],
                           p.adjust.method = "BH", exact = FALSE)
    }, error = function(e) NULL)

    if (!is.null(pw_test)) {
      pw_mat <- pw_test$p.value
      for (ri in seq_len(nrow(pw_mat))) {
        for (ci in seq_len(ncol(pw_mat))) {
          if (!is.na(pw_mat[ri, ci])) {
            nch_a <- colnames(pw_mat)[ci]
            nch_b <- rownames(pw_mat)[ri]

            # Get means for effect size
            vals_a <- vals[grps == nch_a & !is.na(vals)]
            vals_b <- vals[grps == nch_b & !is.na(vals)]

            if (length(vals_a) >= 10 && length(vals_b) >= 10) {
              eff <- (mean(vals_a) - mean(vals_b)) /
                     sqrt((sd(vals_a)^2 + sd(vals_b)^2) / 2)
            } else {
              eff <- NA_real_
            }

            ct_pairwise_list[[length(ct_pairwise_list) + 1]] <- data.frame(
              cell_type = ct,
              score     = sc,
              niche_a   = nch_a,
              niche_b   = nch_b,
              diff      = mean(vals_a) - mean(vals_b),
              cohens_d  = eff,
              p_adj     = pw_mat[ri, ci],
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }

  message(sprintf("  %s: %d neighborhoods tested", ct, length(valid_nbs)))
}

ct_nb_means_df <- do.call(rbind, ct_nb_means_list)
rownames(ct_nb_means_df) <- NULL

ct_kw_df <- do.call(rbind, ct_kw_list)
rownames(ct_kw_df) <- NULL
ct_kw_df$p_adj <- p.adjust(ct_kw_df$p_value, method = "BH")

ct_pairwise_df <- do.call(rbind, ct_pairwise_list)
rownames(ct_pairwise_df) <- NULL

write.csv(ct_nb_means_df,
          file.path(out_13m, "celltype_neighborhood_metabolic_means.csv"),
          row.names = FALSE)
write.csv(ct_kw_df,
          file.path(out_13m, "celltype_kruskal_wallis.csv"),
          row.names = FALSE)
write.csv(ct_pairwise_df,
          file.path(out_13m, "celltype_pairwise_comparisons.csv"),
          row.names = FALSE)

n_sig_kw <- sum(ct_kw_df$p_adj < 0.05, na.rm = TRUE)
n_sig_pw <- sum(ct_pairwise_df$p_adj < 0.05, na.rm = TRUE)
message(sprintf("\n  Kruskal-Wallis significant: %d / %d cell-type x score tests",
                n_sig_kw, nrow(ct_kw_df)))
message(sprintf("  Pairwise significant: %d / %d comparisons",
                n_sig_pw, nrow(ct_pairwise_df)))

# ============================================================================
# STEP 5: Metabolic-immune crosstalk per neighborhood
# ============================================================================

message("\n--- Step 5: Metabolic-immune crosstalk per neighborhood ---")

pw_cols <- grep("^pathway_", colnames(all_meta), value = TRUE)

if (length(pw_cols) > 0) {

  crosstalk_list <- list()

  for (nch in all_niches) {
    nch_idx <- which(all_meta$niche == nch)

    # Subsample for speed
    if (length(nch_idx) > 50000) nch_idx <- sample(nch_idx, 50000)

    for (met_sc in score_cols) {
      for (pw in pw_cols) {
        x <- all_meta[[met_sc]][nch_idx]
        y <- all_meta[[pw]][nch_idx]

        valid <- !is.na(x) & !is.na(y)
        if (sum(valid) < 100) next

        ct <- cor.test(x[valid], y[valid], method = "spearman", exact = FALSE)
        crosstalk_list[[length(crosstalk_list) + 1]] <- data.frame(
          metabolic_score = met_sc,
          pathway         = sub("^pathway_", "", pw),
          neighborhood    = nch,
          rho             = ct$estimate,
          p_value         = ct$p.value,
          n_cells         = sum(valid),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  crosstalk_df <- do.call(rbind, crosstalk_list)
  rownames(crosstalk_df) <- NULL
  crosstalk_df$p_adj <- p.adjust(crosstalk_df$p_value, method = "BH")

  write.csv(crosstalk_df,
            file.path(out_13m, "metabolic_immune_crosstalk_all_neighborhoods.csv"),
            row.names = FALSE)

  message(sprintf("  Metabolic-immune correlations: %d significant / %d total",
                  sum(crosstalk_df$p_adj < 0.05), nrow(crosstalk_df)))
}

# ============================================================================
# STEP 6: Visualizations
# ============================================================================

message("\n--- Step 6: Generating figures ---")

# --- 6a. Neighborhood-level heatmap (mean score per neighborhood) -----------

nb_wide_list <- list()
for (sc in score_cols) {
  sub <- nb_means[nb_means$score == sc, ]
  vals <- sub$mean
  names(vals) <- sub$neighborhood
  nb_wide_list[[sc]] <- vals[all_niches]
}
nb_wide <- do.call(rbind, nb_wide_list)

# Z-score normalize per row (score) for display
nb_z <- t(scale(t(nb_wide)))

pdf(file.path(fig_13m, "neighborhood_metabolic_heatmap.pdf"),
    width = 10, height = 5)
Heatmap(nb_z,
        name = "Z-score",
        col = circlize::colorRamp2(c(-2, 0, 2), c("#87CEFA", "white", "#FF6347")),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 7),
        column_names_rot = 35,
        column_title = "Metabolic scores across all neighborhoods (Z-scored)",
        column_title_gp = gpar(fontsize = 11),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_side = "left")
dev.off()

# Also save the raw version
pdf(file.path(fig_13m, "neighborhood_metabolic_heatmap_raw.pdf"),
    width = 10, height = 5)
Heatmap(nb_wide,
        name = "Mean score",
        col = circlize::colorRamp2(
          c(min(nb_wide, na.rm = TRUE), median(nb_wide, na.rm = TRUE), max(nb_wide, na.rm = TRUE)),
          c("#87CEFA", "white", "#FF6347")),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 7),
        column_names_rot = 35,
        column_title = "Metabolic scores across all neighborhoods (raw means)",
        column_title_gp = gpar(fontsize = 11),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_side = "left")
dev.off()

# --- 6b. Cell-type-matched heatmaps per score (Cohen's d) ------------------

# For each score, make a cell-type x neighborhood heatmap of mean values
for (sc in score_cols) {
  sub <- ct_nb_means_df[ct_nb_means_df$score == sc, ]
  if (nrow(sub) == 0) next

  cts <- unique(sub$cell_type)
  nbs <- unique(sub$neighborhood)

  mat <- matrix(NA, nrow = length(cts), ncol = length(nbs),
                dimnames = list(cts, nbs))
  for (i in seq_len(nrow(sub))) {
    mat[sub$cell_type[i], sub$neighborhood[i]] <- sub$mean[i]
  }

  # Z-score per cell type (row) to show relative differences
  mat_z <- t(scale(t(mat)))

  sc_label <- gsub("_", " ", sc)
  sc_label <- paste0(toupper(substring(sc_label, 1, 1)), substring(sc_label, 2))

  pdf(file.path(fig_13m, paste0("celltype_by_neighborhood_", sc, ".pdf")),
      width = 10, height = 7)
  print(Heatmap(mat_z,
          name = "Z-score",
          col = circlize::colorRamp2(c(-2, 0, 2), c("#87CEFA", "white", "#FF6347")),
          row_names_gp = gpar(fontsize = 9),
          column_names_gp = gpar(fontsize = 7),
          column_names_rot = 35,
          column_title = sprintf("%s score: cell types x neighborhoods (Z-scored within cell type)",
                                 sc_label),
          column_title_gp = gpar(fontsize = 10),
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          row_names_side = "left",
          na_col = "grey90"))
  dev.off()
}

# --- 6c. Box plots: key scores across neighborhoods -----------------------

for (sc in c("glycolysis", "oxphos", "warburg_index", "hypoxia", "mtor_signal")) {
  sc_label <- gsub("_", " ", sc)
  sc_label <- paste0(toupper(substring(sc_label, 1, 1)), substring(sc_label, 2))

  # Order neighborhoods by mean score
  nb_order <- nb_means$neighborhood[nb_means$score == sc]
  nb_order <- nb_order[order(nb_means$mean[nb_means$score == sc])]

  plot_data <- all_meta[, c("niche", sc)]
  colnames(plot_data) <- c("niche", "score_val")
  plot_data$niche <- factor(plot_data$niche, levels = nb_order)
  plot_data <- plot_data[!is.na(plot_data$niche), ]

  # Subsample for plotting speed
  if (nrow(plot_data) > 100000) {
    set.seed(CFG$seed)
    plot_data <- plot_data[sample(nrow(plot_data), 100000), ]
  }

  p <- ggplot(plot_data, aes(x = niche, y = score_val, fill = niche)) +
    geom_boxplot(outlier.size = 0.2, alpha = 0.8, linewidth = 0.3) +
    scale_fill_manual(values = nb_palette, name = "") +
    labs(title = sprintf("%s score across all neighborhoods", sc_label),
         subtitle = "Ordered by mean score",
         x = NULL, y = paste0("Mean ", sc_label, " score")) +
    theme_lab(base_size = 8) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7),
          legend.position = "none")

  ggsave(file.path(fig_13m, paste0("boxplot_all_neighborhoods_", sc, ".pdf")),
         p, width = 10, height = 5)
}

# --- 6d. Cell-type-specific: Macrophage metabolic profile across nbs -------

for (ct in c("Macrophage", "Fibroblast", "T cell")) {
  ct_data <- all_meta[all_meta$cell_type == ct, ]
  nb_counts <- table(ct_data$niche)
  valid_nbs <- names(nb_counts[nb_counts >= 30])

  if (length(valid_nbs) < 3) next

  ct_sub <- ct_data[ct_data$niche %in% valid_nbs, ]

  for (sc in c("glycolysis", "hypoxia", "mtor_signal")) {
    sc_label <- gsub("_", " ", sc)
    sc_label <- paste0(toupper(substring(sc_label, 1, 1)), substring(sc_label, 2))

    # Order by mean
    ct_means <- tapply(ct_sub[[sc]], ct_sub$niche, mean, na.rm = TRUE)
    nb_order <- names(sort(ct_means))

    plot_sub <- ct_sub[, c("niche", sc)]
    colnames(plot_sub) <- c("niche", "score_val")
    plot_sub$niche <- factor(plot_sub$niche, levels = nb_order)

    p <- ggplot(plot_sub, aes(x = niche, y = score_val, fill = niche)) +
      geom_boxplot(outlier.size = 0.2, alpha = 0.8, linewidth = 0.3) +
      scale_fill_manual(values = nb_palette, name = "") +
      labs(title = sprintf("%s: %s score across neighborhoods", ct, sc_label),
           subtitle = sprintf("Cell-type-matched (%s only)", ct),
           x = NULL, y = paste0("Mean ", sc_label, " score")) +
      theme_lab(base_size = 8) +
      theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7),
            legend.position = "none")

    ct_safe <- gsub(" ", "_", tolower(ct))
    ggsave(file.path(fig_13m, paste0("boxplot_", ct_safe, "_", sc, ".pdf")),
           p, width = 10, height = 5)
  }
}

# --- 6e. Metabolic-immune crosstalk heatmap (per neighborhood) -------------

if (exists("crosstalk_df")) {
  # Heatmap of glycolysis-pathway correlations across neighborhoods
  for (met_sc in c("glycolysis", "mtor_signal")) {
    ct_sub <- crosstalk_df[crosstalk_df$metabolic_score == met_sc, ]

    ct_wide <- reshape2::acast(ct_sub, neighborhood ~ pathway, value.var = "rho")

    sc_label <- gsub("_", " ", met_sc)
    sc_label <- paste0(toupper(substring(sc_label, 1, 1)), substring(sc_label, 2))

    pdf(file.path(fig_13m, paste0("crosstalk_", met_sc, "_all_neighborhoods.pdf")),
        width = 11, height = 7)
    print(Heatmap(ct_wide,
            name = "Spearman\nrho",
            col = circlize::colorRamp2(c(-0.4, 0, 0.4),
                                       c("#87CEFA", "white", "#FF6347")),
            row_names_gp = gpar(fontsize = 8),
            column_names_gp = gpar(fontsize = 7),
            column_names_rot = 35,
            column_title = sprintf("%s-immune pathway correlations across neighborhoods",
                                   sc_label),
            column_title_gp = gpar(fontsize = 10),
            cluster_rows = TRUE,
            cluster_columns = TRUE))
    dev.off()
  }
}

# ============================================================================
# STEP 7: Summary table
# ============================================================================

message("\n--- Step 7: Summary ---")

summary_df <- data.frame(
  metric = c(
    "Total cells analyzed",
    "Neighborhoods compared",
    "Metabolic gene sets scored",
    "Metabolic genes on panel",
    "Cell types tested (matched)",
    "KW tests significant (neighborhood-level, padj < 0.05)",
    "KW tests significant (cell-type-matched, padj < 0.05)",
    "Pairwise comparisons significant (padj < 0.05)",
    "Total pairwise comparisons",
    "Metabolic-immune crosstalk correlations significant"
  ),
  value = c(
    nrow(all_meta),
    length(all_niches),
    length(metabolic_sets) + 1,
    length(unique(unlist(metabolic_sets))),
    length(unique(ct_kw_df$cell_type)),
    sum(kw_df$p_value < 0.05, na.rm = TRUE),
    n_sig_kw,
    n_sig_pw,
    nrow(ct_pairwise_df),
    ifelse(exists("crosstalk_df"), sum(crosstalk_df$p_adj < 0.05, na.rm = TRUE), NA)
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_df,
          file.path(out_13m, "analysis_summary.csv"),
          row.names = FALSE)

print(summary_df)

# ============================================================================
# Done
# ============================================================================

message("\n=== Script 13m Complete ===")
message("Outputs saved to: ", out_13m)
message("Figures saved to: ", fig_13m)

log_session()
