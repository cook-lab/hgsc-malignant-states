# ============================================================================
# 10_comprehensive_trajectory.R
# ----------------------------------------------------------------------------
# PURPOSE: Comprehensive secretory trajectory: expand the SecA-dominant -> SecA-mixed -> Intermediate -> SecB succession with adjacent-step and extreme-pole comparisons.
#
# INPUTS:
#   - SFEs (load_sfe) with cell_label + neighborhood
#
# OUTPUTS:
#   - output/13_macrophage_niche/ trajectory DEG + figures
#
# MANUSCRIPT PANEL(S): Backend for Fig 4/5 trajectory framing.
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

message("\n=== Script 13g: Comprehensive Neighborhood Trajectory Analysis ===")

# ============================================================================
# SECTION 0: Configuration
# ============================================================================

# --- Output directories -----------------------------------------------------

out_13g <- file.path(out_dir, "13g_comprehensive_trajectory")
if (!dir.exists(out_13g)) dir.create(out_13g, recursive = TRUE)

fig_13g <- file.path(fig_dir, "13g_comprehensive_trajectory")
if (!dir.exists(fig_13g)) dir.create(fig_13g, recursive = TRUE)

# --- 10 target neighborhoods -----------------------------------------------

# Main secretory trajectory (4)
sec_trajectory <- c("nb_1", "nb_5", "nb_7", "nb_6")

# Comparators (6)
comparators <- c("nb_2", "nb_3", "nb_4", "nb_8", "nb_9", "nb_10")

# All 10 target neighborhoods in display order
target_nbs <- c(sec_trajectory, comparators)

# Labels: short form for plots
nb_labels <- c(
  nb_1  = "SecA-dom",
  nb_5  = "SecA mixed",
  nb_7  = "Intermediate",
  nb_6  = "SecB",
  nb_2  = "Immune",
  nb_3  = "Fibro-stroma",
  nb_4  = "Ciliated-Meso",
  nb_8  = "Meso-niche",
  nb_9  = "Vasc-stromal",
  nb_10 = "Epi-Stroma"
)

# Full labels
nb_labels_full <- c(
  nb_1  = "SecA-dominant epithelium",
  nb_5  = "SecA epithelium (mixed)",
  nb_7  = "Intermediate epithelium",
  nb_6  = "SecB-enriched epithelium",
  nb_2  = "Immune niche",
  nb_3  = "Fibroblast-rich stroma",
  nb_4  = "Ciliated-mesenchymal",
  nb_8  = "Mesothelial niche",
  nb_9  = "Vascular-stromal",
  nb_10 = "Epi-stroma interface"
)

# Colour palette: gold-brown gradient for trajectory + distinct comparators
nb_cols <- c(
  "SecA-dom"       = "#F0C060",
  "SecA mixed"     = "#D9A84C",
  "Intermediate"  = "#C08E48",
  "SecB"           = "#9A7D55",
  "Immune"         = "#6BA3D6",
  "Fibro-stroma"   = "#E07878",
  "Ciliated-Meso"  = "#8BC78B",
  "Meso-niche"     = "#4A9E6E",
  "Vasc-stromal"   = "#6DAE5C",
  "Epi-Stroma"     = "#A97A4A"
)

# --- Immune cell types ------------------------------------------------------

immune_types <- c("Macrophage", "T cell", "B cell", "NK cell",
                  "Plasma cell", "Conventional dendritic cell",
                  "Plasmacytoid dendritic cell", "Mast cell", "Neutrophil")

# Major immune for breakdown (Panel B)
major_immune <- c("Macrophage", "T cell", "NK cell", "B cell")

# --- Immune marker gene lists (from 13f_immune_gene_specific) ---------------

mac_genes <- list(
  `M1-like` = c("CD86", "IRF1", "STAT1", "CXCL9", "CXCL10",
                 "IDO1", "TNF", "ICAM1"),
  `M2/Tissue-remodelling` = c("MRC1", "TREM2", "FOLR1", "CD14", "TGFB1",
                               "VEGFA", "MMP11", "INHBA", "C1QC"),
  `TAM-associated` = c("HAVCR2", "CD274", "LGALS9", "FCGR3A")
)

tc_genes <- list(
  Cytotoxic = c("GZMA", "GZMB", "GZMH", "PRF1", "FGFBP2", "FASLG", "CST7"),
  `Checkpoint/Exhaustion` = c("PDCD1", "HAVCR2", "LAG3", "TIGIT",
                               "CTLA4", "BTLA"),
  Regulatory = c("BATF", "TNFRSF18", "IL2RA", "IL10"),
  `Naive/Memory` = c("IL7R", "SELL", "TCF7", "CCR7", "LTB"),
  Activation = c("CD69", "IFNG", "CD27", "TNF", "ZAP70")
)

nk_genes <- list(
  `CD56dim cytotoxic` = c("FCGR3A", "FGFBP2", "PRF1", "GZMA", "GZMB",
                           "CX3CR1"),
  CD56bright = c("KLRC1", "KLRC2", "NCR1", "SELL"),
  Activation = c("HCST", "KLRD1", "CD69", "IFNG", "FASLG"),
  `Exhaustion/Inhibitory` = c("TIGIT", "HAVCR2", "LILRB1")
)

celltype_genes <- list(
  "Macrophage" = mac_genes,
  "T cell"     = tc_genes,
  "NK cell"    = nk_genes
)

# --- Metabolic gene sets (from 13m) -----------------------------------------

metabolic_sets <- list(
  glycolysis    = c("LDHA", "ALDOA", "ENO1", "PGK1", "SLC2A1",
                    "SLC16A3", "PDK1", "HIF1A", "EPAS1"),
  oxphos        = c("CYC1", "NDUFV1", "IDH1", "IDH2", "CPT1A"),
  hypoxia       = c("HIF1A", "EPAS1", "VEGFA", "SLC2A1", "LDHA", "PDK1",
                    "ENO1", "ALDOA", "PGK1", "SLC16A3", "SERPINE1"),
  antioxidant   = c("NFE2L2", "CAT", "PRDX1", "GPX3", "TXN"),
  mtor_signal   = c("MTOR", "TSC2", "RHEB", "AKT1", "AKT2", "AKT3",
                    "PIK3CA", "PIK3CB", "PIK3R1", "PTEN", "STK11", "RICTOR"),
  aa_catabolism = c("IDO1", "SLC7A2", "IDH1", "IDH2")
)

# --- SFE names --------------------------------------------------------------

sfe_names <- c("sfe_tma", "sfe_OTB_2384", "sfe_OTB_2417", "sfe_OTB_2432",
               "sfe_OTB_2454", "sfe_OTB_2457", "sfe_OTB_2461",
               "sfe_SP24_24824", "sfe_SP24_25573")

# --- Load neighborhood assignments ------------------------------------------

nb_assign <- read.csv(file.path(out_dir, "09_neighborhood",
                                "neighborhood_assignments.csv"),
                      stringsAsFactors = FALSE)
nb_assign$niche_name <- nb_names[nb_assign$neighborhood]
message(sprintf("Loaded %s neighborhood assignments",
                format(nrow(nb_assign), big.mark = ",")))


# ============================================================================
# SECTION 1: Single-pass SFE loading — collect all needed data
# ============================================================================

message("\n--- Section 1: Collecting data from all SFEs (single pass) ---")

# Flatten gene lists for extraction
all_immune_genes <- unique(unlist(celltype_genes))
all_met_genes    <- unique(unlist(metabolic_sets))
all_target_genes <- unique(c(all_immune_genes, all_met_genes))

# Accumulators
comp_list  <- list()   # cell labels for composition (panels A/B)
imm_meta_list  <- list()   # immune cell metadata (panels E/F/G)
imm_expr_list  <- list()   # immune cell expression (panels E/F/G)
met_meta_list  <- list()   # all cell metadata for metabolic scoring (panel D)
met_expr_list  <- list()   # metabolic gene expression (panel D)

for (sname in sfe_names) {

  message("  Loading ", sname, " ...")
  sfe <- load_sfe(sname)

  cd <- as.data.frame(colData(sfe))
  cd$cell_id <- colnames(sfe)

  # Join neighborhood
  nb_match <- nb_assign[match(cd$cell_id, nb_assign$cell_id), ]
  cd$neighborhood <- nb_match$neighborhood

  # Keep cells in 10 target neighborhoods
  keep <- !is.na(cd$neighborhood) & cd$neighborhood %in% target_nbs
  if (sum(keep) < 10) {
    message("    Skipping: only ", sum(keep), " cells in target neighborhoods")
    rm(sfe); gc(verbose = FALSE)
    next
  }

  cd_sub <- cd[keep, ]

  # --- Panel A/B: composition data ---
  comp_list[[sname]] <- data.frame(
    cell_id      = cd_sub$cell_id,
    cell_label   = cd_sub$cell_label,
    neighborhood = cd_sub$neighborhood,
    stringsAsFactors = FALSE
  )

  # --- Panels E/F/G: immune cell expression ---
  imm_types_present <- c("Macrophage", "T cell", "NK cell")
  imm_mask <- keep & cd$cell_label %in% imm_types_present
  if (sum(imm_mask) >= 10) {
    sfe_imm <- sfe[, imm_mask]
    panel_genes <- intersect(all_immune_genes, rownames(sfe_imm))
    lc_imm <- as.matrix(logcounts(sfe_imm[panel_genes, ]))

    imm_meta_list[[sname]] <- data.frame(
      cell_id      = colnames(sfe_imm),
      sample       = sname,
      cell_label   = cd$cell_label[imm_mask],
      neighborhood = cd$neighborhood[imm_mask],
      stringsAsFactors = FALSE
    )
    imm_expr_list[[sname]] <- lc_imm
    rm(sfe_imm, lc_imm)
  }

  # --- Panel D: metabolic gene expression (all cells in target nbs) ---
  sfe_met <- sfe[, keep]
  met_panel <- intersect(all_met_genes, rownames(sfe_met))
  lc_met <- as.matrix(logcounts(sfe_met[met_panel, ]))

  met_meta_list[[sname]] <- data.frame(
    cell_id      = colnames(sfe_met),
    sample       = sname,
    cell_label   = cd_sub$cell_label,
    neighborhood = cd_sub$neighborhood,
    stringsAsFactors = FALSE
  )
  met_expr_list[[sname]] <- lc_met

  # Report counts
  for (nb in target_nbs) {
    n_nb <- sum(cd_sub$neighborhood == nb)
    if (n_nb > 0) {
      message(sprintf("    %s (%s): %s cells",
                      nb, nb_labels[nb], format(n_nb, big.mark = ",")))
    }
  }

  rm(sfe, sfe_met, lc_met, cd, cd_sub); gc(verbose = FALSE)
}

# Combine
comp_df <- do.call(rbind, comp_list)
rownames(comp_df) <- NULL

imm_meta <- do.call(rbind, imm_meta_list)
imm_expr <- do.call(cbind, imm_expr_list)
rownames(imm_meta) <- NULL

met_meta <- do.call(rbind, met_meta_list)
met_expr <- do.call(cbind, met_expr_list)
rownames(met_meta) <- NULL

message(sprintf("\nData collected:"))
message(sprintf("  Composition: %s cells", format(nrow(comp_df), big.mark = ",")))
message(sprintf("  Immune markers: %s cells (%d genes)",
                format(nrow(imm_meta), big.mark = ","), nrow(imm_expr)))
message(sprintf("  Metabolic: %s cells (%d genes)",
                format(nrow(met_meta), big.mark = ","), nrow(met_expr)))

# Add readable labels
imm_meta$nb_label <- factor(nb_labels[imm_meta$neighborhood],
                            levels = nb_labels)
met_meta$nb_label <- factor(nb_labels[met_meta$neighborhood],
                            levels = nb_labels)

rm(comp_list, imm_meta_list, imm_expr_list, met_meta_list, met_expr_list)
gc(verbose = FALSE)


# ============================================================================
# SECTION 2: Cell Type Composition (Panels A & B)
# ============================================================================

message("\n--- Section 2: Cell type composition ---")

# --- Panel A: all cell types, with immune grouped ---

comp_tab <- as.data.frame(table(neighborhood = comp_df$neighborhood,
                                cell_type = comp_df$cell_label))
colnames(comp_tab)[3] <- "count"

comp_tab$neighborhood <- as.character(comp_tab$neighborhood)
comp_tab$cell_type <- as.character(comp_tab$cell_type)
total_per_nb <- tapply(comp_tab$count, comp_tab$neighborhood, sum)
comp_tab$total_in_nb <- total_per_nb[comp_tab$neighborhood]
comp_tab$proportion <- comp_tab$count / comp_tab$total_in_nb
comp_tab$nb_label <- nb_labels[comp_tab$neighborhood]
comp_tab$nb_label_full <- nb_labels_full[comp_tab$neighborhood]

# Grouped version: lump immune types
comp_tab$cell_group <- ifelse(
  comp_tab$cell_type %in% immune_types, "Immune",
  as.character(comp_tab$cell_type)
)

comp_grouped <- aggregate(
  cbind(count, proportion) ~ neighborhood + nb_label + nb_label_full + cell_group,
  data = comp_tab, FUN = sum
)

write.csv(comp_tab,
          file.path(out_13g, "composition_all_celltypes.csv"),
          row.names = FALSE)
write.csv(comp_grouped,
          file.path(out_13g, "composition_grouped.csv"),
          row.names = FALSE)

message("  Composition saved. Cell counts per neighborhood:")
for (nb in target_nbs) {
  message(sprintf("    %s: %s cells",
                  nb_labels[nb], format(total_per_nb[nb], big.mark = ",")))
}

# --- Panel B: immune breakdown ---

imm_comp <- comp_tab[comp_tab$cell_type %in% immune_types, ]

# Total immune per neighborhood
imm_total <- aggregate(count ~ neighborhood, data = imm_comp, FUN = sum)
colnames(imm_total)[2] <- "total_immune"

imm_comp <- merge(imm_comp, imm_total, by = "neighborhood")
imm_comp$pct_of_immune <- 100 * imm_comp$count / imm_comp$total_immune

# Categorise: macrophage, lymphocyte (T/B/NK/Plasma), minor (DC, Mast, Neutrophil)
imm_comp$immune_category <- ifelse(
  imm_comp$cell_type == "Macrophage", "Macrophage",
  ifelse(imm_comp$cell_type %in% c("T cell", "B cell", "NK cell", "Plasma cell"),
         "Lymphocyte",
         "Minor immune")
)

write.csv(imm_comp,
          file.path(out_13g, "composition_immune_breakdown.csv"),
          row.names = FALSE)

message("  Immune breakdown saved.")

rm(comp_df); gc(verbose = FALSE)


# ============================================================================
# SECTION 3: Pathway & Metabolic Means (Panels C & D)
# ============================================================================

message("\n--- Section 3: Pathway and metabolic means ---")

# --- Panel C: Load pre-computed pathway data from 9b ---

pw_data <- read.csv(file.path(out_dir, "9b_scoring",
                               "dotplot_neighborhood_data.csv"),
                    stringsAsFactors = FALSE)

# Target pathways (16)
target_pathways <- c("proliferation", "pi3k_akt_mtor", "notch", "hippo",
                     "angiogenesis", "antigen_presentation", "complement",
                     "nfkb", "type_ii_ifn", "type_i_ifn", "immune_checkpoint",
                     "jak_stat", "cytotoxicity", "chemokine", "tgfb", "emt")

# Subset to target neighborhoods and target pathways
pw_7nb <- pw_data[pw_data$nb_name %in% nb_labels_full &
                  pw_data$pathway %in% target_pathways, ]

# Add short labels
pw_7nb$nb_short <- nb_labels[names(nb_labels_full)[match(pw_7nb$nb_name, nb_labels_full)]]
pw_7nb$nb_short <- factor(pw_7nb$nb_short, levels = nb_labels)

write.csv(pw_7nb, file.path(out_13g, "pathway_means_7nb.csv"),
          row.names = FALSE)
message(sprintf("  Pathways: %d rows (%d pathways × %d neighborhoods)",
                nrow(pw_7nb), length(unique(pw_7nb$pathway)),
                length(unique(pw_7nb$nb_name))))

# --- Panel D: Load pre-computed metabolic data from 13m ---

met_nb_data <- read.csv(file.path(out_dir, "13m_metabolic_niche",
                                   "neighborhood_metabolic_means.csv"),
                        stringsAsFactors = FALSE)

# Subset to target neighborhoods
met_7nb <- met_nb_data[met_nb_data$neighborhood %in% nb_labels_full, ]
met_7nb$nb_short <- nb_labels[names(nb_labels_full)[match(met_7nb$neighborhood, nb_labels_full)]]
met_7nb$nb_short <- factor(met_7nb$nb_short, levels = nb_labels)

write.csv(met_7nb, file.path(out_13g, "metabolic_means_7nb.csv"),
          row.names = FALSE)
message(sprintf("  Metabolic: %d rows (%d scores × %d neighborhoods)",
                nrow(met_7nb), length(unique(met_7nb$score)),
                length(unique(met_7nb$neighborhood))))

# --- Cell-type-specific metabolic data from 13m ---

met_ct_data <- read.csv(file.path(out_dir, "13m_metabolic_niche",
                                   "celltype_neighborhood_metabolic_means.csv"),
                        stringsAsFactors = FALSE)

met_ct_7nb <- met_ct_data[met_ct_data$neighborhood %in% nb_labels_full, ]
met_ct_7nb$nb_short <- nb_labels[names(nb_labels_full)[match(met_ct_7nb$neighborhood, nb_labels_full)]]
met_ct_7nb$nb_short <- factor(met_ct_7nb$nb_short, levels = nb_labels)

write.csv(met_ct_7nb, file.path(out_13g, "metabolic_celltype_means_7nb.csv"),
          row.names = FALSE)
message(sprintf("  Cell-type metabolic: %d rows", nrow(met_ct_7nb)))

# --- Per-cell metabolic scores for Kruskal-Wallis ---

message("  Computing per-cell metabolic scores ...")
panel_genes <- rownames(met_expr)
met_scores <- data.frame(row.names = seq_len(nrow(met_meta)))

for (set_name in names(metabolic_sets)) {
  gs <- intersect(metabolic_sets[[set_name]], panel_genes)
  if (length(gs) >= 2) {
    met_scores[[set_name]] <- colMeans(met_expr[gs, , drop = FALSE])
    message(sprintf("    %s: %d genes", set_name, length(gs)))
  } else {
    met_scores[[set_name]] <- NA_real_
  }
}

# Warburg index
glyc_genes <- intersect(metabolic_sets$glycolysis, panel_genes)
oxph_genes <- intersect(metabolic_sets$oxphos, panel_genes)
met_scores$warburg_index <- (colMeans(met_expr[glyc_genes, , drop = FALSE]) + 0.01) /
                            (colMeans(met_expr[oxph_genes, , drop = FALSE]) + 0.01)

met_meta <- cbind(met_meta, met_scores)
score_cols <- c(names(metabolic_sets), "warburg_index")

# Kruskal-Wallis across all target neighborhoods
met_kw_list <- list()
for (sc in score_cols) {
  vals <- met_meta[[sc]]
  grp  <- met_meta$neighborhood
  valid <- !is.na(vals)
  if (sum(valid) > 100) {
    kw <- kruskal.test(vals[valid] ~ grp[valid])
    met_kw_list[[sc]] <- data.frame(
      score    = sc,
      kw_stat  = kw$statistic,
      kw_p     = kw$p.value,
      stringsAsFactors = FALSE
    )
  }
}
met_kw_df <- do.call(rbind, met_kw_list)
met_kw_df$kw_p_adj <- p.adjust(met_kw_df$kw_p, method = "BH")

write.csv(met_kw_df, file.path(out_13g, "metabolic_kruskalwallis.csv"),
          row.names = FALSE)
message("  Metabolic KW tests saved.")

rm(met_expr, met_scores); gc(verbose = FALSE)


# ============================================================================
# SECTION 4: Immune Marker Expression (Panels E/F/G)
# ============================================================================

message("\n--- Section 4: Immune marker expression ---")

# --- Per gene × cell type × neighborhood: mean, median, % expressing, n ---

marker_means_list <- list()

for (ct in names(celltype_genes)) {
  ct_mask <- imm_meta$cell_label == ct
  n_ct <- sum(ct_mask)
  if (n_ct < 10) next

  gene_groups <- celltype_genes[[ct]]

  for (grp_name in names(gene_groups)) {
    genes <- intersect(gene_groups[[grp_name]], rownames(imm_expr))
    if (length(genes) == 0) next

    for (gene in genes) {
      expr_vec <- imm_expr[gene, ct_mask]
      nb_vec   <- imm_meta$neighborhood[ct_mask]

      for (nb in target_nbs) {
        nb_mask <- nb_vec == nb
        n_cells <- sum(nb_mask)
        if (n_cells < 5) next

        vals <- expr_vec[nb_mask]
        n_pos <- sum(vals > 0)

        marker_means_list[[paste(ct, grp_name, gene, nb)]] <- data.frame(
          cell_type      = ct,
          group          = grp_name,
          gene           = gene,
          neighborhood   = nb,
          nb_label       = nb_labels[nb],
          n_cells        = n_cells,
          mean_expr      = mean(vals),
          median_expr    = median(vals),
          pct_expressing = 100 * n_pos / n_cells,
          mean_expr_pos  = if (n_pos > 0) mean(vals[vals > 0]) else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

marker_means <- do.call(rbind, marker_means_list)
rownames(marker_means) <- NULL

write.csv(marker_means, file.path(out_13g, "immune_marker_means.csv"),
          row.names = FALSE)
message(sprintf("  Marker means: %d rows (%d genes × %d cell types × up to %d nbs)",
                nrow(marker_means),
                length(unique(marker_means$gene)),
                length(unique(marker_means$cell_type)),
                length(target_nbs)))

# --- Pairwise Wilcoxon tests ---

message("  Running pairwise Wilcoxon tests ...")

# Comparisons: adjacent trajectory + SecA vs SecB + each comparator vs SecA + each comparator vs SecB
pairwise_comps <- list(
  list(a = "nb_1",  b = "nb_5",  label = "SecADom_vs_SecAMixed"),
  list(a = "nb_5",  b = "nb_7",  label = "SecAMixed_vs_Trans"),
  list(a = "nb_7",  b = "nb_6",  label = "Trans_vs_SecB"),
  list(a = "nb_1",  b = "nb_6",  label = "SecA_vs_SecB"),
  list(a = "nb_1",  b = "nb_10", label = "SecA_vs_EpiStroma"),
  list(a = "nb_6",  b = "nb_10", label = "SecB_vs_EpiStroma"),
  list(a = "nb_1",  b = "nb_2",  label = "SecA_vs_Immune"),
  list(a = "nb_6",  b = "nb_2",  label = "SecB_vs_Immune"),
  list(a = "nb_1",  b = "nb_3",  label = "SecA_vs_FibroStroma"),
  list(a = "nb_6",  b = "nb_3",  label = "SecB_vs_FibroStroma")
)

stats_list <- list()

for (ct in names(celltype_genes)) {
  ct_mask <- imm_meta$cell_label == ct
  genes_for_ct <- unique(unlist(celltype_genes[[ct]]))
  genes_for_ct <- intersect(genes_for_ct, rownames(imm_expr))

  for (gene in genes_for_ct) {
    expr_vec <- imm_expr[gene, ct_mask]
    nb_vec   <- imm_meta$neighborhood[ct_mask]

    for (comp in pairwise_comps) {
      vals_a <- expr_vec[nb_vec == comp$a]
      vals_b <- expr_vec[nb_vec == comp$b]

      if (length(vals_a) >= 10 && length(vals_b) >= 10) {
        wt <- tryCatch(
          wilcox.test(vals_a, vals_b, exact = FALSE),
          error = function(e) list(p.value = NA_real_)
        )
        pooled_sd <- sqrt((sd(vals_a)^2 + sd(vals_b)^2) / 2)
        d <- if (pooled_sd > 0) (mean(vals_a) - mean(vals_b)) / pooled_sd else 0

        stats_list[[paste(ct, gene, comp$label)]] <- data.frame(
          cell_type  = ct,
          gene       = gene,
          comparison = comp$label,
          mean_a     = mean(vals_a),
          mean_b     = mean(vals_b),
          log2FC     = log2((mean(vals_a) + 0.01) / (mean(vals_b) + 0.01)),
          cohens_d   = d,
          p_value    = wt$p.value,
          n_a        = length(vals_a),
          n_b        = length(vals_b),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

stats_df <- do.call(rbind, stats_list)
rownames(stats_df) <- NULL
stats_df$p_adj <- ave(stats_df$p_value, stats_df$cell_type,
                      FUN = function(x) p.adjust(x, method = "BH"))

# Add gene group labels
group_map <- list()
for (ct in names(celltype_genes)) {
  for (grp in names(celltype_genes[[ct]])) {
    for (g in celltype_genes[[ct]][[grp]]) {
      group_map[[paste(ct, g)]] <- grp
    }
  }
}
stats_df$group <- sapply(paste(stats_df$cell_type, stats_df$gene),
                         function(x) group_map[[x]])

write.csv(stats_df, file.path(out_13g, "immune_marker_stats.csv"),
          row.names = FALSE)

n_sig <- sum(stats_df$p_adj < 0.05, na.rm = TRUE)
message(sprintf("  Marker stats: %d tests, %d significant (padj < 0.05)",
                nrow(stats_df), n_sig))

# --- Proportion expressing analysis ---

message("  Computing proportion expressing ...")

prop_list <- list()

for (ct in names(celltype_genes)) {
  ct_mask <- imm_meta$cell_label == ct
  genes_for_ct <- unique(unlist(celltype_genes[[ct]]))
  genes_for_ct <- intersect(genes_for_ct, rownames(imm_expr))

  for (gene in genes_for_ct) {
    expr_vec <- imm_expr[gene, ct_mask]
    nb_vec   <- imm_meta$neighborhood[ct_mask]

    for (nb in target_nbs) {
      nb_mask <- nb_vec == nb
      n_total <- sum(nb_mask)
      if (n_total < 10) next
      n_pos <- sum(expr_vec[nb_mask] > 0)

      prop_list[[paste(ct, gene, nb)]] <- data.frame(
        cell_type      = ct,
        gene           = gene,
        neighborhood   = nb,
        nb_label       = nb_labels[nb],
        n_total        = n_total,
        n_expressing   = n_pos,
        pct_expressing = 100 * n_pos / n_total,
        mean_expr_pos  = if (n_pos > 0) mean(expr_vec[nb_mask & expr_vec > 0]) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
}

prop_df <- do.call(rbind, prop_list)
rownames(prop_df) <- NULL
prop_df$group <- sapply(paste(prop_df$cell_type, prop_df$gene),
                        function(x) group_map[[x]])

write.csv(prop_df, file.path(out_13g, "immune_marker_proportions.csv"),
          row.names = FALSE)
message(sprintf("  Proportions: %d rows", nrow(prop_df)))


# ============================================================================
# SECTION 5: Summary
# ============================================================================

message("\n--- Section 5: Summary ---")

summary_df <- data.frame(
  neighborhood = target_nbs,
  label        = nb_labels[target_nbs],
  full_label   = nb_labels_full[target_nbs],
  type         = ifelse(target_nbs %in% sec_trajectory, "trajectory", "comparator"),
  n_cells      = total_per_nb[target_nbs],
  stringsAsFactors = FALSE
)

write.csv(summary_df, file.path(out_13g, "analysis_summary.csv"),
          row.names = FALSE)

message("\n=== 13g Complete ===")
message("Output directory: ", out_13g)
message("Files created:")
for (f in list.files(out_13g, pattern = "\\.csv$")) {
  message("  ", f)
}

rm(imm_meta, imm_expr); gc(verbose = FALSE)
