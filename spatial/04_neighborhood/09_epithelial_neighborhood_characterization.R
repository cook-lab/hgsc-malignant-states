# ============================================================================
# 09_epithelial_neighborhood_characterization.R
# ----------------------------------------------------------------------------
# PURPOSE: Characterize the epithelial secretory succession (SecA -> Early Intermediate -> Intermediate -> Late Intermediate -> SecB) across neighborhood membership and nearest-epithelial proximity assignment methods.
#
# INPUTS:
#   - SFEs (load_sfe) with cell_label + neighborhood
#
# OUTPUTS:
#   - output/13_macrophage_niche/ epithelial neighborhood DEG tables + volcanoes
#
# MANUSCRIPT PANEL(S): Backend for Fig 4/5 polarization-trajectory framing.
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

library(circlize)

# ============================================================================
# SECTION 0: Configuration
# ============================================================================

message("\n=== Script 13f: Epithelial Neighborhood Characterization ===")

# --- Output directories -----------------------------------------------------

out_13f <- file.path(out_dir, "13f_epithelial_neighborhood")
if (!dir.exists(out_13f)) dir.create(out_13f, recursive = TRUE)

fig_13f <- file.path(fig_dir, "13f_epithelial_neighborhood")
if (!dir.exists(fig_13f)) dir.create(fig_13f, recursive = TRUE)

# --- Neighborhood definitions -----------------------------------------------

# All epithelial neighborhoods (for composition & pathway analysis)
epi_nbs <- c("nb_1", "nb_5", "nb_7", "nb_6", "nb_10", "nb_4")
epi_nb_labels <- c(
  nb_1  = "SecA-dominant epithelium",
  nb_5  = "SecA epithelium (mixed)",
  nb_7  = "Intermediate epithelium",
  nb_6  = "SecB-enriched epithelium",
  nb_10 = "Epi-stroma interface",
  nb_4  = "Ciliated-mesenchymal"
)

# Secretory trajectory order (4 neighborhoods for pairwise DEG)
sec_trajectory <- c("nb_1", "nb_5", "nb_7", "nb_6")
sec_labels <- epi_nb_labels[sec_trajectory]

# Pairwise comparisons along the trajectory (adjacent + endpoints)
pairwise_comparisons <- list(
  list(a = "nb_1",  b = "nb_6",  label = "SecA_vs_SecB"),
  list(a = "nb_1",  b = "nb_5",  label = "SecA_vs_SecAMixed"),
  list(a = "nb_5",  b = "nb_7",  label = "SecAMixed_vs_Trans"),
  list(a = "nb_7",  b = "nb_6",  label = "Trans_vs_SecB")
)

# --- Cell types --------------------------------------------------------------

immune_types <- c("Macrophage", "T cell", "B cell", "NK cell")

# --- Functional gene sets (from 13/13b/13c/13d) -----------------------------

mac_sets <- list(
  m1_like = c("CD86", "IRF1", "STAT1", "CXCL9", "CXCL10", "CXCL11",
              "IDO1", "NFKB1", "NFKB2", "TNF", "IL15", "ICAM1",
              "COTL1", "TAP1", "CIITA"),
  m2_like = c("MRC1", "C1QA", "C1QB", "C1QC", "TGFB1", "TGFBI",
              "CD14", "TREM2", "INHBA", "VEGFA", "MMP11"),
  tam     = c("HAVCR2", "CD274", "LGALS9", "ADAM17", "ADAM10",
              "CTSS", "ADAMDEC1", "FCGR3A", "SLAMF7", "CCR1")
)

tc_sets <- list(
  cytotoxic  = c("GZMA", "GZMB", "GZMH", "PRF1", "FGFBP2", "CST7",
                 "NKG7", "GNLY", "FASLG", "KLRD1", "KLRB1"),
  treg       = c("BATF", "TIGIT", "TNFRSF18", "IL2RA", "CTLA4",
                 "IL10", "STAT5A", "TGFB1"),
  exhaustion = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4",
                 "BTLA", "CD274", "TOX"),
  naive      = c("IL7R", "SELL", "TCF7", "LTB", "CCR7", "LEF1"),
  activation = c("CD27", "CD28", "ICOS", "CD40LG", "IFNG",
                 "IL2RA", "IL2RG", "TNF", "LCK", "ZAP70")
)

nk_sets <- list(
  cd56dim        = c("FCGR3A", "FGFBP2", "CST7", "KLRF1", "PRF1",
                     "GZMA", "GZMB", "GZMH", "CX3CR1"),
  cd56bright     = c("IGFBP2", "KRT81", "SELL", "NCR1", "CD7",
                     "KLRC1", "KLRC2"),
  nk_activation  = c("NKG7", "GNLY", "HCST", "KLRD1", "CD69",
                      "IFNG", "TNF", "FASLG", "SH2D1A"),
  ifn_responding = c("IFI44L", "IFI6", "ISG20", "MX1", "IFIT1",
                     "IFIT3", "RSAD2", "OAS1"),
  nk_exhaustion  = c("TIGIT", "HAVCR2", "LAG3", "PDCD1",
                     "LILRB1", "LILRB2")
)

bc_sets <- list(
  bc_naive       = c("IGHD", "TCL1A", "IQCG", "SELL", "IL7R"),
  bc_activated   = c("TNFRSF13B", "CD27", "CD40", "CD69", "CD80", "CD86"),
  antigen_pres   = c("CIITA", "TAP1", "TAP2", "TAPBP", "PSME1", "PSME2",
                     "PDIA3", "CALR", "CANX"),
  proliferating  = c("MKI67", "STMN1", "TOP2A", "CCNB1", "CDK1",
                      "CDC20", "PCNA"),
  plasma_diff    = c("XBP1", "MZB1", "JCHAIN", "DERL3", "SSR4", "FKBP11",
                     "SEC11C")
)

celltype_sets <- list(
  "Macrophage" = mac_sets,
  "T cell"     = tc_sets,
  "NK cell"    = nk_sets,
  "B cell"     = bc_sets
)

# --- SFE names ---------------------------------------------------------------

sfe_names <- c("sfe_tma_filtered", sfe_names_wt)

# --- Load neighborhood assignments ------------------------------------------

nb_assign <- read.csv(file.path(out_dir, "09_neighborhood",
                                "neighborhood_assignments.csv"),
                      stringsAsFactors = FALSE)
nb_assign$niche_name <- nb_names[nb_assign$neighborhood]
message(sprintf("Loaded %s neighborhood assignments",
                format(nrow(nb_assign), big.mark = ",")))


# ============================================================================
# SECTION 1: Cell Composition Analysis (all 6 epi neighborhoods)
# ============================================================================

message("\n--- Section 1: Cell type composition per epithelial neighborhood ---")

# Cross-tabulate cell_label × neighborhood from assignments
# Need to join cell_label from SFEs to nb_assign
comp_list <- list()

for (sname in sfe_names) {
  message("  Loading ", sname, " for composition ...")
  sfe <- load_sfe(sname)
  cd <- data.frame(
    cell_id    = colnames(sfe),
    cell_label = colData(sfe)$cell_label,
    stringsAsFactors = FALSE
  )
  # Rename-mismatch fix (idempotent): deposited SFEs still carry the legacy
  # epithelial label. Harmless if absent; keeps composition output/palette
  # consistent with the canonical "Intermediate epithelium" naming.
  cd$cell_label[cd$cell_label == "Transitioning epithelium"] <- "Intermediate epithelium"
  rm(sfe); gc(verbose = FALSE)

  # Join neighborhood
  nb_match <- nb_assign[match(cd$cell_id, nb_assign$cell_id), ]
  cd$neighborhood <- nb_match$neighborhood

  # Keep only cells in epi neighborhoods
  cd <- cd[!is.na(cd$neighborhood) & cd$neighborhood %in% epi_nbs, ]
  if (nrow(cd) > 0) comp_list[[sname]] <- cd
}

comp_df <- do.call(rbind, comp_list)
rownames(comp_df) <- NULL

# Compute counts and proportions
comp_tab <- as.data.frame(table(neighborhood = comp_df$neighborhood,
                                cell_type = comp_df$cell_label))
colnames(comp_tab)[3] <- "count"

# Proportions within each neighborhood
total_per_nb <- tapply(comp_tab$count, comp_tab$neighborhood, sum)
comp_tab$total_in_nb <- total_per_nb[comp_tab$neighborhood]
comp_tab$proportion <- comp_tab$count / comp_tab$total_in_nb

# Add readable names
comp_tab$niche_name <- epi_nb_labels[comp_tab$neighborhood]

write.csv(comp_tab,
          file.path(out_13f, "cell_composition_epi_neighborhoods.csv"),
          row.names = FALSE)

message("  Cell composition saved. Neighborhoods and total cells:")
for (nb in epi_nbs) {
  n <- total_per_nb[nb]
  message(sprintf("    %s (%s): %s cells",
                  nb, epi_nb_labels[nb], format(n, big.mark = ",")))
}

rm(comp_list, comp_df); gc(verbose = FALSE)


# ============================================================================
# SECTION 2: Collect immune cell data for DEG (Method 1: membership)
# ============================================================================

message("\n--- Section 2: Collecting immune cell data (Method 1: membership) ---")

# Accumulate metadata + expression for immune cells in epi neighborhoods
m1_meta_list <- list()
m1_expr_list <- list()

for (sname in sfe_names) {

  message("  Loading ", sname, " ...")
  sfe <- load_sfe(sname)

  cd <- as.data.frame(colData(sfe))
  cd$cell_id <- colnames(sfe)

  # Join neighborhood
  nb_match <- nb_assign[match(cd$cell_id, nb_assign$cell_id), ]
  cd$neighborhood <- nb_match$neighborhood
  cd$niche_name   <- nb_match$niche_name

  # Immune cells in any of the 4 secretory trajectory neighborhoods
  keep <- !is.na(cd$neighborhood) &
          cd$cell_label %in% immune_types &
          cd$neighborhood %in% sec_trajectory

  if (sum(keep) < 10) {
    message("    Skipping: only ", sum(keep), " immune cells in trajectory neighborhoods")
    rm(sfe); gc(verbose = FALSE)
    next
  }

  sfe_sub <- sfe[, keep]
  cd_sub  <- cd[keep, ]
  lc <- as.matrix(logcounts(sfe_sub))

  meta <- data.frame(
    cell_id       = colnames(sfe_sub),
    sample        = sname,
    cell_label    = cd_sub$cell_label,
    neighborhood  = cd_sub$neighborhood,
    niche_name    = cd_sub$niche_name,
    stringsAsFactors = FALSE
  )

  # Pathway scores
  pw_cols <- grep("^pathway_", colnames(cd_sub), value = TRUE)
  if (length(pw_cols) > 0) {
    meta <- cbind(meta, cd_sub[, pw_cols, drop = FALSE])
  }

  # Pseudobulk group_id
  if ("core_id" %in% colnames(cd_sub) && !all(is.na(cd_sub$core_id))) {
    meta$group_id <- paste0(sname, "_core", cd_sub$core_id, "_", cd_sub$neighborhood)
  } else {
    meta$group_id <- paste0(sname, "_", cd_sub$neighborhood)
  }

  m1_meta_list[[sname]] <- meta
  m1_expr_list[[sname]] <- lc

  # Report per-neighborhood counts
  for (nb in sec_trajectory) {
    n_nb <- sum(cd_sub$neighborhood == nb)
    if (n_nb > 0) {
      message(sprintf("    %s in %s: %d immune cells", sname,
                      epi_nb_labels[nb], n_nb))
    }
  }

  rm(sfe, sfe_sub, lc, cd, cd_sub); gc(verbose = FALSE)
}

m1_meta <- do.call(rbind, m1_meta_list)
m1_expr <- do.call(cbind, m1_expr_list)
rownames(m1_meta) <- NULL

message(sprintf("\nMethod 1 total immune cells: %s",
                format(nrow(m1_meta), big.mark = ",")))
for (nb in sec_trajectory) {
  message(sprintf("  %s: %s cells",
                  epi_nb_labels[nb],
                  format(sum(m1_meta$neighborhood == nb), big.mark = ",")))
}


# ============================================================================
# SECTION 2b: Collect immune cell data (Method 2: 13e proximity)
# ============================================================================

message("\n--- Section 2b: Loading 13e proximity-based assignments ---")

rds_path <- file.path(out_dir, "13e_immune_epithelial_interface",
                      "immune_epithelial_distances.rds")

if (file.exists(rds_path)) {
  m2_data <- as.data.frame(readRDS(rds_path))
  message(sprintf("  Loaded %s immune cells from 13e RDS",
                  format(nrow(m2_data), big.mark = ",")))

  # Map nearest_epi (SecA/Trans/SecB) to pseudo-neighborhood for comparison
  # Note: 13e does not distinguish Early Transitioning (nb_10) from Transitioning
  epi_to_nb <- c(SecA = "nb_1", Trans = "nb_7", SecB = "nb_6")
  m2_data$prox_neighborhood <- epi_to_nb[m2_data$nearest_epi]

  message("  Proximity-assigned immune cells per epithelial subtype:")
  for (ep in c("SecA", "Trans", "SecB")) {
    message(sprintf("    %s: %s cells",
                    ep, format(sum(m2_data$nearest_epi == ep), big.mark = ",")))
  }

  m2_available <- TRUE
} else {
  message("  WARNING: 13e RDS not found. Proximity method will be skipped.")
  m2_available <- FALSE
}


# ============================================================================
# SECTION 3: Cell-Matched Pseudobulk DEG Analysis
# ============================================================================

message("\n--- Section 3: Pseudobulk DEG analysis ---")

# Helper function: run pseudobulk DEG for a given cell type + comparison
run_pseudobulk_deg <- function(meta, expr, cell_type, nb_a, nb_b,
                               min_cells_per_group = 20, min_reps = 3) {

  # Subset to cell type + neighborhoods
  keep <- meta$cell_label == cell_type
  if ("neighborhood" %in% colnames(meta)) {
    keep <- keep & meta$neighborhood %in% c(nb_a, nb_b)
  } else if ("prox_neighborhood" %in% colnames(meta)) {
    keep <- keep & meta$prox_neighborhood %in% c(nb_a, nb_b)
  }

  if (sum(keep) < 20) {
    return(data.frame(gene = character(), log2FC = numeric(),
                      mean_a = numeric(), mean_b = numeric(),
                      p_value = numeric(), p_adj = numeric(),
                      sig = logical(), stringsAsFactors = FALSE))
  }

  sub_meta <- meta[keep, ]
  sub_expr <- expr[, keep, drop = FALSE]

  # Determine group assignment
  if ("neighborhood" %in% colnames(sub_meta)) {
    sub_meta$group_nb <- sub_meta$neighborhood
  } else {
    sub_meta$group_nb <- sub_meta$prox_neighborhood
  }

  # Pseudobulk aggregation
  groups <- unique(sub_meta$group_id)
  genes  <- rownames(sub_expr)
  group_sizes <- table(sub_meta$group_id)
  valid_groups <- names(group_sizes[group_sizes >= min_cells_per_group])

  if (length(valid_groups) < 2) return(NULL)

  pb_counts <- matrix(0, nrow = length(genes), ncol = length(valid_groups),
                      dimnames = list(genes, valid_groups))
  for (grp in valid_groups) {
    idx <- which(sub_meta$group_id == grp)
    pb_counts[, grp] <- rowSums(sub_expr[, idx, drop = FALSE])
  }

  pb_meta <- data.frame(
    group_id = valid_groups,
    group_nb = sapply(valid_groups, function(g) sub_meta$group_nb[sub_meta$group_id == g][1]),
    n_cells  = as.integer(group_sizes[valid_groups]),
    stringsAsFactors = FALSE
  )

  # Log2-CPM normalization
  lib_sizes <- colSums(pb_counts)
  pb_expr <- log2(t(t(pb_counts) / lib_sizes) * 1e6 + 1)

  idx_a <- which(pb_meta$group_nb == nb_a)
  idx_b <- which(pb_meta$group_nb == nb_b)

  if (length(idx_a) < min_reps || length(idx_b) < min_reps) return(NULL)

  # Wilcoxon rank-sum per gene
  results <- data.frame(
    gene    = genes,
    log2FC  = NA_real_,
    mean_a  = NA_real_,
    mean_b  = NA_real_,
    p_value = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(genes)) {
    vals_a <- pb_expr[i, idx_a]
    vals_b <- pb_expr[i, idx_b]
    results$mean_a[i] <- mean(vals_a)
    results$mean_b[i] <- mean(vals_b)
    results$log2FC[i] <- mean(vals_a) - mean(vals_b)

    if (length(vals_a) >= 3 && length(vals_b) >= 3) {
      wt <- tryCatch(
        wilcox.test(vals_a, vals_b, exact = FALSE),
        error = function(e) list(p.value = NA_real_)
      )
      results$p_value[i] <- wt$p.value
    }
  }

  results$p_adj <- p.adjust(results$p_value, method = "BH")
  results$sig <- !is.na(results$p_adj) &
                 results$p_adj < 0.05 &
                 abs(results$log2FC) > 0.25
  results <- results[order(-abs(results$log2FC)), ]

  attr(results, "n_reps_a") <- length(idx_a)
  attr(results, "n_reps_b") <- length(idx_b)
  attr(results, "n_cells_a") <- sum(sub_meta$group_nb == nb_a)
  attr(results, "n_cells_b") <- sum(sub_meta$group_nb == nb_b)

  results
}


# --- Run Method 1 DEGs (membership) -----------------------------------------

message("\n  Method 1 (membership) DEG comparisons:")

deg_summary_list <- list()

for (ct in immune_types) {
  for (comp in pairwise_comparisons) {

    label <- paste0(gsub(" ", "", ct), "_", comp$label, "_membership")
    message(sprintf("    %s: %s in %s vs %s",
                    ct, comp$label, epi_nb_labels[comp$a], epi_nb_labels[comp$b]))

    deg <- run_pseudobulk_deg(m1_meta, m1_expr, ct, comp$a, comp$b)

    if (is.null(deg) || nrow(deg) == 0) {
      message("      Insufficient replicates, skipping")
      deg_summary_list[[label]] <- data.frame(
        cell_type = ct, comparison = comp$label, method = "membership",
        n_cells_a = 0, n_cells_b = 0, n_reps_a = 0, n_reps_b = 0,
        n_sig = 0, n_up = 0, n_down = 0,
        stringsAsFactors = FALSE)
      next
    }

    fname <- paste0("deg_", tolower(gsub(" ", "", ct)), "_",
                    comp$label, "_membership.csv")
    write.csv(deg, file.path(out_13f, fname), row.names = FALSE)

    n_sig <- sum(deg$sig, na.rm = TRUE)
    n_up  <- sum(deg$sig & deg$log2FC > 0, na.rm = TRUE)
    n_dn  <- sum(deg$sig & deg$log2FC < 0, na.rm = TRUE)

    message(sprintf("      DEGs: %d total (%d up in %s, %d up in %s), reps: %d vs %d",
                    n_sig, n_up, comp$a, n_dn, comp$b,
                    attr(deg, "n_reps_a"), attr(deg, "n_reps_b")))

    deg_summary_list[[label]] <- data.frame(
      cell_type   = ct,
      comparison  = comp$label,
      method      = "membership",
      n_cells_a   = attr(deg, "n_cells_a"),
      n_cells_b   = attr(deg, "n_cells_b"),
      n_reps_a    = attr(deg, "n_reps_a"),
      n_reps_b    = attr(deg, "n_reps_b"),
      n_sig       = n_sig,
      n_up        = n_up,
      n_down      = n_dn,
      stringsAsFactors = FALSE
    )
  }
}


# --- Run Method 2 DEGs (proximity from 13e) ---------------------------------

if (m2_available) {

  message("\n  Method 2 (proximity) DEG comparisons:")

  # Need to build expression matrix for proximity-assigned cells
  # Re-load SFEs and extract expression for cells in the 13e RDS
  m2_meta_list <- list()
  m2_expr_list <- list()

  for (sname in sfe_names) {
    message("    Loading ", sname, " for proximity DEG ...")
    sfe <- load_sfe(sname)

    cd <- as.data.frame(colData(sfe))
    cd$cell_id <- colnames(sfe)

    # Match cells to 13e data
    m2_match <- m2_data[m2_data$sample_id == sname &
                        m2_data$cell_label %in% immune_types, ]

    if (nrow(m2_match) < 10) {
      rm(sfe); gc(verbose = FALSE)
      next
    }

    # Find matching cells in SFE
    cell_idx <- match(m2_match$cell_id, cd$cell_id)
    cell_idx <- cell_idx[!is.na(cell_idx)]
    m2_match <- m2_match[m2_match$cell_id %in% cd$cell_id[cell_idx], ]

    if (length(cell_idx) < 10) {
      rm(sfe); gc(verbose = FALSE)
      next
    }

    lc <- as.matrix(logcounts(sfe[, cell_idx]))

    meta <- data.frame(
      cell_id            = m2_match$cell_id,
      sample             = sname,
      cell_label         = m2_match$cell_label,
      prox_neighborhood  = m2_match$prox_neighborhood,
      nearest_epi        = m2_match$nearest_epi,
      stringsAsFactors   = FALSE
    )

    # Pseudobulk group_id for proximity method
    if (sname == "sfe_tma_filtered") {
      # Match core_id from original colData
      core_ids <- cd$core_id[cell_idx]
      meta$group_id <- paste0(sname, "_core", core_ids, "_",
                              meta$prox_neighborhood)
    } else {
      meta$group_id <- paste0(sname, "_", meta$prox_neighborhood)
    }

    # Pathway scores from 13e data
    pw_cols <- grep("^pathway_", colnames(m2_match), value = TRUE)
    if (length(pw_cols) > 0) {
      meta <- cbind(meta, m2_match[, pw_cols, drop = FALSE])
    }

    m2_meta_list[[sname]] <- meta
    m2_expr_list[[sname]] <- lc

    rm(sfe, lc, cd); gc(verbose = FALSE)
  }

  m2_meta <- do.call(rbind, m2_meta_list)
  m2_expr <- do.call(cbind, m2_expr_list)
  rownames(m2_meta) <- NULL

  message(sprintf("  Method 2 total immune cells: %s",
                  format(nrow(m2_meta), big.mark = ",")))

  # Proximity comparisons: only 3-way (SecA vs Trans vs SecB)
  # nb_10-specific comparisons not possible
  prox_comparisons <- list(
    list(a = "nb_1", b = "nb_6", label = "SecA_vs_SecB"),
    list(a = "nb_1", b = "nb_7", label = "SecA_vs_Trans"),
    list(a = "nb_7", b = "nb_6", label = "Trans_vs_SecB")
  )

  for (ct in immune_types) {
    for (comp in prox_comparisons) {

      label <- paste0(gsub(" ", "", ct), "_", comp$label, "_proximity")
      message(sprintf("    %s: %s (proximity)", ct, comp$label))

      deg <- run_pseudobulk_deg(m2_meta, m2_expr, ct, comp$a, comp$b)

      if (is.null(deg) || nrow(deg) == 0) {
        message("      Insufficient replicates, skipping")
        deg_summary_list[[label]] <- data.frame(
          cell_type = ct, comparison = comp$label, method = "proximity",
          n_cells_a = 0, n_cells_b = 0, n_reps_a = 0, n_reps_b = 0,
          n_sig = 0, n_up = 0, n_down = 0,
          stringsAsFactors = FALSE)
        next
      }

      fname <- paste0("deg_", tolower(gsub(" ", "", ct)), "_",
                      comp$label, "_proximity.csv")
      write.csv(deg, file.path(out_13f, fname), row.names = FALSE)

      n_sig <- sum(deg$sig, na.rm = TRUE)
      n_up  <- sum(deg$sig & deg$log2FC > 0, na.rm = TRUE)
      n_dn  <- sum(deg$sig & deg$log2FC < 0, na.rm = TRUE)

      message(sprintf("      DEGs: %d total (%d up in %s, %d up in %s)",
                      n_sig, n_up, comp$a, n_dn, comp$b))

      deg_summary_list[[label]] <- data.frame(
        cell_type   = ct,
        comparison  = comp$label,
        method      = "proximity",
        n_cells_a   = attr(deg, "n_cells_a"),
        n_cells_b   = attr(deg, "n_cells_b"),
        n_reps_a    = attr(deg, "n_reps_a"),
        n_reps_b    = attr(deg, "n_reps_b"),
        n_sig       = n_sig,
        n_up        = n_up,
        n_down      = n_dn,
        stringsAsFactors = FALSE
      )
    }
  }
}

# Save DEG summary
deg_summary <- do.call(rbind, deg_summary_list)
rownames(deg_summary) <- NULL
write.csv(deg_summary,
          file.path(out_13f, "deg_summary_all_comparisons.csv"),
          row.names = FALSE)

message("\n  DEG summary saved:")
print(deg_summary[deg_summary$n_sig > 0, ])


# ============================================================================
# SECTION 4: Functional State Scoring
# ============================================================================

message("\n--- Section 4: Functional state scoring ---")

panel_genes <- rownames(m1_expr)

# Helper: compute functional scores for a set of cells
score_functional <- function(expr_mat, cell_meta, cell_type, gene_sets) {
  ct_mask <- cell_meta$cell_label == cell_type
  n_ct <- sum(ct_mask)
  if (n_ct < 10) return(NULL)

  score_list <- list()
  for (sname in names(gene_sets)) {
    gs <- intersect(gene_sets[[sname]], rownames(expr_mat))
    if (length(gs) >= 2) {
      score_list[[sname]] <- colMeans(expr_mat[gs, ct_mask, drop = FALSE])
    } else {
      score_list[[sname]] <- rep(NA_real_, n_ct)
    }
  }
  scores <- as.data.frame(score_list)
  scores$cell_id <- cell_meta$cell_id[ct_mask]
  scores$neighborhood <- if ("neighborhood" %in% colnames(cell_meta)) {
    cell_meta$neighborhood[ct_mask]
  } else {
    cell_meta$prox_neighborhood[ct_mask]
  }
  scores
}

# Helper: pairwise Wilcoxon + Cohen's d for functional scores
test_functional <- function(scores_df, nb_a, nb_b, score_cols) {
  results <- list()
  for (sc in score_cols) {
    vals_a <- scores_df[[sc]][scores_df$neighborhood == nb_a]
    vals_b <- scores_df[[sc]][scores_df$neighborhood == nb_b]
    vals_a <- vals_a[!is.na(vals_a)]
    vals_b <- vals_b[!is.na(vals_b)]

    if (length(vals_a) >= 10 && length(vals_b) >= 10) {
      wt <- tryCatch(
        wilcox.test(vals_a, vals_b, exact = FALSE),
        error = function(e) list(p.value = NA_real_)
      )
      pooled_sd <- sqrt((sd(vals_a)^2 + sd(vals_b)^2) / 2)
      d <- if (pooled_sd > 0) (mean(vals_a) - mean(vals_b)) / pooled_sd else 0

      results[[sc]] <- data.frame(
        score_name = sc,
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
  if (length(results) == 0) return(NULL)
  out <- do.call(rbind, results)
  out$p_adj <- p.adjust(out$p_value, method = "BH")
  out
}

# --- Method 1: functional scores (membership) --------------------------------

func_results_m1 <- list()

for (ct in immune_types) {
  gs <- celltype_sets[[ct]]
  scores <- score_functional(m1_expr, m1_meta, ct, gs)
  if (is.null(scores)) next

  score_cols <- setdiff(colnames(scores), c("cell_id", "neighborhood"))

  for (comp in pairwise_comparisons) {
    res <- test_functional(scores, comp$a, comp$b, score_cols)
    if (!is.null(res)) {
      res$cell_type  <- ct
      res$comparison <- comp$label
      res$method     <- "membership"
      func_results_m1[[paste0(ct, "_", comp$label)]] <- res
    }
  }
}

func_m1_df <- do.call(rbind, func_results_m1)
rownames(func_m1_df) <- NULL

# --- Method 2: functional scores (proximity) ---------------------------------

func_results_m2 <- list()

if (m2_available) {
  for (ct in immune_types) {
    gs <- celltype_sets[[ct]]
    scores <- score_functional(m2_expr, m2_meta, ct, gs)
    if (is.null(scores)) next

    score_cols <- setdiff(colnames(scores), c("cell_id", "neighborhood"))

    for (comp in prox_comparisons) {
      res <- test_functional(scores, comp$a, comp$b, score_cols)
      if (!is.null(res)) {
        res$cell_type  <- ct
        res$comparison <- comp$label
        res$method     <- "proximity"
        func_results_m2[[paste0(ct, "_", comp$label)]] <- res
      }
    }
  }
}

func_m2_df <- if (length(func_results_m2) > 0) {
  do.call(rbind, func_results_m2)
} else {
  data.frame()
}
rownames(func_m2_df) <- NULL

# Combine and save
func_all <- rbind(func_m1_df, func_m2_df)
write.csv(func_all,
          file.path(out_13f, "functional_scores_all_comparisons.csv"),
          row.names = FALSE)

message(sprintf("  Functional scoring: %d tests across %d cell types",
                nrow(func_all), length(unique(func_all$cell_type))))
n_sig_func <- sum(func_all$p_adj < 0.05, na.rm = TRUE)
message(sprintf("  Significant (padj < 0.05): %d / %d", n_sig_func, nrow(func_all)))


# ============================================================================
# SECTION 5: Pathway Enrichment by Epithelial Neighborhood
# ============================================================================

message("\n--- Section 5: Pathway enrichment ---")

pw_cols <- grep("^pathway_", colnames(m1_meta), value = TRUE)

if (length(pw_cols) > 0) {

  # --- 5a: Per cell type, Kruskal-Wallis across all 4 trajectory neighborhoods
  pw_kw_list <- list()
  pw_pairwise_list <- list()

  for (ct in immune_types) {
    ct_mask <- m1_meta$cell_label == ct
    if (sum(ct_mask) < 50) next

    ct_data <- m1_meta[ct_mask, c("neighborhood", pw_cols), drop = FALSE]

    for (pw in pw_cols) {
      vals <- ct_data[[pw]]
      grps <- ct_data$neighborhood

      # Remove NAs and neighborhoods with too few cells
      valid <- !is.na(vals)
      vals <- vals[valid]
      grps <- grps[valid]

      grp_n <- table(grps)
      keep_grps <- names(grp_n[grp_n >= 10])
      keep <- grps %in% keep_grps
      vals <- vals[keep]
      grps <- grps[keep]

      if (length(unique(grps)) < 2) next

      # Kruskal-Wallis
      kw <- kruskal.test(vals ~ grps)
      pw_kw_list[[paste0(ct, "_", pw)]] <- data.frame(
        cell_type = ct,
        pathway   = sub("^pathway_", "", pw),
        kw_stat   = kw$statistic,
        kw_p      = kw$p.value,
        stringsAsFactors = FALSE
      )

      # Pairwise Wilcoxon for trajectory comparisons
      for (comp in pairwise_comparisons) {
        vals_a <- vals[grps == comp$a]
        vals_b <- vals[grps == comp$b]
        if (length(vals_a) >= 10 && length(vals_b) >= 10) {
          wt <- tryCatch(
            wilcox.test(vals_a, vals_b, exact = FALSE),
            error = function(e) list(p.value = NA_real_)
          )
          pooled_sd <- sqrt((sd(vals_a)^2 + sd(vals_b)^2) / 2)
          d <- if (pooled_sd > 0) (mean(vals_a) - mean(vals_b)) / pooled_sd else 0

          pw_pairwise_list[[length(pw_pairwise_list) + 1]] <- data.frame(
            cell_type  = ct,
            pathway    = sub("^pathway_", "", pw),
            comparison = comp$label,
            mean_a     = mean(vals_a),
            mean_b     = mean(vals_b),
            diff       = mean(vals_a) - mean(vals_b),
            cohens_d   = d,
            p_value    = wt$p.value,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  pw_kw_df <- do.call(rbind, pw_kw_list)
  pw_kw_df$kw_p_adj <- p.adjust(pw_kw_df$kw_p, method = "BH")

  pw_pair_df <- do.call(rbind, pw_pairwise_list)
  pw_pair_df$p_adj <- p.adjust(pw_pair_df$p_value, method = "BH")

  write.csv(pw_kw_df,
            file.path(out_13f, "pathway_kruskalwallis_epi_neighborhoods.csv"),
            row.names = FALSE)
  write.csv(pw_pair_df,
            file.path(out_13f, "pathway_pairwise_epi_neighborhoods.csv"),
            row.names = FALSE)

  n_kw_sig <- sum(pw_kw_df$kw_p_adj < 0.05, na.rm = TRUE)
  n_pw_sig <- sum(pw_pair_df$p_adj < 0.05, na.rm = TRUE)
  message(sprintf("  KW tests significant: %d / %d", n_kw_sig, nrow(pw_kw_df)))
  message(sprintf("  Pairwise tests significant: %d / %d", n_pw_sig, nrow(pw_pair_df)))

} else {
  message("  No pathway scores found. Run 9b_pathway_scoring.R first.")
  pw_kw_df <- data.frame()
  pw_pair_df <- data.frame()
}


# ============================================================================
# SECTION 6: Import 13e Nearby Cell Composition
# ============================================================================

message("\n--- Section 6: Importing 13e nearby cell composition ---")

f_comp <- file.path(out_dir, "13e_immune_epithelial_interface",
                    "immune_composition_by_epi_subtype.csv")
f_enrich <- file.path(out_dir, "13e_immune_epithelial_interface",
                      "immune_enrichment_by_epi_subtype.csv")

if (file.exists(f_comp)) {
  nearby_comp <- read.csv(f_comp, stringsAsFactors = FALSE)
  write.csv(nearby_comp,
            file.path(out_13f, "nearby_immune_composition_from_13e.csv"),
            row.names = FALSE)
  message(sprintf("  Imported immune composition: %d rows", nrow(nearby_comp)))
} else {
  message("  WARNING: 13e composition file not found.")
  nearby_comp <- data.frame()
}

if (file.exists(f_enrich)) {
  nearby_enrich <- read.csv(f_enrich, stringsAsFactors = FALSE)
  write.csv(nearby_enrich,
            file.path(out_13f, "nearby_immune_enrichment_from_13e.csv"),
            row.names = FALSE)
  message(sprintf("  Imported immune enrichment: %d rows", nrow(nearby_enrich)))
} else {
  nearby_enrich <- data.frame()
}


# ============================================================================
# SECTION 7: Figures
# ============================================================================

options(ggrastr.default.dpi = 300)

message("\n--- Section 7: Generating figures ---")

# --- 7a. Stacked bar: cell type composition per epi neighborhood --------------

comp_tab <- read.csv(file.path(out_13f, "cell_composition_epi_neighborhoods.csv"),
                     stringsAsFactors = FALSE)

# Order neighborhoods
comp_tab$niche_name <- factor(comp_tab$niche_name,
                              levels = epi_nb_labels[epi_nbs])

p_comp <- ggplot(comp_tab[comp_tab$proportion > 0.005, ],
                 aes(x = niche_name, y = proportion, fill = cell_type)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = ref_palette, name = "Cell type") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(title = "Cell type composition of epithelial neighborhoods",
       x = NULL, y = "Proportion") +
  theme_lab(base_size = 8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.key.size = unit(0.5, "lines"))

ggsave(file.path(fig_13f, "composition_stacked_bar.pdf"),
       p_comp, width = 9, height = 6)


# --- 7b. DEG summary tile plot ------------------------------------------------

if (nrow(deg_summary) > 0 && any(deg_summary$n_sig > 0)) {

  deg_plot <- deg_summary[deg_summary$n_sig > 0, ]
  deg_plot$total_deg <- deg_plot$n_up + deg_plot$n_down
  deg_plot$label <- paste0(deg_plot$n_up, "/", deg_plot$n_down)

  p_deg_tile <- ggplot(deg_plot,
                       aes(x = comparison, y = cell_type,
                           fill = total_deg)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = label), size = 2.5) +
    scale_fill_gradient(low = "lightyellow", high = "#E6A141",
                        name = "Total DEGs") +
    facet_wrap(~ method, ncol = 2) +
    labs(title = "DEG counts per immune cell type and comparison",
         subtitle = "Label: up / down in neighborhood A",
         x = NULL, y = NULL) +
    theme_lab(base_size = 8) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  ggsave(file.path(fig_13f, "deg_summary_tile.pdf"),
         p_deg_tile, width = 10, height = 5)
}


# --- 7c. Volcano: SecA vs SecB (most biologically interesting) ----------------

for (ct in immune_types) {
  ct_short <- tolower(gsub(" ", "", ct))
  deg_file <- file.path(out_13f,
                        paste0("deg_", ct_short, "_SecA_vs_SecB_membership.csv"))

  if (!file.exists(deg_file)) next
  deg <- read.csv(deg_file, stringsAsFactors = FALSE)
  if (nrow(deg) == 0) next

  deg$neg_log10_p <- -log10(pmax(deg$p_adj, 1e-300))
  deg$neg_log10_p[deg$neg_log10_p > 50] <- 50
  deg$category <- ifelse(!deg$sig, "NS",
                   ifelse(deg$log2FC > 0, "Up in SecA-enriched",
                          "Up in SecB-enriched"))

  top_label <- head(deg[deg$sig == TRUE, ], 25)

  p_vol <- ggplot(deg, aes(x = log2FC, y = neg_log10_p, color = category)) +
    geom_point(size = 1.2, alpha = 0.7) +
    scale_color_manual(values = c("NS" = "grey70",
                                  "Up in SecA-enriched" = "#F0C060",
                                  "Up in SecB-enriched" = "#9A7D55"),
                       name = "") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    ggrepel::geom_text_repel(
      data = top_label,
      aes(label = gene, y = pmin(-log10(p_adj), 50)),
      size = 2.2, max.overlaps = 25, color = "black", segment.color = "grey50"
    ) +
    labs(title = sprintf("%s DEGs: SecA-enriched vs SecB-enriched epithelium", ct),
         subtitle = sprintf("Membership method | %d significant",
                            sum(deg$sig, na.rm = TRUE)),
         x = expression(log[2]~fold~change~"(SecA / SecB)"),
         y = expression(-log[10]~adjusted~p)) +
    theme_lab(base_size = 8) +
    theme(legend.position = "top")

  ggsave(file.path(fig_13f, paste0("volcano_", ct_short, "_SecA_vs_SecB.pdf")),
         p_vol, width = 7, height = 6)
}


# --- 7d. Functional score heatmap across trajectory --------------------------

if (nrow(func_all) > 0) {

  # Use membership results for the trajectory heatmap
  func_m1_traj <- func_all[func_all$method == "membership", ]

  if (nrow(func_m1_traj) > 0) {

    # Build matrix: rows = cell_type:score, columns = comparison
    func_m1_traj$row_label <- paste0(func_m1_traj$cell_type, ": ",
                                     func_m1_traj$score_name)

    # Use Cohen's d as the heatmap value
    func_wide <- reshape(func_m1_traj[, c("row_label", "comparison", "cohens_d")],
                         idvar = "row_label", timevar = "comparison",
                         direction = "wide")
    rownames(func_wide) <- func_wide$row_label
    func_wide$row_label <- NULL
    colnames(func_wide) <- gsub("^cohens_d\\.", "", colnames(func_wide))

    # Order columns by trajectory
    traj_order <- c("SecA_vs_EarlyTrans", "EarlyTrans_vs_Trans",
                    "Trans_vs_LateTrans", "LateTrans_vs_SecB", "SecA_vs_SecB")
    traj_order <- intersect(traj_order, colnames(func_wide))
    func_mat <- as.matrix(func_wide[, traj_order, drop = FALSE])

    # Cap extreme values
    func_mat[func_mat > 1] <- 1
    func_mat[func_mat < -1] <- -1

    pdf(file.path(fig_13f, "functional_scores_trajectory_heatmap.pdf"),
        width = 8, height = max(6, nrow(func_mat) * 0.25 + 2))
    ht <- Heatmap(func_mat,
            name = "Cohen's d",
            col = colorRamp2(c(-1, 0, 1), c("#9A7D55", "white", "#F0C060")),
            column_title = "Functional state shifts along secretory trajectory",
            column_title_gp = gpar(fontsize = 10),
            row_names_gp = gpar(fontsize = 7),
            column_names_gp = gpar(fontsize = 8),
            column_names_rot = 35,
            cluster_columns = FALSE,
            clustering_method_rows = "ward.D2",
            rect_gp = gpar(col = "white", lwd = 0.5),
            cell_fun = function(j, i, x, y, width, height, fill) {
              val <- func_mat[i, j]
              if (!is.na(val) && abs(val) > 0.3) {
                grid.text(sprintf("%.2f", val), x, y,
                          gp = gpar(fontsize = 5, col = "black"))
              }
            })
    draw(ht)
    dev.off()
  }
}


# --- 7e. Pathway dotplot: per immune type across neighborhoods ----------------

if (nrow(pw_pair_df) > 0) {

  pw_plot <- pw_pair_df[pw_pair_df$comparison == "SecA_vs_SecB", ]

  if (nrow(pw_plot) > 0) {
    p_pw <- ggplot(pw_plot,
                   aes(x = cohens_d,
                       y = reorder(pathway, cohens_d),
                       size = -log10(pmax(p_adj, 1e-300)),
                       color = cohens_d)) +
      geom_point() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      scale_color_gradient2(low = "#9A7D55", mid = "grey80", high = "#F0C060",
                            midpoint = 0, name = "Cohen's d") +
      scale_size_continuous(name = expression(-log[10]~p[adj]), range = c(1, 5)) +
      facet_wrap(~ cell_type, ncol = 2) +
      labs(title = "Pathway activity: immune cells in SecA vs SecB epithelium",
           subtitle = "Membership method | Positive = higher in SecA-enriched",
           x = "Cohen's d", y = NULL) +
      theme_lab(base_size = 7) +
      theme(strip.text = element_text(size = rel(1.1)))

    ggsave(file.path(fig_13f, "pathway_dotplot_SecA_vs_SecB.pdf"),
           p_pw, width = 10, height = 8)
  }
}


# ============================================================================
# SECTION 8: Summary
# ============================================================================

message("\n--- Section 8: Summary ---")

summary_rows <- list()

# Composition summary
total_cells <- sum(comp_tab$count)
summary_rows[["Total cells in epi neighborhoods"]] <- total_cells

# DEG summary
if (nrow(deg_summary) > 0) {
  summary_rows[["Method 1 comparisons with DEGs"]] <-
    sum(deg_summary$n_sig > 0 & deg_summary$method == "membership")
  summary_rows[["Method 2 comparisons with DEGs"]] <-
    sum(deg_summary$n_sig > 0 & deg_summary$method == "proximity")
  summary_rows[["Total significant DEGs (all comparisons)"]] <-
    sum(deg_summary$n_sig)
}

# Functional scoring summary
if (nrow(func_all) > 0) {
  summary_rows[["Functional score tests"]] <- nrow(func_all)
  summary_rows[["Significant functional scores (padj<0.05)"]] <-
    sum(func_all$p_adj < 0.05, na.rm = TRUE)
}

# Pathway summary
if (nrow(pw_pair_df) > 0) {
  summary_rows[["Pathway pairwise tests"]] <- nrow(pw_pair_df)
  summary_rows[["Significant pathways (padj<0.05)"]] <-
    sum(pw_pair_df$p_adj < 0.05, na.rm = TRUE)
}

summary_df <- data.frame(
  metric = names(summary_rows),
  value  = unlist(summary_rows),
  stringsAsFactors = FALSE
)

write.csv(summary_df,
          file.path(out_13f, "analysis_summary.csv"),
          row.names = FALSE)

print(summary_df)


# ============================================================================
# Done
# ============================================================================

message("\n=== Script 13f Complete ===")
message("Outputs saved to: ", out_13f)
message("Figures saved to: ", fig_13f)
message("\nKey files generated:")
message("  Composition:")
message("    - cell_composition_epi_neighborhoods.csv")
message("  DEG (per cell type × comparison × method):")
message("    - deg_*_membership.csv / deg_*_proximity.csv")
message("    - deg_summary_all_comparisons.csv")
message("  Functional scoring:")
message("    - functional_scores_all_comparisons.csv")
message("  Pathway enrichment:")
message("    - pathway_kruskalwallis_epi_neighborhoods.csv")
message("    - pathway_pairwise_epi_neighborhoods.csv")
message("  13e integration:")
message("    - nearby_immune_composition_from_13e.csv")
message("    - nearby_immune_enrichment_from_13e.csv")
message("  Summary:")
message("    - analysis_summary.csv")
message("  Figures:")
message("    - composition_stacked_bar.pdf")
message("    - deg_summary_tile.pdf")
message("    - volcano_*_SecA_vs_SecB.pdf (per immune type)")
message("    - functional_scores_trajectory_heatmap.pdf")
message("    - pathway_dotplot_SecA_vs_SecB.pdf")

log_session()
