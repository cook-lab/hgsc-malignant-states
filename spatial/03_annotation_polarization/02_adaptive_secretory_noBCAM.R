# ============================================================================
# 02_adaptive_secretory_noBCAM.R — Initial secretory subtyping (noBCAM)
# ============================================================================
# PURPOSE: Assign an initial SecA / Intermediate / SecB subtype to SingleR
#   "Secretory epithelium" cells from a SecB/SecA logcounts-ratio, using the
#   shared noBCAM 7-gene signatures and atlas-aligned thresholds
#   (lower=1.0, upper=2.0). Writes secretory_subtype, cell_label, and per-cell
#   scores into each SFE; this initial cell_label is later refined by the
#   UCell-based 04_reclassification_polarization.R and 05_clean_split_rctd.R.
#
# COHORT PIN: operates over sfe_tma + the 8 published whole tissues.
#
# SIGNATURES: SecA/SecB loaded from shared/signatures.yml (noBCAM 7-gene set)
#   via 00_setup.R (SECA_GENES / SECB_GENES). BCAM is intentionally absent.
#
# NAMING: epithelial label standardized to "Intermediate epithelium"
#   (was "Transitioning epithelium").
#
# INPUTS:
#   - <sfe_dir>/sfe_tma, <sfe_dir>/sfe_<wt>   (annotated SFEs, singler_label)
#
# OUTPUTS:
#   - SFEs updated: secretory_subtype, cell_label, score_SecA, score_SecB,
#       ratio_SecB_SecA
#   - <output_root>/06b_adaptive_secretory_noBCAM/production_summary.csv
#   - <output_root>/06b_adaptive_secretory_noBCAM/reclassification_comparison.csv
#
# MANUSCRIPT PANEL(S): polarization-chain backend (Fig 4–6 composition/ROI).
#
# RUNTIME TIER: moderate
# ============================================================================

source("spatial/00_setup/00_setup.R")

# --- Parameters ---------------------------------------------------------------

secB_genes <- SECB_GENES   # from shared/signatures.yml
secA_genes <- SECA_GENES   # noBCAM 7-gene set from shared/signatures.yml

lower_threshold <- 1.0
upper_threshold <- 2.0
epsilon         <- 0.01

subtype_levels <- c("SecA epithelium",
                     "Intermediate epithelium",
                     "SecB epithelium")

message("\n=== Secretory Epithelium Subtyping (noBCAM, 1.0/2.0) ===")
message("SecB genes (", length(secB_genes), "): ", paste(secB_genes, collapse = ", "))
message("SecA genes (", length(secA_genes), "): ", paste(secA_genes, collapse = ", "))
message("Thresholds: lower = ", lower_threshold, ", upper = ", upper_threshold)

# --- Output directory ---------------------------------------------------------

out_sub <- file.path(out_dir, "06b_adaptive_secretory_noBCAM")
dir.create(out_sub, showWarnings = FALSE, recursive = TRUE)

# --- SFE names ----------------------------------------------------------------

sfe_names <- sfe_names_all   # sfe_tma + published 8 whole tissues (cohort PIN)

# --- Score, classify, and save ------------------------------------------------

summary_list    <- list()
comparison_list <- list()

for (sname in sfe_names) {

  message("\nProcessing ", sname, " ...")
  sfe <- load_sfe(sname)
  n_total <- ncol(sfe)

  is_sec <- sfe$singler_label == "Secretory epithelium"
  is_sec[is.na(is_sec)] <- FALSE
  n_sec <- sum(is_sec)

  # Archive the PRISTINE pre-relabel columns exactly once. Guard on the *_old
  # column already existing so an accidental re-run does not clobber the saved
  # originals with values this script has already overwritten (idempotency).
  if ("secretory_subtype" %in% colnames(colData(sfe)) &&
      !"secretory_subtype_old" %in% colnames(colData(sfe))) {
    sfe$secretory_subtype_old <- sfe$secretory_subtype
  }
  if ("cell_label" %in% colnames(colData(sfe)) &&
      !"cell_label_old" %in% colnames(colData(sfe))) {
    sfe$cell_label_old <- sfe$cell_label
  }

  sfe$score_SecA        <- NA_real_
  sfe$score_SecB        <- NA_real_
  sfe$ratio_SecB_SecA   <- NA_real_
  sfe$secretory_subtype <- NA_character_

  if (n_sec > 0) {
    lc <- as.matrix(logcounts(sfe[c(secB_genes, secA_genes), is_sec]))

    s_SecB <- colMeans(lc[secB_genes, , drop = FALSE])
    s_SecA <- colMeans(lc[secA_genes, , drop = FALSE])
    ratio  <- (s_SecB + epsilon) / (s_SecA + epsilon)

    subtype <- ifelse(ratio < lower_threshold, subtype_levels[1],
               ifelse(ratio <= upper_threshold, subtype_levels[2],
                      subtype_levels[3]))

    sfe$score_SecA[is_sec]        <- s_SecA
    sfe$score_SecB[is_sec]        <- s_SecB
    sfe$ratio_SecB_SecA[is_sec]   <- ratio
    sfe$secretory_subtype[is_sec] <- subtype

    # Build cell_label: subtypes for secretory, singler_label for all others
    sfe$cell_label <- as.character(sfe$singler_label)
    sfe$cell_label[is_sec] <- subtype

    tab <- table(subtype)
    for (lev in subtype_levels) {
      n_lev <- ifelse(lev %in% names(tab), tab[[lev]], 0)
      message("  ", lev, ": ", format(n_lev, big.mark = ","),
              " (", round(100 * n_lev / n_sec, 1), "%)")
    }

    if ("secretory_subtype_old" %in% colnames(colData(sfe))) {
      old_sub <- sfe$secretory_subtype_old[is_sec]
      new_sub <- subtype
      n_changed <- sum(old_sub != new_sub, na.rm = TRUE)
      pct_changed <- round(100 * n_changed / n_sec, 1)
      message("  Reclassified: ", format(n_changed, big.mark = ","),
              " (", pct_changed, "%)")

      conf <- as.data.frame(table(old = old_sub, new = new_sub))
      conf$sample_id <- sname
      comparison_list[[sname]] <- conf
    }

    summary_list[[sname]] <- data.table(
      sample_id       = sname,
      n_total         = n_total,
      n_secretory     = n_sec,
      n_secA          = sum(subtype == subtype_levels[1]),
      n_intermediate  = sum(subtype == subtype_levels[2]),
      n_secB          = sum(subtype == subtype_levels[3]),
      pct_secA         = round(100 * sum(subtype == subtype_levels[1]) / n_sec, 1),
      pct_intermediate = round(100 * sum(subtype == subtype_levels[2]) / n_sec, 1),
      pct_secB         = round(100 * sum(subtype == subtype_levels[3]) / n_sec, 1)
    )
  } else {
    sfe$cell_label <- as.character(sfe$singler_label)
    message("  No secretory cells — skipping subtyping")
    summary_list[[sname]] <- data.table(
      sample_id = sname, n_total = n_total, n_secretory = 0,
      n_secA = 0, n_intermediate = 0, n_secB = 0,
      pct_secA = 0, pct_intermediate = 0, pct_secB = 0)
  }

  # Remove old lasso-based columns if present
  if ("adaptive_label" %in% colnames(colData(sfe))) {
    colData(sfe)$adaptive_label <- NULL
  }
  if ("adaptive_score" %in% colnames(colData(sfe))) {
    colData(sfe)$adaptive_score <- NULL
  }

  message("  Realizing assays...")
  for (a in assayNames(sfe)) {
    assay(sfe, a) <- as(assay(sfe, a), "dgCMatrix")
  }

  message("  Saving...")
  save_sfe(sfe, sname)
  rm(sfe, lc); gc(verbose = FALSE)
}

# --- Summary ------------------------------------------------------------------

summary_df <- rbindlist(summary_list)
write.csv(summary_df, file.path(out_sub, "production_summary.csv"),
          row.names = FALSE)

if (length(comparison_list) > 0) {
  comparison_df <- do.call(rbind, comparison_list)
  write.csv(comparison_df, file.path(out_sub, "reclassification_comparison.csv"),
            row.names = FALSE)
}

message("\n=== Secretory subtyping complete (noBCAM, 1.0/2.0) ===")
message("Total cells:       ", format(sum(summary_df$n_total), big.mark = ","))
message("Total secretory:   ", format(sum(summary_df$n_secretory), big.mark = ","))
message("  SecA:            ", format(sum(summary_df$n_secA), big.mark = ","),
        " (", round(100 * sum(summary_df$n_secA) / sum(summary_df$n_secretory), 1), "%)")
message("  Intermediate:    ", format(sum(summary_df$n_intermediate), big.mark = ","),
        " (", round(100 * sum(summary_df$n_intermediate) / sum(summary_df$n_secretory), 1), "%)")
message("  SecB:            ", format(sum(summary_df$n_secB), big.mark = ","),
        " (", round(100 * sum(summary_df$n_secB) / sum(summary_df$n_secretory), 1), "%)")
message("\nAll SFEs updated with: secretory_subtype, cell_label, score_SecA, score_SecB, ratio_SecB_SecA")
message("Old classification preserved in: secretory_subtype_old, cell_label_old")
message("Summary saved to: ", file.path(out_sub, "production_summary.csv"))

log_session()
