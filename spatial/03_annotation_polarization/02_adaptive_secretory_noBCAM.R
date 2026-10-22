# ============================================================================
# 02_adaptive_secretory_noBCAM.R — Initial SFE cell_label (ratio-based secretory subtypes)
# ============================================================================
# PURPOSE: Write the SFE column cell_label for the first time: the SingleR
#   label, with SingleR "Secretory epithelium" cells given an initial SecA /
#   Intermediate / SecB subtype from the ratio of mean SecB to mean SecA
#   logcounts (ratio < 1.0 SecA, 1.0–2.0 Intermediate, > 2.0 SecB). Also
#   writes secretory_subtype, score_SecA, score_SecB and ratio_SecB_SecA.
#
#   No manuscript panel uses these labels; the manuscript's SecA /
#   Intermediate / SecB calls are the UCell polarization split of the SingleR
#   label tables (03, 04, 06_extract_singler.R). The step is kept because
#   later steps read the column it creates: 04_reclassification_polarization.R
#   re-labels the cells this step marks as secretory (by polarization),
#   05_clean_split_rctd.R overwrites cell_label again, and
#   spatial/02_qc/05b_core_qc.R reads the column.
#
# SAMPLES: every SFE in sfe_names_all.
#
# SIGNATURES: the 7 + 7 SecA / SecB genes of shared/signatures.yml (SECA_GENES /
#   SECB_GENES from 00_setup.R). "noBCAM" in the file name refers to this SecA
#   set, which does not include BCAM.
#
# INPUTS:
#   - <sfe_dir>/sfe_*   (SFEs with singler_label, from 01_annotation_singler.R)
#
# OUTPUTS:
#   - SFEs updated: cell_label, secretory_subtype, score_SecA, score_SecB,
#       ratio_SecB_SecA (on a re-run, the previous secretory_subtype and
#       cell_label are kept once as secretory_subtype_old / cell_label_old)
#
# MANUSCRIPT PANEL(S): none directly (structural step; see PURPOSE).
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

# --- SFE names ----------------------------------------------------------------

sfe_names <- sfe_names_all   # every SFE: TMA A, the cohort.tma_new arrays, WT01–WT08

# --- Score, classify, and save ------------------------------------------------

for (sname in sfe_names) {

  message("\nProcessing ", sname, " ...")
  sfe <- load_sfe(sname)

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
  } else {
    sfe$cell_label <- as.character(sfe$singler_label)
    message("  No secretory cells — skipping subtyping")
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

message("\n=== Secretory subtyping complete (noBCAM, 1.0/2.0) ===")
message("\nAll SFEs updated with: secretory_subtype, cell_label, score_SecA, score_SecB, ratio_SecB_SecA")
message("Old classification preserved in: secretory_subtype_old, cell_label_old")

log_session()
