# ============================================================================
# 03_compute_niche_metabolic_scores.R
# ----------------------------------------------------------------------------
# PURPOSE: Compute per-TME-cell 50um niche metabolic stress (mean epithelial hypoxia+glycolysis within 50um; dbscan frNN), Z-scored within sample + binned.
#
# INPUTS:
#   - 8 WT + sfe_tma_filtered (load_sfe), pathway_hypoxia/glycolysis (9b), 06f override
#
# OUTPUTS:
#   - output/29_macrophage_niche_survival/per_cell_niche_scores.rds
#   - niche_scores_summary_per_sample.csv
#
# MANUSCRIPT PANEL(S): Fig 6B/6E/6F (upstream)
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

suppressPackageStartupMessages({
  library(dbscan); library(data.table)
})

OUT_DIR <- file.path(out_dir, "29_macrophage_niche_survival")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

RADIUS_UM <- 50   # 19e cutoff philosophy

WT_SAMPLES <- sfe_names_wt

EPI_LBLS <- c("SecA epithelium", "Intermediate epithelium", "SecB epithelium")

# ---------------------------------------------------------------------------
# 06f override (same pattern used throughout the project)
# ---------------------------------------------------------------------------
f06f <- file.path(out_dir, "06f_reclassification_polarization",
                  "reclassified_xenium_scores.csv")
stopifnot(file.exists(f06f))
recl <- fread(f06f, select = c("sample", "barcode_orig", "cell_label_06f"))

override_with_06f <- function(sfe, sample_key) {
  sub <- recl[sample == sample_key]
  if (nrow(sub) == 0) return(sfe)
  m <- match(colnames(sfe), sub$barcode_orig)
  hit <- !is.na(m)
  lab <- as.character(sfe$cell_label)
  lab[hit] <- sub$cell_label_06f[m[hit]]
  # Idempotent rename: SFE/06f sources still carry the legacy epithelial label
  # "Transitioning epithelium"; standardize to "Intermediate epithelium" so the
  # downstream EPI_LBLS match (is_epi) does not silently drop the Intermediate epitype.
  lab[lab == "Transitioning epithelium"] <- "Intermediate epithelium"
  sfe$cell_label <- lab
  sfe
}

# ---------------------------------------------------------------------------
# Core per-group compute: for one spatial group (sample for WT, core for TMA),
# find all epi cells within RADIUS_UM of each TME cell and compute the mean
# pathway_hypoxia and pathway_glycolysis.
# ---------------------------------------------------------------------------
compute_niche_for_group <- function(coords, is_epi, is_tme,
                                     pw_hypoxia, pw_glyc) {
  # frNN: for each query point, indices of points within eps
  # Build on epi coordinates, query from TME coordinates
  n_tme <- sum(is_tme)
  if (n_tme == 0 || sum(is_epi) == 0) {
    return(list(niche_hypoxia    = rep(NA_real_, n_tme),
                niche_glycolysis = rep(NA_real_, n_tme),
                n_epi            = rep(0L,        n_tme)))
  }
  epi_coords <- coords[is_epi, , drop = FALSE]
  tme_coords <- coords[is_tme, , drop = FALSE]
  epi_hypo   <- pw_hypoxia[is_epi]
  epi_glyc   <- pw_glyc[is_epi]

  nn <- dbscan::frNN(x = epi_coords, eps = RADIUS_UM, query = tme_coords,
                     sort = FALSE)
  niche_h <- rep(NA_real_, n_tme)
  niche_g <- rep(NA_real_, n_tme)
  n_epi   <- integer(n_tme)
  for (i in seq_len(n_tme)) {
    idx <- nn$id[[i]]
    if (length(idx) > 0) {
      n_epi[i]   <- length(idx)
      niche_h[i] <- mean(epi_hypo[idx], na.rm = TRUE)
      niche_g[i] <- mean(epi_glyc[idx], na.rm = TRUE)
    }
  }
  list(niche_hypoxia = niche_h, niche_glycolysis = niche_g,
       n_epi         = n_epi)
}

# ---------------------------------------------------------------------------
# Process one SFE: applies 06f override, loops groups, returns data.table of
# TME cells with niche scores.
# ---------------------------------------------------------------------------
process_sfe <- function(sfe, cohort, sample_key, group_col = NULL) {
  coords <- SpatialFeatureExperiment::spatialCoords(sfe)
  lbl    <- as.character(sfe$cell_label)
  is_epi <- lbl %in% EPI_LBLS
  is_tme <- !is_epi & lbl != "" & !is.na(lbl)
  pw_h   <- as.numeric(sfe$pathway_hypoxia)
  pw_g   <- as.numeric(sfe$pathway_glycolysis)

  # Group: for WT one group = whole sample; for TMA one group = core_id
  if (is.null(group_col)) {
    group_vec <- rep("all", ncol(sfe))
  } else {
    group_vec <- as.character(colData(sfe)[[group_col]])
  }
  groups <- unique(na.omit(group_vec[group_vec != ""]))

  out_list <- list()
  for (grp in groups) {
    in_grp <- group_vec == grp & !is.na(group_vec)
    epi_g  <- is_epi & in_grp
    tme_g  <- is_tme & in_grp
    if (sum(tme_g) == 0) next
    res <- compute_niche_for_group(coords, epi_g, tme_g, pw_h, pw_g)
    idx_tme <- which(tme_g)
    dt <- data.table(
      cohort       = cohort,
      sample_key   = sample_key,
      group_id     = grp,
      cell_id      = colnames(sfe)[idx_tme],
      cell_label   = lbl[idx_tme],
      x            = coords[idx_tme, 1],
      y            = coords[idx_tme, 2],
      n_epi_within_50um = res$n_epi,
      niche_hypoxia     = res$niche_hypoxia,
      niche_glycolysis  = res$niche_glycolysis
    )
    # Pull patient_id if it exists (TMA)
    if ("patient_id" %in% colnames(colData(sfe))) {
      dt$patient_id <- as.character(colData(sfe)$patient_id[idx_tme])
    } else {
      dt$patient_id <- NA_character_
    }
    out_list[[grp]] <- dt
  }
  rbindlist(out_list)
}

# ---------------------------------------------------------------------------
# Run on 8 WT samples
# ---------------------------------------------------------------------------
wt_list <- list()
for (s in WT_SAMPLES) {
  message("== WT ", s, " ==")
  sfe <- load_sfe(s)
  sfe <- override_with_06f(sfe, s)
  # WT: treat whole SFE as one group (one sample per slide)
  dt <- process_sfe(sfe, cohort = "WT",
                    sample_key = sub("^sfe_", "", s),
                    group_col = NULL)
  wt_list[[s]] <- dt
  message("   TME cells near epi: ", format(nrow(dt), big.mark = ","),
          "   excluded (>50µm): ",
          format(sum(dt$n_epi_within_50um == 0), big.mark = ","))
  rm(sfe); gc(verbose = FALSE)
}
wt <- rbindlist(wt_list)
# Filter to cells WITH ≥1 epi within 50 µm (per plan — cells lacking any
# epi neighbor are excluded, not assigned a zero score)
wt <- wt[n_epi_within_50um >= 1]
message(sprintf("\n[WT] total TME cells with defined niche: %s across %d samples",
                format(nrow(wt), big.mark = ","), length(unique(wt$sample_key))))

# ---------------------------------------------------------------------------
# Run on TMA (grouped by core_id)
# ---------------------------------------------------------------------------
message("\n== TMA sfe_tma_filtered ==")
sfe_t <- load_sfe("sfe_tma_filtered")
sfe_t <- override_with_06f(sfe_t, "sfe_tma")
tma <- process_sfe(sfe_t, cohort = "TMA",
                    sample_key = "sfe_tma_filtered",
                    group_col = "core_id")
tma <- tma[n_epi_within_50um >= 1]
message(sprintf("[TMA] total TME cells with defined niche: %s across %d cores / %d patients",
                format(nrow(tma), big.mark = ","),
                uniqueN(tma$group_id),
                uniqueN(tma$patient_id)))
rm(sfe_t); gc(verbose = FALSE)

# ---------------------------------------------------------------------------
# Composite niche metabolic stress + per-sample z-score + bin
# ---------------------------------------------------------------------------
# Rank-based ntile binning — always produces n_bins even with ties or
# degenerate distributions (TMA cores can have very few cells).
rank_ntile <- function(x, n_bins) {
  if (length(x) < n_bins) return(rep(NA_integer_, length(x)))
  # First-rank: distribute ties evenly; cut on rank gives balanced bins
  r <- rank(x, ties.method = "first")
  pmin(n_bins, as.integer(ceiling(r / length(r) * n_bins)))
}

add_scores <- function(dt, group_col) {
  # Z-score niche_hypoxia and niche_glycolysis within group (sample for WT,
  # core for TMA) to control per-group baseline variation
  # Use ifelse-guarded scale in case sd is 0
  safe_z <- function(v) {
    s <- sd(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) return(rep(0, length(v)))
    (v - mean(v, na.rm = TRUE)) / s
  }
  dt[, niche_hypoxia_z    := safe_z(niche_hypoxia),    by = c(group_col)]
  dt[, niche_glycolysis_z := safe_z(niche_glycolysis), by = c(group_col)]
  dt[, niche_metabolic_stress := (niche_hypoxia_z + niche_glycolysis_z) / 2]
  # Re-Z-score the composite (safer for GLMM/modeling)
  dt[, niche_metabolic_stress_z := safe_z(niche_metabolic_stress),
     by = c(group_col)]
  # Bins within group — rank-based to avoid non-unique quantile breaks
  dt[, stress_decile   := rank_ntile(niche_metabolic_stress, 10),
     by = c(group_col)]
  dt[, stress_quintile := rank_ntile(niche_metabolic_stress,  5),
     by = c(group_col)]
  dt[, stress_tertile  := rank_ntile(niche_metabolic_stress,  3),
     by = c(group_col)]
  dt
}
wt  <- add_scores(wt,  "sample_key")
tma <- add_scores(tma, "group_id")     # TMA z-scored within core

# ---------------------------------------------------------------------------
# Save master cache + summary
# ---------------------------------------------------------------------------
all_cells <- rbindlist(list(wt, tma), fill = TRUE)
OUT_RDS <- file.path(OUT_DIR, "per_cell_niche_scores.rds")
saveRDS(list(wt = wt, tma = tma, all = all_cells), OUT_RDS)
message("\nSaved: ", OUT_RDS)

# Per-sample summary
summ_wt <- wt[, .(n_cells = .N,
                   n_macrophage = sum(cell_label == "Macrophage"),
                   n_tcell      = sum(cell_label == "T cell"),
                   n_nk         = sum(cell_label == "NK cell"),
                   n_bcell      = sum(cell_label == "B cell"),
                   n_plasma     = sum(cell_label == "Plasma cell"),
                   median_n_epi = median(n_epi_within_50um),
                   median_niche_hypoxia   = round(median(niche_hypoxia), 3),
                   median_niche_glycolysis = round(median(niche_glycolysis), 3)),
                by = sample_key]
summ_tma <- tma[, .(n_cells = .N,
                     n_cores = uniqueN(group_id),
                     n_patients = uniqueN(patient_id),
                     n_macrophage = sum(cell_label == "Macrophage"),
                     n_tcell      = sum(cell_label == "T cell"),
                     n_nk         = sum(cell_label == "NK cell"),
                     n_bcell      = sum(cell_label == "B cell"),
                     n_plasma     = sum(cell_label == "Plasma cell")),
                 by = cohort]
fwrite(summ_wt,
       file.path(OUT_DIR, "niche_scores_summary_per_sample_wt.csv"))
fwrite(summ_tma,
       file.path(OUT_DIR, "niche_scores_summary_tma.csv"))

cat("\n=== WT per-sample niche score summary ===\n")
print(summ_wt)
cat("\n=== TMA niche score summary ===\n")
print(summ_tma)
message("\nDone.")
