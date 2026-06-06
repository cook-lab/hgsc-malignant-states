# ============================================================================
# 05_paired_enrichment.R
# ----------------------------------------------------------------------------
# PURPOSE: Paired enrichment of immune populations across low vs high metabolic-stress niches.
#
# INPUTS:
#   - output/29_macrophage_niche_survival/per_cell_niche_scores.rds
#
# OUTPUTS:
#   - output/29_macrophage_niche_survival/ paired enrichment tables
#
# MANUSCRIPT PANEL(S): Fig 6E/6F
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

suppressPackageStartupMessages({
  library(data.table); library(boot)
})

OUT_DIR <- file.path(out_dir, "29_macrophage_niche_survival")
CACHE   <- file.path(OUT_DIR, "per_cell_niche_scores.rds")
stopifnot(file.exists(CACHE))
cached <- readRDS(CACHE)

MIN_CELLS_PER_GROUP <- 30

LYMPH_LBLS <- c("T cell", "NK cell", "B cell", "Plasma cell")

# ---------------------------------------------------------------------------
# Compute enrichment ratio for ONE patient. Takes plain vectors, not .SD.
# ---------------------------------------------------------------------------
enrichment_of_set <- function(is_target, bin_vec) {
  bin_max <- max(bin_vec, na.rm = TRUE)
  top <- bin_vec == bin_max
  bot <- bin_vec == 1L
  n_top <- sum(top); n_bot <- sum(bot)
  n_tar_top <- sum(is_target & top)
  n_tar_bot <- sum(is_target & bot)
  # Laplace smoothing to avoid divide-by-zero
  p_top <- (n_tar_top + 1) / (n_top + 2)
  p_bot <- (n_tar_bot + 1) / (n_bot + 2)
  list(enrichment = p_top / p_bot,
       p_top = p_top, p_bot = p_bot,
       n_target_top = n_tar_top, n_target_bot = n_tar_bot,
       n_top = n_top, n_bot = n_bot)
}

per_patient_LR <- function(dt, group_col, bin_col,
                            lymph_set = LYMPH_LBLS) {
  # Build vectors we need, avoid .SD locking by computing explicitly
  dt[, {
    n_mac   <- sum(cell_label == "Macrophage")
    n_lymph <- sum(cell_label %in% lymph_set)
    if (n_mac < MIN_CELLS_PER_GROUP || n_lymph < MIN_CELLS_PER_GROUP) {
      list(eligible = FALSE,
           n_mac = n_mac, n_lymph = n_lymph,
           E_mac = NA_real_, E_lymph = NA_real_, LR = NA_real_)
    } else {
      is_mac   <- cell_label == "Macrophage"
      is_lymph <- cell_label %in% lymph_set
      bin_vec  <- get(bin_col)
      e_mac    <- enrichment_of_set(is_mac,   bin_vec)
      e_lymph  <- enrichment_of_set(is_lymph, bin_vec)
      list(eligible = TRUE,
           n_mac = n_mac, n_lymph = n_lymph,
           E_mac = e_mac$enrichment,
           E_lymph = e_lymph$enrichment,
           LR = log2(e_mac$enrichment / e_lymph$enrichment))
    }
  }, by = group_col]
}

summarize_LR <- function(lr_dt, label) {
  elig <- lr_dt[eligible == TRUE]
  if (nrow(elig) < 2) {
    return(data.table(cohort = label,
                      n_groups = 0,
                      median_LR = NA_real_,
                      ci95_lo = NA_real_, ci95_hi = NA_real_,
                      p_wilcox = NA_real_,
                      pct_positive = NA_real_))
  }
  # Paired Wilcoxon = one-sample Wilcoxon on LR (because LR is log-ratio of
  # two per-patient enrichments — positive LR = macrophage-enriched)
  wil <- wilcox.test(elig$LR, mu = 0, alternative = "two.sided")
  # Bootstrap CI for median LR
  set.seed(29)
  b <- boot::boot(elig$LR, statistic = function(x, i) median(x[i]), R = 2000)
  ci <- tryCatch(boot::boot.ci(b, type = "perc")$percent[4:5],
                  error = function(e) c(NA, NA))
  data.table(cohort = label,
             n_groups = nrow(elig),
             median_LR = median(elig$LR),
             ci95_lo = ci[1], ci95_hi = ci[2],
             p_wilcox = wil$p.value,
             pct_positive = 100 * mean(elig$LR > 0))
}

# ---------------------------------------------------------------------------
# Per-cohort: WT by sample_key, TMA by patient_id
# For multiple bin definitions AND multiple lymph-subsets
# ---------------------------------------------------------------------------
BINS <- c("stress_decile", "stress_quintile", "stress_tertile")
LYMPH_SETS <- list(
  "T+NK+B+plasma"           = LYMPH_LBLS,
  "T cell only"             = "T cell",
  "NK cell only"            = "NK cell",
  "B cell only"             = "B cell",
  "Plasma cell only"        = "Plasma cell"
)

run_cohort <- function(cells, group_col, cohort_label) {
  summaries <- list()
  per_patient_all <- list()
  for (bn in BINS) {
    for (ls_name in names(LYMPH_SETS)) {
      ls <- LYMPH_SETS[[ls_name]]
      lr_dt <- per_patient_LR(cells, group_col, bn, ls)
      lr_dt[, bin := bn]
      lr_dt[, lymph_set := ls_name]
      per_patient_all[[paste(bn, ls_name, sep = "__")]] <- lr_dt
      sm <- summarize_LR(lr_dt, cohort_label)
      sm[, bin := bn]
      sm[, lymph_set := ls_name]
      summaries[[paste(bn, ls_name, sep = "__")]] <- sm
    }
  }
  list(
    summary = rbindlist(summaries, fill = TRUE),
    per_patient = rbindlist(per_patient_all, fill = TRUE)
  )
}

message("=== WT paired enrichment (by sample_key) ===")
wt_res <- run_cohort(cached$wt, "sample_key", "WT")
setorder(wt_res$summary, bin, lymph_set)
print(wt_res$summary)

message("\n=== TMA paired enrichment (by patient_id) ===")
tma <- cached$tma[!is.na(patient_id) & patient_id != ""]
tma_res <- run_cohort(tma, "patient_id", "TMA")
setorder(tma_res$summary, bin, lymph_set)
print(tma_res$summary)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
all_summary <- rbind(wt_res$summary, tma_res$summary)
all_per_patient <- rbind(wt_res$per_patient, tma_res$per_patient, fill = TRUE)

fwrite(all_summary,
       file.path(OUT_DIR, "paired_enrichment_summary.csv"))
fwrite(wt_res$per_patient,
       file.path(OUT_DIR, "paired_enrichment_per_sample_wt.csv"))
fwrite(tma_res$per_patient,
       file.path(OUT_DIR, "paired_enrichment_per_patient_tma.csv"))
saveRDS(list(wt = wt_res, tma = tma_res),
        file.path(OUT_DIR, "paired_enrichment_results.rds"))
message("\nSaved paired enrichment outputs to ", OUT_DIR)
