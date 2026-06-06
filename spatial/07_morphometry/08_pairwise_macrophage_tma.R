# ============================================================================
# 08_pairwise_macrophage_tma.R
# ----------------------------------------------------------------------------
# PURPOSE: Pairwise TMA macrophage morphometrics across niches (per-patient, subsampled).
#
# INPUTS:
#   - output/34_macrophage_morphometrics/ per-cell rds
#
# OUTPUTS:
#   - output/34_macrophage_morphometrics/ TMA pairwise stats
#
# MANUSCRIPT PANEL(S): Fig 6I/6J (TMA)
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
  library(data.table)
})

OUT_DIR <- file.path(out_dir, "34_macrophage_morphometrics")
cache <- readRDS(file.path(OUT_DIR, "per_cell_macrophage_morphometrics.rds"))
tma <- cache$tma

FEATURES <- c("cell_area","nuc_area","nc_ratio",
                "cell_perimeter","nuc_perimeter",
                "cell_circularity","nuc_circularity",
                "cell_solidity",
                "cell_eccentricity","nuc_eccentricity",
                "nc_centroid_offset")
N_FEATS <- length(FEATURES)
MIN_PER_TERTILE <- 10

cliffs_delta <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 5 || length(b) < 5) return(NA_real_)
  na <- length(a); nb <- length(b)
  ranks_combined <- rank(c(a, b))
  W <- sum(ranks_combined[seq_len(na)]) - na * (na + 1) / 2
  -1 * ((2 * W) / (na * nb) - 1)
}

# ---------------------------------------------------------------------------
# Per-patient tertile assignment (need at least 30 macs with epi neighbours
# for tertiles to be meaningful)
# ---------------------------------------------------------------------------
tma_sub <- tma[!is.na(niche_polarization_mean) & n_total_epi_neighbors >= 5]
# Per-patient eligibility: need at least 30 macs
patient_n <- tma_sub[, .N, by = group_id]
big_patients <- patient_n[N >= 30, group_id]
tma_sub <- tma_sub[group_id %in% big_patients]
cat(sprintf("Patients with >=30 macs (with epi neighbours): %d\n",
             length(big_patients)))

# Per-patient tertiles. If a patient's polarization range is narrow (small
# core, mostly one type), the tertiles may have ties → assign with rank.
tma_sub[, pol_rank := frank(niche_polarization_mean, ties.method = "average"),
          by = group_id]
tma_sub[, pol_pct := pol_rank / max(pol_rank), by = group_id]
tma_sub[, pol_tertile := fcase(
  pol_pct <= 1/3, "low_secA",
  pol_pct >= 2/3, "high_secB",
  default = "mid")]

n_per_t <- tma_sub[, .N, by = .(group_id, pol_tertile)]
n_wide <- dcast(n_per_t, group_id ~ pol_tertile, value.var = "N",
                  fill = 0)
n_wide[, eligible := low_secA >= MIN_PER_TERTILE &
            high_secB >= MIN_PER_TERTILE]
eligible_patients <- n_wide[eligible == TRUE, group_id]
cat(sprintf("Eligible TMA patients (>=%d macs in low AND high tertile): %d / %d\n",
             MIN_PER_TERTILE, length(eligible_patients), nrow(n_wide)))

tma_med <- tma_sub[pol_tertile %in% c("low_secA","high_secB") &
                      group_id %in% eligible_patients,
                     lapply(.SD, median, na.rm = TRUE),
                     by = .(group_id, pol_tertile),
                     .SDcols = FEATURES]

out_rows <- list()
for (feat in FEATURES) {
  per <- dcast(tma_med, group_id ~ pol_tertile, value.var = feat)
  if (!all(c("low_secA","high_secB") %in% colnames(per))) next
  per <- per[!is.na(low_secA) & !is.na(high_secB)]
  if (nrow(per) < 5) next

  delta <- per$high_secB - per$low_secA
  observed_dir <- ifelse(median(delta, na.rm = TRUE) > 0, "up", "down")
  pct_concordant <- 100 * max(mean(delta > 0), mean(delta < 0))
  n_concordant <- ifelse(observed_dir == "up",
                            sum(delta > 0), sum(delta < 0))

  w2 <- wilcox.test(per$low_secA, per$high_secB, paired = TRUE,
                     exact = FALSE)

  a_vals <- tma_sub[group_id %in% eligible_patients & pol_tertile == "low_secA",
                      get(feat)]
  b_vals <- tma_sub[group_id %in% eligible_patients & pol_tertile == "high_secB",
                      get(feat)]
  n_sub <- 30000L
  if (length(a_vals) > n_sub) a_vals <- sample(a_vals, n_sub)
  if (length(b_vals) > n_sub) b_vals <- sample(b_vals, n_sub)
  cd <- tryCatch(cliffs_delta(a_vals, b_vals),
                  error = function(e) NA_real_)

  out_rows[[length(out_rows) + 1L]] <- data.table(
    feature = feat,
    n_patients = nrow(per),
    observed_direction = observed_dir,
    n_concordant = n_concordant,
    pct_concordant = pct_concordant,
    median_delta_per_patient = median(delta, na.rm = TRUE),
    paired_wilcox_p_two_sided = w2$p.value,
    cliffs_delta_pooled = cd
  )
}
result <- rbindlist(out_rows, fill = TRUE)
result[, paired_p_bonf := pmin(paired_wilcox_p_two_sided * N_FEATS, 1)]
result[, validated := pct_concordant >= 70 & paired_p_bonf < 0.05]
setorder(result, feature)

fwrite(result, file.path(OUT_DIR, "tma_pairwise_summary.csv"))

cat("\n=== TMA PAIRWISE SUMMARY (low vs high niche-polarization tertile) ===\n")
print(result[, .(feature, observed_direction,
                   n_patients, pct_concordant = round(pct_concordant, 1),
                   median_delta = round(median_delta_per_patient, 3),
                   p_two = signif(paired_wilcox_p_two_sided, 3),
                   p_bonf = signif(paired_p_bonf, 3),
                   cliffs_delta = round(cliffs_delta_pooled, 3),
                   validated)])

cat(sprintf("\nValidated features in TMA (>=70%% same dir AND Bonferroni p<0.05): %d / %d\n",
             sum(result$validated, na.rm = TRUE), nrow(result)))

message("\nStep 2 complete.")
