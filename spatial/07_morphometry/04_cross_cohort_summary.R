# ============================================================================
# 04_cross_cohort_summary.R
# ----------------------------------------------------------------------------
# PURPOSE: Cross-cohort morphometric summary; monotonic SecA->Intermediate->SecB validation.
#
# INPUTS:
#   - output/33_morphometrics/ WT + TMA pairwise outputs
#
# OUTPUTS:
#   - output/33_morphometrics/cross_cohort_summary.csv
#
# MANUSCRIPT PANEL(S): Fig 5D/5E (summary)
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

OUT_DIR <- file.path(out_dir, "33_morphometrics")

wt  <- fread(file.path(OUT_DIR, "wt_pairwise_summary.csv"))
tma <- fread(file.path(OUT_DIR, "tma_pairwise_summary.csv"))

# Annotate WT with observed direction + n_concordant_either
wt[, observed_direction := ifelse(median_delta_per_sample > 0, "up", "down")]
wt[, n_concordant_either := pmax(n_concordant,
                                    ifelse(observed_direction == "up",
                                            n_concordant,
                                            n_samples - n_concordant))]
# Use sign of median_delta as the canonical "observed" direction count:
wt[, n_concordant_observed := ifelse(observed_direction == "up",
                                        sum(median_delta_per_sample > 0),
                                        sum(median_delta_per_sample < 0)),
     by = .(feature, comparison)]

# For each WT row, the n_concordant_observed is the count of samples agreeing
# with the median direction
wt2 <- wt[, .(feature, comparison, predicted_dir,
               observed_direction,
               n_samples,
               n_concordant_predicted = n_concordant,
               n_concordant_observed = ifelse(observed_direction == predicted_dir,
                                                 n_concordant,
                                                 n_samples - n_concordant),
               wt_median_delta = round(median_delta_per_sample, 3),
               wt_paired_p_two = paired_wilcox_p_two_sided,
               wt_p_bonf_two = paired_p_bonf_two_sided,
               wt_cliffs_delta = round(cliffs_delta_pooled, 3))]

tma2 <- tma[, .(feature, comparison,
                  tma_observed_direction = observed_direction,
                  tma_n_patients = n_patients,
                  tma_pct_concordant = pct_concordant_either,
                  tma_median_delta = round(median_delta_per_patient, 3),
                  tma_paired_p_two = paired_wilcox_p_two_sided,
                  tma_p_bonf_two  = paired_p_bonf_two_sided,
                  tma_cliffs_delta = round(cliffs_delta_pooled, 3))]

merged <- merge(wt2, tma2, by = c("feature","comparison"), all = TRUE)

# Cross-cohort validation: same direction in WT and TMA AND
#   WT n_concordant_observed >= 6/8  AND
#   TMA pct_concordant >= 70%   AND
#   At least one cohort Bonferroni p < 0.05
merged[, same_direction := observed_direction == tma_observed_direction]
merged[, wt_per_sample_OK := n_concordant_observed >= 6]
merged[, tma_concordance_OK := tma_pct_concordant >= 70]
merged[, any_bonf_sig := wt_p_bonf_two < 0.05 | tma_p_bonf_two < 0.05]
merged[, validated_cross_cohort := same_direction & wt_per_sample_OK &
          tma_concordance_OK & any_bonf_sig]

setorder(merged, feature, comparison)
fwrite(merged, file.path(OUT_DIR, "cross_cohort_summary.csv"))

cat("\n=== CROSS-COHORT VALIDATED FINDINGS ===\n")
v <- merged[validated_cross_cohort == TRUE]
print(v[, .(feature, comparison, predicted_dir,
              observed_dir = observed_direction,
              wt_n_concordant_observed = paste0(n_concordant_observed, "/8"),
              tma_pct = tma_pct_concordant,
              wt_p_bonf = signif(wt_p_bonf_two, 3),
              tma_p_bonf = signif(tma_p_bonf_two, 3))])
cat(sprintf("\n%d / 27 (feature × comparison) validated across cohorts\n",
             nrow(v)))

cat("\n=== MONOTONIC VALIDATED FEATURES (all 3 pairwise comparisons validated, same direction) ===\n")
mono_features <- v[, .N, by = feature][N == 3, feature]
cat("Features with monotonic SecA→Trans→SecB validation:",
    paste(mono_features, collapse = ", "), "\n")

cat("\n=== UNEXPECTED REVERSALS (predicted vs observed direction differ) ===\n")
rev <- merged[validated_cross_cohort == TRUE & predicted_dir != observed_direction]
print(rev[, .(feature, comparison, predicted_dir, observed_dir = observed_direction,
                wt_p_bonf = signif(wt_p_bonf_two, 3),
                tma_p_bonf = signif(tma_p_bonf_two, 3))])
cat(sprintf("\n%d / 27 reversals — direction observed is OPPOSITE of histopath textbook prediction\n",
             nrow(rev)))

message("\nStep 3 complete.")
