# ============================================================================
# 09_cross_cohort_summary.R
# ----------------------------------------------------------------------------
# PURPOSE: Cross-cohort macrophage morphometric summary (WT + TMA).
#
# INPUTS:
#   - output/34_macrophage_morphometrics/ WT + TMA pairwise outputs
#
# OUTPUTS:
#   - output/34_macrophage_morphometrics/cross_cohort_summary.csv
#
# MANUSCRIPT PANEL(S): Fig 6I/6J (summary)
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
wt  <- fread(file.path(OUT_DIR, "wt_pairwise_summary.csv"))
tma <- fread(file.path(OUT_DIR, "tma_pairwise_summary.csv"))

setnames(wt,
          c("observed_direction","n_concordant_observed","median_delta_per_sample",
            "paired_wilcox_p_two_sided","paired_p_bonf","cliffs_delta_pooled","validated"),
          c("wt_dir","wt_n_concordant","wt_median_delta",
            "wt_paired_p_two","wt_p_bonf","wt_cliffs_delta","wt_validated"))
setnames(tma,
          c("observed_direction","n_concordant","pct_concordant","median_delta_per_patient",
            "paired_wilcox_p_two_sided","paired_p_bonf","cliffs_delta_pooled","validated"),
          c("tma_dir","tma_n_concordant","tma_pct_concordant","tma_median_delta",
            "tma_paired_p_two","tma_p_bonf","tma_cliffs_delta","tma_validated"))

merged <- merge(wt, tma, by = "feature", all = TRUE)
merged[, same_direction := wt_dir == tma_dir]
# WT n=8 hits Bonferroni floor (× 11 = 0.157) regardless of concordance —
# headline bar: WT ≥ 6/8 same direction AND TMA ≥ 70% concordant AND
# TMA Bonferroni-significant. WT Bonferroni not required.
merged[, validated_cross_cohort := same_direction & wt_n_concordant >= 6 &
          tma_pct_concordant >= 70 & tma_p_bonf < 0.05]
setorder(merged, -wt_n_concordant)

fwrite(merged, file.path(OUT_DIR, "cross_cohort_summary.csv"))

cat("\n=== CROSS-COHORT VALIDATED MACROPHAGE FEATURES ===\n")
v <- merged[validated_cross_cohort == TRUE]
print(v[, .(feature, wt_dir, tma_dir,
              wt_n = paste0(wt_n_concordant, "/", n_samples),
              tma_pct = round(tma_pct_concordant, 1),
              wt_p_bonf = signif(wt_p_bonf, 3),
              tma_p_bonf = signif(tma_p_bonf, 3),
              wt_cd = round(wt_cliffs_delta, 3),
              tma_cd = round(tma_cliffs_delta, 3))])
cat(sprintf("\n%d / %d features fully validated cross-cohort\n",
             nrow(v), nrow(merged)))

cat("\n=== ADDITIONAL DIRECTIONALLY-CONCORDANT FEATURES (8/8 WT but TMA <70% concordant) ===\n")
sug <- merged[same_direction == TRUE & wt_n_concordant >= 8 &
                 validated_cross_cohort == FALSE]
if (nrow(sug) > 0) {
  print(sug[, .(feature, wt_dir, tma_dir,
                  wt_n = paste0(wt_n_concordant, "/", n_samples),
                  tma_pct = round(tma_pct_concordant, 1),
                  wt_p_two = signif(wt_paired_p_two, 3),
                  tma_p_bonf = signif(tma_p_bonf, 3))])
} else {
  cat("  None.\n")
}

message("\nStep 3 complete.")
