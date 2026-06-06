# ============================================================================
# 02_pairwise_wt_morphometrics.R
# ----------------------------------------------------------------------------
# PURPOSE: Pairwise WT morphometric comparisons (SecA->Intermediate, Intermediate->SecB, SecA->SecB) with subsampling + Wilcoxon.
#
# INPUTS:
#   - output/33_morphometrics/ per-cell WT morphometrics
#
# OUTPUTS:
#   - output/33_morphometrics/per_sample_summary_wt.csv + pairwise stats
#
# MANUSCRIPT PANEL(S): Fig 5D/5E (WT)
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

cache <- readRDS(file.path(OUT_DIR, "per_cell_morphometrics.rds"))
wt <- cache$wt

FEATURES <- c("cell_area","nuc_area","nc_ratio",
                "cell_perimeter","nuc_perimeter",
                "cell_circularity","nuc_circularity",
                "cell_solidity",
                "cell_eccentricity","nuc_eccentricity",
                "nc_centroid_offset")

# Predicted direction per feature (positive = increases A→B)
PREDICT <- c(
  cell_area          = "down",
  nuc_area           = "up",
  nc_ratio           = "up",
  cell_perimeter     = "up",
  nuc_perimeter      = "up",
  cell_circularity   = "down",
  nuc_circularity    = "down",
  cell_solidity      = "down",
  cell_eccentricity  = "up",
  nuc_eccentricity   = "up",
  nc_centroid_offset = "down"
)

PAIRWISE <- list(
  list(label = "SecA_to_Trans", a = "SecA epithelium",
       b = "Intermediate epithelium"),
  list(label = "Trans_to_SecB", a = "Intermediate epithelium",
       b = "SecB epithelium"),
  list(label = "SecA_to_SecB", a = "SecA epithelium",
       b = "SecB epithelium")
)

# ---------------------------------------------------------------------------
# Cliff's delta via Mann-Whitney U statistic
#   delta = 2*P(a < b) - 1, computed via combined ranking
# ---------------------------------------------------------------------------
cliffs_delta <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 5 || length(b) < 5) return(NA_real_)
  na <- length(a); nb <- length(b)
  ranks_combined <- rank(c(a, b))
  W <- sum(ranks_combined[seq_len(na)]) - na * (na + 1) / 2
  # Note: convention here gives delta with sign = (b - a) direction
  -1 * ((2 * W) / (na * nb) - 1)
}

# ---------------------------------------------------------------------------
# Per-sample medians per (sample × cell_label × feature)
# ---------------------------------------------------------------------------
wt_med <- wt[cell_label %in% c("SecA epithelium",
                                  "Intermediate epithelium",
                                  "SecB epithelium"),
               lapply(.SD, median, na.rm = TRUE),
               by = .(sample_key, cell_label),
               .SDcols = FEATURES]

# Per-sample n by epitype (eligibility filter)
wt_n <- wt[cell_label %in% c("SecA epithelium",
                                "Intermediate epithelium",
                                "SecB epithelium"),
             .(n = .N), by = .(sample_key, cell_label)]
wt_n_wide <- dcast(wt_n, sample_key ~ cell_label, value.var = "n")
setnames(wt_n_wide,
          c("SecA epithelium","Intermediate epithelium","SecB epithelium"),
          c("n_SecA","n_Int","n_SecB"), skip_absent = TRUE)
wt_n_wide[, eligible := !is.na(n_SecA) & !is.na(n_Int) & !is.na(n_SecB) &
            n_SecA >= 30 & n_Int >= 30 & n_SecB >= 30]

eligible_samples <- wt_n_wide[eligible == TRUE, sample_key]
message(sprintf("Eligible WT samples: %d / 8", length(eligible_samples)))

# ---------------------------------------------------------------------------
# Pairwise tests
# ---------------------------------------------------------------------------
out_rows <- list()
for (feat in FEATURES) {
  for (pp in PAIRWISE) {
    # Per-sample medians (paired)
    per <- dcast(wt_med[cell_label %in% c(pp$a, pp$b) &
                           sample_key %in% eligible_samples],
                  sample_key ~ cell_label, value.var = feat)
    if (!all(c(pp$a, pp$b) %in% colnames(per))) next
    per <- per[!is.na(get(pp$a)) & !is.na(get(pp$b))]
    if (nrow(per) < 3) next
    delta <- per[[pp$b]] - per[[pp$a]]

    # Predicted direction.
    # delta = b - a. If pred = "up" then expect delta > 0 → wilcox.test
    # paired (a, b) with alternative = "less" tests a < b (b > a).
    # If pred = "down" expect delta < 0 → alternative = "greater" (a > b).
    pred <- PREDICT[[feat]]
    if (pred == "up") {
      n_concordant <- sum(delta > 0)
      alt <- "less"            # a < b
    } else {
      n_concordant <- sum(delta < 0)
      alt <- "greater"         # a > b
    }

    w  <- wilcox.test(per[[pp$a]], per[[pp$b]], paired = TRUE,
                       alternative = alt, exact = FALSE)
    w2 <- wilcox.test(per[[pp$a]], per[[pp$b]], paired = TRUE,
                       exact = FALSE)   # two-sided

    # Cliff's δ approximated by per-cell-distribution rank-based effect size:
    # Compute on a random 30k × 30k subsample for tractable runtime.
    a_vals <- wt[sample_key %in% eligible_samples & cell_label == pp$a,
                   get(feat)]
    b_vals <- wt[sample_key %in% eligible_samples & cell_label == pp$b,
                   get(feat)]
    n_sub <- 30000L
    if (length(a_vals) > n_sub) a_vals <- sample(a_vals, n_sub)
    if (length(b_vals) > n_sub) b_vals <- sample(b_vals, n_sub)
    cd <- tryCatch(cliffs_delta(a_vals, b_vals),
                    error = function(e) NA_real_)

    out_rows[[length(out_rows) + 1L]] <- data.table(
      feature       = feat,
      comparison    = pp$label,
      predicted_dir = pred,
      n_samples     = nrow(per),
      n_concordant  = n_concordant,
      median_delta_per_sample = median(delta, na.rm = TRUE),
      paired_wilcox_p_one_sided = w$p.value,
      paired_wilcox_p_two_sided = w2$p.value,
      cliffs_delta_pooled = cd
    )
  }
}
result <- rbindlist(out_rows, fill = TRUE)

# Bonferroni for 27 tests
result[, paired_p_bonf_one_sided := pmin(paired_wilcox_p_one_sided * .N, 1)]
result[, paired_p_bonf_two_sided := pmin(paired_wilcox_p_two_sided * .N, 1)]
result[, validated_per_sample := n_concordant >= 6]
result[, validated_paired_p   := paired_p_bonf_one_sided < 0.05]

setorder(result, feature, comparison)
fwrite(result, file.path(OUT_DIR, "wt_pairwise_summary.csv"))

cat("\n=== WT PAIRWISE SUMMARY ===\n")
print(result[, .(feature, comparison, predicted_dir,
                   n_samples, n_concordant,
                   median_delta = round(median_delta_per_sample, 3),
                   p_one = signif(paired_wilcox_p_one_sided, 3),
                   p_bonf = signif(paired_p_bonf_one_sided, 3),
                   cliffs_delta = round(cliffs_delta_pooled, 3),
                   per_sample_OK = validated_per_sample,
                   bonf_p_OK = validated_paired_p)])

cat("\n=== FEATURES VALIDATED IN WT (per-sample >=6/8 AND Bonferroni p<0.05) ===\n")
ok <- result[validated_per_sample == TRUE & validated_paired_p == TRUE]
print(ok[, .(feature, comparison, predicted_dir, n_concordant,
              p_bonf = signif(paired_p_bonf_one_sided, 3),
              median_delta = round(median_delta_per_sample, 3))])
cat(sprintf("\n%d / 27 (feature × comparison) cells pass both bars in WT\n",
             nrow(ok)))

cat("\n=== MONOTONIC FEATURES (all 3 comparisons concordant + same direction) ===\n")
mono <- result[validated_per_sample == TRUE,
                 .N, by = feature][N == 3, feature]
cat("Features with all 3 pairwise comparisons concordant per-sample:",
    paste(mono, collapse = ", "), "\n")

message("\nStep 1 complete.")
