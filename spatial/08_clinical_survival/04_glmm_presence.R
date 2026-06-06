# ============================================================================
# 04_glmm_presence.R
# ----------------------------------------------------------------------------
# PURPOSE: GLMM of immune-cell presence vs niche metabolic stress.
#
# INPUTS:
#   - output/29_macrophage_niche_survival/per_cell_niche_scores.rds
#
# OUTPUTS:
#   - output/29_macrophage_niche_survival/ GLMM presence results
#
# MANUSCRIPT PANEL(S): Fig 6B
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
  library(lme4); library(data.table)
})

OUT_DIR <- file.path(out_dir, "29_macrophage_niche_survival")
CACHE   <- file.path(OUT_DIR, "per_cell_niche_scores.rds")
stopifnot(file.exists(CACHE))
cached <- readRDS(CACHE)

LYMPH_LBLS <- c("T cell", "NK cell", "B cell", "Plasma cell")

fit_glmm <- function(dt, group_col) {
  dt[, is_mac   := as.integer(cell_label == "Macrophage")]
  dt[, is_lymph := as.integer(cell_label %in% LYMPH_LBLS)]
  dt[, grp := get(group_col)]

  message("  cells: ", format(nrow(dt), big.mark = ","),
          "  groups: ", uniqueN(dt$grp),
          "  macs: ", format(sum(dt$is_mac), big.mark = ","),
          "  lymphs: ", format(sum(dt$is_lymph), big.mark = ","))

  # Fit GLMM for macrophage ~ niche_metabolic_stress_z + (1 + stress | grp)
  message("  fitting GLMM: is_macrophage ~ stress + (1+stress|grp)")
  m_mac <- tryCatch(
    glmer(is_mac ~ niche_metabolic_stress_z +
             (1 + niche_metabolic_stress_z | grp),
          data = dt, family = binomial,
          control = glmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 2e5))),
    error = function(e) { message("  mac model failed: ", e$message); NULL })
  if (is.null(m_mac)) {
    message("  falling back to random intercept only (mac)")
    m_mac <- glmer(is_mac ~ niche_metabolic_stress_z +
                     (1 | grp),
                   data = dt, family = binomial,
                   control = glmerControl(optimizer = "bobyqa",
                                           optCtrl = list(maxfun = 2e5)))
  }

  message("  fitting GLMM: is_lymphocyte ~ stress + (1+stress|grp)")
  m_lymph <- tryCatch(
    glmer(is_lymph ~ niche_metabolic_stress_z +
             (1 + niche_metabolic_stress_z | grp),
          data = dt, family = binomial,
          control = glmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 2e5))),
    error = function(e) { message("  lymph model failed: ", e$message); NULL })
  if (is.null(m_lymph)) {
    message("  falling back to random intercept only (lymph)")
    m_lymph <- glmer(is_lymph ~ niche_metabolic_stress_z + (1 | grp),
                     data = dt, family = binomial,
                     control = glmerControl(optimizer = "bobyqa",
                                             optCtrl = list(maxfun = 2e5)))
  }

  # Extract fixed-effect slopes + SEs
  coef_mac   <- summary(m_mac)$coefficients
  coef_lymph <- summary(m_lymph)$coefficients
  b_mac   <- coef_mac["niche_metabolic_stress_z", "Estimate"]
  se_mac  <- coef_mac["niche_metabolic_stress_z", "Std. Error"]
  p_mac   <- coef_mac["niche_metabolic_stress_z", "Pr(>|z|)"]
  b_lymph <- coef_lymph["niche_metabolic_stress_z", "Estimate"]
  se_lymph <- coef_lymph["niche_metabolic_stress_z", "Std. Error"]
  p_lymph  <- coef_lymph["niche_metabolic_stress_z", "Pr(>|z|)"]

  # Contrast: difference of slopes, z-test (independent model approximation)
  b_diff  <- b_mac - b_lymph
  se_diff <- sqrt(se_mac^2 + se_lymph^2)
  z_diff  <- b_diff / se_diff
  p_diff  <- 2 * pnorm(-abs(z_diff))

  list(
    model_mac = m_mac, model_lymph = m_lymph,
    summary = data.table(
      term      = c("mac_slope", "lymph_slope", "diff_slope"),
      beta      = c(b_mac, b_lymph, b_diff),
      se        = c(se_mac, se_lymph, se_diff),
      z         = c(b_mac / se_mac, b_lymph / se_lymph, z_diff),
      p         = c(p_mac, p_lymph, p_diff),
      # OR per 1 SD increase in stress
      OR_per_SD = c(exp(b_mac), exp(b_lymph), exp(b_diff))
    )
  )
}

# ---- WT ---------------------------------------------------------------------
message("\n=== WT GLMM ===")
res_wt <- fit_glmm(copy(cached$wt), "sample_key")
cat("\n-- WT summary --\n")
print(res_wt$summary)

# ---- TMA --------------------------------------------------------------------
message("\n=== TMA GLMM (patient_id as random) ===")
# For TMA, use patient_id as the random grouping (the unit of generalisation
# is the patient, not the core).
tma <- copy(cached$tma)
tma[, grp_patient := patient_id]
# Keep only cells with a valid patient_id
tma <- tma[!is.na(grp_patient) & grp_patient != ""]
res_tma <- fit_glmm(tma[, .(cell_label,
                              niche_metabolic_stress_z,
                              patient_id)],
                     "patient_id")
cat("\n-- TMA summary --\n")
print(res_tma$summary)

# ---- Sensitivity: bleed-through exclusion -----------------------------------
# TME cells with small n_epi_within_50um are likely ~adjacent to epi (bleed
# risk). Repeat the test restricting to n_epi_within_50um >= median count.
message("\n=== WT GLMM — bleed-through sensitivity (n_epi >= median) ===")
wt_dense <- cached$wt[n_epi_within_50um >= median(n_epi_within_50um)]
res_wt_dense <- fit_glmm(copy(wt_dense), "sample_key")
cat("\n-- WT dense-only summary --\n")
print(res_wt_dense$summary)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
fwrite(res_wt$summary,
       file.path(OUT_DIR, "glmm_summary_wt.csv"))
fwrite(res_tma$summary,
       file.path(OUT_DIR, "glmm_summary_tma.csv"))
fwrite(res_wt_dense$summary,
       file.path(OUT_DIR, "glmm_summary_wt_dense_sensitivity.csv"))

# Save summaries only — full lme4 models are too heavy for robust RDS
# serialisation. Recompute from the cache if models are needed again.
saveRDS(list(wt = res_wt$summary,
              tma = res_tma$summary,
              wt_dense = res_wt_dense$summary),
        file.path(OUT_DIR, "glmm_summaries.rds"))
message("\nSaved GLMM outputs to ", OUT_DIR)
