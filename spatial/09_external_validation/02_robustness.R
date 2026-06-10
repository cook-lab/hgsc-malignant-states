# ============================================================================
# 02_robustness.R
# ----------------------------------------------------------------------------
# PURPOSE: Robustness of the TCGA validation: leave-one-gene-out signatures, GSVA / ssGSEA scoring-method sensitivity, bootstrap stability.
#
# INPUTS:
#   - output/40_tcga_validation/ score tables
#   - SecA/SecB signatures from shared/signatures.yml
#
# OUTPUTS:
#   - output/40_tcga_validation/ robustness tables
#
# MANUSCRIPT PANEL(S): Robustness for Fig 7E/F/G validation.
# RUNTIME TIER: heavy
#
# Migrated from 2026_final_xenium_analysis/scripts/. Analytical logic preserved;
# paths routed through central config, epithelial label "Transitioning" ->
# "Intermediate", SecA/SecB from shared/signatures.yml.
#
# SEEDING: a global set.seed(CFG$seed) is set at startup, BUT the bootstrap
# block below intentionally re-seeds with a FIXED LITERAL seed,
# set.seed(20260508), immediately before resampling. That literal seed is what
# produced the published bootstrap 95% CIs (and the downstream c-index bootstrap
# inherits the same RNG stream), so it is preserved verbatim — do NOT switch it
# to CFG$seed (that would change the published CIs). See flagged_for_user.
# ============================================================================

# --- Config + shared setup (replaces hardcoded /Volumes/CookLab/Sarah paths) ---
here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
source(file.path(here, "..", "..", "config", "config.R"))   # CFG, cfg_obj, cfg_path
source(file.path(here, "..", "00_setup", "00_setup.R"))      # load_sfe, save_sfe, theme_lab, nb_names, palettes
set.seed(CFG$seed)

# SecA/SecB signatures loaded from the shared source of truth (noBCAM 7-gene set)
.sigs      <- yaml::read_yaml(file.path(here, "..", "..", "shared", "signatures.yml"))
secA_genes <- .sigs$SecA
secB_genes <- .sigs$SecB
sig_list   <- list(SecA = secA_genes, SecB = secB_genes)

suppressPackageStartupMessages({
  library(UCell)
  library(survival)
  library(survRM2)
  library(singscore)
  library(GSVA)
})

step_dir <- file.path(out_dir, "40_tcga_validation")
run_log  <- file.path(step_dir, "run.log")
sink_log <- function(msg) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg, "\n",
                              file = run_log, append = TRUE)
sink_log("Phase 40b started")

phase40 <- readRDS(file.path(step_dir, "phase40_primary.rds"))
cohort  <- phase40$cohort
fam     <- phase40$fam

primary_scores <- c("polar_UCell_z", "log2_ratio_UCell_z",
                    "prop_SecB_of_sec_z", "log2_ratio_csx_z")
endpoints      <- c("OS5", "PFS5")

# --- Helper: Cox HR for a score with covariates -----------------------------
fit_cox_hr <- function(df, score, time_col, event_col, covars = c("Age", "stage_grp")) {
  rhs <- paste(c(score, covars), collapse = " + ")
  f <- as.formula(sprintf("Surv(%s, %s) ~ %s", time_col, event_col, rhs))
  df_fit <- df[!is.na(get(score)) & !is.na(get(time_col)) & !is.na(get(event_col))]
  for (cv in covars) df_fit <- df_fit[!is.na(get(cv))]
  if (nrow(df_fit) < 30 || sum(df_fit[[event_col]]) < 10) return(NULL)
  fit <- tryCatch(coxph(f, data = df_fit), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  list(fit = fit, n = nrow(df_fit), df = df_fit)
}

# --- Bootstrap 95% CI -------------------------------------------------------
sink_log("Bootstrap 95% CIs (1000 resamples)")
# Fixed literal seed — reproduces the PUBLISHED bootstrap CIs. Preserved verbatim
# (NOT CFG$seed) by design; see header. Changing it would change published numbers.
set.seed(20260508)
B <- 1000

boot_results <- list()
for (ep in endpoints) {
  tcol <- paste0(ep, "_months"); ecol <- paste0(ep, "_event")
  for (sc in primary_scores) {
    df_sc <- if (sc %in% c("prop_SecB_of_sec_z", "log2_ratio_csx_z"))
      cohort[csx_pvalue < 0.05 & !is.na(csx_pvalue)] else cohort
    res <- fit_cox_hr(df_sc, sc, tcol, ecol)
    if (is.null(res)) next
    df_fit <- res$df
    n_fit <- nrow(df_fit)
    boot_HR <- numeric(B)
    for (b in 1:B) {
      idx <- sample.int(n_fit, n_fit, replace = TRUE)
      r2 <- fit_cox_hr(df_fit[idx], sc, tcol, ecol)
      boot_HR[b] <- if (is.null(r2)) NA_real_ else exp(coef(r2$fit)[[sc]])
    }
    ci <- quantile(boot_HR, c(0.025, 0.975), na.rm = TRUE)
    boot_results[[length(boot_results) + 1]] <- data.table(
      endpoint = ep, score = sc,
      HR_point = exp(coef(res$fit)[sc]),
      HR_boot_lo = ci[1], HR_boot_hi = ci[2],
      boot_n_valid = sum(!is.na(boot_HR)))
  }
}
boot_dt <- rbindlist(boot_results)
fwrite(boot_dt, file.path(step_dir, "robustness_bootstrap.csv"))

# --- Schoenfeld PH check + RMST fallback -----------------------------------
sink_log("Schoenfeld PH check + RMST fallback")
ph_results <- list()
rmst_results <- list()
for (ep in endpoints) {
  tcol <- paste0(ep, "_months"); ecol <- paste0(ep, "_event")
  for (sc in primary_scores) {
    df_sc <- if (sc %in% c("prop_SecB_of_sec_z", "log2_ratio_csx_z"))
      cohort[csx_pvalue < 0.05 & !is.na(csx_pvalue)] else cohort
    res <- fit_cox_hr(df_sc, sc, tcol, ecol)
    if (is.null(res)) next
    zph <- tryCatch(cox.zph(res$fit), error = function(e) NULL)
    p_score <- if (!is.null(zph)) zph$table[sc, "p"] else NA_real_
    p_global <- if (!is.null(zph)) zph$table["GLOBAL", "p"] else NA_real_
    ph_results[[length(ph_results) + 1]] <- data.table(
      endpoint = ep, score = sc, ph_p_score = p_score, ph_p_global = p_global,
      ph_violated = !is.na(p_score) & p_score < 0.05)

    # RMST at 5y by score-tertile (low vs high)
    df_fit <- res$df
    df_fit[, score_tert := cut(get(sc), quantile(get(sc), c(0, 1/3, 2/3, 1),
                               na.rm = TRUE), include.lowest = TRUE,
                               labels = c("low", "mid", "high"))]
    df_lh <- df_fit[score_tert %in% c("low", "high")]
    if (nrow(df_lh) >= 50) {
      arm <- as.integer(df_lh$score_tert == "high")
      rm <- tryCatch(rmst2(time = df_lh[[tcol]], status = df_lh[[ecol]],
                           arm = arm, tau = 60),
                     error = function(e) NULL)
      if (!is.null(rm)) {
        diff_row <- rm$unadjusted.result[1, ]
        rmst_results[[length(rmst_results) + 1]] <- data.table(
          endpoint = ep, score = sc,
          rmst_diff_high_minus_low = diff_row["Est."],
          rmst_p = diff_row["p"])
      }
    }
  }
}
ph_dt   <- rbindlist(ph_results)
rmst_dt <- rbindlist(rmst_results, fill = TRUE)
fwrite(ph_dt,   file.path(step_dir, "robustness_ph.csv"))
fwrite(rmst_dt, file.path(step_dir, "robustness_rmst.csv"))

# --- C-index: age+stage vs age+stage+score ---------------------------------
sink_log("C-index comparison age+stage vs age+stage+score (1000 bootstraps)")
cindex_results <- list()
for (ep in endpoints) {
  tcol <- paste0(ep, "_months"); ecol <- paste0(ep, "_event")
  for (sc in primary_scores) {
    df_sc <- if (sc %in% c("prop_SecB_of_sec_z", "log2_ratio_csx_z"))
      cohort[csx_pvalue < 0.05 & !is.na(csx_pvalue)] else cohort
    df_sc <- df_sc[!is.na(get(sc)) & !is.na(Age) & !is.na(stage_grp) &
                    !is.na(get(tcol)) & !is.na(get(ecol))]
    if (nrow(df_sc) < 50) next

    f_base <- as.formula(sprintf("Surv(%s, %s) ~ Age + stage_grp", tcol, ecol))
    f_full <- as.formula(sprintf("Surv(%s, %s) ~ Age + stage_grp + %s",
                                 tcol, ecol, sc))
    fit_b <- coxph(f_base, data = df_sc)
    fit_f <- coxph(f_full, data = df_sc)
    c_b <- summary(fit_b)$concordance[1]
    c_f <- summary(fit_f)$concordance[1]

    # Bootstrap delta-c
    n_fit <- nrow(df_sc)
    delta_boot <- numeric(B)
    for (b in 1:B) {
      idx <- sample.int(n_fit, n_fit, replace = TRUE)
      df_b <- df_sc[idx]
      fit_b_b <- tryCatch(coxph(f_base, data = df_b), error = function(e) NULL)
      fit_f_b <- tryCatch(coxph(f_full, data = df_b), error = function(e) NULL)
      delta_boot[b] <- if (is.null(fit_b_b) || is.null(fit_f_b)) NA_real_
                       else summary(fit_f_b)$concordance[1] -
                            summary(fit_b_b)$concordance[1]
    }
    p_boot <- 2 * min(mean(delta_boot <= 0, na.rm = TRUE),
                       mean(delta_boot >= 0, na.rm = TRUE))
    cindex_results[[length(cindex_results) + 1]] <- data.table(
      endpoint = ep, score = sc, n = nrow(df_sc),
      c_base = c_b, c_full = c_f, delta_c = c_f - c_b,
      delta_c_boot_lo = quantile(delta_boot, 0.025, na.rm = TRUE),
      delta_c_boot_hi = quantile(delta_boot, 0.975, na.rm = TRUE),
      delta_c_boot_p  = p_boot)
  }
}
cindex_dt <- rbindlist(cindex_results, fill = TRUE)
fwrite(cindex_dt, file.path(step_dir, "robustness_cindex.csv"))

# --- Leave-one-gene-out (LOGO) ---------------------------------------------
sink_log("Leave-one-gene-out for SecA/SecB on 5-yr OS")
all_sig <- c(secA_genes, secB_genes)

# Need the original log-expression matrix on locked aliquots
tcga_rna <- fread(file.path(data_dir, "TCGA_data", "tcga_data_plotting.csv"))
expr_cols <- setdiff(names(tcga_rna), c("V1", "hgnc_symbol", "X"))
expr_mat <- as.matrix(tcga_rna[, ..expr_cols]); rownames(expr_mat) <- tcga_rna$hgnc_symbol
expr_mat <- expr_mat[!is.na(rownames(expr_mat)) & nzchar(rownames(expr_mat)), , drop = FALSE]
if (any(duplicated(rownames(expr_mat)))) {
  ord <- order(rowMeans(expr_mat, na.rm = TRUE), decreasing = TRUE)
  expr_mat <- expr_mat[ord, , drop = FALSE]
  expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), , drop = FALSE]
}
expr_locked <- log2(expr_mat[, cohort$aliquot_id_csv, drop = FALSE] + 1)

# Compute UCell once per gene-drop scenario, then run Cox for BOTH endpoints
score_with_drop <- function(drop_gene) {
  if (is.null(drop_gene)) {
    a <- secA_genes; b <- secB_genes
  } else {
    a <- if (drop_gene %in% secA_genes) setdiff(secA_genes, drop_gene) else secA_genes
    b <- if (drop_gene %in% secB_genes) setdiff(secB_genes, drop_gene) else secB_genes
  }
  uc <- ScoreSignatures_UCell(expr_locked,
            features = list(SecA = a, SecB = b),
            ncores = 1, chunk.size = 500)
  z <- scale(uc[, "SecB_UCell"] - uc[, "SecA_UCell"])[, 1]
  z
}

logo_results <- list()
# Full signature
cohort[, polar_logo_z := score_with_drop(NULL)]
for (ep in endpoints) {
  tcol <- paste0(ep, "_months"); ecol <- paste0(ep, "_event")
  r <- fit_cox_hr(cohort, "polar_logo_z", tcol, ecol)
  if (is.null(r)) next
  logo_results[[length(logo_results) + 1]] <- data.table(
    endpoint = ep, dropped_gene = "(none)", set = "(full)",
    HR = exp(coef(r$fit)["polar_logo_z"]),
    p  = summary(r$fit)$coefficients["polar_logo_z", "Pr(>|z|)"])
}
for (g in all_sig) {
  cohort[, polar_logo_z := score_with_drop(g)]
  set_lbl <- if (g %in% secA_genes) "SecA" else "SecB"
  for (ep in endpoints) {
    tcol <- paste0(ep, "_months"); ecol <- paste0(ep, "_event")
    r <- fit_cox_hr(cohort, "polar_logo_z", tcol, ecol)
    if (is.null(r)) next
    logo_results[[length(logo_results) + 1]] <- data.table(
      endpoint = ep, dropped_gene = g, set = set_lbl,
      HR = exp(coef(r$fit)["polar_logo_z"]),
      p  = summary(r$fit)$coefficients["polar_logo_z", "Pr(>|z|)"])
  }
}
logo_dt <- rbindlist(logo_results)
fwrite(logo_dt, file.path(step_dir, "robustness_logo.csv"))

# --- Scoring-method sensitivity --------------------------------------------
sink_log("Scoring-method sensitivity (singscore + GSVA ssGSEA + GSVA gsva)")
# singscore
rs <- rankGenes(expr_locked)
sec_a_ss <- simpleScore(rs, upSet = secA_genes, knownDirection = FALSE)$TotalScore
sec_b_ss <- simpleScore(rs, upSet = secB_genes, knownDirection = FALSE)$TotalScore
cohort[, polar_singscore := sec_b_ss - sec_a_ss]
cohort[, polar_singscore_z := scale(polar_singscore)[, 1]]

# GSVA — gsva method
gsva_out <- tryCatch(
  gsva(expr_locked, sig_list,
       method = "gsva", kcdf = "Gaussian", verbose = FALSE),
  error = function(e) {
    par <- gsvaParam(exprData = expr_locked, geneSets = sig_list,
                     kcdf = "Gaussian")
    gsva(par, verbose = FALSE)
  })
cohort[, polar_gsva := gsva_out["SecB", ] - gsva_out["SecA", ]]
cohort[, polar_gsva_z := scale(polar_gsva)[, 1]]

# ssGSEA
ssgsea_out <- tryCatch(
  gsva(expr_locked, sig_list, method = "ssgsea", verbose = FALSE),
  error = function(e) {
    par <- ssgseaParam(exprData = expr_locked, geneSets = sig_list)
    gsva(par, verbose = FALSE)
  })
cohort[, polar_ssgsea := ssgsea_out["SecB", ] - ssgsea_out["SecA", ]]
cohort[, polar_ssgsea_z := scale(polar_ssgsea)[, 1]]

scorer_methods <- c(UCell = "polar_UCell_z",
                    singscore = "polar_singscore_z",
                    GSVA = "polar_gsva_z",
                    ssGSEA = "polar_ssgsea_z")
method_results <- list()
for (ep in endpoints) {
  tcol <- paste0(ep, "_months"); ecol <- paste0(ep, "_event")
  for (m in names(scorer_methods)) {
    sc <- scorer_methods[[m]]
    r <- fit_cox_hr(cohort, sc, tcol, ecol)
    if (is.null(r)) next
    s <- summary(r$fit)
    method_results[[length(method_results) + 1]] <- data.table(
      endpoint = ep, method = m, score = sc,
      n = r$n, HR = exp(coef(r$fit)[sc]),
      HR_lo = s$conf.int[sc, "lower .95"],
      HR_hi = s$conf.int[sc, "upper .95"],
      p = s$coefficients[sc, "Pr(>|z|)"])
  }
}
method_dt <- rbindlist(method_results, fill = TRUE)
fwrite(method_dt, file.path(step_dir, "robustness_methods.csv"))

# --- Save combined RDS -----------------------------------------------------
saveRDS(list(bootstrap = boot_dt, ph = ph_dt, rmst = rmst_dt,
             cindex = cindex_dt, logo = logo_dt, methods = method_dt),
        file.path(step_dir, "phase40_robustness.rds"))

cat("\n=== ROBUSTNESS — Bootstrap CIs (age+stage primary model) ===\n")
print(boot_dt[, .(endpoint, score,
                  HR = round(HR_point, 3),
                  CI_boot = sprintf("[%.2f, %.2f]", HR_boot_lo, HR_boot_hi))])

cat("\n=== ROBUSTNESS — PH check (Schoenfeld) ===\n")
print(ph_dt[, .(endpoint, score,
                ph_p_score = signif(ph_p_score, 3),
                ph_p_global = signif(ph_p_global, 3),
                ph_violated)])

cat("\n=== ROBUSTNESS — c-index Δ ===\n")
print(cindex_dt[, .(endpoint, score, n,
                    c_base = round(c_base, 3),
                    c_full = round(c_full, 3),
                    dC = round(delta_c, 4),
                    boot_p = signif(delta_c_boot_p, 3))])

cat("\n=== ROBUSTNESS — LOGO summary (5-yr OS, polar UCell) ===\n")
print(logo_dt[endpoint == "OS5",
              .(dropped_gene, set, HR = round(HR, 3), p = signif(p, 3))])

cat("\n=== ROBUSTNESS — Scoring-method sensitivity (polar, both endpoints) ===\n")
print(method_dt[, .(endpoint, method, n,
                    HR = round(HR, 3),
                    CI = sprintf("[%.2f, %.2f]", HR_lo, HR_hi),
                    p = signif(p, 3))])

sink_log("Phase 40b robustness complete")
