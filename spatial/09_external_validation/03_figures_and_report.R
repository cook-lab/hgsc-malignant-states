# ============================================================================
# 03_figures_and_report.R
# ----------------------------------------------------------------------------
# PURPOSE: Figures and report for the TCGA external validation (discovery vs validation effect sizes, scoring-method sensitivity).
#
# INPUTS:
#   - output/40_tcga_validation/ Cox + robustness tables
#
# OUTPUTS:
#   - output/40_tcga_validation/figures/ + report
#
# MANUSCRIPT PANEL(S): Supports Fig 7E/F/G validation narrative.
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
  library(survival)
  library(survminer)
  library(maxstat)
})

step_dir <- file.path(out_dir, "40_tcga_validation")
fig_dir40 <- file.path(step_dir, "figures")
dir.create(fig_dir40, showWarnings = FALSE, recursive = TRUE)

run_log <- file.path(step_dir, "run.log")
sink_log <- function(msg) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg, "\n",
                              file = run_log, append = TRUE)
sink_log("Phase 40c figures started")

phase40 <- readRDS(file.path(step_dir, "phase40_primary.rds"))
robust  <- if (file.exists(file.path(step_dir, "phase40_robustness.rds")))
              readRDS(file.path(step_dir, "phase40_robustness.rds")) else NULL
cohort  <- phase40$cohort
fam     <- phase40$fam
primary <- phase40$primary_dt

# Discovery effects
disc <- fread(file.path(out_dir, "10_clinical_v2", "cox_univariate_results.csv"))
disc[, endpoint_short := fifelse(endpoint == "OS_5yr", "OS5", "PFS5")]

# --- Fig40A: CONSORT-style cohort flow ---------------------------------------
flow <- fread(file.path(step_dir, "cohort_flowchart.csv"))
flow[, step := factor(step, levels = step)]

p_consort <- ggplot(flow, aes(x = 1, y = -as.integer(step))) +
  geom_tile(aes(width = 0.6, height = 0.7), fill = "#DDD5CA", color = "black") +
  geom_text(aes(label = sprintf("%s\nn = %d", step, n)), size = 2.6) +
  scale_y_continuous(expand = c(0.05, 0.05)) +
  scale_x_continuous(expand = c(0.5, 0.5)) +
  theme_void() +
  ggtitle("Phase 40 — TCGA-OV cohort lock") +
  theme(plot.title = element_text(size = 10, hjust = 0.5))
ggsave(file.path(fig_dir40, "fig40A_consort_flow.png"),
       p_consort, width = 5, height = 0.5 + 0.7 * nrow(flow), units = "in",
       dpi = 300)
ggsave(file.path(fig_dir40, "fig40A_consort_flow.pdf"),
       p_consort, width = 5, height = 0.5 + 0.7 * nrow(flow))

# --- Fig40B: side-by-side forest plot ----------------------------------------
# Map TCGA score names to discovery feature names where comparable
disc_map <- data.table(
  tcga_score = c("polar_UCell_z", "log2_ratio_UCell_z",
                 "prop_SecB_of_sec_z", "log2_ratio_csx_z"),
  disc_feature = c("polar_mean", "log2_ratio",
                   "prop_SecB_of_sec", "log2_ratio"),
  pretty = c("Polarization (SecB - SecA, UCell)",
             "log2(SecA / SecB), UCell",
             "Epi_SecB / (Epi_SecA + Epi_SecB), CIBERSORTx",
             "log2(Epi_SecA / Epi_SecB), CIBERSORTx"))

# Build a single long forest table: discovery (Xenium) + TCGA primary (age+stage)
# + TCGA kitchen-sink
tcga_age_stage <- primary[model == "age_stage",
                          .(score, endpoint, HR, HR_lo, HR_hi)]
tcga_age_stage[, source := "TCGA (age+stage)"]
tcga_ks <- primary[model == "kitchen_sink",
                   .(score, endpoint, HR, HR_lo, HR_hi)]
tcga_ks[, source := "TCGA (kitchen-sink)"]
tcga_long <- rbind(tcga_age_stage, tcga_ks)
tcga_long[, score_pretty := disc_map$pretty[match(score, disc_map$tcga_score)]]
tcga_long[, endpoint_pretty := fifelse(endpoint == "OS5",
                                       "5-yr OS", "5-yr PFS")]

disc_pull <- function(feat, ep) {
  r <- disc[feature == feat & endpoint_short == ep]
  if (!nrow(r)) return(NULL)
  r[1, .(HR, HR_lo = HR_lower, HR_hi = HR_upper, p_value)]
}
disc_rows <- list()
for (i in seq_len(nrow(disc_map))) {
  for (ep in c("OS5", "PFS5")) {
    rr <- disc_pull(disc_map$disc_feature[i], ep)
    if (is.null(rr)) next
    disc_rows[[length(disc_rows) + 1]] <- data.table(
      score = disc_map$tcga_score[i],
      endpoint = ep,
      HR = rr$HR, HR_lo = rr$HR_lo, HR_hi = rr$HR_hi,
      source = "Xenium (discovery)",
      score_pretty = disc_map$pretty[i],
      endpoint_pretty = fifelse(ep == "OS5", "5-yr OS", "5-yr PFS"))
  }
}
disc_long <- rbindlist(disc_rows, fill = TRUE)

forest_dt <- rbind(disc_long, tcga_long, fill = TRUE)
forest_dt[, source := factor(source,
                             levels = c("Xenium (discovery)",
                                        "TCGA (age+stage)",
                                        "TCGA (kitchen-sink)"))]
forest_dt[, score_pretty := factor(score_pretty, levels = disc_map$pretty)]

cohort_pal <- c("Xenium (discovery)" = "#9A7D55",
                "TCGA (age+stage)"   = "#56AFC4",
                "TCGA (kitchen-sink)" = "#5665B6")

p_forest <- ggplot(forest_dt,
                   aes(y = source, x = HR, xmin = HR_lo, xmax = HR_hi,
                       color = source)) +
  geom_vline(xintercept = 1, color = "grey60", linetype = "dashed") +
  geom_pointrange(size = 0.4, fatten = 2) +
  facet_grid(score_pretty ~ endpoint_pretty, scales = "free_y", switch = "y") +
  scale_x_log10(breaks = c(0.5, 0.7, 1, 1.5, 2, 3),
                labels = c("0.5", "0.7", "1", "1.5", "2", "3")) +
  scale_color_manual(values = cohort_pal) +
  labs(x = "HR (log scale)", y = NULL,
       title = "Phase 40 — SecA/SecB polarization signature: discovery vs TCGA validation",
       color = NULL) +
  theme_lab(base_size = 9) +
  theme(legend.position = "bottom",
        strip.text.y.left = element_text(angle = 0, hjust = 1),
        strip.placement = "outside")
ggsave(file.path(fig_dir40, "fig40B_forest_xenium_vs_tcga.png"),
       p_forest, width = 9, height = 6, dpi = 300)
ggsave(file.path(fig_dir40, "fig40B_forest_xenium_vs_tcga.pdf"),
       p_forest, width = 9, height = 6)

# --- Fig40C: KM with optimal cutpoint (Lausen-Schumacher) --------------------
make_km <- function(score_col, endpoint_lbl) {
  tcol <- paste0(endpoint_lbl, "_months")
  ecol <- paste0(endpoint_lbl, "_event")
  df <- as.data.frame(cohort[!is.na(get(score_col)) & !is.na(get(tcol)) & !is.na(get(ecol))])
  df$.time   <- df[[tcol]]
  df$.event  <- df[[ecol]]
  df$.score  <- df[[score_col]]
  ms <- maxstat.test(Surv(.time, .event) ~ .score,
                     data = df, smethod = "LogRank",
                     pmethod = "Lau94", abseps = 0.01)
  cut_val <- as.numeric(ms$estimate)
  df$group <- factor(ifelse(df$.score >= cut_val, "high", "low"),
                     levels = c("low", "high"))
  fit <- survfit(Surv(.time, .event) ~ group, data = df)
  surv_diff <- survdiff(Surv(.time, .event) ~ group, data = df)
  p_lr <- 1 - pchisq(surv_diff$chisq, df = 1)

  ggsurvplot(fit, data = df, risk.table = TRUE, pval = TRUE,
             conf.int = FALSE,
             palette = c("#E6A141", "#9A7D55"),
             legend.labs = c(sprintf("low (n=%d)", sum(df$group == "low")),
                             sprintf("high (n=%d)", sum(df$group == "high"))),
             xlab = sprintf("Months (%s, capped at 60)", endpoint_lbl),
             ylab = "Survival probability",
             title = sprintf("%s by polar_UCell (cutpoint=%.3f, Lau94 p_LR=%.3g)",
                             endpoint_lbl, cut_val, p_lr),
             xlim = c(0, 60), break.time.by = 12,
             ggtheme = theme_lab(base_size = 9))
}

km_os  <- make_km("polar_UCell_z", "OS5")
km_pfs <- make_km("polar_UCell_z", "PFS5")
save_km <- function(km, path) {
  png(path, width = 6 * 300, height = 6 * 300, res = 300)
  print(km)
  dev.off()
}
save_km(km_os,  file.path(fig_dir40, "fig40C_KM_polar_OS5.png"))
save_km(km_pfs, file.path(fig_dir40, "fig40C_KM_polar_PFS5.png"))

# --- Fig40D: c-index Δ ------------------------------------------------------
if (!is.null(robust)) {
  cdt <- robust$cindex
  cdt[, score_pretty := disc_map$pretty[match(score, disc_map$tcga_score)]]
  cdt[, endpoint_pretty := fifelse(endpoint == "OS5", "5-yr OS", "5-yr PFS")]
  p_cidx <- ggplot(cdt, aes(x = score_pretty, y = delta_c,
                            ymin = delta_c_boot_lo, ymax = delta_c_boot_hi)) +
    geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
    geom_pointrange(color = "#5665B6") +
    coord_flip() +
    facet_wrap(~ endpoint_pretty, ncol = 2) +
    labs(x = NULL, y = "Δ c-index (full vs age+stage), with 95% bootstrap CI",
         title = "Added prognostic discrimination from polarization signature") +
    theme_lab(base_size = 9)
  ggsave(file.path(fig_dir40, "fig40D_cindex_delta.png"),
         p_cidx, width = 8, height = 4, dpi = 300)
  ggsave(file.path(fig_dir40, "fig40D_cindex_delta.pdf"),
         p_cidx, width = 8, height = 4)

  # LOGO summary figure
  logo <- robust$logo[endpoint == "OS5"]
  logo[, label := ifelse(dropped_gene == "(none)", "FULL (14 genes)", dropped_gene)]
  logo_full_HR <- logo[dropped_gene == "(none)", HR]
  p_logo <- ggplot(logo, aes(x = reorder(label, HR), y = HR,
                             color = set, fill = set)) +
    geom_hline(yintercept = 1, color = "grey60", linetype = "dashed") +
    geom_hline(yintercept = logo_full_HR, color = "black", linetype = "dotted") +
    geom_point(size = 3, shape = 21) +
    coord_flip() +
    scale_color_manual(values = c("SecA" = "#E6A141", "SecB" = "#9A7D55",
                                   "(full)" = "black")) +
    scale_fill_manual(values = c("SecA" = "#E6A141", "SecB" = "#9A7D55",
                                  "(full)" = "black")) +
    labs(x = NULL, y = "HR for polar_UCell (5-yr OS, age+stage adjusted)",
         title = "Leave-one-gene-out sensitivity",
         subtitle = sprintf("Dotted line = full-signature HR (%.3f)", logo_full_HR)) +
    theme_lab(base_size = 9)
  ggsave(file.path(fig_dir40, "fig40E_LOGO.png"),
         p_logo, width = 7, height = 5, dpi = 300)
  ggsave(file.path(fig_dir40, "fig40E_LOGO.pdf"),
         p_logo, width = 7, height = 5)

  # Scoring-method sensitivity figure
  meth <- robust$methods
  meth[, endpoint_pretty := fifelse(endpoint == "OS5", "5-yr OS", "5-yr PFS")]
  p_meth <- ggplot(meth, aes(x = method, y = HR, ymin = HR_lo, ymax = HR_hi,
                             color = method)) +
    geom_hline(yintercept = 1, color = "grey60", linetype = "dashed") +
    geom_pointrange(size = 0.5) +
    facet_wrap(~ endpoint_pretty, ncol = 2) +
    scale_y_log10() +
    scale_color_manual(values = c(UCell = "#5665B6", singscore = "#56AFC4",
                                   GSVA = "#8FBC8F", ssGSEA = "#E6A141")) +
    labs(x = NULL, y = "HR per SD (log scale)",
         title = "Scoring-method sensitivity: polarization (SecB - SecA)") +
    theme_lab(base_size = 9) +
    theme(legend.position = "none")
  ggsave(file.path(fig_dir40, "fig40F_scoring_methods.png"),
         p_meth, width = 7, height = 4, dpi = 300)
}

# --- Numerical results table -------------------------------------------------
out_tab <- primary[, .(endpoint, score, model,
                       n, n_event,
                       HR = round(HR, 3),
                       CI_lo = round(HR_lo, 3),
                       CI_hi = round(HR_hi, 3),
                       p = signif(p, 3),
                       p_BH = signif(p_bh, 3),
                       p_Bonf = signif(p_bonferroni, 3),
                       expected_direction, direction_match,
                       ci_excludes_1, replicated)]
fwrite(out_tab, file.path(step_dir, "results_table.csv"))

# --- Replication call summary -----------------------------------------------
n_sign_concordant <- sum(fam$direction_match, na.rm = TRUE)
n_total           <- nrow(fam)
sign_test_p       <- pbinom(n_sign_concordant - 1, n_total, 0.5,
                            lower.tail = FALSE)
n_unadj_sig_dir   <- sum(fam$direction_match & fam$p < 0.05, na.rm = TRUE)
n_replicated      <- sum(fam$replicated, na.rm = TRUE)

summary_lines <- c(
  sprintf("Phase 40 — TCGA-OV external validation summary"),
  sprintf("Cohort: n=%d patients, %d OS5 events, %d PFS5 events",
          nrow(cohort), sum(cohort$OS5_event), sum(cohort$PFS5_event)),
  sprintf("Primary tests (model = age+stage): %d", n_total),
  sprintf("HR direction concordant with discovery: %d/%d (sign-test p=%.4g)",
          n_sign_concordant, n_total, sign_test_p),
  sprintf("Direction-concordant AND unadjusted p<0.05: %d/%d", n_unadj_sig_dir, n_total),
  sprintf("Replicated (direction + Bonferroni p<0.05 + CI excludes 1): %d/%d",
          n_replicated, n_total),
  sprintf("Replication call: %s", phase40$call))
writeLines(summary_lines, file.path(step_dir, "replication_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
sink_log("Phase 40c figures complete")
