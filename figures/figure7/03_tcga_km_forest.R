#!/usr/bin/env Rscript
# ============================================================================
# Figure 7E,7F,7G — TCGA-OV KM curves + stepwise Cox forest
# ----------------------------------------------------------------------------
# PURPOSE
#   TCGA-OV survival validation of the secretory polarization axis:
#     KM curves, OS (7E) and PFS (7F), polarization tertiles (SecB-like "High"
#       vs SecA-like "Low"; mid tertile dropped), 5-year censored, log-rank.
#     Stepwise Cox forest (7G): unadjusted -> + epithelial fraction ->
#       + stage + age; polarization inverted so HR>1 = SecB-enriched = worse.
#       (Platinum intentionally omitted to retain the full cohort.)
#
# INPUTS
#   data_root/2026_final_atlas/output/22_tcga_deconvolution/
#     22d_signature_scores.csv  (per-patient scores + clinical + polar_tertile)
#   Shared helpers: config/config.R, spatial/00_setup/00_setup.R (ref_palette, theme_lab).
#
# OUTPUTS
#   figures_dir/figure7/tcga_km_os.{svg,pdf}            (7E)
#   figures_dir/figure7/tcga_km_pfs.{svg,pdf}           (7F)
#   figures_dir/figure7/tcga_forest_stepwise.{svg,pdf}  (7G)
#   figures_dir/figure7/tcga_forest_stepwise_data.csv
#
# MANUSCRIPT PANEL(S): Fig 7E (KM OS), Fig 7F (KM PFS), Fig 7G (stepwise Cox forest)
# RUNTIME TIER: fast (per-patient table; in-script Cox)
#
# NOTE: BayesPrism multivariate covariate set is a documented discrepancy
#   (docs/REPRODUCIBILITY.md); the in-script Cox is migrated faithfully.
# ============================================================================

.here     <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))
source(file.path(.here, "..", "..", "spatial", "00_setup", "00_setup.R"))

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(patchwork)
})

FIG_DIR <- cfg_path("figures_dir", "figure7")

FA <- 6; FK <- 5.5; FN <- 5

# =============================================================================
# 1. LOAD DATA
# =============================================================================
message("=== Loading atlas 22d TCGA data ===")

scores <- read.csv(cfg_path("data_root", "2026_final_atlas", "output",
                            "22_tcga_deconvolution", "22d_signature_scores.csv"),
                   stringsAsFactors = FALSE)

message("  Patients with scores: ", nrow(scores))
message("  Tertile distribution:")
print(table(scores$polar_tertile, useNA = "ifany"))

# Cap at 5 years
scores$os_time_5y  <- pmin(scores$os_months, 60)
scores$os_event_5y <- ifelse(scores$os_event == 1 & scores$os_months <= 60, 1, 0)
scores$pfs_time_5y <- pmin(scores$pfs_months, 60)
scores$pfs_event_5y <- ifelse(scores$pfs_event == 1 & scores$pfs_months <= 60, 1, 0)

# Binary: drop mid tertile
scores_binary <- scores[scores$polar_tertile %in% c("SecA-like", "SecB-like"), ]
scores_binary$group <- factor(
  ifelse(scores_binary$polar_tertile == "SecB-like", "High", "Low"),
  levels = c("Low", "High")
)

message("  Binary KM cohort: n = ", nrow(scores_binary),
        " (Low=", sum(scores_binary$group == "Low"),
        ", High=", sum(scores_binary$group == "High"), ")")

# =============================================================================
# 2. KM CURVES (Fig 7E / 7F)
# =============================================================================
message("\nPlotting KM curves...")

km_pal <- c(unname(ref_palette["SecA epithelium"]),
            unname(ref_palette["SecB epithelium"]))

make_km <- function(dat, time_col, event_col, endpoint_label) {
  dat$.time  <- dat[[time_col]]
  dat$.event <- dat[[event_col]]
  fit <- survfit(Surv(.time, .event) ~ group, data = dat)

  n_low  <- sum(dat$group == "Low")
  n_high <- sum(dat$group == "High")

  sd <- survdiff(Surv(.time, .event) ~ group, data = dat)
  p_lr <- 1 - pchisq(sd$chisq, 1)

  ggsurvplot(fit, data = dat,
             palette = km_pal,
             legend.labs = c(sprintf("SecA-dominant (n=%d)", n_low),
                             sprintf("SecB-enriched (n=%d)", n_high)),
             legend.title = "",
             risk.table = FALSE,
             conf.int = FALSE,
             pval = FALSE,
             xlab = "Months",
             ylab = "Survival probability",
             xlim = c(0, 60),
             break.time.by = 12,
             ggtheme = theme_lab(base_size = 6) +
               theme(legend.position = c(0.65, 0.92),
                     legend.text = element_text(size = FN),
                     legend.key.size = unit(0.3, "cm"),
                     legend.background = element_blank(),
                     plot.margin = margin(4, 6, 4, 4)),
             surv.median.line = "none")$plot +
    annotate("text", x = 2, y = 0.05,
             label = sprintf("p = %.3f", p_lr),
             hjust = 0, size = FN / .pt, color = "grey30") +
    labs(subtitle = paste0("5-yr ", endpoint_label))
}

p_km_os  <- make_km(scores_binary, "os_time_5y",  "os_event_5y",  "OS")
p_km_pfs <- make_km(scores_binary, "pfs_time_5y", "pfs_event_5y", "PFS")

# =============================================================================
# 3. STEPWISE COX (5-year, no platinum)
# =============================================================================
message("\nComputing stepwise Cox models...")

scores$polar_inv <- -scores$polarization

run_cox <- function(dat, time_col, event_col, formula_str) {
  dat$.time  <- dat[[time_col]]
  dat$.event <- dat[[event_col]]
  f <- as.formula(paste0("Surv(.time, .event) ~ ", formula_str))
  fit <- coxph(f, data = dat)
  s <- summary(fit)
  coef_row <- s$conf.int["polar_inv", , drop = FALSE]
  p_row    <- s$coefficients["polar_inv", "Pr(>|z|)"]
  data.frame(
    HR       = coef_row[1, "exp(coef)"],
    HR_lower = coef_row[1, "lower .95"],
    HR_upper = coef_row[1, "upper .95"],
    p_value  = p_row,
    n        = s$n,
    events   = s$nevent,
    stringsAsFactors = FALSE
  )
}

forest_list <- list()
for (ep in c("OS", "PFS")) {
  tc <- ifelse(ep == "OS", "os_time_5y", "pfs_time_5y")
  ec <- ifelse(ep == "OS", "os_event_5y", "pfs_event_5y")

  r1 <- run_cox(scores, tc, ec, "polar_inv")
  r1$endpoint <- ep; r1$model <- "Polarization score"
  forest_list[[length(forest_list) + 1]] <- r1

  r2 <- run_cox(scores, tc, ec, "polar_inv + epi_fraction")
  r2$endpoint <- ep; r2$model <- "+ Epithelial fraction"
  forest_list[[length(forest_list) + 1]] <- r2

  r3 <- run_cox(scores, tc, ec, "polar_inv + epi_fraction + stage_coded + age")
  r3$endpoint <- ep; r3$model <- "+ Stage + Age"
  forest_list[[length(forest_list) + 1]] <- r3
}

forest_df <- do.call(rbind, forest_list)
rownames(forest_df) <- NULL

message("\nStepwise Cox (SecB-enriched direction, HR>1 = worse):")
print(forest_df[, c("endpoint", "model", "n", "events", "HR", "HR_lower", "HR_upper", "p_value")])

# =============================================================================
# 4. FOREST PLOT (Fig 7G) — side-by-side OS | PFS
# =============================================================================
message("\nPlotting forest...")

model_levels <- c("Polarization score",
                   "+ Epithelial fraction",
                   "+ Stage + Age")

forest_df$model <- factor(forest_df$model, levels = model_levels)
forest_df$hr_text <- sprintf("%.2f (%.2f-%.2f)",
                             forest_df$HR, forest_df$HR_lower, forest_df$HR_upper)
forest_df$p_text  <- ifelse(forest_df$p_value < 0.001, "<.001",
                            sub("^0\\.", ".", formatC(forest_df$p_value, format = "f", digits = 3)))

secb_col <- unname(ref_palette["SecB epithelium"])
forest_df$pt_color <- secb_col

x_hr_pos <- 4.0
x_p_pos  <- 12.0
n_feat   <- length(model_levels)

plot_tcga_forest <- function(fd, endpoint_filter, endpoint_label) {
  fd <- fd[fd$endpoint == endpoint_filter, ]

  p <- ggplot(fd, aes(x = HR, y = model))

  p <- p +
    annotate("rect",
             xmin = 0.5, xmax = 25,
             ymin = 0.55, ymax = 1.45,
             fill = secb_col, alpha = 0.08)

  p <- p +
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "grey50", linewidth = 0.3) +
    geom_errorbar(aes(xmin = HR_lower, xmax = HR_upper, color = pt_color),
                  width = 0.2, linewidth = 0.35, show.legend = FALSE,
                  orientation = "y") +
    geom_point(aes(fill = pt_color, color = pt_color), size = 2.2, shape = 21,
               stroke = 0.5, alpha = 0.7, show.legend = FALSE) +
    geom_point(aes(color = pt_color), size = 2.2, shape = 1,
               stroke = 0.5, show.legend = FALSE) +
    scale_color_identity() +
    scale_fill_identity() +
    geom_text(aes(x = x_hr_pos, label = hr_text),
              hjust = 0, size = FN / .pt, color = "black") +
    geom_text(aes(x = x_p_pos, label = p_text),
              hjust = 0, size = FN / .pt, color = "black") +
    annotate("text", x = x_hr_pos, y = n_feat + 0.8,
             label = "HR (95% CI)", hjust = 0,
             size = FN / .pt, fontface = "bold") +
    annotate("text", x = x_p_pos, y = n_feat + 0.8,
             label = "p", hjust = 0,
             size = FN / .pt, fontface = "bold") +
    scale_x_log10(breaks = c(0.5, 1, 2, 3),
                  labels = c("0.5", "1.0", "2.0", "3.0"),
                  limits = c(0.5, 3.5)) +
    scale_y_discrete(drop = FALSE) +
    coord_cartesian(clip = "off") +
    labs(x = paste0("Hazard Ratio - ", endpoint_label), y = NULL) +
    theme_lab(base_size = 6) +
    theme(
      plot.title         = element_blank(),
      axis.text.y        = element_text(size = FK),
      axis.text.x        = element_text(size = FK),
      axis.title.x       = element_text(size = FA),
      panel.grid.major.y = element_line(color = "grey93", linewidth = 0.15),
      plot.margin        = margin(10, 60, 4, 4)
    )

  p
}

p_forest_os  <- plot_tcga_forest(forest_df, "OS",  "5-yr OS")
p_forest_pfs <- plot_tcga_forest(forest_df, "PFS", "5-yr PFS")

p_forest <- p_forest_os / p_forest_pfs

# =============================================================================
# 5. SAVE
# =============================================================================
for (ext in c("svg", "pdf")) {
  ggsave(file.path(FIG_DIR, paste0("tcga_km_os.", ext)),
         p_km_os, width = 3.5, height = 3.5, bg = "white")
  ggsave(file.path(FIG_DIR, paste0("tcga_km_pfs.", ext)),
         p_km_pfs, width = 3.5, height = 3.5, bg = "white")
}

for (ext in c("svg", "pdf")) {
  ggsave(file.path(FIG_DIR, paste0("tcga_forest_stepwise.", ext)),
         p_forest, width = 4, height = 2, bg = "white")
}

write.csv(forest_df, file.path(FIG_DIR, "tcga_forest_stepwise_data.csv"),
          row.names = FALSE)

message("\n  Saved: tcga_km_os.svg")
message("  Saved: tcga_km_pfs.svg")
message("  Saved: tcga_forest_stepwise.svg")
message("  Saved: forest data CSV")
message("\nDone.")
