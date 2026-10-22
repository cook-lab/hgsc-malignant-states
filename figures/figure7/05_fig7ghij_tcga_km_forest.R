#!/usr/bin/env Rscript
# ============================================================================
# Figure 7g-7j — TCGA-OV KM curves + stepwise Cox forests (one combined SVG)
# ----------------------------------------------------------------------------
# PURPOSE
#   TCGA-OV survival validation of the secretory polarization axis:
#     KM curves, OS (7g) and PFS (7h), polarization tertiles (SecB-like "High"
#       vs SecA-like "Low"; mid tertile dropped), 5-year censored, log-rank.
#     Stepwise Cox forests (7i OS, 7j PFS): unadjusted -> + epithelial fraction ->
#       + stage + age; polarization inverted so HR>1 = SecB-enriched = worse.
#       (Platinum intentionally omitted to retain the full cohort.)
#
# INPUTS
#   data_root/2026_final_atlas/output/22_tcga_deconvolution/
#     22d_signature_scores.csv  (per-patient scores + clinical + polar_tertile)
#   Shared helpers: config/config.R, spatial/00_setup/00_setup.R (ref_palette, theme_lab).
#
# OUTPUTS
#   figures_dir/figure7/tcga_survival_row_EFG.{svg,png}  (Fig 7g-j; the stem keeps earlier letters)
#   figures_dir/figure7/tcga_forest_stepwise_data.csv     (HR, CI, p, BH q, n, events per model and endpoint)
#
# MANUSCRIPT PANEL(S): Fig 7g (KM OS), 7h (KM PFS), 7i/7j (stepwise Cox forests, OS/PFS)
# RUNTIME TIER: fast (per-patient table; in-script Cox)
#
# NOTE: the Cox models are fitted here, from the scores, clinical covariates and BayesPrism epithelial
#   fraction in 22d_signature_scores.csv.
# ============================================================================

.here     <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))
source(file.path(.here, "..", "..", "spatial", "00_setup", "00_setup.R"))
source(file.path(.here, "_km_forest_helpers.R"))   # km_panel(), forest_theme(), save_survival_row()

suppressPackageStartupMessages({
  library(survival)
  library(patchwork)
})

FIG_DIR <- cfg_path("figures_dir", "figure7")

# Figure 7 house typography (shared with the helpers): 8pt labels, 6pt ticks.
FA <- F7_FA; FK <- F7_FK; FN <- F7_FA; FL <- F7_FL

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
# 2. KM CURVES (Fig 7g / 7h)
# =============================================================================
message("\nPlotting KM curves...")

km_pal <- c(unname(ref_palette["SecA epithelium"]),
            unname(ref_palette["SecB epithelium"]))

# Native KM panel (geom_step via km_panel(); branding matches Figs 1-3). Second
# panel drops the y axis — it shares the first panel's 0-1 survival scale.
make_km <- function(dat, time_col, event_col, endpoint_label, show_y = TRUE) {
  dat$.time  <- dat[[time_col]]
  dat$.event <- dat[[event_col]]
  fit <- survfit(Surv(.time, .event) ~ group, data = dat)   # strata: Low, High

  # Key n = patients analysed for this endpoint (survfit drops missing time/event), not the
  # tertile size: OS 104 / 106, PFS 104 / 105 of the 108 / 108 tertile patients.
  ok     <- !is.na(dat$.time) & !is.na(dat$.event)
  n_low  <- sum(dat$group[ok] == "Low")
  n_high <- sum(dat$group[ok] == "High")

  sd   <- survdiff(Surv(.time, .event) ~ group, data = dat)
  p_lr <- 1 - pchisq(sd$chisq, 1)
  p_text <- ifelse(p_lr < 0.001, "p < 0.001", sprintf("p = %.3f", p_lr))

  km_panel(fit,
           cols = km_pal,                                   # SecA gold, SecB brown
           labs = c(sprintf("SecA-dominant (n=%d)", n_low),
                    sprintf("SecB-enriched (n=%d)", n_high)),
           p_text = p_text,
           subtitle = NULL,   # OS/PFS labelled in Affinity; endpoint also on the forest x-axis
           show_y = show_y)
}

p_km_os  <- make_km(scores_binary, "os_time_5y",  "os_event_5y",  "OS",  show_y = TRUE)
p_km_pfs <- make_km(scores_binary, "pfs_time_5y", "pfs_event_5y", "PFS", show_y = FALSE)

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
# 4. FOREST PLOTS (Fig 7i OS, 7j PFS) — stacked OS over PFS
# =============================================================================
message("\nPlotting forest...")

model_levels <- c("Polarization score",
                   "+ Epithelial fraction",
                   "+ Stage + Age")

forest_df$model <- factor(forest_df$model, levels = model_levels)
# HR value only (95% CI dropped from the text; error bars show the interval).
forest_df$hr_text <- sprintf("%.2f", forest_df$HR)
forest_df$p_text  <- ifelse(forest_df$p_value < 0.001, "<.001",
                            sub("^0\\.", ".", formatC(forest_df$p_value, format = "f", digits = 3)))
# q = Benjamini-Hochberg FDR across the 6 tests shown (3 models x OS/PFS), as in Fig 7c,d and 7m,n
forest_df$q_value <- p.adjust(forest_df$p_value, method = "BH")
forest_df$q_text  <- ifelse(forest_df$q_value < 0.001, "<.001",
                            sub("^0\\.", ".", formatC(forest_df$q_value, format = "f", digits = 3)))

secb_col <- unname(ref_palette["SecB epithelium"])
forest_df$pt_color <- secb_col

# Tight HR / p / q columns just right of the forest (no CI text) so the data region fills the panel.
G <- forest_geometry(lo = 0.45, data_hi = 3.4)   # shared Fig 7 forest layout (_km_forest_helpers.R)
x_hr_pos <- G$x[["HR"]]
x_p_pos  <- G$x[["p"]]
x_q_pos  <- G$x[["q"]]
band_hi  <- G$limits[2] * (G$limits[2] / G$limits[1])^0.05   # highlight band to the panel's right edge
ypad     <- forest_ylab_pad(model_levels)
n_feat   <- length(model_levels)

plot_tcga_forest <- function(fd, endpoint_filter, endpoint_label) {
  fd <- fd[fd$endpoint == endpoint_filter, ]

  p <- ggplot(fd, aes(x = HR, y = model))

  p <- p +
    annotate("rect",
             xmin = 0.5, xmax = band_hi,
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
              hjust = 0, size = FA / .pt, color = "black") +
    geom_text(aes(x = x_p_pos, label = p_text),
              hjust = 0, size = FA / .pt, color = "black") +
    geom_text(aes(x = x_q_pos, label = q_text),
              hjust = 0, size = FA / .pt, color = "black") +
    annotate("text", x = x_hr_pos, y = n_feat + 0.8,
             label = "HR", hjust = 0,
             size = FA / .pt, fontface = "bold") +
    annotate("text", x = x_p_pos, y = n_feat + 0.8,
             label = "p", hjust = 0,
             size = FA / .pt, fontface = "bold") +
    annotate("text", x = x_q_pos, y = n_feat + 0.8,
             label = "q", hjust = 0,
             size = FA / .pt, fontface = "bold") +
    # Scale limit keeps the HR/p/q text columns (right of the forest) on-panel;
    # breaks emphasise the 0.5-3.0 range.
    scale_x_log10(breaks = c(0.5, 1, 2, 3),
                  labels = c("0.5", "1.0", "2.0", "3.0"),
                  limits = G$limits, oob = scales::oob_keep) +
    scale_y_discrete(drop = FALSE) +
    coord_cartesian(clip = "off") +
    labs(x = paste0("Hazard Ratio - ", endpoint_label), y = NULL) +
    forest_theme() + theme(plot.margin = margin(10, 14, 4, 4 + ypad))

  p
}

p_forest_os  <- plot_tcga_forest(forest_df, "OS",  "5-yr OS")
p_forest_pfs <- plot_tcga_forest(forest_df, "PFS", "5-yr PFS")

# Forest for the right column = OS rows over PFS rows (one visual column).
p_forest <- p_forest_os / p_forest_pfs

# =============================================================================
# 5. SAVE — ONE combined SVG for the whole survival row (7g | 7h | 7i/7j).
#    No sub-panel letters (added in Affinity). Same call builds the TMA (a-d) and ICON7 (k-n)
#    rows, so the survival rows are pixel-for-pixel consistent.
# =============================================================================
save_survival_row(p_km_os, p_km_pfs, p_forest,
                  stem = "tcga_survival_row_EFG", fig_dir = FIG_DIR,
                  width_in = 7.2, height_in = 2.6, widths = c(1, 1, 1.7))

write.csv(forest_df, file.path(FIG_DIR, "tcga_forest_stepwise_data.csv"),
          row.names = FALSE)

message("  Saved: forest data CSV")
message("\nDone.")
