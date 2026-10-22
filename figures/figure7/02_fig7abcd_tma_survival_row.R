#!/usr/bin/env Rscript
# ============================================================================
# Figure 7a-7d — Xenium TMA survival row (one combined SVG), same layout as the TCGA and ICON7 rows
# ----------------------------------------------------------------------------
# PURPOSE
#   Primary-endpoint view of secretory polarization in our TMA cohort (patients with primary,
#   treatment-naive HGSC represented by two tumour cores; 01_fig7abcd_tma_cox_km_tables.R):
#     7a  5-yr OS KM, SecA-dominant (lowest polarization tertile) vs SecB-enriched (highest tertile)
#         patients of the analysed cohort, middle tertile dropped, log-rank p: the design of the
#         TCGA row (7g,h).
#     7b  5-yr PFS KM, same split
#     7c  OS forest of the three secretory measures (SecA proportion, SecB proportion, polarization):
#         multivariable Cox HR per +1 SD (adjusted for stage, residual disease and age), p and q
#     7d  PFS forest, same measures; q = Benjamini-Hochberg FDR across the 6 tests of 7c + 7d
#   The HRs are the models in cox_2core_results.csv / Table S11 (the three secretory measures).
#   KM and Cox use the same complete-case cohort (final TMA cores: OS n = 104, PFS n = 103);
#   the KM compares its outer tertiles (OS 35 vs 35, PFS 35 vs 34; tertiles of each endpoint's cohort).
#   Same layout as the TCGA (05) and ICON7 (06) rows: km_panel() / forest_theme() / save_survival_row()
#   plus the shared forest geometry (forest_geometry(), forest_ylab_pad()), group names, p format and
#   stacked forests; the OS/PFS labels over the KM panels are added in Affinity, as for the TCGA row.
#
# INPUTS
#   figures_dir/figure7/tma_2core_patients.csv   (01_fig7abcd_tma_cox_km_tables.R; per-patient measures + clinical)
#   figures_dir/figure7/cox_2core_results.csv    (01_fig7abcd_tma_cox_km_tables.R; adjusted HRs)
#   Shared helpers: config/config.R, spatial/00_setup/00_setup.R, _km_forest_helpers.R
#
# OUTPUTS
#   figures_dir/figure7/xenium_tma_survival_row_ABC.{svg,png}   (Fig 7a-d; the stem keeps earlier letters)
#   figures_dir/figure7/xenium_tma_km_logrank.csv        (per endpoint: cohort n, tertile cuts, n and events per group,
#                                                         5-yr KM survival, log-rank p; also the 3-tertile log-rank
#                                                         and trend p, source data only)
#   figures_dir/figure7/xenium_tma_forest_fdr.csv        (6 rows: HR, CI, p, BH q)
#   (the two CSVs are copied to final_panels as Fig7_ABCD_km_logrank.csv / Fig7_ABCD_forest_fdr.csv, which Table S11
#    checks against)
#
# MANUSCRIPT PANEL(S): Fig 7a, 7b, 7c, 7d
# RUNTIME TIER: fast
# ============================================================================
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))
source(file.path(.here, "..", "..", "spatial", "00_setup", "00_setup.R"))
source(file.path(.here, "_km_forest_helpers.R"))
suppressPackageStartupMessages({ library(data.table); library(survival); library(patchwork) })

FIG_DIR <- cfg_path("figures_dir", "figure7")
FA <- F7_FA
COL_SECA <- unname(ref_palette["SecA epithelium"])
COL_SECB <- unname(ref_palette["SecB epithelium"])
COL_POL  <- "#4D4D4D"
fmt_p <- function(x) ifelse(x < 0.001, "<.001", sub("^0\\.", ".", formatC(x, format = "f", digits = 3)))

pt <- fread(file.path(FIG_DIR, "tma_2core_patients.csv"))
cox <- fread(file.path(FIG_DIR, "cox_2core_results.csv"))

# ---- 7a / 7b: KM of the outer polarization tertiles of each endpoint's analysed (complete-case) cohort ----
km_rows <- list()
make_km <- function(ep, show_y) {
  tcol <- if (ep == "OS") "os_time" else "pfs_time"; ecol <- if (ep == "OS") "os_event" else "pfs_event"
  d <- pt[is.finite(get(tcol)) & get(tcol) > 0 & !is.na(get(ecol)) & !is.na(stage_late) &
          debulk %in% c("optimal", "suboptimal") & !is.na(age) & is.finite(polar_mean)]   # = the Cox cohort
  n_cohort <- nrow(d)
  cut3 <- quantile(d$polar_mean, c(1/3, 2/3))
  # source data only: the test with all three tertiles kept (3-group log-rank; trend = Cox on the tertile index)
  k <- cut(d$polar_mean, c(-Inf, cut3, Inf))
  sv <- Surv(d[[tcol]], d[[ecol]])
  p_3 <- pchisq(survdiff(sv ~ k)$chisq, 2, lower.tail = FALSE)
  p_trend <- summary(coxph(sv ~ as.numeric(k)))$coefficients[1, 5]
  d <- d[polar_mean <= cut3[1] | polar_mean > cut3[2]]          # middle tertile dropped
  d[, grp := factor(ifelse(polar_mean > cut3[2], "SecB", "SecA"), levels = c("SecA", "SecB"))]
  d[, `:=`(.time = get(tcol), .event = get(ecol))]
  fit <- survfit(Surv(.time, .event) ~ grp, data = d)
  sd <- survdiff(Surv(.time, .event) ~ grp, data = d); p_lr <- pchisq(sd$chisq, 1, lower.tail = FALSE)
  s60 <- summary(fit, times = 60, extend = TRUE)   # 5-yr KM estimate per group (strata order SecA, SecB)
  km_rows[[ep]] <<- data.table(endpoint = paste0(ep, "_5yr"), n_cohort = n_cohort, n = nrow(d), events = sum(d$.event),
                               cut_low = unname(cut3[1]), cut_high = unname(cut3[2]), n_SecA = sum(d$grp == "SecA"), n_SecB = sum(d$grp == "SecB"),
                               events_SecA = sum(d$.event[d$grp == "SecA"]), events_SecB = sum(d$.event[d$grp == "SecB"]),
                               surv5y_SecA = s60$surv[1], surv5y_SecA_lower = s60$lower[1], surv5y_SecA_upper = s60$upper[1],
                               surv5y_SecB = s60$surv[2], surv5y_SecB_lower = s60$lower[2], surv5y_SecB_upper = s60$upper[2],
                               logrank_p = p_lr, logrank_p_3groups = p_3, trend_p_3groups = p_trend)
  km_panel(fit, cols = c(COL_SECA, COL_SECB),
           labs = c(sprintf("SecA-dominant (n=%d)", sum(d$grp == "SecA")), sprintf("SecB-enriched (n=%d)", sum(d$grp == "SecB"))),
           p_text = ifelse(p_lr < 0.001, "p < 0.001", sprintf("p = %.3f", p_lr)),
           subtitle = NULL, xlim = c(0, 60), xbreak = 12, show_y = show_y)   # OS/PFS labelled in Affinity
}
p_A <- make_km("OS", show_y = TRUE)
p_B <- make_km("PFS", show_y = FALSE)
km_tab <- rbindlist(km_rows)
fwrite(km_tab, file.path(FIG_DIR, "xenium_tma_km_logrank.csv"))
print(km_tab)

# ---- 7c / 7d: forests of the three secretory measures (OS over PFS), BH FDR across the 6 tests ----
FE <- c(prop_SecA_of_sec = "SecA proportion", prop_SecB_of_sec = "SecB proportion", polar_mean = "Polarization")
COL_FE <- c(`SecA proportion` = COL_SECA, `SecB proportion` = COL_SECB, Polarization = COL_POL)
fd <- cox[feature %in% names(FE)]
fd[, `:=`(measure = factor(FE[feature], levels = rev(FE)), ep = sub("_5yr", "", endpoint))]
fd[, q_fdr := p.adjust(p_value, "BH")]
fwrite(fd[, .(endpoint, measure, n, n_events, HR, HR_lower, HR_upper, p_value, q_fdr)],
       file.path(FIG_DIR, "xenium_tma_forest_fdr.csv"))
print(fd[, .(endpoint, measure, n, HR = round(HR, 2), lo = round(HR_lower, 2), hi = round(HR_upper, 2),
             p = signif(p_value, 3), q = signif(q_fdr, 3))])
fd[, `:=`(hr_text = sprintf("%.2f", HR), p_text = fmt_p(p_value), q_text = fmt_p(q_fdr),
          col = COL_FE[as.character(measure)])]
G <- forest_geometry(lo = 0.45, data_hi = 2.2)              # shared Fig 7 forest layout (_km_forest_helpers.R)
x_hr <- G$x[["HR"]]; x_p <- G$x[["p"]]; x_q <- G$x[["q"]]
ypad <- forest_ylab_pad(FE)
plot_forest <- function(ep_f) {
  d <- fd[ep == ep_f]
  ggplot(d, aes(x = HR, y = measure)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_errorbar(aes(xmin = HR_lower, xmax = HR_upper, colour = col), width = 0.2, linewidth = 0.35,
                  orientation = "y") +
    geom_point(aes(fill = col, colour = col), size = 2.2, shape = 21, stroke = 0.5, alpha = 0.7) +
    geom_point(aes(colour = col), size = 2.2, shape = 1, stroke = 0.5) +
    scale_colour_identity() + scale_fill_identity() +
    geom_text(aes(x = x_hr, label = hr_text), hjust = 0, size = FA / .pt, colour = "black") +
    geom_text(aes(x = x_p, label = p_text), hjust = 0, size = FA / .pt, colour = "black") +
    geom_text(aes(x = x_q, label = q_text), hjust = 0, size = FA / .pt, colour = "black") +
    annotate("text", x = c(x_hr, x_p, x_q), y = length(FE) + 0.8, label = c("HR", "p", "q"),
             hjust = 0, size = FA / .pt, fontface = "bold") +
    scale_x_log10(breaks = c(0.5, 1, 2), labels = c("0.5", "1.0", "2.0"), limits = G$limits,
                  oob = scales::oob_keep) +
    coord_cartesian(clip = "off") +
    labs(x = sprintf("Adjusted HR per SD - 5-yr %s (n=%d)", ep_f, d$n[1]), y = NULL) +
    forest_theme() + theme(plot.margin = margin(10, 14, 4, 4 + ypad))
}
p_C <- plot_forest("OS") / plot_forest("PFS")   # 7c over 7d

save_survival_row(p_A, p_B, p_C, stem = "xenium_tma_survival_row_ABC", fig_dir = FIG_DIR,
                  width_in = 7.2, height_in = 2.6, widths = c(1, 1, 1.7))
message("Done.")
