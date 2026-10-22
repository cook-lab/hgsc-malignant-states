#!/usr/bin/env Rscript
# ============================================================================
# Figure 7k-7n — ICON7 bevacizumab survival row (one combined SVG)
# ----------------------------------------------------------------------------
# PURPOSE
#   The bevacizumab PFS benefit is concentrated in, and time-limited within,
#   SecB-polarized tumours (phase III ICON7 / GSE140082):
#     7k  SecA-polarized: PFS KM, Standard vs Bevacizumab (bev HR ~0.80, ns)
#     7l  SecB-polarized: PFS KM, Standard vs Bevacizumab (bev HR ~0.36, p<0.001)
#     7m/7n  time-varying bev HR (15-mo period split): SecB-polarized (m) over SecA-polarized (n)
#   Built with the SAME km_panel()/forest_theme()/save_survival_row() helpers as
#   the TMA (7a-7d) and TCGA (7g-7j) rows, so the survival rows are pixel-for-pixel
#   consistent.
#
#   Colours: k/l arms match Fig 7f proteins — Standard = KRT7 blue,
#   Bevacizumab = KRT19 pink; m/n groups = atlas SecA gold / SecB brown.
#
# INPUTS
#   cfg_obj("icon7_cohort")  -> ICON7/data/processed/cohort_filtered.tsv
#   Shared helpers: config/config.R, spatial/00_setup/00_setup.R, _km_forest_helpers.R
#
# OUTPUTS
#   figures_dir/figure7/icon7_survival_row_HIJ.{svg,png}  (Fig 7k-n; the stem keeps earlier letters)
#   figures_dir/figure7/icon7_km_pfs_within_group_hr.csv  (7k,l: adjusted bevacizumab HR per group; Table S12 checks it)
#   figures_dir/figure7/icon7_timevarying_forest_15mo_data.csv (7m,n; Table S12 checks it)
#
# MANUSCRIPT PANEL(S): Fig 7k, 7l, 7m, 7n
# RUNTIME TIER: fast (per-patient table; in-script Cox)
# ============================================================================

.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))
source(file.path(.here, "..", "..", "spatial", "00_setup", "00_setup.R"))
source(file.path(.here, "_km_forest_helpers.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(patchwork)
})
set.seed(1)

FIG_DIR <- cfg_path("figures_dir", "figure7")
FA <- F7_FA; FK <- F7_FK; FN <- F7_FA

# Arm palette (k/l) = Fig 7f proteins.  Group palette (m/n) = atlas SecA/SecB.
COL_STD  <- "#5665B6"                                   # KRT7 blue  -> Standard
COL_BEV  <- "#D14E6C"                                   # KRT19 pink -> Bevacizumab
COL_SECA <- unname(ref_palette["SecA epithelium"])      # gold
COL_SECB <- unname(ref_palette["SecB epithelium"])      # brown
CUT_DAYS <- round(15 * 30.44)                           # 15-month split = 456 d

# =============================================================================
# 1. LOAD + PREPARE (FIGO III/IV, complete cases, median split of polarization)
# =============================================================================
cohort_file <- cfg_obj("icon7_cohort")
if (!file.exists(cohort_file)) {
  stop("ICON7 cohort table not found: ", cohort_file,
       "\n  Point DATA_ROOT at the deposited bundle.")
}
d <- fread(cohort_file)
d$treatment <- factor(d$treatment, levels = c("standard", "bevacizumab"))
d <- d[figo_stage %in% c("III", "IV")]
d$figo_stage <- factor(d$figo_stage, levels = c("III", "IV"))
d$debulking  <- factor(ifelse(d$debulking_status == "OPTIMAL", "OPTIMAL", "SUBOPT_INOP"),
                       levels = c("OPTIMAL", "SUBOPT_INOP"))
d$age_num <- as.numeric(d$age)

mvars <- c("final_pfstm", "final_pfsid", "polarization_ucell",
           "treatment", "figo_stage", "debulking", "age_num")
ana <- d[complete.cases(d[, ..mvars])]
med <- median(ana$polarization_ucell)
ana[, pol_bin := factor(ifelse(polarization_ucell > med, "SecB", "SecA"),
                        levels = c("SecA", "SecB"))]
ana$arm <- factor(ana$treatment, labels = c("Standard", "Bevacizumab"))
message(sprintf("ICON7 cohort: n = %d (FIGO III/IV); SecA n = %d, SecB n = %d",
                nrow(ana), sum(ana$pol_bin == "SecA"), sum(ana$pol_bin == "SecB")))

COVARS <- "figo_stage + debulking + age_num"

# --- within-group adjusted bevacizumab HR, events and log-rank p (CSV) --------
within_group_hr <- function(group) {
  s <- ana[pol_bin == group]
  f <- coxph(as.formula(sprintf("Surv(final_pfstm, final_pfsid) ~ treatment + %s", COVARS)),
             data = s)
  ci <- summary(f)$conf.int["treatmentbevacizumab", ]
  p  <- summary(f)$coefficients["treatmentbevacizumab", "Pr(>|z|)"]
  lr <- survdiff(Surv(final_pfstm, final_pfsid) ~ treatment, data = s)
  data.frame(group = group, n = nrow(s), events = sum(s$final_pfsid),
             HR = ci[1], lo = ci[3], hi = ci[4], adj_p = p,
             logrank_p = pchisq(lr$chisq, 1, lower.tail = FALSE), row.names = NULL)
}
hr_tab <- rbind(within_group_hr("SecA"), within_group_hr("SecB"))
write.csv(hr_tab, file.path(FIG_DIR, "icon7_km_pfs_within_group_hr.csv"), row.names = FALSE)
message("\nWithin-group PFS bevacizumab HR (adjusted):"); print(hr_tab)

# =============================================================================
# 2. KM PANELS 7k / 7l (native km_panel; PFS in months, arm palette)
# =============================================================================
# Clip the x-axis to the actual ICON7 PFS follow-up (curves end well before the
# TCGA row's 60-mo / 5-yr window), rounded up to a clean 12-mo tick. Both KM
# panels share this xmax so k and l stay directly comparable.
icon7_xmax <- ceiling(max(ana$final_pfstm / 30.44, na.rm = TRUE) / 6) * 6
message(sprintf("ICON7 KM x-axis clipped to %d mo (max follow-up %.1f mo)",
                icon7_xmax, max(ana$final_pfstm / 30.44, na.rm = TRUE)))

make_icon7_km <- function(group, subtitle, show_y = TRUE) {
  dat <- as.data.frame(ana[pol_bin == group])
  dat$.time  <- dat$final_pfstm / 30.44
  dat$.event <- dat$final_pfsid
  fit <- survfit(Surv(.time, .event) ~ arm, data = dat)   # strata: Standard, Bevacizumab
  n_std <- sum(dat$arm == "Standard"); n_bev <- sum(dat$arm == "Bevacizumab")
  sd <- survdiff(Surv(.time, .event) ~ arm, data = dat)
  p_lr <- 1 - pchisq(sd$chisq, 1)
  km_panel(fit,
           cols = c(COL_STD, COL_BEV),
           labs = c(sprintf("Standard (n=%d)", n_std),
                    sprintf("Bevacizumab (n=%d)", n_bev)),
           p_text = ifelse(p_lr < 0.001, "p < 0.001", sprintf("p = %.3f", p_lr)),
           subtitle = subtitle, show_y = show_y,
           xlim = c(0, icon7_xmax), xbreak = 12)
}
p_H <- make_icon7_km("SecA", NULL, show_y = TRUE)   # 7k (SecA-polarized; label in Affinity)
p_I <- make_icon7_km("SecB", NULL, show_y = FALSE)  # 7l (SecB-polarized; label in Affinity)

# =============================================================================
# 3. TIME-VARYING bev HR at the 15-month cut, per group (7m, 7n)
# =============================================================================
timevarying_at <- function(group, cut_days) {
  s <- as.data.frame(ana[pol_bin == group])
  s$.t <- s$final_pfstm; s$.e <- s$final_pfsid
  sp <- survSplit(Surv(.t, .e) ~ ., data = s, cut = cut_days, episode = "period")
  sp$period <- factor(sp$period, levels = 1:2, labels = c("early", "late"))
  m <- coxph(Surv(tstart, .t, .e) ~ treatment * period + figo_stage +
               debulking + age_num, data = sp)
  cf <- summary(m)$coefficients; v <- vcov(m)
  ti <- "treatmentbevacizumab"; ii <- "treatmentbevacizumab:periodlate"
  le_lhr <- cf[ti, "coef"]; le_se <- cf[ti, "se(coef)"]
  gt_lhr <- cf[ti, "coef"] + cf[ii, "coef"]
  gt_se  <- sqrt(v[ti, ti] + v[ii, ii] + 2 * v[ti, ii])
  cut_mo <- round(cut_days / 30.44, 0)
  data.table(group = group,
             period = c(sprintf("0-%d mo", cut_mo), sprintf(">%d mo", cut_mo)),
             HR = exp(c(le_lhr, gt_lhr)),
             lo = exp(c(le_lhr - 1.96 * le_se, gt_lhr - 1.96 * gt_se)),
             hi = exp(c(le_lhr + 1.96 * le_se, gt_lhr + 1.96 * gt_se)),
             p  = c(cf[ti, "Pr(>|z|)"], 2 * pnorm(-abs(gt_lhr / gt_se))),
             period_interaction_p = cf[ii, "Pr(>|z|)"])
}
fdat <- rbindlist(lapply(c("SecA", "SecB"), timevarying_at, cut_days = CUT_DAYS))
# Benjamini-Hochberg FDR: q across the 4 bevacizumab HRs shown (as Fig 7c,d and 7i,j); interaction q across the 2 groups
fdat[, q := p.adjust(p, "BH")]
ig <- unique(fdat[, .(group, period_interaction_p)])[, period_interaction_q := p.adjust(period_interaction_p, "BH")]
fdat[ig, period_interaction_q := i.period_interaction_q, on = "group"]
write.csv(fdat, file.path(FIG_DIR, "icon7_timevarying_forest_15mo_data.csv"), row.names = FALSE)
message("\n15-month time-varying bevacizumab HR (PFS):")
print(fdat[, .(group, period, HR = round(HR, 3), lo = round(lo, 3), hi = round(hi, 3),
               p = signif(p, 3), q = signif(q, 3), period_interaction_p = signif(period_interaction_p, 3),
               period_interaction_q = signif(period_interaction_q, 3))])

lab_lo <- round(CUT_DAYS / 30.44)
PERIODS <- c(sprintf("0-%d mo", lab_lo), sprintf(">%d mo", lab_lo))
fdat[, period_f := factor(period, levels = rev(PERIODS))]  # 0-15 mo on top
fdat[, early := grepl("^0-", period)]
fdat[, hr_text := sprintf("%.2f", HR)]                    # HR only (CI dropped)
fdat[, p_text  := ifelse(p < 0.001, "<.001",
                         sub("^0\\.", ".", formatC(p, format = "f", digits = 3)))]
fdat[, q_text  := ifelse(q < 0.001, "<.001",
                         sub("^0\\.", ".", formatC(q, format = "f", digits = 3)))]

# Same layout as the TMA (a-d) and TCGA (g-j) rows: one forest per KM panel, stacked (SecB = m over
# SecA = n), shared forest geometry; circle = 0-15 mo, triangle = >15 mo; q across all 4 estimates.
G <- forest_geometry(lo = 0.08, data_hi = 3)              # shared Fig 7 forest layout (_km_forest_helpers.R)
x_hr_pos <- G$x[["HR"]]; x_p_pos <- G$x[["p"]]; x_q_pos <- G$x[["q"]]
ypad <- forest_ylab_pad(PERIODS)
plot_icon7_forest <- function(g) {
  d <- fdat[group == g]
  col <- c(SecA = COL_SECA, SecB = COL_SECB)[[g]]
  ggplot(d, aes(x = HR, y = period_f)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_errorbar(aes(xmin = pmax(lo, G$limits[1]), xmax = hi), width = 0.2, linewidth = 0.35,
                  color = col, orientation = "y") +
    geom_point(aes(shape = early), fill = col, color = col, size = 2.2, stroke = 0.5, alpha = 0.7) +
    geom_point(aes(shape = early), fill = NA, color = col, size = 2.2, stroke = 0.5) +
    geom_text(aes(x = x_hr_pos, label = hr_text), hjust = 0, size = FA / .pt, color = "black") +
    geom_text(aes(x = x_p_pos, label = p_text),  hjust = 0, size = FA / .pt, color = "black") +
    geom_text(aes(x = x_q_pos, label = q_text),  hjust = 0, size = FA / .pt, color = "black") +
    annotate("text", x = c(x_hr_pos, x_p_pos, x_q_pos), y = length(PERIODS) + 0.8,
             label = c("HR", "p", "q"), hjust = 0, size = FA / .pt, fontface = "bold") +
    scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 24), guide = "none") +
    scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2),
                  labels = c("0.1", "0.25", "0.5", "1.0", "2.0"),
                  limits = G$limits, oob = scales::oob_keep) +
    coord_cartesian(clip = "off") +
    labs(x = sprintf("Bevacizumab HR (PFS) - %s-polarized", g), y = NULL) +
    forest_theme() + theme(plot.margin = margin(10, 14, 4, 4 + ypad))
}
p_J <- plot_icon7_forest("SecB") / plot_icon7_forest("SecA")   # 7m (SecB) over 7n (SecA)

# =============================================================================
# 4. SAVE — ONE combined SVG (7k | 7l | 7m/7n), same footprint as the TMA and TCGA rows.
# =============================================================================
save_survival_row(p_H, p_I, p_J,
                  stem = "icon7_survival_row_HIJ", fig_dir = FIG_DIR,
                  width_in = 7.2, height_in = 2.6, widths = c(1, 1, 1.7))

message("\nDone.")
