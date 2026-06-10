#!/usr/bin/env Rscript
# ============================================================================
# ICON7 — SecB polarization x bevacizumab: per-arm slope reversal   [PLACEHOLDER]
# ----------------------------------------------------------------------------
# *** PLACEHOLDER PANEL — final figure/panel letter NOT yet assigned. ***
# Slated as an "expanded Figure 7" (the bevacizumab external validation, sibling
# to Fig 7E-G TCGA-OV). Rename this directory + script and set MANUSCRIPT
# PANEL(S) once the panel is decided. See ./README.md and the source module
# ../../2026_final_xenium_analysis/davids side quests/ICON7/ (+ its report.html).
#
# PURPOSE
#   The canonical ICON7 result: baseline SecB-SecA polarization is prognostic in
#   the chemotherapy-only arm and FLAT in the bevacizumab arm -> the
#   polarization->outcome slope reverses sign between arms (the cleanest
#   signature of a true treatment-modifier interaction). This script reproduces
#   that headline panel as a per-arm "HR per 1-SD" forest, computing the Cox
#   models in-script from the per-patient ICON7 cohort table (mirrors the
#   in-script Cox in figure7/03_tcga_km_forest.R). The full polished narrative
#   panels (KM extremes, median-rescue, time-varying, composite) are produced by
#   the source backend script 04_figures.R (Section A).
#
# INPUTS
#   cfg_obj("icon7_cohort")  -> data/processed/cohort_filtered.tsv
#     per-patient: polarization_ucell, treatment, final_{pfstm,pfsid,ostm,osid},
#     figo_stage, debulking_status, age, t1_cluster_name
#   Shared helpers: config/config.R, spatial/00_setup/00_setup.R (theme_lab, ref_palette)
#
# OUTPUTS
#   figures_dir/figure_icon7_bevacizumab/icon7_per_arm_slope_forest.{svg,pdf}
#   figures_dir/figure_icon7_bevacizumab/icon7_per_arm_slope_forest_data.csv
#
# MANUSCRIPT PANEL(S): TBD (expanded Fig 7)
# RUNTIME TIER: fast (per-patient table; in-script Cox)
# ============================================================================

.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))
source(file.path(.here, "..", "..", "spatial", "00_setup", "00_setup.R"))

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
})

FIG_DIR <- cfg_path("figures_dir", "figure_icon7_bevacizumab")
FA <- 6; FK <- 5.5; FN <- 5

# =============================================================================
# 1. LOAD per-patient ICON7 cohort
# =============================================================================
cohort_file <- cfg_obj("icon7_cohort")
if (!file.exists(cohort_file)) {
  stop("ICON7 cohort table not found: ", cohort_file,
       "\n  Run the ICON7 backend (scripts 01-04) or point DATA_ROOT at the ",
       "deposited bundle. See ./README.md.")
}
d <- read.delim(cohort_file, stringsAsFactors = FALSE)
message("ICON7 cohort rows: ", nrow(d))

# FIGO III/IV survival cohort (n=191); the n=21 FIGO I/II subset has 0 OS events.
d <- d[d$figo_stage %in% c("III", "IV"), ]
d$figo_stage <- factor(d$figo_stage, levels = c("III", "IV"))
d$treatment  <- factor(d$treatment, levels = c("standard", "bevacizumab"))
d$debulking  <- factor(ifelse(d$debulking_status == "OPTIMAL", "OPTIMAL", "SUBOPT_INOP"),
                       levels = c("OPTIMAL", "SUBOPT_INOP"))
d$age_num    <- as.numeric(d$age)
# Polarization standardized across the survival cohort -> "per 1-SD" HRs.
d$polarization_z <- as.numeric(scale(d$polarization_ucell))
message("FIGO III/IV n = ", nrow(d),
        " (standard=", sum(d$treatment == "standard"),
        ", bev=", sum(d$treatment == "bevacizumab"), ")")

# =============================================================================
# 2. PER-ARM CONTINUOUS Cox: HR per 1-SD of SecB polarization, by arm
# =============================================================================
cox_per_sd <- function(dat, time_col, event_col) {
  dat$.time <- dat[[time_col]]; dat$.event <- dat[[event_col]]
  fit <- coxph(Surv(.time, .event) ~ polarization_z + figo_stage + debulking + age_num,
               data = dat)
  s  <- summary(fit)
  ci <- s$conf.int["polarization_z", , drop = FALSE]
  data.frame(HR = ci[1, "exp(coef)"], HR_lower = ci[1, "lower .95"],
             HR_upper = ci[1, "upper .95"],
             p_value = s$coefficients["polarization_z", "Pr(>|z|)"],
             n = s$n, events = s$nevent, stringsAsFactors = FALSE)
}

rows <- list()
for (ep in c("OS", "PFS")) {
  tc <- if (ep == "OS") "final_ostm"  else "final_pfstm"
  ec <- if (ep == "OS") "final_osid"  else "final_pfsid"
  for (arm in c("standard", "bevacizumab")) {
    r <- cox_per_sd(d[d$treatment == arm, ], tc, ec)
    r$endpoint <- ep
    r$arm <- if (arm == "standard") "Chemo only" else "Chemo + bevacizumab"
    rows[[length(rows) + 1]] <- r
  }
}
forest_df <- do.call(rbind, rows)
forest_df$arm      <- factor(forest_df$arm,
                             levels = c("Chemo + bevacizumab", "Chemo only"))
forest_df$endpoint <- factor(forest_df$endpoint, levels = c("OS", "PFS"))
forest_df$hr_text  <- sprintf("%.2f (%.2f-%.2f)",
                              forest_df$HR, forest_df$HR_lower, forest_df$HR_upper)
forest_df$p_text   <- ifelse(forest_df$p_value < 0.001, "<.001",
                             sub("^0\\.", ".", formatC(forest_df$p_value, format = "f", digits = 3)))
message("\nPer-arm HR per 1-SD SecB polarization (HR>1 = more SecB, worse outcome):")
print(forest_df[, c("endpoint", "arm", "n", "events", "HR", "HR_lower", "HR_upper", "p_value")])

# =============================================================================
# 3. FOREST: the slope reversal (chemo prognostic, bev flat)
# =============================================================================
arm_pal <- c("Chemo only" = unname(ref_palette["SecB epithelium"]),
             "Chemo + bevacizumab" = "grey55")

p <- ggplot(forest_df, aes(x = HR, y = arm, color = arm)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_errorbar(aes(xmin = HR_lower, xmax = HR_upper), width = 0.18,
                linewidth = 0.4, orientation = "y") +
  geom_point(size = 2.2) +
  geom_text(aes(x = 2.6, label = hr_text), hjust = 0, size = FN / .pt, color = "black") +
  facet_wrap(~ endpoint, ncol = 1) +
  scale_color_manual(values = arm_pal, guide = "none") +
  scale_x_log10(breaks = c(0.5, 1, 1.5, 2), limits = c(0.5, 4.5)) +
  coord_cartesian(clip = "off") +
  labs(x = "HR per 1-SD SecB polarization", y = NULL,
       title = "ICON7: SecB prognostic under chemo, abolished by bevacizumab") +
  theme_lab(base_size = 6) +
  theme(plot.title = element_text(size = FA, face = "bold"),
        axis.text.y = element_text(size = FK),
        plot.margin = margin(8, 70, 6, 6))

# =============================================================================
# 4. SAVE
# =============================================================================
for (ext in c("svg", "pdf")) {
  ggsave(file.path(FIG_DIR, paste0("icon7_per_arm_slope_forest.", ext)),
         p, width = 4, height = 2.6, bg = "white")
}
write.csv(forest_df[, c("endpoint", "arm", "n", "events",
                        "HR", "HR_lower", "HR_upper", "p_value")],
          file.path(FIG_DIR, "icon7_per_arm_slope_forest_data.csv"), row.names = FALSE)
message("\n  Saved: icon7_per_arm_slope_forest.{svg,pdf} (+ _data.csv)")
message("Done. [PLACEHOLDER — assign panel letter, then rename per ./README.md]")
