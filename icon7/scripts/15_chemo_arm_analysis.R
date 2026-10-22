# Trimmed copy (2026-10-21) of 2026_final_xenium_analysis/davids side quests/ICON7/scripts/15_chemo_arm_analysis.R,
# the script that ran. Kept: section B.2, which writes 17_chemo_arm_continuous_polarization.tsv and
# 17b_continuous_polarization_per_arm.tsv (Table S12 model A). Removed: sections B.1, B.3 and B.4 and their plots.
# The untrimmed script is in the lab project folder.
# ============================================================================
# 15_chemo_arm_analysis.R
#
# Task B: What does baseline polarization say about outcomes in the
# chemotherapy-only arm, independent of bevacizumab?
#
# B.2 Continuous polarization Cox in standard arm
#
# Outputs:
#   results/tables/17_chemo_arm_continuous_polarization.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
})

ROOT <- "/Users/dpcook/My Drive/Lab/Manuscripts/hgsc_epitypes/supplemental_analysis/ICON7"
IN_FILE <- file.path(ROOT, "data/processed/cohort_filtered.tsv")
TAB_DIR <- file.path(ROOT, "results/tables")
FIG_DIR <- file.path(ROOT, "results/figs")

d <- fread(IN_FILE)
d$treatment <- factor(d$treatment, levels = c("standard", "bevacizumab"))
d <- d[figo_stage %in% c("III", "IV")]
d$figo_stage <- factor(d$figo_stage, levels = c("III", "IV"))
d$debulking <- factor(ifelse(d$debulking_status == "OPTIMAL",
                             "OPTIMAL", "SUBOPT_INOP"),
                      levels = c("OPTIMAL", "SUBOPT_INOP"))
d$age_num <- as.numeric(d$age)

# Compute tertile cuts on the FULL FIGO III/IV cohort (so that the same cuts
# are used as the main analysis) — i.e. quantiles of polarization over n=191.
breaks <- quantile(d$polarization_ucell, c(0, 1/3, 2/3, 1), na.rm = TRUE)
d$pol_tertile <- cut(d$polarization_ucell, breaks = breaks,
                     include.lowest = TRUE,
                     labels = c("Low", "Mid", "High"))
d$pol_tertile <- factor(d$pol_tertile, levels = c("Low", "Mid", "High"))

mvars <- c("final_pfstm", "final_pfsid", "final_ostm", "final_osid",
           "polarization_ucell", "pol_tertile", "treatment", "figo_stage",
           "debulking", "age_num")
ana <- d[complete.cases(d[, ..mvars])]
message(sprintf("Full FIGO III/IV cohort n = %d", nrow(ana)))

# Standardized polarization (z-scored within full cohort) so the
# continuous-polarization HR is interpretable per 1 SD
ana$polarization_z <- as.numeric(scale(ana$polarization_ucell))

# Standard-arm subset
std <- ana[treatment == "standard"]
message(sprintf("Standard-arm n = %d (FIGO III/IV)", nrow(std)))
print(table(std$pol_tertile, useNA = "ifany"))

bev <- ana[treatment == "bevacizumab"]
message(sprintf("Bev-arm n = %d (FIGO III/IV)", nrow(bev)))

cox_tidy <- function(m) {
  s <- summary(m)
  coef <- s$coefficients
  ci <- s$conf.int
  data.table(
    term = rownames(coef),
    HR = coef[, "exp(coef)"],
    HR_lower = ci[, "lower .95"],
    HR_upper = ci[, "upper .95"],
    se = coef[, "se(coef)"],
    z = coef[, "z"],
    p_value = coef[, "Pr(>|z|)"],
    n = s$n,
    nevent = s$nevent
  )
}

# ============================================================================
# B.2 — Continuous polarization Cox in standard arm only
# ============================================================================

m_pfs_cont_std <- coxph(Surv(final_pfstm, final_pfsid) ~ polarization_z +
                          figo_stage + debulking + age_num, data = std)
m_os_cont_std <- coxph(Surv(final_ostm, final_osid) ~ polarization_z +
                         figo_stage + debulking + age_num, data = std)

t_cont_pfs <- cox_tidy(m_pfs_cont_std); t_cont_pfs$endpoint <- "PFS"
t_cont_os  <- cox_tidy(m_os_cont_std);  t_cont_os$endpoint <- "OS"
cont_tab <- rbind(t_cont_pfs, t_cont_os)
# Annotate which row is the polarization_z row
cont_tab$is_polarization <- cont_tab$term == "polarization_z"
fwrite(cont_tab, file.path(TAB_DIR, "17_chemo_arm_continuous_polarization.tsv"),
       sep = "\t")
print(cont_tab[is_polarization == TRUE])

pol_pfs_HR <- cont_tab[is_polarization == TRUE & endpoint == "PFS", HR]
pol_pfs_p  <- cont_tab[is_polarization == TRUE & endpoint == "PFS", p_value]
pol_os_HR  <- cont_tab[is_polarization == TRUE & endpoint == "OS", HR]
pol_os_p   <- cont_tab[is_polarization == TRUE & endpoint == "OS", p_value]

message(sprintf("Standard-arm continuous polarization HR per 1-SD (PFS) = %.3f, p = %.4g",
                pol_pfs_HR, pol_pfs_p))
message(sprintf("Standard-arm continuous polarization HR per 1-SD (OS) = %.3f, p = %.4g",
                pol_os_HR, pol_os_p))

# Compare to bev arm (for symmetry)
m_pfs_cont_bev <- coxph(Surv(final_pfstm, final_pfsid) ~ polarization_z +
                          figo_stage + debulking + age_num, data = bev)
m_os_cont_bev <- coxph(Surv(final_ostm, final_osid) ~ polarization_z +
                         figo_stage + debulking + age_num, data = bev)
bev_pfs_HR <- summary(m_pfs_cont_bev)$conf.int["polarization_z", "exp(coef)"]
bev_pfs_p  <- summary(m_pfs_cont_bev)$coefficients["polarization_z", "Pr(>|z|)"]
bev_os_HR  <- summary(m_os_cont_bev)$conf.int["polarization_z", "exp(coef)"]
bev_os_p   <- summary(m_os_cont_bev)$coefficients["polarization_z", "Pr(>|z|)"]

# Append the bev-arm rows for symmetric reporting
arm_compare <- data.table(
  arm = c("standard", "standard", "bevacizumab", "bevacizumab"),
  endpoint = c("PFS", "OS", "PFS", "OS"),
  n = c(m_pfs_cont_std$n, m_os_cont_std$n,
        m_pfs_cont_bev$n, m_os_cont_bev$n),
  events = c(m_pfs_cont_std$nevent, m_os_cont_std$nevent,
             m_pfs_cont_bev$nevent, m_os_cont_bev$nevent),
  polarization_HR_per_1SD = c(pol_pfs_HR, pol_os_HR, bev_pfs_HR, bev_os_HR),
  polarization_HR_lower = c(
    summary(m_pfs_cont_std)$conf.int["polarization_z", "lower .95"],
    summary(m_os_cont_std)$conf.int["polarization_z", "lower .95"],
    summary(m_pfs_cont_bev)$conf.int["polarization_z", "lower .95"],
    summary(m_os_cont_bev)$conf.int["polarization_z", "lower .95"]),
  polarization_HR_upper = c(
    summary(m_pfs_cont_std)$conf.int["polarization_z", "upper .95"],
    summary(m_os_cont_std)$conf.int["polarization_z", "upper .95"],
    summary(m_pfs_cont_bev)$conf.int["polarization_z", "upper .95"],
    summary(m_os_cont_bev)$conf.int["polarization_z", "upper .95"]),
  polarization_p = c(pol_pfs_p, pol_os_p, bev_pfs_p, bev_os_p)
)
fwrite(arm_compare,
       file.path(TAB_DIR, "17b_continuous_polarization_per_arm.tsv"),
       sep = "\t")
print(arm_compare)
