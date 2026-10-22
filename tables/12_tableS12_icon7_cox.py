#!/usr/bin/env python3
"""Supplementary Table 12 -- ICON7 (GSE140082) Cox models: secretory polarization and bevacizumab benefit
(Figure 7k-n; Results, ICON7).

Regenerates the table end-to-end from the per-patient ICON7 cohort (cfg obj "icon7_cohort"):
  A. Chemotherapy-alone (standard-arm) prognostic model: continuous polarization (z, per +1 SD) + FIGO stage +
     debulking + age, OS and PFS; also the bevacizumab arm. Faithful extraction of section B.2 of the ICON7
     module script 15_chemo_arm_analysis.R (copy in icon7/scripts/; source of the text's OS HR 1.38);
     checked against its output results/tables/17_chemo_arm_continuous_polarization.tsv.
  B. Within-group bevacizumab PFS effect (SecA- / SecB-polarized, median split of polarization_ucell) adjusted for
     FIGO stage + debulking + age, KM medians and log-rank p: figures/figure7/06_fig7klmn_icon7_survival_row.R section 1-2.
  C. 15-month time-varying bevacizumab effect (treatment x period, survSplit at 456 d) adjusted for the same
     covariates: 06_fig7klmn_icon7_survival_row.R section 3.
  D. Treatment x polarization-group interaction for PFS (median split; same covariates) in all 191 patients, with the
     likelihood-ratio test against the model without the interaction (the test the Results and Methods report).
  B and C are checked against the Figure 7 generator outputs icon7_km_pfs_within_group_hr.csv and
  icon7_timevarying_forest_15mo_data.csv. Model code, cohort filters and covariates are unchanged; the only
  additions are exporting all coefficients, KM medians and group sizes, and the Benjamini-Hochberg q values
  (model C: across the 4 bevacizumab HRs of Figure 7m,n, and the period interaction across the 2 groups), which
  are checked against the q columns the generator writes.

Intermediates: $OUTPUT_ROOT/supp_tables/ST12/ (inputs read-only).
Outputs: SuppTable12_icon7_polarization_bevacizumab_cox.xlsx (+ one CSV per sheet) via supp_xlsx.write_table.
Manuscript item(s): Table S12; Results ICON7 numbers (OS HR 1.38; interaction HR 0.44, LRT p 0.025).
"""
import os, sys, subprocess
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE); sys.path.insert(0, REPO)
from supp_xlsx import write_table  # noqa: E402
from config import config as cfg  # noqa: E402

COHORT = cfg.obj("icon7_cohort")
ICON7_TAB = str(cfg.DATA_ROOT / cfg.CFG["objects"]["icon7_results_dir"])
FIG7_OUT = os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "figures/figure7")   # where 06_fig7klmn_icon7_survival_row.R wrote its CSVs
WORK = os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "supp_tables/ST12")
os.makedirs(WORK, exist_ok=True)

R_CODE = r'''
suppressPackageStartupMessages({library(data.table); library(survival)})
OUT <- Sys.getenv("ST12_WORK"); d0 <- fread(Sys.getenv("ST12_COHORT"))
CUT_DAYS <- round(15 * 30.44)   # 456 d, as 06_fig7klmn_icon7_survival_row.R
tidy <- function(m, model, population, endpoint) {
  s <- summary(m); cf <- s$coefficients; ci <- s$conf.int
  data.table(model = model, population = population, endpoint = endpoint, term = rownames(cf),
             coef = cf[, "coef"], se = cf[, "se(coef)"], HR = ci[, "exp(coef)"],
             HR_lower = ci[, "lower .95"], HR_upper = ci[, "upper .95"], z = cf[, "z"],
             p_value = cf[, "Pr(>|z|)"], n = s$n, events = s$nevent)
}
km_med <- function(dat, tcol, ecol) {      # median (months) per arm with 95% CI
  f <- survfit(as.formula(sprintf("Surv(%s/30.44, %s) ~ treatment", tcol, ecol)), data = dat)
  q <- summary(f)$table
  data.table(arm = sub("treatment=", "", rownames(q)), n = q[, "records"], events = q[, "events"],
             median_mo = q[, "median"], median_lo = q[, "0.95LCL"], median_hi = q[, "0.95UCL"])
}
full <- list(); rows <- list(); meds <- list()

## ---------- prep (identical to 06_fig7klmn_icon7_survival_row.R section 1) ----------
d <- copy(d0)
d$treatment <- factor(d$treatment, levels = c("standard", "bevacizumab"))
d <- d[figo_stage %in% c("III", "IV")]
d$figo_stage <- factor(d$figo_stage, levels = c("III", "IV"))
d$debulking  <- factor(ifelse(d$debulking_status == "OPTIMAL", "OPTIMAL", "SUBOPT_INOP"), levels = c("OPTIMAL", "SUBOPT_INOP"))
d$age_num <- as.numeric(d$age)
mvars <- c("final_pfstm", "final_pfsid", "polarization_ucell", "treatment", "figo_stage", "debulking", "age_num")
ana <- d[complete.cases(d[, ..mvars])]
med <- median(ana$polarization_ucell)
ana[, pol_bin := factor(ifelse(polarization_ucell > med, "SecB", "SecA"), levels = c("SecA", "SecB"))]
COVARS <- "figo_stage + debulking + age_num"
cat(sprintf("ICON7 FIGO III/IV n=%d; SecA %d SecB %d; median cut %.6f\n", nrow(ana), sum(ana$pol_bin=="SecA"), sum(ana$pol_bin=="SecB"), med))

## ---------- A. standard / bev arm continuous polarization (15_chemo_arm_analysis.R B.2) ----------
brk <- quantile(d$polarization_ucell, c(0, 1/3, 2/3, 1), na.rm = TRUE)
d$pol_tertile <- factor(cut(d$polarization_ucell, breaks = brk, include.lowest = TRUE, labels = c("Low","Mid","High")), levels = c("Low","Mid","High"))
mv2 <- c("final_pfstm","final_pfsid","final_ostm","final_osid","polarization_ucell","pol_tertile","treatment","figo_stage","debulking","age_num")
ana2 <- d[complete.cases(d[, ..mv2])]
ana2$polarization_z <- as.numeric(scale(ana2$polarization_ucell))
cat(sprintf("chemo-arm module cohort n=%d; polarization SD=%.6f\n", nrow(ana2), sd(ana2$polarization_ucell)))
for (arm in c("standard", "bevacizumab")) for (ep in c("OS", "PFS")) {
  s <- ana2[treatment == arm]
  tc <- if (ep == "OS") "final_ostm" else "final_pfstm"; ec <- if (ep == "OS") "final_osid" else "final_pfsid"
  m <- coxph(as.formula(sprintf("Surv(%s, %s) ~ polarization_z + %s", tc, ec, COVARS)), data = s)
  t <- tidy(m, "A_arm_continuous_polarization", paste0(arm, " arm"), ep); full[[length(full)+1]] <- t
  r <- t[term == "polarization_z"]
  rows[[length(rows)+1]] <- r[, .(id = paste0("A_", arm, "_", ep), model = "A_arm_continuous_polarization",
      population = paste0(arm, " arm"), endpoint = ep, term, HR, HR_lower, HR_upper, p_value, n, events,
      polarization_SD = sd(ana2$polarization_ucell))]
}

## ---------- B. within-group adjusted bev PFS HR + KM (06_fig7klmn_icon7_survival_row.R) ----------
for (g in c("SecA", "SecB")) {
  s <- ana[pol_bin == g]
  m <- coxph(as.formula(sprintf("Surv(final_pfstm, final_pfsid) ~ treatment + %s", COVARS)), data = s)
  t <- tidy(m, "B_within_group_bev_PFS", paste0(g, "-polarized"), "PFS"); full[[length(full)+1]] <- t
  lr <- survdiff(Surv(final_pfstm, final_pfsid) ~ treatment, data = s)
  mm <- km_med(s, "final_pfstm", "final_pfsid"); mm[, `:=`(group = g, endpoint = "PFS")]; meds[[length(meds)+1]] <- mm
  mo <- km_med(s, "final_ostm", "final_osid"); mo[, `:=`(group = g, endpoint = "OS")]; meds[[length(meds)+1]] <- mo
  r <- t[term == "treatmentbevacizumab"]
  rows[[length(rows)+1]] <- r[, .(id = paste0("B_", g), model = "B_within_group_bev_PFS", population = paste0(g, "-polarized"),
      endpoint = "PFS", term, HR, HR_lower, HR_upper, p_value, n, events,
      logrank_p = pchisq(lr$chisq, 1, lower.tail = FALSE),
      median_standard_mo = mm[arm == "standard", median_mo], median_bev_mo = mm[arm == "bevacizumab", median_mo])]
}

## ---------- C. time-varying (06_fig7klmn_icon7_survival_row.R section 3) ----------
for (g in c("SecA", "SecB")) {
  s <- as.data.frame(ana[pol_bin == g]); s$.t <- s$final_pfstm; s$.e <- s$final_pfsid
  sp <- survSplit(Surv(.t, .e) ~ ., data = s, cut = CUT_DAYS, episode = "period")
  sp$period <- factor(sp$period, levels = 1:2, labels = c("early", "late"))
  m <- coxph(Surv(tstart, .t, .e) ~ treatment * period + figo_stage + debulking + age_num, data = sp)
  t <- tidy(m, "C_timevarying_bev_PFS_15mo", paste0(g, "-polarized"), "PFS"); t[, n := nrow(s)]; full[[length(full)+1]] <- t
  cf <- summary(m)$coefficients; v <- vcov(m)
  ti <- "treatmentbevacizumab"; ii <- "treatmentbevacizumab:periodlate"
  gl <- cf[ti, "coef"] + cf[ii, "coef"]; gse <- sqrt(v[ti, ti] + v[ii, ii] + 2 * v[ti, ii])
  ev_e <- sum(sp$.e[sp$period == "early"]); ev_l <- sum(sp$.e[sp$period == "late"])
  at_risk_late <- sum(sp$period == "late")
  rows[[length(rows)+1]] <- data.table(id = paste0("C_", g, c("_early", "_late")), model = "C_timevarying_bev_PFS_15mo",
      population = paste0(g, "-polarized"), endpoint = "PFS", term = c("bevacizumab, 0-15 mo", "bevacizumab, >15 mo"),
      HR = exp(c(cf[ti, "coef"], gl)), HR_lower = exp(c(cf[ti, "coef"] - 1.96 * cf[ti, "se(coef)"], gl - 1.96 * gse)),
      HR_upper = exp(c(cf[ti, "coef"] + 1.96 * cf[ti, "se(coef)"], gl + 1.96 * gse)),
      p_value = c(cf[ti, "Pr(>|z|)"], 2 * pnorm(-abs(gl / gse))), n = c(nrow(s), at_risk_late), events = c(ev_e, ev_l),
      period_interaction_p = cf[ii, "Pr(>|z|)"])
}

## ---------- D. treatment x polarization-group interaction, PFS (Results ICON7; Methods Survival) ----------
mi <- coxph(as.formula(sprintf("Surv(final_pfstm, final_pfsid) ~ treatment * pol_bin + %s", COVARS)), data = ana)
m0 <- coxph(as.formula(sprintf("Surv(final_pfstm, final_pfsid) ~ treatment + pol_bin + %s", COVARS)), data = ana)
t <- tidy(mi, "D_interaction_bev_x_polarization_PFS", "all FIGO III/IV patients", "PFS"); full[[length(full)+1]] <- t
r <- t[term == "treatmentbevacizumab:pol_binSecB"]
rows[[length(rows)+1]] <- r[, .(id = "D_interaction", model = "D_interaction_bev_x_polarization_PFS",
    population = "all FIGO III/IV patients", endpoint = "PFS", term, HR, HR_lower, HR_upper, p_value, n, events,
    interaction_LRT_p = anova(m0, mi)[2, "Pr(>|Chi|)"])]

## ---------- group sizes ----------
gs <- ana[, .(n = .N, PFS_events = sum(final_pfsid), OS_available = sum(!is.na(final_osid) & !is.na(final_ostm)),
              OS_events = sum(final_osid, na.rm = TRUE), FIGO_IV = sum(figo_stage == "IV"),
              suboptimal_or_inoperable = sum(debulking == "SUBOPT_INOP"), age_median = median(age_num),
              polarization_min = min(polarization_ucell), polarization_max = max(polarization_ucell)),
          by = .(group = pol_bin, arm = treatment)][order(group, arm)]
fwrite(gs, file.path(OUT, "st12_group_sizes.csv"))
fwrite(data.table(median_cut = med, n = nrow(ana), n_module = nrow(ana2)), file.path(OUT, "st12_meta.csv"))
fwrite(rbindlist(rows, fill = TRUE), file.path(OUT, "st12_models.csv"))
fwrite(rbindlist(full, fill = TRUE), file.path(OUT, "st12_full_coefficients.csv"))
fwrite(rbindlist(meds), file.path(OUT, "st12_km_medians.csv"))
cat("ST12_R_DONE\n")
'''
rfile = os.path.join(WORK, "st12_models.R")
with open(rfile, "w") as fh:
    fh.write(R_CODE)
subprocess.run(["Rscript", rfile], check=True, env={**os.environ, "ST12_WORK": WORK, "ST12_COHORT": COHORT})

M = pd.read_csv(f"{WORK}/st12_models.csv")
FULL = pd.read_csv(f"{WORK}/st12_full_coefficients.csv")
GS = pd.read_csv(f"{WORK}/st12_group_sizes.csv")
MED = pd.read_csv(f"{WORK}/st12_km_medians.csv")
META = pd.read_csv(f"{WORK}/st12_meta.csv").iloc[0]
assert int(META.n) == int(META.n_module) == 191, META

# ---- reproduction checks against the saved generator outputs ------------------------------------------
def close(a, b, tol=1e-8): return abs(float(a) - float(b)) < tol
t17 = pd.read_csv(f"{ICON7_TAB}/17_chemo_arm_continuous_polarization.tsv", sep="\t")
for ep in ["OS", "PFS"]:
    ref = t17[(t17.endpoint == ep)]
    got = FULL[(FULL.model == "A_arm_continuous_polarization") & (FULL.population == "standard arm") & (FULL.endpoint == ep)]
    for tm in ref.term:
        assert close(ref[ref.term == tm].HR.iloc[0], got[got.term == tm].HR.iloc[0]), (ep, tm)
t17b = pd.read_csv(f"{ICON7_TAB}/17b_continuous_polarization_per_arm.tsv", sep="\t")
for r in t17b.itertuples():
    g = M[M.id == f"A_{r.arm}_{r.endpoint}"].iloc[0]
    assert close(r.polarization_HR_per_1SD, g.HR) and close(r.polarization_p, g.p_value) and r.n == g.n, r
wg = pd.read_csv(f"{FIG7_OUT}/icon7_km_pfs_within_group_hr.csv")
for r in wg.itertuples():
    g = M[M.id == f"B_{r.group}"].iloc[0]
    assert close(r.HR, g.HR) and close(r.adj_p, g.p_value) and close(r.logrank_p, g.logrank_p) and r.n == g.n and r.events == g.events, r
# Benjamini-Hochberg (as R p.adjust "BH"): Figure 7m,n family = the 4 bevacizumab HRs of model C; the period
# interaction family = the 2 polarization groups
def bh(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p)[::-1]
    q = np.empty(m); q[o] = np.minimum.accumulate(p[o] * m / np.arange(m, 0, -1))
    return np.minimum(q, 1)
Mi = M.set_index("id")
C_IDS = [i for i in M.id if i.startswith("C_")]
assert len(C_IDS) == 4, C_IDS
Q_C = dict(zip(C_IDS, bh(Mi.loc[C_IDS, "p_value"].values)))
Q_INT = dict(zip(["SecA", "SecB"], bh([Mi.loc[f"C_{g}_early", "period_interaction_p"] for g in ["SecA", "SecB"]])))
tv = pd.read_csv(f"{FIG7_OUT}/icon7_timevarying_forest_15mo_data.csv")
for r in tv.itertuples():
    key = f"C_{r.group}_{'early' if r.period.startswith('0-') else 'late'}"
    g = M[M.id == key].iloc[0]
    assert close(r.HR, g.HR) and close(r.p, g.p_value) and close(r.period_interaction_p, g.period_interaction_p), r
    assert close(r.q, Q_C[key]) and close(r.period_interaction_q, Q_INT[r.group]), r
print("reproduction checks passed: 17_chemo_arm_continuous_polarization.tsv, 17b, icon7_km_pfs_within_group_hr.csv, "
      "icon7_timevarying_forest_15mo_data.csv")
D = Mi.loc["D_interaction"]
assert int(D.n) == 191, D
print(f"model D interaction (PFS): HR {D.HR:.4f} ({D.HR_lower:.4f}-{D.HR_upper:.4f}), Wald p {D.p_value:.5f}, "
      f"LRT p {D.interaction_LRT_p:.5f}, n {int(D.n)}, events {int(D.events)}")

# ---- format -------------------------------------------------------------------------------------------
COV = "FIGO stage (IV vs III), debulking (suboptimal/inoperable vs optimal), age (years)"
DESC = {
    "A_standard_OS": ("Results (ICON7)", "Carboplatin-paclitaxel alone (standard arm)", "Polarization score (continuous, per +1 SD): higher = SecB-leaning",
                      "OS HR 1.38 (1.02-1.86), p = 0.038"),
    "A_standard_PFS": ("context (not in text)", "Carboplatin-paclitaxel alone (standard arm)", "Polarization score (continuous, per +1 SD)", ""),
    "A_bevacizumab_OS": ("context (not in text)", "Carboplatin-paclitaxel + bevacizumab arm", "Polarization score (continuous, per +1 SD)", ""),
    "A_bevacizumab_PFS": ("context (not in text)", "Carboplatin-paclitaxel + bevacizumab arm", "Polarization score (continuous, per +1 SD)", ""),
    "B_SecA": ("Fig 7k", "SecA-polarized (polarization <= cohort median)", "Bevacizumab vs standard",
               "PFS HR 0.80 (0.48-1.35), p = 0.41; median 17.9 vs 21.6 mo"),
    "B_SecB": ("Fig 7l", "SecB-polarized (polarization > cohort median)", "Bevacizumab vs standard",
               "PFS HR 0.36 (0.22-0.59), p < 0.001; median 11.1 vs 18.7 mo"),
    "C_SecA_early": ("Fig 7n", "SecA-polarized", "Bevacizumab vs standard, 0-15 months", "interaction p = 0.19"),
    "C_SecA_late": ("Fig 7n", "SecA-polarized", "Bevacizumab vs standard, >15 months", "interaction p = 0.19"),
    "C_SecB_early": ("Fig 7m", "SecB-polarized", "Bevacizumab vs standard, 0-15 months", "HR 0.23 (0.12-0.44), p < 0.001"),
    "C_SecB_late": ("Fig 7m", "SecB-polarized", "Bevacizumab vs standard, >15 months", "HR 0.87, p = 0.77; interaction p = 0.022"),
    "D_interaction": ("Results (ICON7) and Methods (Survival)", "All FIGO III/IV patients",
                      "Bevacizumab x SecB-polarized (treatment x polarization-group interaction)",
                      "median split HR 0.44 (0.22-0.90), likelihood-ratio p = 0.025"),
}
MODEL = {"A_arm_continuous_polarization": "Cox: Surv ~ polarization_z + covariates (within arm)",
         "B_within_group_bev_PFS": "Cox: Surv ~ treatment + covariates (within polarization group)",
         "C_timevarying_bev_PFS_15mo": "Cox, counting process split at 15 mo (456 d): Surv(tstart, t, event) ~ treatment * period + covariates",
         "D_interaction_bev_x_polarization_PFS": "Cox: Surv ~ treatment * polarization group + covariates; likelihood-ratio test against the model without the interaction"}
def f3(v): return None if pd.isna(v) else round(float(v), 3)
def pg(v): return None if pd.isna(v) else float(f"{float(v):.3g}")
out = []
for r in M.itertuples():
    src, pop, contrast, txt = DESC[r.id]
    out.append({"row_id": r.id, "reported_in": src, "population": pop, "endpoint": r.endpoint,
                "contrast": contrast, "model": MODEL[r.model], "covariates": COV,
                "HR": f3(r.HR), "CI95": f"{r.HR_lower:.2f}-{r.HR_upper:.2f}", "HR_lower": f3(r.HR_lower), "HR_upper": f3(r.HR_upper),
                "p_value": pg(r.p_value), "n": int(r.n), "events": int(r.events),
                "median_PFS_standard_mo": f3(getattr(r, "median_standard_mo")), "median_PFS_bevacizumab_mo": f3(getattr(r, "median_bev_mo")),
                "logrank_p": pg(getattr(r, "logrank_p")), "period_interaction_p": pg(getattr(r, "period_interaction_p")),
                "FDR_BH_timevarying": pg(Q_C.get(r.id, np.nan)),
                "period_interaction_FDR_BH": pg(Q_INT[r.id.split("_")[1]]) if r.id.startswith("C_") else None,
                "interaction_LRT_p": pg(getattr(r, "interaction_LRT_p")),
                "manuscript_text_value": txt})
models = pd.DataFrame(out)
models.loc[models.row_id.str.startswith("C_") & models.row_id.str.endswith("_late"), "n"] = None
models["n_note"] = np.where(models.row_id.str.startswith("C_"),
                            np.where(models.row_id.str.endswith("_late"),
                                     "n = patients still at risk after 15 mo (see n_at_risk_after_15mo); events = events after 15 mo",
                                     "n = all patients in group; events = events within 0-15 mo"), "")
late_n = M.set_index("id").n
models["n_at_risk_after_15mo"] = [int(late_n[i]) if i.endswith("_late") else None for i in models.row_id]

TERM = {"polarization_z": "polarization (per +1 SD)",
        "treatmentbevacizumab": "bevacizumab vs standard (0-15 mo in model C; in SecA-polarized tumours in model D)",
        "pol_binSecB": "SecB-polarized vs SecA-polarized (standard arm, model D)",
        "treatmentbevacizumab:pol_binSecB": "bevacizumab x SecB-polarized interaction (model D)",
        "periodlate": "period >15 mo vs 0-15 mo (not estimable: a pure function of time, absorbed by the baseline hazard; NA by design)", "treatmentbevacizumab:periodlate": "bevacizumab x period (>15 mo) interaction",
        "figo_stageIV": "FIGO stage IV vs III", "debulkingSUBOPT_INOP": "debulking suboptimal/inoperable vs optimal",
        "age_num": "age (per year)"}
fc = FULL.copy()
fc.insert(4, "term_label", fc.term.map(TERM))
for c in ["coef", "se", "HR", "HR_lower", "HR_upper", "z"]: fc[c] = fc[c].round(4)
fc["p_value"] = fc.p_value.map(pg)
fc["model"] = fc.model.map(lambda m: f"{m}: {MODEL[m]}")

gs = GS.copy()
gs["group"] = gs.group + "-polarized"
for c in ["age_median", "polarization_min", "polarization_max"]: gs[c] = gs[c].round(4)
med = MED.copy(); med["group"] = med.group + "-polarized"
for c in ["median_mo", "median_lo", "median_hi"]: med[c] = med[c].round(2)
med = med[["group", "endpoint", "arm", "n", "events", "median_mo", "median_lo", "median_hi"]].rename(
    columns={"median_mo": "KM_median_mo", "median_lo": "KM_median_95CI_lo", "median_hi": "KM_median_95CI_hi"})

readme = [
    "Cohort: ICON7 phase III trial (GSE140082; carboplatin-paclitaxel +/- bevacizumab), pre-treatment tumour expression profiles, FIGO stage III/IV, "
    f"complete data for PFS, polarization, treatment, FIGO stage, debulking and age: n = {int(META.n)} (standard 95, bevacizumab 96).",
    f"Polarization: per-sample polarization_ucell = z-scored SecB UCell score minus z-scored SecA UCell score (z across the GSE140082 samples; precomputed in the ICON7 cohort table). Groups: SecB-polarized = polarization > cohort median "
    f"({META.median_cut:.4f}); SecA-polarized = <= median (SecA n = 96, SecB n = 95).",
    "Covariates in every model: " + COV + ". Survival times in days in the models; medians reported in months (days / 30.44).",
    "Model A (Results, ICON7: 'each 1 SD increase in polarization predicted worse OS under standard chemotherapy'): polarization as a CONTINUOUS variable, z-scored across all 191 "
    "patients (HR per +1 SD), fitted within the standard arm (n = 95). This is not a SecB-vs-SecA group contrast. PFS and bevacizumab-arm fits shown for context.",
    "Model B (Figure 7k,l): bevacizumab vs standard within each polarization group, PFS; adjusted HR (Wald p) plus unadjusted log-rank p (shown on the KM panels) and KM medians.",
    "Model C (Figure 7m,n): time-varying bevacizumab effect; follow-up split at 15 months (456 days; spanning the ~12 months of bevacizumab treatment and its wash-out); the >15-month HR and CI are the linear "
    "combination treatment + treatment:period with its covariance; period_interaction_p = Wald p of treatment:period.",
    "Model D (Results, ICON7; Methods, Survival): treatment x polarization-group interaction for PFS in all 191 patients, same covariates. HR = the "
    "interaction term (the bevacizumab HR in SecB-polarized tumours divided by that in SecA-polarized tumours); interaction_LRT_p = likelihood-ratio "
    "test against the model without the interaction (the p reported in the text); p_value = its Wald p.",
    "FDR_BH_timevarying: Benjamini-Hochberg q across the 4 model C bevacizumab HRs shown in Figure 7m,n; period_interaction_FDR_BH: "
    "Benjamini-Hochberg q of the period interaction across the 2 polarization groups. Blank for models A, B and D.",
    "Full_coefficients: in model C the period main effect is NA by design (period is a deterministic function of follow-up time and is absorbed by the Cox baseline hazard); "
    "only treatment and treatment x period are identifiable. KM_medians: blank = median not reached (or its upper confidence limit not reached).",
    "Sheets: Models = one row per reported contrast; Full_coefficients = every term of every model; KM_medians = PFS and OS medians per group x arm; Group_sizes = n, events and covariates per group x arm.",
    "Sources: Model A = 'davids side quests/ICON7/scripts/15_chemo_arm_analysis.R' section B.2 (verified identical to its results/tables/17_chemo_arm_continuous_polarization.tsv "
    "and 17b); Models B and C = figures/figure7/06_fig7klmn_icon7_survival_row.R (verified identical to icon7_km_pfs_within_group_hr.csv and icon7_timevarying_forest_15mo_data.csv); "
    "Model D = this builder, on the cohort, median split and covariates of 06_fig7klmn_icon7_survival_row.R. "
    "Generated by tables/12_tableS12_icon7_cox.py. R survival::coxph (Efron ties).",
]
write_table(12, "icon7_polarization_bevacizumab_cox",
            title="Supplementary Table 12. Cox regression of secretory polarization and bevacizumab benefit in the ICON7 trial (GSE140082; FIGO III/IV, n = 191; Figure 7k-n).",
            readme=readme,
            sheets=[("Models", models), ("Full_coefficients", fc), ("KM_medians", med), ("Group_sizes", gs)])
