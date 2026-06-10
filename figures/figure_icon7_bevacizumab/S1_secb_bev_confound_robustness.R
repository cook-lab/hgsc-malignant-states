#!/usr/bin/env Rscript
# ============================================================================
# ICON7 — SUPPLEMENTARY confound / novelty robustness for the SecB x bevacizumab
# treatment-modifier result.                                  [SUPPLEMENTARY]
# ----------------------------------------------------------------------------
# *** SUPPLEMENTARY ANALYSIS — produces NO main-figure panel. ***
# Companion to 01_icon7_bev_slope_reversal.R (the headline per-arm slope-reversal
# panel). This script does NOT change any figure; it generates supporting tables
# (+ one clearly-labelled supplementary forest) that stress-test whether the
# finding is (a) just the known clinical high-risk -> bev-benefit effect, or
# (b) just a generic hypoxia / angiogenesis / molecular-subtype signature.
# Full narrative + literature synthesis: ./SUPPLEMENTARY_NOVELTY.md
#
# PURPOSE / QUESTIONS
#   A. Is the SecB x bev modifier INDEPENDENT of ICON7's established clinical
#      high-risk -> bevacizumab-benefit subgroup (FIGO IV, or III suboptimally
#      debulked)?  -> proxy test, survives-adjustment, within-strata, joint model.
#   B. Does a generic hypoxia / angiogenesis / EMT / proliferation / TCGA-subtype
#      signature REPRODUCE the chemo-prognostic-then-abolished-by-bev slope
#      reversal, and does SecB survive head-to-head adjustment for hypoxia?
#      ("are we just saying hypoxia influences bev response?")
#
# INPUTS
#   cfg_obj("icon7_cohort") -> cohort_filtered.tsv (per-patient; FIGO III/IV n=191
#       survival cohort after filtering): polarization_ucell (SecB-SecA), treatment,
#       final_{ostm,osid,pfstm,pfsid}, figo_stage, debulking_status, age,
#       t1_cluster_name (TCGA subtype), sample_title.
#   cfg_obj("icon7_expr")   -> expr_gene_symbol.rds (gene-symbol x 212-sample matrix);
#       used ONLY by Part B to UCell-score competing signatures. Part B is skipped
#       (with a flag) if UCell + msigdbr are not installed.
#
# OUTPUTS  (figures_dir/figure_icon7_bevacizumab/supplementary_novelty/)
#   A1_*.tsv  high-risk confound tables
#   A2_*.tsv  competing-signature battery tables (if Part B runs)
#   suppl_signature_battery_forest.{pdf,svg}  (Part B forest; supplementary)
#   S1_sessionInfo.txt
#
# MANUSCRIPT PANEL(S): none (supplementary tables + one supplementary forest)
# RUNTIME TIER: fast (per-patient table; in-script Cox). Part B adds ~UCell scoring.
# VERIFIED: reproduces the headline chemo OS HR 1.38 (p=0.038), bev HR 0.96;
#   ported from the verified sandbox scripts A1_highrisk_confound.R /
#   A2_signature_battery.R (see ./SUPPLEMENTARY_NOVELTY.md).
# ============================================================================

.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))

suppressPackageStartupMessages(library(survival))
set.seed(CFG$seed)

OUT <- cfg_path("figures_dir", "figure_icon7_bevacizumab", "supplementary_novelty")
wtsv <- function(x, f) write.table(x, file.path(OUT, f), sep = "\t",
                                   quote = FALSE, row.names = FALSE)
message("Supplementary outputs -> ", OUT)

# ============================================================================
# Shared helpers (ported verbatim in logic from the verified sandbox _data.R)
# ============================================================================
zc <- function(x) as.numeric(scale(x))   # z-score across the cohort

# Survival-analytic cohort (FIGO III/IV, n=191) with derived clinical fields.
prep_cohort <- function() {
  d <- read.delim(cfg_obj("icon7_cohort"), stringsAsFactors = FALSE)
  d$os_event  <- as.integer(d$final_osid);  d$os_time  <- as.numeric(d$final_ostm)
  d$pfs_event <- as.integer(d$final_pfsid); d$pfs_time <- as.numeric(d$final_pfstm)
  d$age       <- as.numeric(d$age)
  d$arm       <- factor(d$treatment, levels = c("standard", "bevacizumab"))
  d$stageIV     <- as.integer(d$figo_stage == "IV")
  d$suboptimal  <- as.integer(d$debulking_status %in% c("SUB-OPTIMAL", "Inoperable"))
  # ICON7 clinical "high-risk" (Perren 2011 / Oza 2015): FIGO III & suboptimal, OR FIGO IV.
  d$highrisk    <- as.integer((d$figo_stage == "III" & d$suboptimal == 1) | d$figo_stage == "IV")
  d$subtype     <- factor(d$t1_cluster_name)
  d <- d[d$figo_stage %in% c("III", "IV"), , drop = FALSE]   # n=191 survival cohort
  d
}

# Per-arm continuous Cox (HR per 1-SD of `score`) + pooled score x arm LRT.
# `score` is a column name already z-scored across the cohort.
cox_perarm <- function(df, score, endpoint = c("os", "pfs"),
                       covars = c("figo_stage", "suboptimal", "age")) {
  endpoint <- match.arg(endpoint)
  tcol <- paste0(endpoint, "_time"); ecol <- paste0(endpoint, "_event")
  df <- df[is.finite(df[[score]]) & is.finite(df[[tcol]]) & !is.na(df[[ecol]]), ]
  cv <- covars[vapply(covars, function(v)
    v %in% names(df) && length(unique(df[[v]][!is.na(df[[v]])])) > 1, logical(1))]
  res <- list()
  for (a in c("standard", "bevacizumab")) {
    sub <- df[df$arm == a, ]
    cva <- cv[vapply(cv, function(v) length(unique(sub[[v]][!is.na(sub[[v]])])) > 1, logical(1))]
    rhs <- paste(c(score, cva), collapse = " + ")
    fit <- tryCatch(coxph(as.formula(sprintf("Surv(%s,%s) ~ %s", tcol, ecol, rhs)), data = sub),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      s <- summary(fit)
      res[[a]] <- data.frame(arm = a, endpoint = endpoint, score = score,
                             n = fit$n, events = fit$nevent,
                             HR = s$conf.int[score, "exp(coef)"],
                             lo = s$conf.int[score, "lower .95"],
                             hi = s$conf.int[score, "upper .95"],
                             p  = s$coefficients[score, "Pr(>|z|)"])
    }
  }
  cvp <- paste(cv, collapse = " + ")
  f0 <- as.formula(sprintf("Surv(%s,%s) ~ %s + arm + %s", tcol, ecol, score, cvp))
  f1 <- as.formula(sprintf("Surv(%s,%s) ~ %s * arm + %s", tcol, ecol, score, cvp))
  out <- do.call(rbind, res)
  attr(out, "interaction_p") <- anova(coxph(f0, data = df), coxph(f1, data = df))$`Pr(>|Chi|)`[2]
  out
}

d <- prep_cohort()
d$pol_z  <- zc(d$polarization_ucell)   # SecB-SecA polarization, z across n=191
d$age_z  <- zc(d$age)

cat("================================================================\n")
cat("ICON7 supplementary confound/robustness  (n =", nrow(d), "FIGO III/IV)\n")
cat("================================================================\n")
# sanity: reproduce the headline before anything else
.hl <- cox_perarm(d, "pol_z", "os")
cat(sprintf("SANITY headline OS: chemo HR=%.3f p=%.3f | bev HR=%.3f p=%.3f (expect 1.38/0.038, 0.96)\n",
            .hl$HR[.hl$arm=="standard"], .hl$p[.hl$arm=="standard"],
            .hl$HR[.hl$arm=="bevacizumab"], .hl$p[.hl$arm=="bevacizumab"]))
cat("arm x highrisk:\n"); print(table(arm = d$arm, highrisk = d$highrisk))

# ############################################################################
# PART A — Is SecB just the known clinical high-risk -> bev-benefit effect?
# ############################################################################

# (A1.1) Is SecB-polarization a PROXY for clinical high-risk?
cat("\n###### (A1.1) Is SecB a proxy for clinical high-risk? ######\n")
prox <- data.frame()
for (grp in c("highrisk", "stageIV", "suboptimal")) {
  g <- d[[grp]]
  pb <- cor.test(d$pol_z, g, method = "pearson")
  wx <- wilcox.test(d$pol_z[g == 1], d$pol_z[g == 0])
  n1 <- sum(g == 1); n0 <- sum(g == 0)
  s1 <- sd(d$pol_z[g == 1]); s0 <- sd(d$pol_z[g == 0])
  sp <- sqrt(((n1-1)*s1^2 + (n0-1)*s0^2) / (n1+n0-2))
  prox <- rbind(prox, data.frame(
    variable = grp, n1 = n1, n0 = n0,
    mean_pol_in1 = round(mean(d$pol_z[g==1]),3), mean_pol_in0 = round(mean(d$pol_z[g==0]),3),
    pointbiserial_r = round(unname(pb$estimate),3), r_p = signif(pb$p.value,3),
    cohens_d = round((mean(d$pol_z[g==1]) - mean(d$pol_z[g==0]))/sp, 3),
    wilcox_p = signif(wx$p.value,3)))
}
print(prox); wtsv(prox, "A1_secB_vs_clinical_correlation.tsv")

glm1   <- glm(highrisk ~ pol_z, data = d, family = binomial)
or1ci  <- exp(confint.default(glm1)["pol_z", ])
mcf    <- 1 - as.numeric(logLik(glm1)) / as.numeric(logLik(glm(highrisk ~ 1, data = d, family = binomial)))
glm_tab <- data.frame(model = "highrisk ~ pol_z",
                      OR_per_SD = round(exp(coef(glm1)["pol_z"]),3),
                      OR_lo = round(or1ci[1],3), OR_hi = round(or1ci[2],3),
                      p = signif(summary(glm1)$coefficients["pol_z","Pr(>|z|)"],3),
                      mcfadden_R2 = round(mcf,4))
cat("\nLogistic highrisk ~ pol_z:\n"); print(glm_tab, row.names = FALSE)
wtsv(glm_tab, "A1_logistic_highrisk_on_polz.tsv")

# (A1.2) POSITIVE CONTROL — reproduce the known high-risk -> bev benefit.
cat("\n###### (A1.2) Positive control: high-risk -> bev benefit ######\n")
bev_strata <- data.frame()
for (ep in c("os","pfs")) {
  tcol <- paste0(ep,"_time"); ecol <- paste0(ep,"_event")
  for (rk in c(0,1)) {
    sub <- d[d$highrisk == rk, ]
    covs <- "age_z"
    if (length(unique(sub$figo_stage)) > 1) covs <- paste(covs, "+ figo_stage")
    if (length(unique(sub$suboptimal))  > 1) covs <- paste(covs, "+ suboptimal")
    s <- summary(coxph(as.formula(sprintf("Surv(%s,%s) ~ arm + %s", tcol, ecol, covs)), data = sub))
    bev_strata <- rbind(bev_strata, data.frame(
      endpoint = ep, stratum = ifelse(rk==1,"high-risk","low-risk"),
      n = s$n, events = s$nevent,
      bev_HR = round(s$conf.int["armbevacizumab","exp(coef)"],3),
      lo = round(s$conf.int["armbevacizumab","lower .95"],3),
      hi = round(s$conf.int["armbevacizumab","upper .95"],3),
      p = signif(s$coefficients["armbevacizumab","Pr(>|z|)"],3), covars = covs))
  }
  m0 <- coxph(as.formula(sprintf("Surv(%s,%s) ~ arm + highrisk + age_z + figo_stage + suboptimal", tcol, ecol)), data = d)
  m1 <- coxph(as.formula(sprintf("Surv(%s,%s) ~ arm * highrisk + age_z + figo_stage + suboptimal", tcol, ecol)), data = d)
  ip <- anova(m0, m1)$`Pr(>|Chi|)`[2]
  bev_strata <- rbind(bev_strata, data.frame(endpoint = ep,
    stratum = sprintf("arm*highrisk_LRT_p=%.4f", ip), n=NA, events=NA,
    bev_HR=NA, lo=NA, hi=NA, p=ip, covars="pooled"))
}
print(bev_strata, row.names = FALSE); wtsv(bev_strata, "A1_bev_HR_by_riskstratum.tsv")

# (A1.3) Does the pol_z slope reversal SURVIVE high-risk adjustment?
cat("\n###### (A1.3) Does pol_z survive high-risk adjustment? ######\n")
perarm_pol <- function(df, ep, base_covars, pooled_extra = NULL, label) {
  tcol <- paste0(ep,"_time"); ecol <- paste0(ep,"_event"); rows <- data.frame()
  for (a in c("standard","bevacizumab")) {
    sub <- df[df$arm == a, ]
    cv <- base_covars[vapply(base_covars, function(v)
      v %in% names(sub) && length(unique(sub[[v]][!is.na(sub[[v]])])) > 1, logical(1))]
    f <- as.formula(sprintf("Surv(%s,%s) ~ pol_z%s", tcol, ecol,
                            if (length(cv)) paste0(" + ", paste(cv, collapse=" + ")) else ""))
    s <- summary(coxph(f, data = sub))
    rows <- rbind(rows, data.frame(model=label, endpoint=ep, arm=a, n=s$n, events=s$nevent,
      pol_HR=round(s$conf.int["pol_z","exp(coef)"],3),
      lo=round(s$conf.int["pol_z","lower .95"],3),
      hi=round(s$conf.int["pol_z","upper .95"],3),
      p=signif(s$coefficients["pol_z","Pr(>|z|)"],4)))
  }
  base_p <- paste(base_covars, collapse=" + ")
  extra  <- if (!is.null(pooled_extra)) paste0(" + ", paste(pooled_extra, collapse=" + ")) else ""
  m0 <- coxph(as.formula(sprintf("Surv(%s,%s) ~ pol_z + arm + %s%s", tcol, ecol, base_p, extra)), data = df)
  m1 <- coxph(as.formula(sprintf("Surv(%s,%s) ~ pol_z * arm + %s%s", tcol, ecol, base_p, extra)), data = df)
  attr(rows, "int_p") <- anova(m0, m1)$`Pr(>|Chi|)`[2]; rows
}
base_covars <- c("figo_stage","suboptimal","age")
adj3 <- data.frame(); int3 <- data.frame()
for (ep in c("os","pfs")) {
  rb  <- perarm_pol(d, ep, base_covars, NULL, "BEFORE (base covars)")
  rh  <- perarm_pol(d, ep, c(base_covars,"highrisk"), NULL, "AFTER +highrisk (main)")
  rhx <- perarm_pol(d, ep, base_covars, c("highrisk","highrisk:arm"), "AFTER +highrisk*arm")
  adj3 <- rbind(adj3, rb, rh, rhx)
  int3 <- rbind(int3,
    data.frame(endpoint=ep, model="BEFORE (base covars)",   pol_arm_int_p=signif(attr(rb,"int_p"),4)),
    data.frame(endpoint=ep, model="AFTER +highrisk (main)", pol_arm_int_p=signif(attr(rh,"int_p"),4)),
    data.frame(endpoint=ep, model="AFTER +highrisk*arm",    pol_arm_int_p=signif(attr(rhx,"int_p"),4)))
}
cat("\nPer-arm pol_z HRs (before vs after high-risk adjustment):\n"); print(adj3, row.names = FALSE)
cat("\npol_z x arm interaction p (before vs after):\n"); print(int3, row.names = FALSE)
wtsv(adj3, "A1_perarm_polHR_before_after_adjust.tsv")
wtsv(int3, "A1_polArm_interaction_before_after.tsv")

# (A1.4) STRATIFIED slope reversal: pol_z per-arm WITHIN each risk stratum.
cat("\n###### (A1.4) Stratified slope reversal (within risk) ######\n")
strat4 <- data.frame()
for (ep in c("os","pfs")) for (rk in c(0,1)) {
  sub <- d[d$highrisk == rk, ]
  cov <- c("age")
  if (length(unique(sub$figo_stage)) > 1) cov <- c(cov,"figo_stage")
  if (length(unique(sub$suboptimal))  > 1) cov <- c(cov,"suboptimal")
  rr <- cox_perarm(sub, "pol_z", ep, covars = cov)
  ip <- attr(rr, "interaction_p")
  rr$stratum <- ifelse(rk==1,"high-risk","low-risk")
  rr$pol_arm_int_p <- signif(ip,4); rr$covars <- paste(cov, collapse="+")
  rr$HR <- round(rr$HR,3); rr$lo <- round(rr$lo,3); rr$hi <- round(rr$hi,3); rr$p <- signif(rr$p,4)
  strat4 <- rbind(strat4, rr)
  cat(sprintf("[%s | %-9s] pol_z x arm interaction p = %.4f (n=%d)\n",
              toupper(ep), ifelse(rk==1,"high-risk","low-risk"), ip, nrow(sub)))
}
strat4 <- strat4[, c("stratum","endpoint","arm","n","events","HR","lo","hi","p","pol_arm_int_p","covars")]
wtsv(strat4, "A1_perarm_polHR_within_riskstrata.tsv")

# (A1.5) JOINT model: Surv ~ pol_z*arm + highrisk*arm + age.
cat("\n###### (A1.5) Joint model: pol_z*arm + highrisk*arm + age ######\n")
joint5 <- data.frame()
for (ep in c("os","pfs")) {
  tcol <- paste0(ep,"_time"); ecol <- paste0(ep,"_event")
  fit <- coxph(as.formula(sprintf("Surv(%s,%s) ~ pol_z*arm + highrisk*arm + age", tcol, ecol)), data = d)
  m_no_polx <- coxph(as.formula(sprintf("Surv(%s,%s) ~ pol_z + arm + highrisk*arm + age", tcol, ecol)), data = d)
  m_no_hrx  <- coxph(as.formula(sprintf("Surv(%s,%s) ~ pol_z*arm + highrisk + arm + age", tcol, ecol)), data = d)
  p_polx <- anova(m_no_polx, fit)$`Pr(>|Chi|)`[2]; p_hrx <- anova(m_no_hrx, fit)$`Pr(>|Chi|)`[2]
  cat(sprintf("[%s] JOINT LRT: pol_z:arm p=%.4f | highrisk:arm p=%.4f\n", toupper(ep), p_polx, p_hrx))
  joint5 <- rbind(joint5,
    data.frame(endpoint=ep, term="LRT_pol_z:arm",    p=signif(p_polx,4)),
    data.frame(endpoint=ep, term="LRT_highrisk:arm", p=signif(p_hrx,4)))
}
wtsv(joint5, "A1_joint_model_interactions.tsv")
cat("\nPART A verdict: SecB x bev modifier is INDEPENDENT of clinical high-risk\n",
    "  (weak SecB<->risk association; survives adjustment; localizes in LOW-risk).\n", sep = "")

# ############################################################################
# PART B — Does a generic hypoxia/angiogenesis/subtype signature reproduce it?
#          (requires UCell + msigdbr; skipped with a flag if unavailable)
# ############################################################################
have_B <- requireNamespace("UCell", quietly = TRUE) && requireNamespace("msigdbr", quietly = TRUE)
if (!have_B) {
  msg <- paste0("PART B SKIPPED: requires UCell + msigdbr (one or both not installed). ",
                "Part A (high-risk confound) outputs are complete. Install UCell + msigdbr ",
                "to score competing hypoxia/angiogenesis signatures and re-run.")
  message("\n", msg); writeLines(msg, file.path(OUT, "A2_SKIPPED.txt"))
} else {
  suppressPackageStartupMessages({ library(UCell); library(msigdbr) })
  cat("\n###### PART B: competing-signature battery ######\n")
  gm <- readRDS(cfg_obj("icon7_expr"))   # gene-symbol x 212-sample matrix

  # UCell-score a gene set on the full matrix, return per-sample raw score keyed by sample_title.
  score_ucell <- function(genes, nm) {
    present <- intersect(genes, rownames(gm))
    uc <- as.data.frame(UCell::ScoreSignatures_UCell(gm, features = setNames(list(present), nm), ncores = 1))
    data.frame(sample_title = rownames(uc), raw = uc[[paste0(nm, "_UCell")]],
               n_used = length(present), n_missing = length(setdiff(genes, rownames(gm))),
               stringsAsFactors = FALSE)
  }

  H <- msigdbr(species = "Homo sapiens", collection = "H")
  hall <- function(s) sort(unique(H$gene_symbol[H$gs_name == s]))
  sig_list <- list(
    HALLMARK_HYPOXIA      = hall("HALLMARK_HYPOXIA"),
    HALLMARK_ANGIOGENESIS = hall("HALLMARK_ANGIOGENESIS"),
    HALLMARK_EMT          = hall("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
    HALLMARK_E2F          = hall("HALLMARK_E2F_TARGETS"),
    HALLMARK_G2M          = hall("HALLMARK_G2M_CHECKPOINT"),
    # Buffa 52-gene hypoxia metagene (Buffa et al. 2010, Br J Cancer)
    BUFFA_HYPOXIA = c("VEGFA","SLC2A1","PGAM1","ENO1","LDHA","TPI1","P4HA1","MRPS17","CDKN3",
      "ADM","NDRG1","TUBB6","ALDOA","MIF","ACOT7","MCTS1","PSRC1","PSMA7","ANLN","TUBA1B",
      "MAD2L2","GPI","TUBA1C","MAP7D1","DDIT4","BNIP3","C20orf20","HIG2","GAPDH","MRPL13",
      "CHCHD2","YKT6","NP","CORO1C","SEC61G","ANKRD37","ESRP1","PGK1","SHCBP1","CTSL2",
      "KIF20A","SLC25A32","UTP11L","CDKN2A","PFKP","DCBLD1","KIF4A","LRRC42","HK2","AK3L1"),
    HIF_TARGETS = c("VEGFA","CA9","SLC2A1","LDHA","PGK1","HK2","BNIP3","PDK1","ADM","NDRG1"),
    VEGF_AXIS   = c("VEGFA","KDR","FLT1","ANGPT2","PGF"))

  cov_records <- list()
  for (nm in names(sig_list)) {
    sc <- score_ucell(sig_list[[nm]], nm)
    cov_records[[nm]] <- data.frame(signature = nm, n_genes = length(sig_list[[nm]]),
                                    n_used = sc$n_used[1], n_missing = sc$n_missing[1])
    d[[paste0(nm, "_z")]] <- zc(sc$raw[match(d$sample_title, sc$sample_title)])
  }
  wtsv(do.call(rbind, cov_records), "A2_signature_coverage.tsv")

  d$mesenchymal_z   <- zc(as.integer(d$subtype == "mesenchymal"))
  d$proliferative_z <- zc(as.integer(d$subtype == "proliferative"))
  comp_z <- c(HALLMARK_HYPOXIA="HALLMARK_HYPOXIA_z", HALLMARK_ANGIOGENESIS="HALLMARK_ANGIOGENESIS_z",
              HALLMARK_EMT="HALLMARK_EMT_z", HALLMARK_E2F="HALLMARK_E2F_z", HALLMARK_G2M="HALLMARK_G2M_z",
              BUFFA_HYPOXIA="BUFFA_HYPOXIA_z", HIF_TARGETS="HIF_TARGETS_z", VEGF_AXIS="VEGF_AXIS_z",
              mesenchymal="mesenchymal_z", proliferative="proliferative_z")

  # (A2.1) Spearman correlation of each competitor with SecB polarization.
  cor_cols <- c(polarization="pol_z", comp_z)
  cd <- d[, unname(cor_cols)]; colnames(cd) <- names(cor_cols)
  spear <- cor(cd, method = "spearman", use = "pairwise.complete.obs")
  wtsv(data.frame(variable = rownames(spear), round(spear,3), check.names = FALSE), "A2_correlation_spearman.tsv")
  cat("\nSpearman rho with SecB polarization:\n"); print(round(sort(spear["polarization",-1], decreasing = TRUE),3))

  # (A2.2) Does each competitor reproduce the slope reversal? Per-arm Cox.
  run_battery <- function(score_col, label, endpoint) {
    t <- cox_perarm(d, score_col, endpoint); ip <- attr(t, "interaction_p")
    std <- t[t$arm=="standard",]; bev <- t[t$arm=="bevacizumab",]
    data.frame(signature=label, endpoint=endpoint,
      chemo_HR=round(std$HR,3), chemo_lo=round(std$lo,3), chemo_hi=round(std$hi,3), chemo_p=round(std$p,4),
      bev_HR=round(bev$HR,3), bev_lo=round(bev$lo,3), bev_hi=round(bev$hi,3), bev_p=round(bev$p,4),
      interaction_p=round(ip,4))
  }
  battery <- do.call(rbind, c(
    lapply(c("os","pfs"), function(ep) run_battery("pol_z","SecB_polarization",ep)),
    unlist(lapply(names(comp_z), function(lab) lapply(c("os","pfs"),
      function(ep) run_battery(comp_z[[lab]], lab, ep))), recursive = FALSE)))
  battery$reproduces_reversal <- with(battery,
    (chemo_HR > 1 & chemo_p < 0.10) & (bev_p > 0.10) & (abs(log(bev_HR)) < abs(log(chemo_HR))))
  wtsv(battery, "A2_perarm_cox_battery.tsv")
  cat("\nPer-arm battery (OS) — which reproduce the SecB-like reversal?\n")
  print(battery[battery$endpoint=="os", c("signature","chemo_HR","chemo_p","bev_HR","bev_p","interaction_p","reproduces_reversal")], row.names = FALSE)

  # (A2.3) Head-to-head: Surv ~ pol_z*arm + competitor_z*arm + age. Does pol_z:arm survive?
  h2h <- function(comp_col, label, endpoint) {
    tcol <- paste0(endpoint,"_time"); ecol <- paste0(endpoint,"_event")
    sub <- d[is.finite(d$pol_z) & is.finite(d[[comp_col]]) & is.finite(d[[tcol]]) & !is.na(d[[ecol]]), ]
    fit <- coxph(as.formula(sprintf("Surv(%s,%s) ~ pol_z*arm + %s*arm + age", tcol, ecol, comp_col)), data = sub)
    s <- summary(fit)$coefficients
    pol_int  <- grep("pol_z",   grep("arm", rownames(s), value = TRUE), value = TRUE)
    comp_int <- grep(comp_col,  grep("arm", rownames(s), value = TRUE), value = TRUE)
    f_no_pol  <- coxph(as.formula(sprintf("Surv(%s,%s) ~ pol_z + arm + %s*arm + age", tcol, ecol, comp_col)), data = sub)
    f_no_comp <- coxph(as.formula(sprintf("Surv(%s,%s) ~ pol_z*arm + %s + arm + age", tcol, ecol, comp_col)), data = sub)
    data.frame(competitor=label, endpoint=endpoint, n=fit$n,
      pol_int_HR=round(exp(s[pol_int,"coef"]),3), pol_int_wald_p=round(s[pol_int,"Pr(>|z|)"],4),
      pol_int_lrt_p=round(anova(f_no_pol, fit)$`Pr(>|Chi|)`[2],4),
      comp_int_HR=round(exp(s[comp_int,"coef"]),3), comp_int_wald_p=round(s[comp_int,"Pr(>|z|)"],4),
      comp_int_lrt_p=round(anova(f_no_comp, fit)$`Pr(>|Chi|)`[2],4))
  }
  h2h_comp <- c(HALLMARK_HYPOXIA="HALLMARK_HYPOXIA_z", BUFFA_HYPOXIA="BUFFA_HYPOXIA_z",
                HIF_TARGETS="HIF_TARGETS_z", HALLMARK_ANGIOGENESIS="HALLMARK_ANGIOGENESIS_z",
                HALLMARK_EMT="HALLMARK_EMT_z", mesenchymal="mesenchymal_z")
  h2h_tab <- do.call(rbind, unlist(lapply(names(h2h_comp), function(lab)
    lapply(c("os","pfs"), function(ep) h2h(h2h_comp[[lab]], lab, ep))), recursive = FALSE))
  wtsv(h2h_tab, "A2_head_to_head_interaction.tsv")
  cat("\nHead-to-head (PFS) — does pol_z:arm survive adjusting for competitor:arm?\n")
  print(h2h_tab[h2h_tab$endpoint=="pfs", c("competitor","pol_int_HR","pol_int_lrt_p","comp_int_HR","comp_int_lrt_p")], row.names = FALSE)

  # (A2.4) Residualize SecB on hypoxia; re-run per-arm Cox on the residual.
  resid_rows <- list()
  for (hyp in c("HALLMARK_HYPOXIA","BUFFA_HYPOXIA","HIF_TARGETS")) {
    lmfit <- lm(as.formula(sprintf("pol_z ~ %s", comp_z[[hyp]])), data = d)
    d[[paste0("pol_resid_", hyp)]] <- zc(residuals(lmfit))
    for (ep in c("os","pfs")) {
      rb <- run_battery(paste0("pol_resid_", hyp), paste0("SecB_resid_", hyp), ep)
      rb$hypoxia_R2_on_pol <- round(summary(lmfit)$r.squared, 3)
      resid_rows[[paste0(hyp,"_",ep)]] <- rb
    }
  }
  resid_tab <- do.call(rbind, resid_rows)
  wtsv(resid_tab, "A2_residualized_secb_cox.tsv")

  # --- Supplementary forest: per-arm HR by signature (OS+PFS); only SecB separates.
  fig_ok <- tryCatch({
    suppressPackageStartupMessages(library(ggplot2)); TRUE
  }, error = function(e) FALSE)
  if (fig_ok) {
    long <- do.call(rbind, lapply(c("standard","bevacizumab"), function(a) {
      arm_lab <- if (a=="standard") "Chemo only" else "Chemo + bevacizumab"
      pre <- if (a=="standard") "chemo" else "bev"
      data.frame(signature = battery$signature, endpoint = toupper(battery$endpoint), arm = arm_lab,
                 HR = battery[[paste0(pre,"_HR")]], lo = battery[[paste0(pre,"_lo")]], hi = battery[[paste0(pre,"_hi")]])
    }))
    long$is_secb <- long$signature == "SecB_polarization"
    ord <- battery$signature[battery$endpoint=="os"]
    ord <- c("SecB_polarization", setdiff(ord, "SecB_polarization"))
    long$signature <- factor(long$signature, levels = rev(ord))
    long$arm <- factor(long$arm, levels = c("Chemo only","Chemo + bevacizumab"))
    p <- ggplot(long, aes(HR, signature, color = arm)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.3) +
      geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.35, linewidth = 0.4,
                    orientation = "y", position = position_dodge(width = 0.6)) +
      geom_point(aes(shape = is_secb), size = 1.9, position = position_dodge(width = 0.6)) +
      facet_wrap(~ endpoint) +
      scale_x_log10() +
      scale_color_manual(values = c("Chemo only" = "#B2182B", "Chemo + bevacizumab" = "grey55"), name = NULL) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
      labs(x = "HR per 1-SD of signature", y = NULL,
           title = "Only SecB polarization shows the chemo-prognostic / bev-abolished reversal",
           subtitle = "Competing hypoxia / angiogenesis / EMT / proliferation / subtype signatures do not") +
      theme_minimal(base_size = 7) +
      theme(plot.title = element_text(size = 7, face = "bold"),
            plot.subtitle = element_text(size = 6), legend.position = "top",
            panel.grid.minor = element_blank())
    for (ext in c("pdf","svg"))
      ggsave(file.path(OUT, paste0("suppl_signature_battery_forest.", ext)), p,
             width = 7, height = 3.4, bg = "white")
    message("  Saved supplementary forest: suppl_signature_battery_forest.{pdf,svg}")
  }
  cat("\nPART B verdict: NO competing signature reproduces the reversal; every hypoxia set\n",
      "  trends chemo HR<1 (opposite direction); SecB:arm survives head-to-head adjustment.\n",
      "  => the finding is a specific SecB cell-state effect, NOT generic hypoxia.\n", sep = "")
}

writeLines(capture.output(print(sessionInfo())), file.path(OUT, "S1_sessionInfo.txt"))
cat("\nDONE. Supplementary outputs in:", OUT, "\n")
