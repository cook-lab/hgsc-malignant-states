# ============================================================================
# 01_tcga_external_validation.R
# ----------------------------------------------------------------------------
# PURPOSE: Pre-registered external validation of the SecA/SecB polarization signature in TCGA-OV: cohort lock, UCell + CIBERSORTx score construction, primary Cox tests x adjustment models (BH + Bonferroni).
#
# INPUTS:
#   - data/TCGA_data/ (TPM matrix, external clinical, CIBERSORTx)
#   - SecA/SecB signatures from shared/signatures.yml
#
# OUTPUTS:
#   - output/40_tcga_validation/ Cox tables + run log
#
# MANUSCRIPT PANEL(S): Cross-validates atlas 22d (Fig 7E/F/G).
# RUNTIME TIER: heavy
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

# SecA/SecB signatures loaded from the shared source of truth (noBCAM 7-gene set)
.sigs      <- yaml::read_yaml(file.path(here, "..", "..", "shared", "signatures.yml"))
secA_genes <- .sigs$SecA
secB_genes <- .sigs$SecB
sig_list   <- list(SecA = secA_genes, SecB = secB_genes)

suppressPackageStartupMessages({
  library(UCell)
  library(survival)
  library(readxl)
  library(stringr)
})

# --- Paths -------------------------------------------------------------------

tcga_dir   <- file.path(data_dir, "TCGA_data")
ext_dir    <- file.path(tcga_dir, "external_clinical")
step_dir   <- file.path(out_dir, "40_tcga_validation")
fig40_dir  <- file.path(step_dir, "figures")
dir.create(step_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig40_dir, showWarnings = FALSE, recursive = TRUE)

run_log <- file.path(step_dir, "run.log")
sink_log <- function(msg) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg, "\n",
                              file = run_log, append = TRUE)
file.create(run_log, showWarnings = FALSE)
sink_log("Phase 40 started")

# --- Frozen signature definition (no-BCAM production set) --------------------

# --- Load TCGA expression matrix --------------------------------------------

sink_log("Loading TCGA RNA matrix")
tcga_rna <- fread(file.path(tcga_dir, "tcga_data_plotting.csv"))

# First 3 cols: row index, hgnc_symbol, ensembl
expr_cols <- setdiff(names(tcga_rna), c("V1", "hgnc_symbol", "X"))
gene_sym  <- tcga_rna$hgnc_symbol
expr_mat  <- as.matrix(tcga_rna[, ..expr_cols])
rownames(expr_mat) <- gene_sym

# Drop genes with NA / empty symbol; collapse duplicate symbols by max-mean
keep_g <- !is.na(rownames(expr_mat)) & nzchar(rownames(expr_mat))
expr_mat <- expr_mat[keep_g, , drop = FALSE]
if (any(duplicated(rownames(expr_mat)))) {
  rmean <- rowMeans(expr_mat, na.rm = TRUE)
  ord   <- order(rmean, decreasing = TRUE)
  expr_mat <- expr_mat[ord, , drop = FALSE]
  expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), , drop = FALSE]
}
storage.mode(expr_mat) <- "numeric"
sink_log(sprintf("Loaded expression matrix: %d genes x %d aliquots",
                 nrow(expr_mat), ncol(expr_mat)))

stopifnot(all(c(secA_genes, secB_genes) %in% rownames(expr_mat)))

# Aliquot barcodes -> patient barcodes
parse_patient <- function(x) {
  x <- gsub("\\.", "-", x)
  substr(x, 1, 12)
}
aliquot_bc  <- gsub("\\.", "-", colnames(expr_mat))
patient_bc  <- parse_patient(colnames(expr_mat))

# Aliquot type code (.01 = primary tumour)
aliquot_type <- substr(sapply(strsplit(aliquot_bc, "-"), `[`, 4), 1, 2)

aliq_dt <- data.table(
  aliquot_id_csv = colnames(expr_mat),
  aliquot_bc     = aliquot_bc,
  patient_bc     = patient_bc,
  aliquot_type   = aliquot_type
)
sink_log(sprintf("Aliquot type distribution: %s",
                 paste(names(table(aliquot_type)), table(aliquot_type),
                       sep = "=", collapse = ", ")))

# --- Load TCGA clinical / mutation file --------------------------------------

sink_log("Loading TCGA clinical (tcga_meta.csv)")
clin <- fread(file.path(tcga_dir, "tcga_meta.csv"),
              select = c("PATIENT.NUMBER", "Age", "Grade", "Clinical.Stage",
                         "Therapy.Response..Coded.", "Platinum",
                         "OS", "OS..Months.", "PFS", "PFS..Months.", "TP53.x"))
setnames(clin, "PATIENT.NUMBER", "patient_bc")
clin[, patient_bc := substr(gsub("-", "-", patient_bc), 1, 12)]

# --- Load CIBERSORTx ---------------------------------------------------------

sink_log("Loading CIBERSORTx Job21 results")
csx <- fread(file.path(tcga_dir, "CIBERSORTx_Job21_Results.csv"))
setnames(csx, "Mixture", "aliquot_id_csv")
csx[, aliquot_bc := gsub("\\.", "-", aliquot_id_csv)]
csx[, patient_bc := parse_patient(aliquot_bc)]

# --- Load ABSOLUTE purity ----------------------------------------------------

sink_log("Loading ABSOLUTE purity (TCGA_mastercalls.abs_tables_JSedit.fixed.txt)")
abs_purity <- fread(file.path(ext_dir, "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"))
# Columns: array, sample, call status, purity, ploidy, ...
setnames(abs_purity, "array", "patient_aliquot")  # array col is TCGA-XX-XXXX-01
abs_purity[, patient_bc := substr(patient_aliquot, 1, 12)]
abs_purity[, aliquot_short := substr(patient_aliquot, 14, 15)]  # "01"
# Keep .01 only (primary tumour)
abs_purity_01 <- abs_purity[aliquot_short == "01"]
abs_purity_pt <- abs_purity_01[, .(absolute_purity = max(purity, na.rm = TRUE)),
                               by = patient_bc]
abs_purity_pt[is.infinite(absolute_purity), absolute_purity := NA_real_]
sink_log(sprintf("ABSOLUTE purity: %d patients (any cancer type)", nrow(abs_purity_pt)))

# --- Load TCGA-CDR (residual disease + better OS/PFS) ------------------------

sink_log("Loading TCGA-CDR (Liu 2018) for OV residual disease")
cdr <- as.data.table(read_xlsx(file.path(ext_dir, "TCGA-CDR-SupplementalTableS1.xlsx"),
                               sheet = "TCGA-CDR"))
cdr_ov <- cdr[type == "OV"]
sink_log(sprintf("CDR OV: %d patients", nrow(cdr_ov)))
# Residual disease coding (TCGA-CDR field is `residual_tumor` in some releases;
# fall back to `tumor_residual_disease` if needed)
if ("residual_tumor" %in% names(cdr_ov)) {
  cdr_ov[, residual_disease := residual_tumor]
} else if ("tumor_residual_disease" %in% names(cdr_ov)) {
  cdr_ov[, residual_disease := tumor_residual_disease]
} else {
  sink_log("WARNING: CDR has no residual_tumor / tumor_residual_disease column; residual disease will be NA")
  cdr_ov[, residual_disease := NA_character_]
}
cdr_ov_keep <- cdr_ov[, .(patient_bc = bcr_patient_barcode,
                          residual_disease)]

# --- Cohort lock with CONSORT-style flowchart --------------------------------

sink_log("Locking cohort")
flow <- data.table(step = character(), n = integer(), description = character())
add_flow <- function(label, n, desc = "") {
  flow <<- rbind(flow, data.table(step = label, n = n, description = desc))
}

# Start from RNA aliquots
add_flow("Step 0: RNA aliquots in tcga_data_plotting.csv", nrow(aliq_dt))

# Filter 1: primary tumour (.01) only
aliq_dt <- aliq_dt[aliquot_type == "01"]
add_flow("Step 1: .01 aliquots (primary tumour)", nrow(aliq_dt))

# Collapse to patient (one aliquot per patient — take first if duplicates)
setorder(aliq_dt, patient_bc, aliquot_id_csv)
aliq_dt <- aliq_dt[, .SD[1], by = patient_bc]
add_flow("Step 1b: collapse to one aliquot per patient (first .01 listed)",
         nrow(aliq_dt))

# Merge clinical
patient_dt <- merge(aliq_dt, clin, by = "patient_bc", all.x = TRUE)

# Filter 2: high-grade serous (G2 / G3)
patient_dt[, Grade_clean := toupper(gsub('"', '', Grade))]
patient_dt <- patient_dt[Grade_clean %in% c("G2", "G3")]
add_flow("Step 2: Grade in {G2, G3}", nrow(patient_dt),
         "Drops G1, GB, GX, missing")

# Recode OS / PFS event from text
recode_event <- function(x) {
  x <- as.character(x)
  ifelse(grepl("DECEASED|PROGRESSION", x, ignore.case = TRUE), 1L,
         ifelse(grepl("LIVING|CENSORED|DISEASEFREE|DISEASE FREE", x, ignore.case = TRUE), 0L,
                NA_integer_))
}
patient_dt[, OS_event  := recode_event(OS)]
patient_dt[, PFS_event := recode_event(PFS)]
patient_dt[, OS_months  := as.numeric(`OS..Months.`)]
patient_dt[, PFS_months := as.numeric(`PFS..Months.`)]

# Filter 3: valid OS and PFS time
patient_dt <- patient_dt[!is.na(OS_months) & !is.na(PFS_months) &
                          !is.na(OS_event)  & !is.na(PFS_event)]
add_flow("Step 3: valid OS_months/PFS_months/events", nrow(patient_dt))

# Merge ABSOLUTE purity
patient_dt <- merge(patient_dt, abs_purity_pt, by = "patient_bc", all.x = TRUE)

# Filter 4: ABSOLUTE purity >= 0.4
n_pre_purity <- nrow(patient_dt)
n_missing_purity <- sum(is.na(patient_dt$absolute_purity))
patient_dt <- patient_dt[!is.na(absolute_purity) & absolute_purity >= 0.4]
add_flow("Step 4: ABSOLUTE purity >= 0.4", nrow(patient_dt),
         sprintf("Pre-filter %d, missing purity %d, dropped %d",
                 n_pre_purity, n_missing_purity,
                 n_pre_purity - n_missing_purity - nrow(patient_dt)))

# Merge residual disease (kitchen-sink covariate; do NOT filter on it)
patient_dt <- merge(patient_dt, cdr_ov_keep, by = "patient_bc", all.x = TRUE)

# Standardize stage: collapse to {II_or_III, IV} (II rare in TCGA-OV)
patient_dt[, stage_clean := toupper(gsub('"', '', Clinical.Stage))]
patient_dt[, stage_iv := fifelse(grepl("STAGE IV", stage_clean), "IV",
                          fifelse(grepl("STAGE III", stage_clean), "III",
                          fifelse(grepl("STAGE II", stage_clean), "II",
                          fifelse(grepl("STAGE I", stage_clean), "I",
                                  NA_character_))))]
patient_dt[, stage_grp := fifelse(stage_iv == "IV", "IV", "II_III")]
patient_dt[, stage_grp := factor(stage_grp, levels = c("II_III", "IV"))]

# Standardize platinum
patient_dt[, platinum_clean := toupper(gsub('"', '', Platinum))]
patient_dt[, platinum_resist := fifelse(grepl("RESISTANT|REFRACTORY", platinum_clean), 1L,
                                 fifelse(grepl("SENSITIVE", platinum_clean), 0L,
                                         NA_integer_))]
patient_dt[, platinum_resist := factor(platinum_resist, levels = c(0L, 1L),
                                       labels = c("sensitive", "resistant"))]

# Residual disease as ordered factor: no_residual < <=1cm < >1cm
patient_dt[, residual_clean := toupper(gsub('"', '', residual_disease))]
patient_dt[, residual_grp := fifelse(grepl("NO MACRO|NO RESIDUAL", residual_clean), "no_residual",
                              fifelse(grepl("1-10|<=10|<10|<=1\\s?CM|MICRO", residual_clean), "le_10mm",
                              fifelse(grepl(">10|>1\\s?CM", residual_clean), "gt_10mm",
                                      NA_character_)))]
patient_dt[, residual_grp := factor(residual_grp,
                                    levels = c("no_residual", "le_10mm", "gt_10mm"),
                                    ordered = TRUE)]

# TP53 binary (NA-aware): 1 if mutation called, 0 if explicitly WT, NA otherwise
patient_dt[, tp53_status := fifelse(!is.na(TP53.x) & TP53.x != "" & TP53.x != "NA", 1L, NA_integer_)]
# (NA stays NA — not the same as WT for this dataset)

# Grade as factor
patient_dt[, grade_grp := factor(Grade_clean, levels = c("G2", "G3"))]

# Cap survival at 5 years for primary tests
patient_dt[, OS5_months  := pmin(OS_months, 60)]
patient_dt[, OS5_event   := fifelse(OS_event == 1L & OS_months <= 60, 1L, 0L)]
patient_dt[, PFS5_months := pmin(PFS_months, 60)]
patient_dt[, PFS5_event  := fifelse(PFS_event == 1L & PFS_months <= 60, 1L, 0L)]

sink_log(sprintf("Locked cohort: n=%d patients (OS5 events=%d, PFS5 events=%d)",
                 nrow(patient_dt), sum(patient_dt$OS5_event),
                 sum(patient_dt$PFS5_event)))

fwrite(flow, file.path(step_dir, "cohort_flowchart.csv"))

# --- UCell scoring on log2(TPM+1) (only on cohort-locked aliquots) -----------

# Subset expression matrix to locked aliquots (sample columns) — but we need
# original aliquot_id_csv (with .) to subset
locked_aliquot_ids <- patient_dt$aliquot_id_csv
expr_locked <- expr_mat[, locked_aliquot_ids, drop = FALSE]

# log2(x + 1) — UCell is rank-based so only monotone matters, but log helps
# the auxiliary GSVA / ssGSEA / singscore sensitivity scorers downstream.
expr_log <- log2(expr_locked + 1)

sink_log("Computing UCell SecA/SecB scores")
uc <- ScoreSignatures_UCell(matrix = expr_log,
                            features = sig_list,
                            ncores = 1,
                            chunk.size = 500)
uc <- as.data.frame(uc)
uc$aliquot_id_csv <- rownames(uc)
setnames(uc, c("SecA_UCell", "SecB_UCell", "aliquot_id_csv"))

scores_dt <- as.data.table(uc)

# UCell-derived composite scores
eps <- 1e-6
scores_dt[, polar_UCell      := SecB_UCell - SecA_UCell]
scores_dt[, log2_ratio_UCell := log2((SecA_UCell + eps) / (SecB_UCell + eps))]

# Merge CIBERSORTx-derived scores
csx_keep <- csx[, .(aliquot_id_csv, Epi_SecA, Epi_SecB,
                    csx_pvalue = `P-value`)]
csx_keep[, prop_SecB_of_sec := Epi_SecB / pmax(Epi_SecA + Epi_SecB, 1e-9)]
csx_keep[, log2_ratio_csx   := log2((Epi_SecA + eps) / (Epi_SecB + eps))]
scores_dt <- merge(scores_dt, csx_keep, by = "aliquot_id_csv", all.x = TRUE)

# Z-standardize all 4 primary scores (within locked TCGA cohort)
zsc <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
scores_dt[, polar_UCell_z       := zsc(polar_UCell)]
scores_dt[, log2_ratio_UCell_z  := zsc(log2_ratio_UCell)]
scores_dt[, prop_SecB_of_sec_z  := zsc(prop_SecB_of_sec)]
scores_dt[, log2_ratio_csx_z    := zsc(log2_ratio_csx)]

patient_dt <- merge(patient_dt, scores_dt, by = "aliquot_id_csv", all.x = TRUE)

fwrite(scores_dt, file.path(step_dir, "scores_per_sample.csv"))
saveRDS(patient_dt, file.path(step_dir, "tcga_cohort_locked.rds"))
sink_log(sprintf("Saved scores_per_sample.csv (n=%d)", nrow(scores_dt)))

# --- Step 3: Primary Cox tests -----------------------------------------------

sink_log("Fitting primary Cox tests (8 = 4 scores x 2 endpoints)")

primary_scores <- c("polar_UCell_z", "log2_ratio_UCell_z",
                    "prop_SecB_of_sec_z", "log2_ratio_csx_z")
expected_dir   <- c(polar_UCell_z = "+",   # SecB-pole = worse = HR>1
                    log2_ratio_UCell_z = "-",   # SecA num: HR<1 = better SecA
                    prop_SecB_of_sec_z = "+",
                    log2_ratio_csx_z = "-")
endpoints      <- c("OS5", "PFS5")

# CIBERSORTx primary tests use the cohort additionally filtered by csx_pvalue<0.05
csx_filter <- patient_dt$csx_pvalue < 0.05 & !is.na(patient_dt$csx_pvalue)
sink_log(sprintf("CIBERSORTx P<0.05 sub-cohort: n=%d / %d",
                 sum(csx_filter), nrow(patient_dt)))

# Univariate covariate Cox to drive data-driven selection (NEVER conditioning on score)
# Excluded:
#   - platinum_resist: defined post-hoc by progression timing -> conditioning on outcome
#   - residual_grp: TCGA-CDR has residual_tumor 100% NA for OV (would need GDC XML or
#     Bell 2011 supplement to recover; documented limitation)
covariate_candidates <- c("Age", "stage_grp", "grade_grp",
                          "absolute_purity", "tp53_status")

uni_cov_results <- list()
for (ep in endpoints) {
  tcol <- paste0(ep, "_months"); ecol <- paste0(ep, "_event")
  for (cv in covariate_candidates) {
    df_cv <- patient_dt[!is.na(get(cv)) & !is.na(get(tcol)) & !is.na(get(ecol))]
    if (nrow(df_cv) < 30) next
    f <- as.formula(sprintf("Surv(%s, %s) ~ %s", tcol, ecol, cv))
    fit <- tryCatch(coxph(f, data = df_cv), error = function(e) NULL)
    if (is.null(fit)) next
    s <- summary(fit)
    p <- s$logtest["pvalue"]
    uni_cov_results[[length(uni_cov_results) + 1]] <- data.table(
      endpoint = ep, covariate = cv, n = nrow(df_cv),
      lr_p = unname(p))
  }
}
uni_cov_dt <- rbindlist(uni_cov_results)
fwrite(uni_cov_dt, file.path(step_dir, "univariate_covariate_screen.csv"))

# Data-driven covariate set (per endpoint): force {Age, stage_grp}, then add any
# candidate with univariate LR p < 0.10 in this cohort.
build_data_driven <- function(ep) {
  forced <- c("Age", "stage_grp")
  cand <- uni_cov_dt[endpoint == ep & lr_p < 0.10 & !covariate %in% forced,
                     covariate]
  union(forced, cand)
}

# Cox fit helper
fit_cox <- function(df, score, time_col, event_col, covars = character()) {
  rhs <- if (length(covars)) paste(c(score, covars), collapse = " + ")
         else score
  f <- as.formula(sprintf("Surv(%s, %s) ~ %s", time_col, event_col, rhs))
  df_fit <- df[!is.na(get(score)) & !is.na(get(time_col)) & !is.na(get(event_col))]
  for (cv in covars) df_fit <- df_fit[!is.na(get(cv))]
  if (nrow(df_fit) < 30 || sum(df_fit[[event_col]]) < 10) return(NULL)
  fit <- tryCatch(coxph(f, data = df_fit), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  cf <- s$coefficients[score, , drop = FALSE]
  ci <- s$conf.int[score, , drop = FALSE]
  data.table(score = score,
             n = fit$n, n_event = sum(df_fit[[event_col]]),
             HR = ci[1], HR_lo = ci[3], HR_hi = ci[4],
             p = cf[, "Pr(>|z|)"],
             logtest_p = unname(s$logtest["pvalue"]),
             covariates = paste(covars, collapse = "+"))
}

primary_results <- list()
for (ep in endpoints) {
  tcol <- paste0(ep, "_months"); ecol <- paste0(ep, "_event")
  dd_covs <- build_data_driven(ep)
  ks_covs <- c("Age", "stage_grp", "grade_grp",
               "absolute_purity", "tp53_status")
  for (sc in primary_scores) {
    df_sc <- patient_dt
    if (sc %in% c("prop_SecB_of_sec_z", "log2_ratio_csx_z")) df_sc <- df_sc[csx_filter]
    for (model_lbl in c("crude", "age_stage", "data_driven", "kitchen_sink")) {
      covars <- switch(model_lbl,
                       crude        = character(),
                       age_stage    = c("Age", "stage_grp"),
                       data_driven  = dd_covs,
                       kitchen_sink = ks_covs)
      r <- fit_cox(df_sc, sc, tcol, ecol, covars)
      if (is.null(r)) next
      r[, endpoint := ep]
      r[, model    := model_lbl]
      r[, expected_direction := expected_dir[[sc]]]
      r[, direction_match := fifelse(
           expected_direction == "+", HR > 1,
           HR < 1)]
      primary_results[[length(primary_results) + 1]] <- r
    }
  }
}
primary_dt <- rbindlist(primary_results, fill = TRUE)
setcolorder(primary_dt, c("endpoint", "score", "model", "n", "n_event",
                          "HR", "HR_lo", "HR_hi", "p", "logtest_p",
                          "expected_direction", "direction_match",
                          "covariates"))

# Multiplicity: applied within the family of 8 PRIMARY tests (model = age_stage)
fam <- primary_dt[model == "age_stage"]
fam[, p_bh         := p.adjust(p, method = "BH")]
fam[, p_bonferroni := p.adjust(p, method = "bonferroni")]
fam[, ci_excludes_1 := (HR_lo > 1) | (HR_hi < 1)]
fam[, replicated   := direction_match & p_bonferroni < 0.05 & ci_excludes_1]

# Merge multiplicity columns back to full primary table
primary_dt <- merge(primary_dt,
                    fam[, .(endpoint, score, p_bh, p_bonferroni,
                            ci_excludes_1, replicated)],
                    by = c("endpoint", "score"), all.x = TRUE)

fwrite(primary_dt, file.path(step_dir, "cox_primary_results.csv"))
sink_log(sprintf("Wrote cox_primary_results.csv (%d rows)", nrow(primary_dt)))

# Print summary to console
cat("\n=== PRIMARY Cox: model = age_stage (the pre-specified inferential model) ===\n")
print(fam[, .(endpoint, score, n, n_event,
              HR = round(HR, 3), CI = sprintf("[%.2f, %.2f]", HR_lo, HR_hi),
              p = signif(p, 3),
              p_BH = signif(p_bh, 3), p_Bonf = signif(p_bonferroni, 3),
              replicated)])

# Replication call
n_rep <- sum(fam$replicated, na.rm = TRUE)
n_dir_unadj <- sum(fam$direction_match & fam$p < 0.05, na.rm = TRUE)
call_str <- if (n_rep >= 1) {
  "REPLICATED"
} else if (n_dir_unadj >= 1) {
  "PARTIALLY REPLICATED"
} else {
  "NOT REPLICATED"
}

sink_log(sprintf("Replication call: %s (rep=%d, dir+unadj.p<0.05=%d, n_tests=%d)",
                 call_str, n_rep, n_dir_unadj, nrow(fam)))
cat("\n=== Replication call:", call_str, "===\n")

saveRDS(list(primary_dt = primary_dt, fam = fam, call = call_str,
             cohort = patient_dt),
        file.path(step_dir, "phase40_primary.rds"))

sink_log("Phase 40 primary tests complete")
