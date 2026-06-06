# ============================================================================
# 02_clinical_comprehensive.R
# ----------------------------------------------------------------------------
# PURPOSE: Comprehensive clinical analysis: composite per-patient features (log2 SecA:SecB ratio, proportions) and extended survival reporting.
#
# INPUTS:
#   - output/10_clinical_v2/per_patient_features_v2.csv + clinical_data_clean.csv
#
# OUTPUTS:
#   - output/10_clinical_v2/ comprehensive survival tables + report
#
# MANUSCRIPT PANEL(S): Fig 7A, Fig 7B (supporting).
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

# Survival analysis packages
if (!requireNamespace("survival", quietly = TRUE))
  install.packages("survival", repos = "https://cloud.r-project.org")
if (!requireNamespace("survminer", quietly = TRUE))
  install.packages("survminer", repos = "https://cloud.r-project.org")

library(survival)
library(survminer)
library(broom)
library(scales)
library(ggrepel)
library(circlize)

# --- Output directory -------------------------------------------------------

comp_dir <- file.path(out_dir, "10_clinical_comprehensive")
fig_comp <- file.path(comp_dir, "figures")
for (d in c(comp_dir, fig_comp)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# --- Constants / palettes ---------------------------------------------------

epi_types     <- c("Ciliated epithelium", "SecA epithelium",
                   "Intermediate epithelium", "SecB epithelium")
immune_types  <- c("T cell", "NK cell", "B cell", "Plasma cell",
                   "Macrophage", "Conventional dendritic cell",
                   "Plasmacytoid dendritic cell", "Neutrophil", "Mast cell")
stromal_types <- c("Fibroblast", "Smooth muscle", "Mesothelial",
                   "Pericyte", "Endothelial")
celltype_order <- c(epi_types, stromal_types, immune_types)

km_pal    <- c("Low" = "#0072B2", "High" = "#D55E00")
chemo_pal <- c("chemosensitive" = "#56B4E9", "chemoresistant" = "#D55E00")
debulk_pal <- c("optimal" = "#56B4E9", "suboptimal" = "#D55E00")
sec_pal   <- c("SecA" = "#E6A141", "Intermediate" = "#C08E48",
               "SecB" = "#9A7D55")
strat_pal <- c("High ratio / Sensitive"  = "#0072B2",
               "High ratio / Resistant"  = "#56B4E9",
               "Low ratio / Sensitive"   = "#D55E00",
               "Low ratio / Resistant"   = "#E69F00")
feat_type_pal <- c("Cell density" = "#0072B2", "Neighborhood" = "#009E73",
                   "Composite" = "#D55E00", "Pathway" = "#CC79A7",
                   "UCell" = "#E69F00")

# Residual disease ordered levels
residual_order  <- c("R0", "optimal_microscopic", "optimal_le1cm",
                     "suboptimal_1to2cm", "suboptimal_gt2cm",
                     "suboptimal_miliary")
residual_labels <- c("R0", "Microscopic", "<=1 cm", "1-2 cm",
                     ">2 cm", "Miliary")


# ============================================================================
# PART 1: LOAD DATA & FEATURE ENGINEERING
# ============================================================================

message("\n=== Part 1: Feature Engineering ===")

per_patient <- fread(file.path(out_dir, "10_clinical_v2",
                               "per_patient_features_v2.csv"))

# --- Derived features -------------------------------------------------------

eps_ratio <- 0.1

# Ratio (may already exist; recompute for safety)
per_patient[, ratio_SecA_SecB := (dens_SecA_epithelium + eps_ratio) /
                                  (dens_SecB_epithelium + eps_ratio)]
per_patient[, log2_ratio := log2(ratio_SecA_SecB)]

# Total secretory density
per_patient[, total_secretory := dens_SecA_epithelium + dens_SecB_epithelium +
                                  dens_Intermediate_epithelium]

# Secretory subtype proportions (use precomputed if present, else compute)
if (!"prop_SecA" %in% names(per_patient)) {
  per_patient[, prop_SecA  := dens_SecA_epithelium / (total_secretory + eps_ratio)]
  per_patient[, prop_SecB  := dens_SecB_epithelium / (total_secretory + eps_ratio)]
  per_patient[, prop_Int := dens_Intermediate_epithelium /
                                (total_secretory + eps_ratio)]
} else {
  # Use existing columns (prop_SecA_of_sec, etc.)
  if ("prop_SecA_of_sec" %in% names(per_patient)) {
    per_patient[, prop_SecA  := prop_SecA_of_sec]
    per_patient[, prop_SecB  := prop_SecB_of_sec]
    per_patient[, prop_Int := prop_Int_of_sec]
  }
}

# Clean empty strings to NA for key clinical columns
clin_cols <- c("residual_binary", "chemo_status_6months", "chemo_status_1year",
               "residual_category", "treatment_status")
for (v in intersect(clin_cols, names(per_patient))) {
  per_patient[get(v) == "", (v) := NA_character_]
}

# Stage binary
per_patient[, stage_binary := fifelse(stage <= 2, "I-II", "III-IV")]

# Debulking binary outcome for logistic
per_patient[!is.na(residual_binary),
            debulk_optimal := as.integer(residual_binary == "optimal")]

# Factor residual_category
per_patient[residual_category %in% residual_order,
            resid_cat_ordered := factor(residual_category,
                                        levels = residual_order,
                                        labels = residual_labels)]

# --- Feature lists ----------------------------------------------------------

dens_features   <- grep("^dens_", names(per_patient), value = TRUE)
nb_features     <- grep("^prop_nb_", names(per_patient), value = TRUE)
pathway_features <- grep("^pathway_", names(per_patient), value = TRUE)
ucell_features  <- grep("^polar_", names(per_patient), value = TRUE)

composite_features <- c("log2_ratio", "prop_SecA", "prop_SecB", "prop_Int",
                        "total_secretory")
# Only keep composite features that actually exist
composite_features <- intersect(composite_features, names(per_patient))

lineage_features <- intersect(c("prop_epi", "prop_immune", "prop_stromal"),
                              names(per_patient))

all_features <- c(dens_features, nb_features, composite_features)

message("Cohort: ", nrow(per_patient), " patients")
message("Density features: ", length(dens_features))
message("Neighborhood features: ", length(nb_features))
message("Pathway features: ", length(pathway_features))
message("UCell polarization features: ", length(ucell_features))
message("Clinical columns: chemo_response_6mo N=",
        sum(!is.na(per_patient$chemo_response_6mo)),
        ", chemo_response_12mo N=",
        sum(!is.na(per_patient$chemo_response_12mo)),
        ", residual_binary N=",
        sum(!is.na(per_patient$residual_binary)))


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# --- Feature label helper ---------------------------------------------------
make_feature_label <- function(feat) {
  lab <- gsub("^dens_|^prop_", "", feat)
  lab <- gsub("_", " ", lab)
  if (grepl("^prop_nb_", feat)) {
    nb_key <- gsub("^prop_", "", feat)
    if (nb_key %in% names(nb_names)) lab <- nb_names[nb_key]
  }
  if (feat == "log2_ratio") lab <- "log2(SecA:SecB)"
  if (feat == "prop_SecA")  lab <- "SecA proportion"
  if (feat == "prop_SecB")  lab <- "SecB proportion"
  if (feat == "prop_Int") lab <- "Trans proportion"
  if (feat == "total_secretory") lab <- "Total secretory density"
  lab
}

# --- Feature type classifier ------------------------------------------------
get_feature_type <- function(feat) {
  ifelse(grepl("^dens_", feat), "Cell density",
  ifelse(grepl("^prop_nb_", feat), "Neighborhood",
  ifelse(grepl("^pathway_", feat), "Pathway",
  ifelse(grepl("^polar_", feat), "UCell",
         "Composite"))))
}

# --- Univariate logistic regression -----------------------------------------
run_logistic <- function(dat, features, outcome_var, endpoint_label = "") {
  results <- lapply(features, function(feat) {
    df <- dat[!is.na(get(feat)) & !is.na(get(outcome_var))]
    n_total <- nrow(df)
    if (n_total < 10) {
      return(data.table(endpoint = endpoint_label, feature = feat,
                        n = n_total, n_pos = NA_integer_, n_neg = NA_integer_,
                        OR = NA_real_, OR_lower = NA_real_, OR_upper = NA_real_,
                        p_value = NA_real_))
    }
    n_pos <- sum(df[[outcome_var]], na.rm = TRUE)
    n_neg <- n_total - n_pos
    n_min <- min(n_pos, n_neg)
    if (n_min < 5) {
      return(data.table(endpoint = endpoint_label, feature = feat,
                        n = n_total, n_pos = n_pos, n_neg = n_neg,
                        OR = NA_real_, OR_lower = NA_real_, OR_upper = NA_real_,
                        p_value = NA_real_))
    }
    df[, x_scaled := as.numeric(scale(get(feat)))]
    fit <- tryCatch(
      glm(get(outcome_var) ~ x_scaled, data = df, family = binomial),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(data.table(endpoint = endpoint_label, feature = feat,
                        n = n_total, n_pos = n_pos, n_neg = n_neg,
                        OR = NA_real_, OR_lower = NA_real_, OR_upper = NA_real_,
                        p_value = NA_real_))
    }
    ci <- tryCatch(confint(fit), error = function(e) matrix(NA, 2, 2))
    coef_row <- coef(summary(fit))[2, ]
    data.table(
      endpoint    = endpoint_label,
      feature     = feat,
      n           = n_total,
      n_pos       = n_pos,
      n_neg       = n_neg,
      OR          = exp(coef(fit)[2]),
      OR_lower    = exp(ci[2, 1]),
      OR_upper    = exp(ci[2, 2]),
      p_value     = coef_row[4]
    )
  })
  rbindlist(results)
}

# --- Univariate Cox PH ------------------------------------------------------
run_cox_univariate <- function(dat, features, time_var, event_var,
                               endpoint_label = "") {
  results <- lapply(features, function(feat) {
    df <- dat[!is.na(get(feat)) & !is.na(get(time_var)) & !is.na(get(event_var))]
    n_total <- nrow(df)
    if (n_total < 15) {
      return(data.table(endpoint = endpoint_label, feature = feat,
                        n = n_total, n_events = NA_integer_,
                        HR = NA_real_, HR_lower = NA_real_, HR_upper = NA_real_,
                        p_value = NA_real_, concordance = NA_real_))
    }
    df[, x_scaled := as.numeric(scale(get(feat)))]
    fit <- tryCatch(
      coxph(Surv(get(time_var), get(event_var)) ~ x_scaled, data = df),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(data.table(endpoint = endpoint_label, feature = feat,
                        n = n_total, n_events = NA_integer_,
                        HR = NA_real_, HR_lower = NA_real_, HR_upper = NA_real_,
                        p_value = NA_real_, concordance = NA_real_))
    }
    s <- summary(fit)
    data.table(
      endpoint    = endpoint_label,
      feature     = feat,
      n           = s$n,
      n_events    = s$nevent,
      HR          = s$conf.int[1, 1],
      HR_lower    = s$conf.int[1, 3],
      HR_upper    = s$conf.int[1, 4],
      p_value     = s$coefficients[1, 5],
      concordance = s$concordance[1]
    )
  })
  rbindlist(results)
}

# --- Within-stratum Wilcoxon tests ------------------------------------------
run_within_stratum_wilcox <- function(data, features, stratum_col = "ratio_group",
                                      outcome_col = "chemo_group") {
  strata <- unique(data[[stratum_col]])
  results <- rbindlist(lapply(strata, function(s) {
    sub <- data[get(stratum_col) == s]
    rbindlist(lapply(features, function(feat) {
      vals_s <- sub[get(outcome_col) == "Sensitive"][[feat]]
      vals_r <- sub[get(outcome_col) == "Resistant"][[feat]]
      vals_s <- vals_s[!is.na(vals_s)]
      vals_r <- vals_r[!is.na(vals_r)]
      if (length(vals_s) < 3 || length(vals_r) < 3) {
        return(data.table(stratum = s, feature = feat,
                          n_sensitive = length(vals_s),
                          n_resistant = length(vals_r),
                          median_sensitive = NA_real_,
                          median_resistant = NA_real_,
                          hodges_lehmann = NA_real_,
                          p_value = NA_real_))
      }
      wt <- tryCatch(wilcox.test(vals_s, vals_r, conf.int = TRUE),
                     error = function(e) NULL)
      data.table(
        stratum          = s,
        feature          = feat,
        n_sensitive      = length(vals_s),
        n_resistant      = length(vals_r),
        median_sensitive = round(median(vals_s), 4),
        median_resistant = round(median(vals_r), 4),
        hodges_lehmann   = if (!is.null(wt)) round(wt$estimate, 4) else NA_real_,
        p_value          = if (!is.null(wt)) wt$p.value else NA_real_
      )
    }))
  }))
  results[, q_value := p.adjust(p_value, method = "BH"), by = stratum]
  results
}

# --- GLM summary extractor --------------------------------------------------
extract_glm <- function(fit, model_name) {
  if (is.null(fit)) {
    return(data.table(Model = model_name, Term = "fit_failed",
                      OR = NA, CI_lower = NA, CI_upper = NA,
                      p_value = NA, AIC = NA))
  }
  tidy_fit <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  tidy_fit <- tidy_fit[tidy_fit$term != "(Intercept)", ]
  data.table(
    Model    = model_name,
    Term     = tidy_fit$term,
    OR       = round(tidy_fit$estimate, 3),
    CI_lower = round(tidy_fit$conf.low, 3),
    CI_upper = round(tidy_fit$conf.high, 3),
    p_value  = signif(tidy_fit$p.value, 3),
    AIC      = round(AIC(fit), 1)
  )
}

# --- Save figure helper (ragg) ----------------------------------------------
save_fig <- function(filename, plot_obj = last_plot(), width = 10, height = 7,
                     dpi = 150) {
  fp <- file.path(fig_comp, filename)
  ragg::agg_png(fp, width = width, height = height, units = "in", res = dpi)
  print(plot_obj)
  invisible(dev.off())
  message("  Saved: ", fp)
  fp
}


# ============================================================================
# PART 2: SURVIVAL ANALYSIS
# ============================================================================

message("\n=== Part 2: Survival Analysis ===")

# --- 2a. Univariate Cox PH for all features --------------------------------

cox_os <- run_cox_univariate(per_patient, all_features,
                             "os_time_5y", "os_event_5y", "OS")
cox_pfs <- run_cox_univariate(per_patient, all_features,
                              "pfs_time_5y", "pfs_event_5y", "PFS")
cox_all <- rbind(cox_os, cox_pfs)
cox_all[, q_value := p.adjust(p_value, method = "BH"), by = endpoint]
cox_all[, feature_label := sapply(feature, make_feature_label)]
cox_all[, feature_type := sapply(feature, get_feature_type)]

fwrite(cox_all, file.path(comp_dir, "cox_univariate_comprehensive.csv"))
message("  Univariate Cox: ", nrow(cox_all[!is.na(HR)]), " valid models")

# --- 2b. Forest plot of HRs ------------------------------------------------

p_forest_cox <- {
  fdt <- cox_all[!is.na(HR) & is.finite(HR_lower) & is.finite(HR_upper)]
  fdt[, feature_label := reorder(feature_label, HR)]
  ggplot(fdt, aes(x = HR, y = feature_label, color = feature_type)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
    geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper),
                   height = 0.3, linewidth = 0.4) +
    geom_point(size = 2) +
    scale_x_log10() +
    scale_color_manual(values = feat_type_pal, name = "Feature type") +
    facet_wrap(~ endpoint) +
    labs(x = "Hazard ratio (per SD, log scale)", y = NULL,
         title = "Univariate Cox PH: all features vs OS and PFS",
         subtitle = "HR > 1 = worse prognosis; HR < 1 = better prognosis") +
    theme_lab(base_size = 9) +
    theme(strip.text = element_text(size = rel(1.1)))
}
save_fig("part2_forest_cox_univariate.png", p_forest_cox, width = 14, height = 10)

# --- 2c. KM curves for SecB density and log2 ratio -------------------------

km_features <- list(
  list(feat = "dens_SecB_epithelium", label = "SecB density"),
  list(feat = "log2_ratio", label = "log2(SecA:SecB) ratio")
)

km_plots <- list()
for (km_item in km_features) {
  feat  <- km_item$feat
  label <- km_item$label
  for (ep in c("os", "pfs")) {
    time_var  <- paste0(ep, "_time_5y")
    event_var <- paste0(ep, "_event_5y")
    ep_label  <- toupper(ep)

    df <- per_patient[!is.na(get(feat)) & !is.na(get(time_var)) &
                      !is.na(get(event_var))]
    med_val <- median(df[[feat]], na.rm = TRUE)
    df[, km_group := fifelse(get(feat) > med_val, "High", "Low")]

    fit_km <- survfit(Surv(get(time_var), get(event_var)) ~ km_group,
                      data = df)

    fname <- sprintf("part2_km_%s_%s.png",
                     gsub("[^a-z0-9]", "_", tolower(label)), ep)
    ragg::agg_png(file.path(fig_comp, fname), width = 10, height = 7,
                  units = "in", res = 150)
    print(ggsurvplot(fit_km, data = df,
      pval = TRUE, risk.table = TRUE,
      palette = unname(km_pal),
      xlab = "Time (months)", ylab = "Survival probability",
      title = sprintf("5-yr %s by %s (median split)", ep_label, label),
      legend.labs = c("High", "Low"),
      ggtheme = theme_lab(base_size = 9),
      risk.table.height = 0.28, fontsize = 2.5))
    invisible(dev.off())
    message("  Saved: ", fname)
  }
}

# --- 2d. Multivariate Cox: SecA + SecB + immune proportion ------------------

for (ep in c("os", "pfs")) {
  time_var  <- paste0(ep, "_time_5y")
  event_var <- paste0(ep, "_event_5y")

  df_mv <- per_patient[!is.na(dens_SecA_epithelium) &
                        !is.na(dens_SecB_epithelium) &
                        !is.na(get(time_var)) & !is.na(get(event_var))]

  # If prop_immune exists, include it
  if ("prop_immune" %in% names(df_mv) &&
      sum(!is.na(df_mv$prop_immune)) > 15) {
    df_mv <- df_mv[!is.na(prop_immune)]
    df_mv[, z_secA := as.numeric(scale(dens_SecA_epithelium))]
    df_mv[, z_secB := as.numeric(scale(dens_SecB_epithelium))]
    df_mv[, z_imm  := as.numeric(scale(prop_immune))]

    fit_mv <- tryCatch(
      coxph(Surv(get(time_var), get(event_var)) ~ z_secA + z_secB + z_imm,
            data = df_mv),
      error = function(e) NULL
    )
  } else {
    df_mv[, z_secA := as.numeric(scale(dens_SecA_epithelium))]
    df_mv[, z_secB := as.numeric(scale(dens_SecB_epithelium))]
    fit_mv <- tryCatch(
      coxph(Surv(get(time_var), get(event_var)) ~ z_secA + z_secB,
            data = df_mv),
      error = function(e) NULL
    )
  }

  if (!is.null(fit_mv)) {
    message("  Multivariate Cox (", toupper(ep), "): n=", fit_mv$n,
            ", events=", fit_mv$nevent)
    s_mv <- summary(fit_mv)
    print(round(s_mv$conf.int, 3))
  }
}


# ============================================================================
# PART 3: CHEMO SENSITIVITY (replicating 10b)
# ============================================================================

message("\n=== Part 3: Chemo Sensitivity ===")

# --- 3a. Univariate logistic: all features vs chemo response ----------------

logistic_6mo  <- run_logistic(per_patient, all_features,
                              "chemo_response_6mo", "Chemo_6mo")
logistic_12mo <- run_logistic(per_patient, all_features,
                              "chemo_response_12mo", "Chemo_12mo")
logistic_all  <- rbind(logistic_6mo, logistic_12mo)
logistic_all[, q_value := p.adjust(p_value, method = "BH"), by = endpoint]
logistic_all[, feature_label := sapply(feature, make_feature_label)]
logistic_all[, feature_type := sapply(feature, get_feature_type)]

fwrite(logistic_all, file.path(comp_dir, "chemo_logistic_univariate.csv"))
message("  Logistic univariate: ",
        nrow(logistic_all[!is.na(OR)]), " valid models")

# --- 3b. Forest plot — logistic ---------------------------------------------

p_forest_logistic <- {
  fdt <- logistic_all[!is.na(OR)]
  fdt[, feature_label := reorder(feature_label, OR)]
  ggplot(fdt, aes(x = OR, y = feature_label, color = feature_type)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
    geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper),
                   height = 0.3, linewidth = 0.4) +
    geom_point(size = 2) +
    scale_x_log10() +
    scale_color_manual(values = feat_type_pal, name = "Feature type") +
    facet_wrap(~ endpoint, scales = "free_y") +
    labs(x = "Odds ratio (per SD increase, log scale)", y = NULL,
         color = "Feature type",
         title = "Univariate logistic regression: chemo sensitivity",
         subtitle = paste0("OR > 1 = higher feature -> more likely sensitive; ",
                           "OR < 1 = more likely resistant")) +
    theme_lab(base_size = 9) +
    theme(strip.text = element_text(size = rel(1.1)))
}
save_fig("part3_forest_logistic_chemo.png", p_forest_logistic,
         width = 14, height = 10)

# --- 3c. Opposing effects: SecA vs SecB direction check ---------------------

dir_feats <- c("dens_SecA_epithelium", "dens_SecB_epithelium",
               "prop_SecA", "prop_SecB", "log2_ratio")
dir_labs  <- c("SecA density", "SecB density",
               "SecA proportion", "SecB proportion", "log2(SecA:SecB)")

dir_results <- rbindlist(lapply(seq_along(dir_feats), function(i) {
  feat <- dir_feats[i]; lab <- dir_labs[i]
  rbindlist(lapply(c("chemo_response_6mo", "chemo_response_12mo"), function(cv) {
    tp <- ifelse(grepl("6", cv), "6mo", "12mo")
    df <- per_patient[!is.na(get(cv)) & !is.na(get(feat))]
    if (nrow(df) < 10) return(NULL)
    df[, x_scaled := as.numeric(scale(get(feat)))]
    fit <- tryCatch(glm(get(cv) ~ x_scaled, data = df, family = binomial),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    ci <- tryCatch(confint(fit), error = function(e) matrix(NA, 2, 2))
    data.table(Feature = lab, Endpoint = tp,
               OR = exp(coef(fit)[2]),
               OR_lower = exp(ci[2, 1]),
               OR_upper = exp(ci[2, 2]),
               p_value = coef(summary(fit))[2, 4],
               direction = ifelse(exp(coef(fit)[2]) > 1,
                                  "-> sensitive", "-> resistant"))
  }))
}))

p_opposing <- {
  dir_results[, Feature := factor(Feature, levels = rev(dir_labs))]
  ggplot(dir_results, aes(x = OR, y = Feature, color = direction)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
    geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper),
                   height = 0.3, linewidth = 0.5) +
    geom_point(size = 3) +
    scale_x_log10() +
    scale_color_manual(values = c("-> sensitive" = "#56B4E9",
                                  "-> resistant" = "#D55E00"),
                       name = "Direction") +
    facet_wrap(~ Endpoint) +
    labs(x = "Odds ratio (per SD, log scale)", y = NULL,
         title = "Do SecA and SecB pull in opposite directions?",
         subtitle = paste0("If SecA -> sensitive and SecB -> resistant, ",
                           "the ratio would cancel out")) +
    theme_lab(base_size = 10)
}
save_fig("part3_opposing_effects.png", p_opposing, width = 12, height = 5)

# --- 3d. Core test: SecB at matched SecA levels (12-month) ------------------

dat_12 <- per_patient[!is.na(chemo_response_12mo)]
dat_12[, y := chemo_response_12mo]
dat_12[, z_secA  := as.numeric(scale(dens_SecA_epithelium))]
dat_12[, z_secB  := as.numeric(scale(dens_SecB_epithelium))]
dat_12[, z_ratio := as.numeric(scale(log2_ratio))]
dat_12[, z_propA := as.numeric(scale(prop_SecA))]
dat_12[, z_propB := as.numeric(scale(prop_SecB))]

core_models <- list(
  mA = tryCatch(glm(y ~ z_secB, data = dat_12, family = binomial),
                error = function(e) NULL),
  mB = tryCatch(glm(y ~ z_secA, data = dat_12, family = binomial),
                error = function(e) NULL),
  mC = tryCatch(glm(y ~ z_secB + z_secA, data = dat_12, family = binomial),
                error = function(e) NULL),
  mD = tryCatch(glm(y ~ z_ratio + z_secA, data = dat_12, family = binomial),
                error = function(e) NULL),
  mE = tryCatch(glm(y ~ z_propB + z_propA, data = dat_12, family = binomial),
                error = function(e) NULL),
  mF = tryCatch(glm(y ~ z_ratio, data = dat_12, family = binomial),
                error = function(e) NULL)
)

core_results <- rbind(
  extract_glm(core_models$mA, "A: SecB alone"),
  extract_glm(core_models$mB, "B: SecA alone"),
  extract_glm(core_models$mC, "C: SecB + SecA (key test)"),
  extract_glm(core_models$mD, "D: ratio + SecA"),
  extract_glm(core_models$mE, "E: prop_SecB + prop_SecA"),
  extract_glm(core_models$mF, "F: ratio alone")
)
fwrite(core_results, file.path(comp_dir, "chemo_secB_confound_models.csv"))
message("  Core confound models saved")

# VIF check
cor_AB <- cor(dat_12$z_secA, dat_12$z_secB, use = "complete.obs")
message("  Pearson r(SecA, SecB) = ", round(cor_AB, 3),
        ", VIF = ", round(1 / (1 - cor_AB^2), 2))

# --- 3e. SecA-residualized SecB analysis ------------------------------------

dat_resid <- per_patient[!is.na(chemo_response_12mo) &
                          !is.na(dens_SecB_epithelium) &
                          !is.na(dens_SecA_epithelium)]
fit_secB_on_secA <- lm(dens_SecB_epithelium ~ dens_SecA_epithelium,
                        data = dat_resid)
dat_resid[, secB_resid := residuals(fit_secB_on_secA)]
dat_resid[, z_secB_resid := as.numeric(scale(secB_resid))]

wt_secB_resid <- tryCatch(
  wilcox.test(secB_resid ~ chemo_status_1year, data = dat_resid),
  error = function(e) NULL
)
fit_resid_logistic <- tryCatch(
  glm(chemo_response_12mo ~ z_secB_resid, data = dat_resid, family = binomial),
  error = function(e) NULL
)

if (!is.null(fit_resid_logistic)) {
  resid_or <- exp(coef(fit_resid_logistic)[2])
  resid_ci <- tryCatch(exp(confint(fit_resid_logistic)[2, ]),
                        error = function(e) c(NA, NA))
  resid_p <- coef(summary(fit_resid_logistic))[2, 4]
  message("  SecB residualized on SecA: OR = ", round(resid_or, 3),
          ", p = ", signif(resid_p, 3))
}

p_secB_resid <- {
  p1 <- ggplot(dat_resid, aes(x = dens_SecA_epithelium,
                                y = dens_SecB_epithelium,
                                color = chemo_status_1year)) +
    geom_point(size = 2.5, alpha = 0.6) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                color = "grey50", linewidth = 0.5, linetype = "dashed") +
    scale_color_manual(values = chemo_pal, name = "12mo status") +
    labs(x = expression("SecA density (cells/mm"^2*")"),
         y = expression("SecB density (cells/mm"^2*")"),
         title = "SecA vs SecB with regression line") +
    theme_lab(base_size = 9)

  p2 <- ggplot(dat_resid, aes(x = chemo_status_1year, y = secB_resid,
                                fill = chemo_status_1year)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, linewidth = 0.3) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
    scale_fill_manual(values = chemo_pal, guide = "none") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    labs(x = "12-month chemo status",
         y = "SecB residual\n(SecA regressed out)",
         title = "Excess SecB (controlling for SecA)") +
    theme_lab(base_size = 9)

  p1 | p2
}
save_fig("part3_secB_residualized.png", p_secB_resid, width = 12, height = 5)

# --- 3f. Interaction: ratio x chemo response on survival --------------------

int_results <- list()
for (ep in c("os", "pfs")) {
  time_var  <- paste0(ep, "_time_5y")
  event_var <- paste0(ep, "_event_5y")
  ep_label  <- toupper(ep)

  df <- per_patient[!is.na(chemo_response_12mo) & !is.na(log2_ratio) &
                    !is.na(get(time_var)) & !is.na(get(event_var))]
  if (nrow(df) < 15) next

  fit_int <- tryCatch(
    coxph(as.formula(paste("Surv(", time_var, ",", event_var,
                           ") ~ scale(log2_ratio) * factor(chemo_response_12mo)")),
          data = df),
    error = function(e) NULL
  )
  if (!is.null(fit_int)) {
    s <- summary(fit_int)
    int_results[[ep_label]] <- data.table(
      Endpoint = ep_label,
      Term = rownames(s$coefficients),
      HR = round(s$conf.int[, 1], 3),
      CI = sprintf("[%.2f-%.2f]", s$conf.int[, 3], s$conf.int[, 4]),
      p_value = signif(s$coefficients[, 5], 3)
    )
  }
}
if (length(int_results) > 0) {
  int_dt <- rbindlist(int_results)
  fwrite(int_dt, file.path(comp_dir, "chemo_interaction_cox.csv"))
  message("  Interaction models saved")
}

# --- 3g. Stratified KM: 2x2 (ratio x chemo) --------------------------------

df_int <- per_patient[!is.na(chemo_response_12mo) & !is.na(log2_ratio)]
med_ratio <- median(df_int$log2_ratio, na.rm = TRUE)
df_int[, ratio_group := fifelse(log2_ratio > med_ratio,
                                "High ratio", "Low ratio")]
df_int[, chemo_group := fifelse(chemo_response_12mo == 1,
                                "Sensitive", "Resistant")]
df_int[, strat_group := factor(
  paste(ratio_group, chemo_group, sep = " / "),
  levels = names(strat_pal)
)]

for (ep in c("os", "pfs")) {
  time_var  <- paste0(ep, "_time_5y")
  event_var <- paste0(ep, "_event_5y")

  df_km <- data.frame(time = df_int[[time_var]],
                       event = df_int[[event_var]],
                       group = df_int$strat_group)
  df_km <- df_km[complete.cases(df_km), ]
  if (nrow(df_km) < 10) next

  fit_km <- survfit(Surv(time, event) ~ group, data = df_km)
  fname <- sprintf("part3_km_2x2_%s.png", ep)
  ragg::agg_png(file.path(fig_comp, fname), width = 12, height = 7,
                units = "in", res = 150)
  print(ggsurvplot(fit_km, data = df_km,
    pval = TRUE, risk.table = TRUE,
    palette = unname(strat_pal),
    xlab = "Time (months)", ylab = "Survival probability",
    title = sprintf("5-yr %s: SecA:SecB ratio x 12-month chemo response",
                    toupper(ep)),
    legend.labs = names(strat_pal),
    ggtheme = theme_lab(base_size = 8),
    risk.table.height = 0.3, fontsize = 2.5))
  invisible(dev.off())
  message("  Saved: ", fname)
}


# ============================================================================
# PART 4: DEBULKING COMPOSITION (replicating 10c)
# ============================================================================

message("\n=== Part 4: Debulking Composition ===")

# --- 4a. Wilcoxon tests per feature ----------------------------------------

debulk_dat <- per_patient[!is.na(residual_binary)]
n_debulk <- nrow(debulk_dat)
message("  Debulking cohort: ", n_debulk, " patients (",
        sum(debulk_dat$residual_binary == "optimal"), " optimal, ",
        sum(debulk_dat$residual_binary == "suboptimal"), " suboptimal)")

if (n_debulk >= 10) {

  wilcox_debulk <- rbindlist(lapply(all_features, function(feat) {
    df <- debulk_dat[!is.na(get(feat))]
    n_opt <- sum(df$residual_binary == "optimal")
    n_sub <- sum(df$residual_binary == "suboptimal")
    wt <- tryCatch(wilcox.test(get(feat) ~ residual_binary, data = df),
                   error = function(e) NULL)
    data.table(
      feature = feat,
      label = make_feature_label(feat),
      feature_type = get_feature_type(feat),
      median_optimal = round(df[residual_binary == "optimal",
                                 median(get(feat), na.rm = TRUE)], 4),
      median_suboptimal = round(df[residual_binary == "suboptimal",
                                    median(get(feat), na.rm = TRUE)], 4),
      n_optimal = n_opt,
      n_suboptimal = n_sub,
      p_value = if (!is.null(wt)) wt$p.value else NA_real_
    )
  }))
  wilcox_debulk[, q_value := p.adjust(p_value, method = "BH")]
  wilcox_debulk <- wilcox_debulk[order(p_value)]

  fwrite(wilcox_debulk, file.path(comp_dir, "debulking_wilcoxon_tests.csv"))
  message("  Features with p < 0.05: ",
          sum(wilcox_debulk$p_value < 0.05, na.rm = TRUE))

  # --- 4b. Boxplots: secretory subtypes by debulking -------------------------

  sec_feats <- c("dens_SecA_epithelium", "dens_SecB_epithelium",
                 "dens_Intermediate_epithelium", "dens_Ciliated_epithelium",
                 "log2_ratio")
  sec_labs  <- c("SecA density", "SecB density", "Transitioning density",
                 "Ciliated density", "log2(SecA:SecB)")

  plots_debulk <- lapply(seq_along(sec_feats), function(i) {
    feat <- sec_feats[i]; lab <- sec_labs[i]
    df <- debulk_dat[!is.na(get(feat)), .(val = get(feat),
                                           group = residual_binary)]
    wt <- tryCatch(wilcox.test(val ~ group, data = df)$p.value,
                   error = function(e) NA)
    ggplot(df, aes(x = group, y = val, fill = group)) +
      geom_boxplot(outlier.shape = NA, linewidth = 0.3, alpha = 0.7) +
      geom_jitter(width = 0.15, size = 1.2, alpha = 0.5) +
      scale_fill_manual(values = debulk_pal, guide = "none") +
      labs(x = NULL, y = lab,
           subtitle = sprintf("Wilcoxon p = %s",
                              ifelse(is.na(wt), "NA", signif(wt, 3)))) +
      theme_lab(base_size = 9)
  })

  p_debulk_box <- (plots_debulk[[1]] | plots_debulk[[2]] |
                    plots_debulk[[3]]) /
                   (plots_debulk[[4]] | plots_debulk[[5]] | plot_spacer()) +
    plot_annotation(
      title = "Secretory subtype composition by debulking outcome"
    )
  save_fig("part4_debulking_boxplots.png", p_debulk_box, width = 12, height = 8)

  # --- 4c. Ternary plot (SecA, SecB, Transitioning) --------------------------

  tern_dt <- debulk_dat[total_secretory > 0,
    .(patient_id, residual_binary, prop_SecA, prop_SecB, prop_Int)]
  tern_dt[, row_sum := prop_SecA + prop_SecB + prop_Int]
  tern_dt[row_sum > 0, `:=`(
    a = prop_SecA / row_sum,
    b = prop_SecB / row_sum,
    cc = prop_Int / row_sum
  )]
  tern_dt[, x := b + 0.5 * cc]
  tern_dt[, y := (sqrt(3) / 2) * cc]

  tri <- data.frame(x = c(0, 1, 0.5, 0), y = c(0, 0, sqrt(3)/2, 0))

  p_ternary <- ggplot(tern_dt[!is.na(x)], aes(x = x, y = y,
                                                 color = residual_binary)) +
    geom_polygon(data = tri, aes(x = x, y = y),
                 fill = NA, color = "grey40", linewidth = 0.5,
                 inherit.aes = FALSE) +
    geom_point(size = 2.5, alpha = 0.6) +
    scale_color_manual(values = debulk_pal, name = "Debulking") +
    annotate("text", x = -0.05, y = -0.03, label = "SecA",
             size = 3.5, fontface = "bold") +
    annotate("text", x = 1.05, y = -0.03, label = "SecB",
             size = 3.5, fontface = "bold") +
    annotate("text", x = 0.5, y = sqrt(3)/2 + 0.04, label = "Trans",
             size = 3.5, fontface = "bold") +
    coord_fixed() +
    labs(title = "Secretory subtype composition (ternary projection)",
         subtitle = "Each point = one patient") +
    theme_lab(base_size = 10) +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          axis.line = element_blank(), axis.ticks = element_blank())
  save_fig("part4_ternary_debulking.png", p_ternary, width = 8, height = 7)

  # --- 4d. Logistic: debulking outcome, stage-adjusted ----------------------

  df_deb_mv <- debulk_dat[!is.na(debulk_optimal) & !is.na(stage)]
  df_deb_mv[, z_ratio := as.numeric(scale(log2_ratio))]
  df_deb_mv[, z_secA  := as.numeric(scale(dens_SecA_epithelium))]
  df_deb_mv[, z_secB  := as.numeric(scale(dens_SecB_epithelium))]

  deb_models <- list(
    m1 = tryCatch(glm(debulk_optimal ~ z_ratio + stage, data = df_deb_mv,
                       family = binomial), error = function(e) NULL),
    m2 = tryCatch(glm(debulk_optimal ~ z_secA + stage, data = df_deb_mv,
                       family = binomial), error = function(e) NULL),
    m3 = tryCatch(glm(debulk_optimal ~ z_secB + stage, data = df_deb_mv,
                       family = binomial), error = function(e) NULL)
  )

  deb_results <- rbind(
    extract_glm(deb_models$m1, "ratio + stage"),
    extract_glm(deb_models$m2, "SecA + stage"),
    extract_glm(deb_models$m3, "SecB + stage")
  )
  fwrite(deb_results, file.path(comp_dir, "debulking_logistic_models.csv"))
  message("  Debulking logistic models saved")

} else {
  message("  Skipping debulking analysis: insufficient data (N = ", n_debulk, ")")
}


# ============================================================================
# PART 5: STRATIFIED TME ANALYSIS (replicating 10d)
# ============================================================================

message("\n=== Part 5: Stratified TME Analysis ===")

# --- 5a. Define strata: median ratio x chemo response ----------------------

dat_strat <- per_patient[!is.na(chemo_response_12mo) & !is.na(log2_ratio)]
med_ratio <- median(dat_strat$log2_ratio, na.rm = TRUE)
dat_strat[, ratio_group := fifelse(log2_ratio > med_ratio,
                                   "High ratio", "Low ratio")]
dat_strat[, chemo_group := fifelse(chemo_response_12mo == 1,
                                   "Sensitive", "Resistant")]
dat_strat[, strat_group := factor(
  paste(ratio_group, chemo_group, sep = " / "),
  levels = names(strat_pal)
)]

message("  Median log2(SecA:SecB) ratio: ", round(med_ratio, 3))
message("  Strata counts:")
print(table(dat_strat$ratio_group, dat_strat$chemo_group))

# --- 5b. Within-stratum Wilcoxon: cell type densities + composites ----------

dens_composite <- c(dens_features, composite_features)
wilcox_dens <- run_within_stratum_wilcox(dat_strat, dens_composite)
wilcox_dens[, feature_label := sapply(feature, make_feature_label)]
wilcox_dens[, feature_type := sapply(feature, get_feature_type)]

# --- 5c. Within-stratum Wilcoxon: neighborhood proportions -----------------

wilcox_nb <- run_within_stratum_wilcox(dat_strat, nb_features)
wilcox_nb[, feature_label := sapply(feature, make_feature_label)]
wilcox_nb[, feature_type := "Neighborhood"]

# --- 5d. Within-stratum Wilcoxon: pathway scores ---------------------------

wilcox_pw <- data.table()
if (length(pathway_features) > 0) {
  wilcox_pw <- run_within_stratum_wilcox(dat_strat, pathway_features)
  wilcox_pw[, feature_label := sapply(feature, make_feature_label)]
  wilcox_pw[, feature_type := "Pathway"]
}

# --- 5e. Within-stratum Wilcoxon: UCell polarization metrics ----------------

wilcox_uc <- data.table()
if (length(ucell_features) > 0) {
  wilcox_uc <- run_within_stratum_wilcox(dat_strat, ucell_features)
  wilcox_uc[, feature_label := sapply(feature, make_feature_label)]
  wilcox_uc[, feature_type := "UCell"]
}

# Combine all
wilcox_strat_all <- rbind(wilcox_dens, wilcox_nb, wilcox_pw, wilcox_uc,
                           fill = TRUE)
fwrite(wilcox_strat_all,
       file.path(comp_dir, "stratified_wilcoxon_all_features.csv"))
message("  Stratified Wilcoxon tests: ",
        nrow(wilcox_strat_all[!is.na(p_value)]), " valid tests")

# --- 5f. Top hits comparison: boxplots --------------------------------------

top_feats_all <- wilcox_strat_all[!is.na(p_value),
  .(min_p = min(p_value, na.rm = TRUE)), by = feature]
top_feats_all <- top_feats_all[order(min_p)][1:min(6, .N)]$feature

if (length(top_feats_all) > 0) {
  top_plots <- lapply(top_feats_all, function(feat) {
    lab <- make_feature_label(feat)
    p_high <- wilcox_strat_all[feature == feat & stratum == "High ratio",
                                p_value]
    p_low  <- wilcox_strat_all[feature == feat & stratum == "Low ratio",
                                p_value]
    if (length(p_high) == 0) p_high <- NA
    if (length(p_low) == 0) p_low <- NA

    sub <- dat_strat[, .(value = get(feat), chemo = chemo_group, ratio_group)]
    ggplot(sub, aes(x = chemo, y = value, fill = chemo)) +
      geom_boxplot(outlier.shape = NA, linewidth = 0.3, alpha = 0.7) +
      geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
      scale_fill_manual(values = c("Sensitive" = "#56B4E9",
                                    "Resistant" = "#D55E00"),
                        guide = "none") +
      facet_wrap(~ ratio_group) +
      labs(x = NULL, y = lab,
           subtitle = sprintf("High p=%.3g | Low p=%.3g",
                              ifelse(is.na(p_high), NA, p_high),
                              ifelse(is.na(p_low), NA, p_low))) +
      theme_lab(base_size = 8)
  })

  p_top_hits <- wrap_plots(top_plots, ncol = 3) +
    plot_annotation(
      title = paste0("Top discriminating features: Sensitive vs Resistant ",
                     "within each ratio stratum"),
      subtitle = "Selected by minimum p-value across either stratum"
    )
  save_fig("part5_top_hits_stratified.png", p_top_hits, width = 14, height = 10)
}

# --- 5g. Effect size concordance across strata ------------------------------

hl_wide <- dcast(wilcox_strat_all[!is.na(hodges_lehmann)],
                  feature + feature_label + feature_type ~ stratum,
                  value.var = "hodges_lehmann")

# Safe column rename (the columns depend on actual stratum names)
strat_cols <- setdiff(names(hl_wide), c("feature", "feature_label",
                                         "feature_type"))
if (length(strat_cols) == 2) {
  setnames(hl_wide, strat_cols, c("HL_1", "HL_2"))

  pvals <- dcast(wilcox_strat_all[!is.na(p_value)],
                  feature ~ stratum, value.var = "p_value")
  setnames(pvals, strat_cols, c("p_1", "p_2"))
  hl_wide <- merge(hl_wide, pvals[, .(feature, p_1, p_2)], by = "feature")
  hl_wide[, min_p := pmin(p_1, p_2, na.rm = TRUE)]
  hl_wide[, neg_log10_p := -log10(min_p)]

  rho <- cor(hl_wide$HL_1, hl_wide$HL_2, use = "complete.obs",
             method = "spearman")
  message("  Effect concordance Spearman rho = ", round(rho, 3))

  p_concordance <- ggplot(hl_wide[!is.na(HL_1) & !is.na(HL_2)],
                           aes(x = HL_1, y = HL_2, color = feature_type)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(aes(size = neg_log10_p), alpha = 0.7) +
    geom_text_repel(aes(label = feature_label), size = 2.5,
                    max.overlaps = 15,
                    segment.color = "grey70", segment.size = 0.3) +
    scale_color_manual(values = feat_type_pal, name = "Feature type") +
    scale_size_continuous(range = c(1.5, 5), name = "-log10(min p)") +
    labs(x = sprintf("H-L shift (%s: Sensitive - Resistant)", strat_cols[1]),
         y = sprintf("H-L shift (%s: Sensitive - Resistant)", strat_cols[2]),
         title = "Effect size concordance across strata",
         subtitle = sprintf("Spearman rho = %.3f", rho)) +
    theme_lab(base_size = 9) +
    coord_fixed()
  save_fig("part5_effect_concordance.png", p_concordance, width = 10, height = 10)
}

# --- 5h. Within-stratum survival KMs ----------------------------------------

dat_high <- dat_strat[ratio_group == "High ratio"]
dat_low  <- dat_strat[ratio_group == "Low ratio"]

for (ep in c("os", "pfs")) {
  time_var  <- paste0(ep, "_time_5y")
  event_var <- paste0(ep, "_event_5y")

  for (grp in list(list(d = dat_high, name = "high_ratio"),
                    list(d = dat_low,  name = "low_ratio"))) {
    df_km <- data.frame(time = grp$d[[time_var]],
                         event = grp$d[[event_var]],
                         group = grp$d$chemo_group)
    df_km <- df_km[complete.cases(df_km), ]
    if (nrow(df_km) < 6) next

    fit_km <- survfit(Surv(time, event) ~ group, data = df_km)
    fname <- sprintf("part5_km_%s_%s.png", grp$name, ep)
    ragg::agg_png(file.path(fig_comp, fname), width = 10, height = 7,
                  units = "in", res = 150)
    print(ggsurvplot(fit_km, data = df_km,
      pval = TRUE, risk.table = TRUE,
      palette = c("#D55E00", "#56B4E9"),
      xlab = "Time (months)", ylab = "Survival probability",
      title = sprintf("5-yr %s within %s stratum",
                      toupper(ep), toupper(gsub("_", " ", grp$name))),
      legend.labs = c("Resistant", "Sensitive"),
      ggtheme = theme_lab(base_size = 9),
      risk.table.height = 0.28, fontsize = 2.5))
    invisible(dev.off())
    message("  Saved: ", fname)
  }
}


# ============================================================================
# PART 6: UCell-SPECIFIC ANALYSES
# ============================================================================

message("\n=== Part 6: UCell Polarization Analyses ===")

ucell_cols <- intersect(c("polar_mean", "polar_median", "polar_sd",
                           "polar_iqr", "polar_n"), names(per_patient))

if (length(ucell_cols) >= 2) {

  # --- 6a. Polarization variability as clinical predictor -------------------

  # Cox PH: polar_sd and polar_iqr vs OS/PFS
  ucell_cox <- list()
  for (uc in ucell_cols) {
    for (ep in c("os", "pfs")) {
      time_var  <- paste0(ep, "_time_5y")
      event_var <- paste0(ep, "_event_5y")
      df <- per_patient[!is.na(get(uc)) & !is.na(get(time_var)) &
                        !is.na(get(event_var))]
      if (nrow(df) < 15) next
      df[, x_scaled := as.numeric(scale(get(uc)))]
      fit <- tryCatch(
        coxph(Surv(get(time_var), get(event_var)) ~ x_scaled, data = df),
        error = function(e) NULL
      )
      if (!is.null(fit)) {
        s <- summary(fit)
        ucell_cox[[paste(uc, ep)]] <- data.table(
          feature = uc, endpoint = toupper(ep),
          n = s$n, n_events = s$nevent,
          HR = round(s$conf.int[1, 1], 3),
          HR_lower = round(s$conf.int[1, 3], 3),
          HR_upper = round(s$conf.int[1, 4], 3),
          p_value = signif(s$coefficients[1, 5], 3),
          concordance = round(s$concordance[1], 3)
        )
      }
    }
  }

  if (length(ucell_cox) > 0) {
    ucell_cox_dt <- rbindlist(ucell_cox)
    fwrite(ucell_cox_dt, file.path(comp_dir, "ucell_survival_cox.csv"))
    message("  UCell survival models: ", nrow(ucell_cox_dt))
  }

  # --- 6b. Mean polarization by clinical group ------------------------------

  ucell_clinical <- list()

  # By chemo response
  for (cvar in c("chemo_status_6months", "chemo_status_1year")) {
    if (!cvar %in% names(per_patient)) next
    for (uc in ucell_cols) {
      df <- per_patient[!is.na(get(cvar)) & !is.na(get(uc))]
      if (nrow(df) < 10) next
      wt <- tryCatch(wilcox.test(get(uc) ~ get(cvar), data = df),
                     error = function(e) NULL)
      ucell_clinical[[paste(uc, cvar)]] <- data.table(
        feature = uc, grouping = cvar,
        group1 = "chemosensitive", group2 = "chemoresistant",
        median_g1 = df[get(cvar) == "chemosensitive",
                       median(get(uc), na.rm = TRUE)],
        median_g2 = df[get(cvar) == "chemoresistant",
                       median(get(uc), na.rm = TRUE)],
        p_value = if (!is.null(wt)) wt$p.value else NA_real_
      )
    }
  }

  # By debulking
  if ("residual_binary" %in% names(per_patient)) {
    for (uc in ucell_cols) {
      df <- per_patient[!is.na(residual_binary) & !is.na(get(uc))]
      if (nrow(df) < 10) next
      wt <- tryCatch(wilcox.test(get(uc) ~ residual_binary, data = df),
                     error = function(e) NULL)
      ucell_clinical[[paste(uc, "debulk")]] <- data.table(
        feature = uc, grouping = "residual_binary",
        group1 = "optimal", group2 = "suboptimal",
        median_g1 = df[residual_binary == "optimal",
                       median(get(uc), na.rm = TRUE)],
        median_g2 = df[residual_binary == "suboptimal",
                       median(get(uc), na.rm = TRUE)],
        p_value = if (!is.null(wt)) wt$p.value else NA_real_
      )
    }
  }

  if (length(ucell_clinical) > 0) {
    ucell_clinical_dt <- rbindlist(ucell_clinical)
    ucell_clinical_dt[, q_value := p.adjust(p_value, method = "BH")]
    fwrite(ucell_clinical_dt,
           file.path(comp_dir, "ucell_clinical_associations.csv"))
    message("  UCell clinical associations: ", nrow(ucell_clinical_dt))
  }

  # --- 6c. Boxplot: polar_sd by chemo response and debulking ----------------

  uc_box_plots <- list()
  for (uc in intersect(c("polar_sd", "polar_iqr", "polar_mean"), ucell_cols)) {
    uc_label <- gsub("_", " ", uc)

    # Chemo 12mo
    df_c <- per_patient[!is.na(chemo_status_1year) & !is.na(get(uc))]
    if (nrow(df_c) >= 10) {
      wt_c <- tryCatch(
        wilcox.test(get(uc) ~ chemo_status_1year, data = df_c)$p.value,
        error = function(e) NA
      )
      uc_box_plots[[paste0(uc, "_chemo")]] <- ggplot(
        df_c, aes(x = chemo_status_1year, y = get(uc),
                   fill = chemo_status_1year)) +
        geom_boxplot(outlier.shape = NA, linewidth = 0.3, alpha = 0.7) +
        geom_jitter(width = 0.15, size = 1.2, alpha = 0.5) +
        scale_fill_manual(values = chemo_pal, guide = "none") +
        labs(x = NULL, y = uc_label,
             title = paste(uc_label, "by 12mo chemo"),
             subtitle = sprintf("p = %s",
                                ifelse(is.na(wt_c), "NA", signif(wt_c, 3)))) +
        theme_lab(base_size = 9)
    }

    # Debulking
    df_d <- per_patient[!is.na(residual_binary) & !is.na(get(uc))]
    if (nrow(df_d) >= 10) {
      wt_d <- tryCatch(
        wilcox.test(get(uc) ~ residual_binary, data = df_d)$p.value,
        error = function(e) NA
      )
      uc_box_plots[[paste0(uc, "_debulk")]] <- ggplot(
        df_d, aes(x = residual_binary, y = get(uc),
                   fill = residual_binary)) +
        geom_boxplot(outlier.shape = NA, linewidth = 0.3, alpha = 0.7) +
        geom_jitter(width = 0.15, size = 1.2, alpha = 0.5) +
        scale_fill_manual(values = debulk_pal, guide = "none") +
        labs(x = NULL, y = uc_label,
             title = paste(uc_label, "by debulking"),
             subtitle = sprintf("p = %s",
                                ifelse(is.na(wt_d), "NA", signif(wt_d, 3)))) +
        theme_lab(base_size = 9)
    }
  }

  if (length(uc_box_plots) >= 2) {
    p_ucell_box <- wrap_plots(uc_box_plots, ncol = 3) +
      plot_annotation(title = "UCell polarization metrics by clinical group")
    save_fig("part6_ucell_clinical_boxplots.png", p_ucell_box,
             width = 14, height = 8)
  }

  # --- 6d. Correlation: polarization metrics vs survival time ---------------

  p_ucell_scatter_list <- list()
  for (uc in intersect(c("polar_sd", "polar_mean"), ucell_cols)) {
    for (ep in c("os", "pfs")) {
      time_var  <- paste0(ep, "_time_5y")
      event_var <- paste0(ep, "_event_5y")
      df <- per_patient[!is.na(get(uc)) & !is.na(get(time_var))]
      if (nrow(df) < 15) next

      rho_val <- cor(df[[uc]], df[[time_var]], use = "complete.obs",
                      method = "spearman")

      p_ucell_scatter_list[[paste(uc, ep)]] <- ggplot(
        df, aes(x = get(uc), y = get(time_var))) +
        geom_point(size = 1.5, alpha = 0.5) +
        geom_smooth(method = "lm", se = TRUE, color = "#0072B2",
                    linewidth = 0.6) +
        labs(x = gsub("_", " ", uc),
             y = paste0(toupper(ep), " time (months)"),
             subtitle = sprintf("rho = %.3f", rho_val)) +
        theme_lab(base_size = 9)
    }
  }

  if (length(p_ucell_scatter_list) >= 2) {
    p_ucell_scatter <- wrap_plots(p_ucell_scatter_list, ncol = 2) +
      plot_annotation(
        title = "UCell polarization vs survival time (Spearman correlation)"
      )
    save_fig("part6_ucell_survival_scatter.png", p_ucell_scatter,
             width = 12, height = 8)
  }

} else {
  message("  Skipping UCell analyses: insufficient polarization columns")
}


# ============================================================================
# PART 7: HTML REPORT
# ============================================================================

message("\n=== Part 7: Generating HTML Report ===")

report_path <- file.path(comp_dir, "10_clinical_comprehensive_report.html")

# --- Collect all figure paths -----------------------------------------------
fig_files <- sort(list.files(fig_comp, pattern = "\\.png$", full.names = TRUE))

# --- Base64 encode images ---------------------------------------------------
encode_img <- function(path) {
  if (!file.exists(path)) return("")
  raw <- readBin(path, "raw", file.info(path)$size)
  b64 <- base64enc::base64encode(raw)
  paste0("data:image/png;base64,", b64)
}

# Check for base64enc
if (!requireNamespace("base64enc", quietly = TRUE))
  install.packages("base64enc", repos = "https://cloud.r-project.org")

# --- Build HTML sections ----------------------------------------------------

html_header <- '<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>10 Clinical Comprehensive Analysis</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
         sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px;
         line-height: 1.6; color: #333; }
  h1 { color: #2c3e50; border-bottom: 2px solid #2c3e50; padding-bottom: 10px; }
  h2 { color: #34495e; border-bottom: 1px solid #bdc3c7; padding-bottom: 6px;
       margin-top: 40px; }
  h3 { color: #4a6785; }
  img { max-width: 100%%; height: auto; border: 1px solid #eee;
        margin: 10px 0; }
  table { border-collapse: collapse; width: 100%%; margin: 15px 0;
          font-size: 0.9em; }
  th { background-color: #2c3e50; color: white; padding: 8px 12px;
       text-align: left; }
  td { border-bottom: 1px solid #ddd; padding: 6px 12px; }
  tr:nth-child(even) { background-color: #f8f9fa; }
  .callout { padding: 12px 16px; margin: 15px 0; border-left: 4px solid;
             border-radius: 0 4px 4px 0; }
  .callout-warn { background-color: #fff3e0; border-color: #FF9800; }
  .callout-info { background-color: #e3f2fd; border-color: #2196F3; }
  .callout-ok   { background-color: #e8f5e9; border-color: #4CAF50; }
  .callout-alert { background-color: #fce4ec; border-color: #E91E63; }
  .section { margin-bottom: 30px; }
  code { background: #f4f4f4; padding: 2px 6px; border-radius: 3px;
         font-size: 0.9em; }
  .timestamp { color: #888; font-size: 0.85em; }
</style>
</head>
<body>
<h1>10 Clinical Comprehensive Analysis</h1>
<p class="timestamp">Generated: %s</p>
<p>Consolidation of Phases 10b (chemo sensitivity), 10c (debulking composition),
   and 10d (ratio-chemo stratified TME profiling), updated with reclassified
   SecA/SecB/Transitioning annotations (v2).</p>
<p>Cohort: <strong>%d patients</strong> from
   <code>per_patient_features_v2.csv</code>.</p>
<hr>
'

html_body <- sprintf(html_header, format(Sys.time(), "%Y-%m-%d %H:%M"),
                      nrow(per_patient))

# --- Part 2 section --------------------------------------------------------
html_body <- paste0(html_body, '
<div class="section">
<h2>Part 2: Survival Analysis</h2>
<p>Univariate Cox proportional hazards models for all density, neighborhood,
   and composite features vs 5-year OS and PFS.</p>
')

# Survival summary table
cox_top <- cox_all[!is.na(HR)][order(p_value)][1:min(15, .N)]
html_body <- paste0(html_body, '<h3>Top 15 survival associations (by p-value)</h3>')
html_body <- paste0(html_body, '<table><tr><th>Endpoint</th><th>Feature</th>',
                     '<th>HR</th><th>95% CI</th><th>p</th><th>FDR</th></tr>')
for (i in seq_len(nrow(cox_top))) {
  row <- cox_top[i]
  html_body <- paste0(html_body, sprintf(
    '<tr><td>%s</td><td>%s</td><td>%.3f</td><td>[%.2f-%.2f]</td><td>%s</td><td>%s</td></tr>',
    row$endpoint, row$feature_label, row$HR, row$HR_lower, row$HR_upper,
    signif(row$p_value, 3), signif(row$q_value, 3)
  ))
}
html_body <- paste0(html_body, '</table>')

# Insert forest plot
fp <- grep("forest_cox", fig_files, value = TRUE)
if (length(fp) > 0) {
  html_body <- paste0(html_body, sprintf(
    '<img src="%s" alt="Cox forest plot">', encode_img(fp[1])))
}

# Insert KM plots
km_files <- grep("part2_km_", fig_files, value = TRUE)
for (kf in km_files) {
  html_body <- paste0(html_body, sprintf(
    '<img src="%s" alt="KM plot" style="max-width:80%%">', encode_img(kf)))
}

html_body <- paste0(html_body, '</div>')

# --- Part 3 section --------------------------------------------------------
html_body <- paste0(html_body, '
<div class="section">
<h2>Part 3: Chemo Sensitivity</h2>
<p><strong>Central hypothesis</strong>: SecA (proliferative) and SecB (aggressive)
   exert opposing effects on apparent chemo response. At matched SecA levels,
   SecB predicts chemo resistance.</p>
')

# Opposing effects
fp_opp <- grep("opposing_effects", fig_files, value = TRUE)
if (length(fp_opp) > 0) {
  html_body <- paste0(html_body,
    '<h3>Opposing effects check</h3>',
    sprintf('<img src="%s" alt="Opposing effects">', encode_img(fp_opp[1])))
}

# Core confound models table
if (exists("core_results") && nrow(core_results) > 0) {
  html_body <- paste0(html_body, '<h3>Core confound test: SecB at matched SecA</h3>')
  html_body <- paste0(html_body,
    '<table><tr><th>Model</th><th>Term</th><th>OR</th><th>95% CI</th>',
    '<th>p</th><th>AIC</th></tr>')
  for (i in seq_len(nrow(core_results))) {
    row <- core_results[i]
    html_body <- paste0(html_body, sprintf(
      '<tr><td>%s</td><td>%s</td><td>%s</td><td>[%s-%s]</td><td>%s</td><td>%s</td></tr>',
      row$Model, row$Term, row$OR, row$CI_lower, row$CI_upper,
      row$p_value, row$AIC
    ))
  }
  html_body <- paste0(html_body, '</table>')
}

# SecB residualized plot
fp_resid <- grep("secB_residualized", fig_files, value = TRUE)
if (length(fp_resid) > 0) {
  html_body <- paste0(html_body,
    '<h3>SecA-residualized SecB</h3>',
    sprintf('<img src="%s" alt="SecB residualized">', encode_img(fp_resid[1])))
}

# Logistic forest plot
fp_log <- grep("forest_logistic_chemo", fig_files, value = TRUE)
if (length(fp_log) > 0) {
  html_body <- paste0(html_body,
    '<h3>Univariate logistic: chemo response</h3>',
    sprintf('<img src="%s" alt="Logistic forest">', encode_img(fp_log[1])))
}

# 2x2 KM plots
km2x2 <- grep("part3_km_2x2", fig_files, value = TRUE)
for (kf in km2x2) {
  html_body <- paste0(html_body,
    sprintf('<img src="%s" alt="2x2 KM" style="max-width:80%%">',
            encode_img(kf)))
}

html_body <- paste0(html_body, '</div>')

# --- Part 4 section --------------------------------------------------------
html_body <- paste0(html_body, '
<div class="section">
<h2>Part 4: Debulking Composition</h2>
<p>Comparison of secretory subtype composition between optimal vs suboptimal
   surgical debulking.</p>
')

fp_deb <- grep("debulking_boxplots", fig_files, value = TRUE)
if (length(fp_deb) > 0) {
  html_body <- paste0(html_body,
    sprintf('<img src="%s" alt="Debulking boxplots">', encode_img(fp_deb[1])))
}

fp_tern <- grep("ternary_debulking", fig_files, value = TRUE)
if (length(fp_tern) > 0) {
  html_body <- paste0(html_body,
    sprintf('<img src="%s" alt="Ternary plot">', encode_img(fp_tern[1])))
}

# Debulking Wilcoxon top hits
if (exists("wilcox_debulk") && nrow(wilcox_debulk) > 0) {
  deb_top <- wilcox_debulk[!is.na(p_value)][order(p_value)][1:min(10, .N)]
  html_body <- paste0(html_body,
    '<h3>Top 10 Wilcoxon tests (optimal vs suboptimal)</h3>')
  html_body <- paste0(html_body,
    '<table><tr><th>Feature</th><th>Med(opt)</th><th>Med(sub)</th>',
    '<th>p</th><th>FDR</th></tr>')
  for (i in seq_len(nrow(deb_top))) {
    row <- deb_top[i]
    html_body <- paste0(html_body, sprintf(
      '<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>',
      row$label, row$median_optimal, row$median_suboptimal,
      signif(row$p_value, 3), signif(row$q_value, 3)
    ))
  }
  html_body <- paste0(html_body, '</table>')
}

html_body <- paste0(html_body, '</div>')

# --- Part 5 section --------------------------------------------------------
html_body <- paste0(html_body, '
<div class="section">
<h2>Part 5: Stratified TME Analysis</h2>
<p>Patients split by median log2(SecA:SecB) ratio into High (SecA-dominant) vs
   Low (SecB-dominant), crossed with 12-month chemo response. Within each ratio
   stratum: what TME features distinguish resistant from sensitive?</p>
')

html_body <- paste0(html_body, sprintf(
  '<div class="callout callout-info"><strong>Median log2 ratio</strong>: %.3f.
   High ratio = SecA-dominant; Low ratio = SecB-dominant.</div>',
  med_ratio))

fp_top <- grep("top_hits_stratified", fig_files, value = TRUE)
if (length(fp_top) > 0) {
  html_body <- paste0(html_body,
    '<h3>Top discriminating features</h3>',
    sprintf('<img src="%s" alt="Top hits">', encode_img(fp_top[1])))
}

fp_conc <- grep("effect_concordance", fig_files, value = TRUE)
if (length(fp_conc) > 0) {
  html_body <- paste0(html_body,
    '<h3>Effect concordance across strata</h3>',
    sprintf('<img src="%s" alt="Concordance">', encode_img(fp_conc[1])))
}

# Within-stratum KM plots
km_strat <- grep("part5_km_", fig_files, value = TRUE)
for (kf in km_strat) {
  html_body <- paste0(html_body,
    sprintf('<img src="%s" alt="Stratum KM" style="max-width:80%%">',
            encode_img(kf)))
}

html_body <- paste0(html_body, '</div>')

# --- Part 6 section --------------------------------------------------------
html_body <- paste0(html_body, '
<div class="section">
<h2>Part 6: UCell Polarization Analyses</h2>
<p>SecA-SecB polarization variability (SD, IQR) and mean polarization
   examined as clinical predictors.</p>
')

fp_uc_box <- grep("ucell_clinical_boxplots", fig_files, value = TRUE)
if (length(fp_uc_box) > 0) {
  html_body <- paste0(html_body,
    '<h3>Polarization by clinical group</h3>',
    sprintf('<img src="%s" alt="UCell boxplots">', encode_img(fp_uc_box[1])))
}

fp_uc_sc <- grep("ucell_survival_scatter", fig_files, value = TRUE)
if (length(fp_uc_sc) > 0) {
  html_body <- paste0(html_body,
    '<h3>Polarization vs survival</h3>',
    sprintf('<img src="%s" alt="UCell scatter">', encode_img(fp_uc_sc[1])))
}

# UCell survival table
ucell_cox_path <- file.path(comp_dir, "ucell_survival_cox.csv")
if (file.exists(ucell_cox_path)) {
  uc_surv <- fread(ucell_cox_path)
  html_body <- paste0(html_body,
    '<h3>UCell features: Cox PH survival</h3>')
  html_body <- paste0(html_body,
    '<table><tr><th>Feature</th><th>Endpoint</th><th>HR</th>',
    '<th>95% CI</th><th>p</th><th>C-index</th></tr>')
  for (i in seq_len(nrow(uc_surv))) {
    row <- uc_surv[i]
    html_body <- paste0(html_body, sprintf(
      '<tr><td>%s</td><td>%s</td><td>%s</td><td>[%s-%s]</td><td>%s</td><td>%s</td></tr>',
      row$feature, row$endpoint, row$HR, row$HR_lower, row$HR_upper,
      row$p_value, row$concordance
    ))
  }
  html_body <- paste0(html_body, '</table>')
}

html_body <- paste0(html_body, '</div>')

# --- Footer -----------------------------------------------------------------
html_body <- paste0(html_body, '
<hr>
<div class="section">
<h2>Methods Summary</h2>
<ul>
  <li>Input: <code>per_patient_features_v2.csv</code> (reclassified SecA/SecB/Transitioning annotations)</li>
  <li>Survival: Cox PH, 5-year OS and PFS, Kaplan-Meier with log-rank</li>
  <li>Chemo sensitivity: logistic regression (OR per SD), residualization approach</li>
  <li>Debulking: Wilcoxon rank-sum, stage-adjusted logistic regression</li>
  <li>Stratification: median log2(SecA:SecB) ratio split, within-stratum Wilcoxon</li>
  <li>Multiple testing: Benjamini-Hochberg FDR correction</li>
  <li>Neighborhoods: 10-neighborhood assignments (from <code>neighborhood_assignments.csv</code>)</li>
</ul>
</div>

<p class="timestamp">Analysis script: <code>scripts/10_clinical_comprehensive.R</code></p>
</body>
</html>')

# --- Write HTML file --------------------------------------------------------
writeLines(html_body, report_path)
message("Report saved: ", report_path)


# ============================================================================
# SESSION
# ============================================================================

message("\n=== Done ===")
message("Output directory: ", comp_dir)
message("Figures: ", length(fig_files), " PNGs")
message("Report: ", report_path)
log_session()
