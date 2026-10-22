# ============================================================================
# 05b_core_qc.R — TMA core QC: annotation-level outlier detection
# ----------------------------------------------------------------------------
# PURPOSE: Flag non-representative tumour cores of a TMA array and write the
#   per-core QC table that spatial/02_qc/06b_final_tma_cores.py reads to build
#   the final TMA core set. Flags:
#     - flag_small: < 200 cells
#     - flag_low_epi: < 10% epithelial
#     - flag_low_confidence: median SingleR score < p5 of the cores in the run
#       OR > 10% pruned
#     - flag_composition_outlier: Mahalanobis distance on CLR-transformed 18-type
#       composition exceeds the chi-sq p<0.01 threshold
#     - flag_discordant_pair: within-patient core pair with JSD > p90 OR
#       Spearman composition correlation < 0.5
#   recommendation = exclude (>= 2 flags) / review (1) / pass (0).
#   Fully deterministic (quantiles / Mahalanobis / JSD — no RNG).
#
# INPUTS:
#   - the SFE(s) named in CORE_QC_TMAS (load_sfe; default sfe_tma), carrying
#     singler_label/score/pruned, polarization_UCell, core_id, sample_type,
#     patient_id and the cell_label column (read, then replaced by the SingleR
#     labels when CORE_QC_LABELS=singler)
#
# FINAL RUN (2026-10-19, author): one run per array with the same method and SingleR labels
#   (CORE_QC_LABELS=singler, the default): CORE_QC_TMAS=sfe_tma | sfe_COEUR_TMA_G | sfe_COEUR_TMA_H,
#   CORE_QC_OUT_TAG=_singler_<array>. 06b_final_tma_cores.py combines the three into the one TMA core set
#   used by every manuscript panel.
#
# OUTPUTS (<output_root>/07_core_qc/):
#   - core_qc_summary<CORE_QC_OUT_TAG>.csv  (per-core metrics, flags, recommendation;
#       final run: core_qc_summary_singler_{tma,COEUR_TMA_G,COEUR_TMA_H}.csv)
#
# MANUSCRIPT PANEL(S): none directly — the core-QC flags behind the final TMA
#   core set (TMA arms of Figs 4–7, Figures S12–S15, Tables S8 and S11). Runs
#   after spatial/03_annotation_polarization/ (needs the SingleR labels and the
#   polarization scores).
# RUNTIME TIER: moderate
#
# Paths come from the central config. The epithelial label "Transitioning" is
# standardized to "Intermediate" on read, before any composition counting.
# ============================================================================

# --- Config + shared setup ---------------------------------------------------
here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
source(file.path(here, "..", "..", "config", "config.R"))   # CFG, cfg_obj, cfg_path
source(file.path(here, "..", "00_setup", "00_setup.R"))      # load_sfe, data.table, log_session
set.seed(CFG$seed)

qc_dir <- file.path(out_dir, "07_core_qc")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# --- Lineage groupings (cell_label; "Intermediate" was "Transitioning") ------
epi_types <- c("Ciliated epithelium", "SecA epithelium",
               "Intermediate epithelium", "SecB epithelium")
secretory_subtypes <- c("SecA epithelium", "Intermediate epithelium", "SecB epithelium")
immune_types <- c("T cell", "NK cell", "B cell", "Plasma cell", "Macrophage",
                  "Conventional dendritic cell", "Plasmacytoid dendritic cell",
                  "Neutrophil", "Mast cell")
stromal_types  <- c("Fibroblast", "Smooth muscle")
vascular_types <- c("Pericyte", "Endothelial")
celltype_order <- c("Ciliated epithelium", "SecA epithelium",
                    "Intermediate epithelium", "SecB epithelium",
                    "Mesothelial", "Fibroblast", "Smooth muscle", "Pericyte",
                    "Endothelial", "T cell", "NK cell", "B cell", "Plasma cell",
                    "Macrophage", "Conventional dendritic cell",
                    "Plasmacytoid dendritic cell", "Neutrophil", "Mast cell")

# Jensen-Shannon divergence between two probability vectors
jsd <- function(p, q) {
  eps <- 1e-10
  p <- p + eps; p <- p / sum(p)
  q <- q + eps; q <- q / sum(q)
  m <- 0.5 * (p + q)
  0.5 * sum(p * log2(p / m)) + 0.5 * sum(q * log2(q / m))
}

# --- 1. Load TMA data --------------------------------------------------------
# CORE_QC_TMAS names the array(s) to QC (comma-separated; several are pooled
# into one QC, "v2"). The final run used one array per run. Default = sfe_tma.
TMA_SFES <- {
  .v <- Sys.getenv("CORE_QC_TMAS", "")
  if (nzchar(.v)) trimws(strsplit(.v, ",")[[1]]) else "sfe_tma"
}
V2 <- length(TMA_SFES) > 1L
cat("Core QC TMAs:", paste(TMA_SFES, collapse = ", "), if (V2) "  [v2]" else "", "\n")

load_one_tma <- function(nm) {
  sfe <- load_sfe(nm)
  dt <- data.table(
    tma            = nm,
    cell_id        = colnames(sfe),
    cell_label     = as.character(sfe$cell_label),
    singler_label  = sfe$singler_label,
    singler_score  = sfe$singler_score,
    singler_pruned = sfe$singler_pruned,
    core_id_raw    = as.character(sfe$core_id),
    sample_type    = as.character(sfe$sample_type),
    patient_id_raw = as.character(sfe$patient_id)
  )
  # CORE_QC_LABELS=singler (the default, used in the final run): composition flags use the SingleR labels,
  # with Secretory epithelium split into SecA / Intermediate / SecB at the frozen polarization thresholds
  # (as spatial/03_annotation_polarization/06_extract_singler.R). Any other value keeps the SFE cell_label.
  if (Sys.getenv("CORE_QC_LABELS", "singler") == "singler") {
    pol <- if ("polarization_UCell" %in% colnames(colData(sfe))) sfe$polarization_UCell else NA_real_
    lab <- as.character(sfe$singler_label)
    lab[lab %in% "Transitioning epithelium"] <- "Intermediate epithelium"
    is_sec <- !is.na(lab) & lab == "Secretory epithelium"
    split <- ifelse(is.na(pol), NA_character_,
                    ifelse(pol >= 0.210501, "SecB epithelium",
                           ifelse(pol <= 0.109429, "SecA epithelium", "Intermediate epithelium")))
    lab[is_sec] <- split[is_sec]
    dt[, cell_label := lab]
  }
  rm(sfe); gc(verbose = FALSE)
  dt
}
tma_dt <- rbindlist(lapply(TMA_SFES, load_one_tma))
cat("Composition labels:", Sys.getenv("CORE_QC_LABELS", "singler"), "\n")

# Namespace core_id / patient_id by TMA so numeric IDs never collide across TMAs
# (raw values kept for reporting). The logic below keys on the namespaced
# core_id / patient_id.
tma_dt[, core_id    := fifelse(is.na(core_id_raw),    NA_character_, paste(tma, core_id_raw,    sep = "::"))]
tma_dt[, patient_id := fifelse(is.na(patient_id_raw), NA_character_, paste(tma, patient_id_raw, sep = "::"))]

# Standardize the legacy epithelial label BEFORE any composition counting so the
# epithelial / secretory percentages count the same cells under the new name.
tma_dt[cell_label == "Transitioning epithelium", cell_label := "Intermediate epithelium"]

# On-core tumour cells only
tma_dt <- tma_dt[!is.na(core_id_raw) & core_id_raw != "Off core"]
tumour_dt <- tma_dt[sample_type == "tumour"]
cat("Tumour cells:", format(nrow(tumour_dt), big.mark = ","),
    "| cores:", uniqueN(tumour_dt$core_id),
    "| patients:", uniqueN(tumour_dt$patient_id),
    "| TMAs:", uniqueN(tumour_dt$tma), "\n")

# --- 2. Per-core composition (18 types) + summary metrics --------------------
core_comp   <- tumour_dt[, .N, by = .(core_id, cell_label)]
core_totals <- tumour_dt[, .(n_cells = .N), by = core_id]
core_comp   <- merge(core_comp, core_totals, by = "core_id")
core_comp[, pct := 100 * N / n_cells]

# Fill zeros for missing types
all_combos <- CJ(core_id = unique(core_comp$core_id), cell_label = celltype_order)
core_comp  <- merge(all_combos, core_comp, by = c("core_id", "cell_label"), all.x = TRUE)
core_comp[is.na(N), `:=`(N = 0L, pct = 0)]
totals_map <- core_totals[, setNames(n_cells, core_id)]
core_comp[is.na(n_cells), n_cells := totals_map[core_id]]
core_comp[, pct := 100 * N / n_cells]

core_meta <- unique(tumour_dt[, .(core_id, patient_id, tma, core_id_raw, patient_id_raw)])
core_summary <- tumour_dt[, .(
  n_cells              = .N,
  pct_epithelial       = 100 * sum(cell_label %in% epi_types) / .N,
  pct_secretory        = 100 * sum(cell_label %in% secretory_subtypes) / .N,
  pct_immune           = 100 * sum(cell_label %in% immune_types) / .N,
  pct_stromal          = 100 * sum(cell_label %in% stromal_types) / .N,
  pct_vascular         = 100 * sum(cell_label %in% vascular_types) / .N,
  pct_mesothelial      = 100 * sum(cell_label == "Mesothelial") / .N,
  median_singler_score = median(singler_score, na.rm = TRUE),
  pct_pruned           = 100 * sum(is.na(singler_pruned)) / .N
), by = core_id]

# Secretory subtype proportions (among secretory cells only)
sec_comp <- tumour_dt[cell_label %in% secretory_subtypes,
                      .(n_sec = .N,
                        pct_progenitor   = 100 * sum(cell_label == secretory_subtypes[1]) / .N,
                        pct_intermediate = 100 * sum(cell_label == secretory_subtypes[2]) / .N,
                        pct_adaptive     = 100 * sum(cell_label == secretory_subtypes[3]) / .N),
                      by = core_id]
core_summary <- merge(core_summary, core_meta, by = "core_id")
core_summary <- merge(core_summary, sec_comp, by = "core_id", all.x = TRUE)
cat("Per-core summary computed for", nrow(core_summary), "tumour cores\n")

# --- 3. Flag problem cores ---------------------------------------------------
core_summary[, flag_small := n_cells < 200]
core_summary[, flag_low_epi := pct_epithelial < 10]
score_p5 <- quantile(core_summary$median_singler_score, 0.05)
core_summary[, flag_low_confidence := median_singler_score < score_p5 | pct_pruned > 10]

# Composition outlier (Mahalanobis on CLR-transformed composition)
comp_wide <- dcast(core_comp[core_id %in% core_summary$core_id],
                   core_id ~ cell_label, value.var = "pct", fill = 0)
comp_mat <- as.matrix(comp_wide[, -1])
rownames(comp_mat) <- comp_wide$core_id

clr_transform <- function(mat) {
  mat_pseudo <- mat + 0.5          # pseudocount for zeros
  log_mat <- log(mat_pseudo)
  log_mat - rowMeans(log_mat)
}
comp_clr <- clr_transform(comp_mat)
center  <- colMeans(comp_clr)
cov_mat <- cov(comp_clr)
maha_dist <- tryCatch(
  mahalanobis(comp_clr, center, cov_mat),
  error = function(e) {
    cov_reg <- cov_mat + diag(0.01, ncol(cov_mat))   # regularized fallback
    mahalanobis(comp_clr, center, cov_reg)
  })
core_summary[, maha_dist := maha_dist[match(core_id, names(maha_dist))]]
maha_threshold <- qchisq(0.99, df = ncol(comp_clr) - 1)
core_summary[, flag_composition_outlier := maha_dist > maha_threshold]
cat("Mahalanobis threshold (chi-sq p<0.01, df=", ncol(comp_clr) - 1, "):",
    round(maha_threshold, 1), "| outliers:",
    sum(core_summary$flag_composition_outlier), "\n")

flag_cols <- c("flag_small", "flag_low_epi", "flag_low_confidence",
               "flag_composition_outlier")
core_summary[, n_flags := rowSums(.SD), .SDcols = flag_cols]

# --- 4. Within-patient concordance -> discordant-pair flag -------------------
cores_per_patient <- core_summary[, .(n_cores = .N), by = patient_id]
multi_core_patients <- cores_per_patient[n_cores >= 2, patient_id]

concordance_list <- list()
# Guard: only cores that carry a composition row (a core of entirely unassigned
# cells has none) can be compared pairwise.
.n_no_comp <- sum(!core_summary$core_id %in% rownames(comp_mat))
if (.n_no_comp > 0)
  cat("  Note:", .n_no_comp, "core(s) lack a composition row (all-unassigned); excluded from concordance.\n")
for (pid in multi_core_patients) {
  pid_cores <- core_summary[patient_id == pid, core_id]
  pid_cores <- pid_cores[pid_cores %in% rownames(comp_mat)]
  if (length(pid_cores) < 2) next
  for (i in 1:(length(pid_cores) - 1)) {
    for (j in (i + 1):length(pid_cores)) {
      c1 <- pid_cores[i]; c2 <- pid_cores[j]
      p1 <- comp_mat[c1, ] / 100
      p2 <- comp_mat[c2, ] / 100
      jsd_val <- jsd(p1, p2)
      cor_val <- cor(p1, p2, method = "spearman")
      s1 <- core_summary[core_id == c1]; s2 <- core_summary[core_id == c2]
      delta_epi <- abs(s1$pct_epithelial - s2$pct_epithelial)
      delta_adaptive <- abs(ifelse(is.na(s1$pct_adaptive), 0, s1$pct_adaptive) -
                            ifelse(is.na(s2$pct_adaptive), 0, s2$pct_adaptive))
      concordance_list[[paste(c1, c2, sep = "_")]] <- data.table(
        tma = s1$tma, patient_id = pid, core_1 = c1, core_2 = c2,
        n_cells_1 = s1$n_cells, n_cells_2 = s2$n_cells,
        jsd = jsd_val, spearman_cor = cor_val,
        delta_epi = delta_epi, delta_adaptive = delta_adaptive,
        pct_epi_1 = s1$pct_epithelial, pct_epi_2 = s2$pct_epithelial)
    }
  }
}
concordance <- rbindlist(concordance_list)

jsd_p90 <- quantile(concordance$jsd, 0.90)
concordance[, flag_discordant := jsd > jsd_p90 | spearman_cor < 0.5]
discordant_cores <- unique(c(concordance[flag_discordant == TRUE, core_1],
                             concordance[flag_discordant == TRUE, core_2]))
core_summary[, flag_discordant_pair := core_id %in% discordant_cores]

# Recompute n_flags with the discordance flag
flag_cols_all <- c(flag_cols, "flag_discordant_pair")
core_summary[, n_flags := rowSums(.SD), .SDcols = flag_cols_all]

# --- 5. Recommendation -------------------------------------------------------
core_summary[, recommendation := fifelse(
  n_flags >= 2, "exclude",
  fifelse(n_flags == 1, "review", "pass"))]

cat("\nRecommendations:\n"); print(table(core_summary$recommendation))
cat("Flag breakdown:\n")
for (f in flag_cols_all) cat(sprintf("  %s: %d cores\n", f, sum(core_summary[[f]])))

# --- 6. Save outputs ---------------------------------------------------------
# Output suffix. CORE_QC_OUT_TAG names the per-array report (final run:
# "_singler_<array>"); without it a pooled run writes "_v2" and a single-array
# run no suffix.
suffix <- {
  .t <- Sys.getenv("CORE_QC_OUT_TAG", "")
  if (nzchar(.t)) .t else if (V2) "_v2" else ""
}
front <- c("tma", "core_id_raw", "patient_id_raw", "core_id", "patient_id")
setcolorder(core_summary, c(intersect(front, names(core_summary)),
                            setdiff(names(core_summary), front)))
write.csv(core_summary[order(tma, core_id_raw)],
          file.path(qc_dir, paste0("core_qc_summary", suffix, ".csv")), row.names = FALSE)
cat("  ->", qc_dir, "\n")
log_session()
