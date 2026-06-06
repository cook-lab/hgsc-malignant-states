# ============================================================================
# 06_functional_survival.R
# ----------------------------------------------------------------------------
# PURPOSE: Functional survival association of niche metabolic stress / immune exhaustion.
#
# INPUTS:
#   - output/29_macrophage_niche_survival/ niche scores + clinical
#
# OUTPUTS:
#   - output/29_macrophage_niche_survival/ survival results
#
# MANUSCRIPT PANEL(S): Fig 6 survival support
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

suppressPackageStartupMessages({
  library(data.table); library(SummarizedExperiment)
})

OUT_DIR <- file.path(out_dir, "29_macrophage_niche_survival")
NICHE_CACHE <- file.path(OUT_DIR, "per_cell_niche_scores.rds")
stopifnot(file.exists(NICHE_CACHE))

MIN_CELLS_PER_GROUP <- 15

WT_SAMPLES <- c("sfe_OTB_2384", "sfe_OTB_2417", "sfe_OTB_2432",
                "sfe_OTB_2454", "sfe_OTB_2457", "sfe_OTB_2461",
                "sfe_SP24_24824", "sfe_SP24_25573")

APOPTOSIS_GENES <- c("CASP3", "CASP8", "BAX", "FAS")
PROLIF_GENES    <- c("MKI67", "TOP2A", "CCNB1", "CDC20", "CDK1",
                      "PCNA", "STMN1", "TUBB")
FUNC_GENES <- list(
  mac_M1           = c("CD86","IRF1","STAT1","CXCL9","CXCL10","CXCL11",
                         "IDO1","NFKB1","NFKB2","TNF","IL15","ICAM1",
                         "COTL1","TAP1","CIITA"),
  mac_M2           = c("MRC1","C1QA","C1QB","C1QC","TGFB1","TGFBI",
                         "CD14","TREM2","INHBA","VEGFA","MMP11"),
  mac_TAM          = c("HAVCR2","CD274","LGALS9","ADAM17","ADAM10",
                         "CTSS","ADAMDEC1","FCGR3A","SLAMF7","CCR1"),
  tcell_cytotoxic  = c("GZMA","GZMB","GZMH","PRF1","FGFBP2","CST7",
                         "NKG7","GNLY","FASLG","KLRD1","KLRB1"),
  tcell_exhaustion = c("PDCD1","HAVCR2","LAG3","TIGIT","CTLA4",
                         "BTLA","CD274","TOX"),
  nk_activation    = c("NKG7","GNLY","HCST","KLRD1","CD69",
                         "IFNG","TNF","FASLG","SH2D1A"),
  nk_exhaustion    = c("TIGIT","HAVCR2","LAG3","PDCD1",
                         "LILRB1","LILRB2"),
  bcell_activated  = c("TNFRSF13B","CD27","CD40","CD69","CD80","CD86"),
  plasma_diff      = c("XBP1","MZB1","JCHAIN","DERL3","SSR4","FKBP11","SEC11C")
)

ALL_SCORING_GENES <- unique(c(APOPTOSIS_GENES, PROLIF_GENES,
                               unlist(FUNC_GENES)))
message("Total unique scoring genes needed: ", length(ALL_SCORING_GENES))

CELL_TYPE_SCORES <- list(
  "Macrophage"  = c("apoptosis_score", "proliferation_score",
                     "func_mac_M1", "func_mac_M2", "func_mac_TAM"),
  "T cell"      = c("apoptosis_score", "proliferation_score",
                     "func_tcell_cytotoxic", "func_tcell_exhaustion"),
  "NK cell"     = c("apoptosis_score", "proliferation_score",
                     "func_nk_activation", "func_nk_exhaustion"),
  "B cell"      = c("apoptosis_score", "proliferation_score",
                     "func_bcell_activated"),
  "Plasma cell" = c("apoptosis_score", "proliferation_score",
                     "func_plasma_diff")
)

# ---------------------------------------------------------------------------
# Helper: extract scores for cells of interest from a single SFE
# ---------------------------------------------------------------------------
compute_score <- function(mat_rows, genes) {
  g_ok <- intersect(genes, rownames(mat_rows))
  if (length(g_ok) == 0) return(rep(NA_real_, ncol(mat_rows)))
  m <- as.matrix(mat_rows[g_ok, , drop = FALSE])
  colMeans(m, na.rm = TRUE)
}

extract_scores_for_sfe <- function(sfe, cell_ids_wanted) {
  lc <- assay(sfe, "logcounts")
  # only rows that are scoring genes
  avail_genes <- intersect(ALL_SCORING_GENES, rownames(lc))
  lc_sub <- lc[avail_genes, , drop = FALSE]
  # subset columns to wanted cells (in SFE column order, matching those in cell_ids_wanted)
  keep <- match(cell_ids_wanted, colnames(sfe))
  keep_ok <- !is.na(keep)
  if (!any(keep_ok)) return(NULL)
  lc_sub <- lc_sub[, keep[keep_ok], drop = FALSE]
  # Compute scores
  out <- data.table(cell_id = cell_ids_wanted[keep_ok])
  out[, apoptosis_score    := compute_score(lc_sub, APOPTOSIS_GENES)]
  out[, proliferation_score := compute_score(lc_sub, PROLIF_GENES)]
  for (sg in names(FUNC_GENES)) {
    col <- paste0("func_", sg)
    out[[col]] <- compute_score(lc_sub, FUNC_GENES[[sg]])
  }
  out
}

# ---------------------------------------------------------------------------
# Load niche cache
# ---------------------------------------------------------------------------
message("[cache] loading niche scores...")
niche <- readRDS(NICHE_CACHE)

# ---------------------------------------------------------------------------
# WT: loop samples, extract scores, merge with niche data
# ---------------------------------------------------------------------------
wt_list <- list()
for (s in WT_SAMPLES) {
  samp <- sub("^sfe_", "", s)
  message("== WT scoring ", samp, " ==")
  sfe <- load_sfe(s)
  sample_cells <- niche$wt[sample_key == samp]
  scores <- extract_scores_for_sfe(sfe, sample_cells$cell_id)
  if (is.null(scores)) { rm(sfe); next }
  merged <- merge(sample_cells, scores, by = "cell_id")
  wt_list[[samp]] <- merged
  message("   merged cells: ", format(nrow(merged), big.mark = ","))
  rm(sfe); gc(verbose = FALSE)
}
wt <- rbindlist(wt_list, fill = TRUE)

# ---------------------------------------------------------------------------
# TMA
# ---------------------------------------------------------------------------
message("\n== TMA scoring ==")
sfe_t <- load_sfe("sfe_tma_filtered")
tma_scores <- extract_scores_for_sfe(sfe_t, niche$tma$cell_id)
tma <- merge(niche$tma, tma_scores, by = "cell_id")
message("   TMA merged cells: ", format(nrow(tma), big.mark = ","))
rm(sfe_t); gc(verbose = FALSE)

# ---------------------------------------------------------------------------
# Per-group top-vs-bottom decile delta, per cell type, per score
# ---------------------------------------------------------------------------
compute_deltas <- function(dt, group_col) {
  out <- list()
  for (ct in names(CELL_TYPE_SCORES)) {
    scores <- CELL_TYPE_SCORES[[ct]]
    ct_dt <- dt[cell_label == ct]
    if (nrow(ct_dt) == 0) next
    grp_counts <- ct_dt[, .(
      n_top = sum(stress_decile == 10, na.rm = TRUE),
      n_bot = sum(stress_decile == 1, na.rm = TRUE)
    ), by = group_col]
    elig_grp <- grp_counts[n_top >= MIN_CELLS_PER_GROUP &
                             n_bot >= MIN_CELLS_PER_GROUP][[group_col]]
    ct_elig <- ct_dt[get(group_col) %in% elig_grp]
    if (nrow(ct_elig) == 0) next
    for (sc in scores) {
      d <- ct_elig[, .(
        med_top = median(get(sc)[stress_decile == 10], na.rm = TRUE),
        med_bot = median(get(sc)[stress_decile == 1],  na.rm = TRUE)
      ), by = group_col]
      d[, delta := med_top - med_bot]
      d <- d[is.finite(delta)]
      if (nrow(d) < 2) next
      w <- wilcox.test(d$delta, mu = 0, alternative = "two.sided")
      out[[paste(ct, sc, sep = "__")]] <- data.table(
        cell_type = ct, score = sc,
        n_groups = nrow(d),
        median_delta = median(d$delta),
        pct_positive = 100 * mean(d$delta > 0),
        p_wilcox = w$p.value
      )
    }
  }
  if (length(out) == 0) return(data.table())
  rbindlist(out, fill = TRUE)
}

message("\n=== WT functional survival deltas ===")
wt_res <- compute_deltas(wt, "sample_key")
if (nrow(wt_res) > 0) setorder(wt_res, cell_type, score)
print(wt_res)

message("\n=== TMA functional survival deltas (by patient) ===")
tma_res <- compute_deltas(tma[!is.na(patient_id) & patient_id != ""],
                           "patient_id")
if (nrow(tma_res) > 0) setorder(tma_res, cell_type, score)
print(tma_res)

wt_res[,  cohort := "WT"]
tma_res[, cohort := "TMA"]
all_res <- rbind(wt_res, tma_res, fill = TRUE)
fwrite(all_res, file.path(OUT_DIR, "functional_survival_summary.csv"))
message("\nSaved: ", file.path(OUT_DIR, "functional_survival_summary.csv"))

# Per-patient deltas for figures
compute_per_group <- function(dt, group_col) {
  out <- list()
  for (ct in names(CELL_TYPE_SCORES)) {
    ct_dt <- dt[cell_label == ct]
    if (nrow(ct_dt) == 0) next
    for (sc in CELL_TYPE_SCORES[[ct]]) {
      d <- ct_dt[, .(
        n_top = sum(stress_decile == 10, na.rm = TRUE),
        n_bot = sum(stress_decile == 1, na.rm = TRUE),
        med_top = median(get(sc)[stress_decile == 10], na.rm = TRUE),
        med_bot = median(get(sc)[stress_decile == 1],  na.rm = TRUE)
      ), by = group_col]
      d[, delta := med_top - med_bot]
      d[, cell_type := ct]
      d[, score := sc]
      out[[paste(ct, sc, sep = "__")]] <- d
    }
  }
  rbindlist(out, fill = TRUE)
}
wt_per  <- compute_per_group(wt, "sample_key")
tma_per <- compute_per_group(tma[!is.na(patient_id) & patient_id != ""],
                               "patient_id")
wt_per[, cohort := "WT"]; tma_per[, cohort := "TMA"]
all_per <- rbind(wt_per, tma_per, fill = TRUE)
fwrite(all_per, file.path(OUT_DIR, "functional_survival_per_patient.csv"))
message("Saved: ", file.path(OUT_DIR, "functional_survival_per_patient.csv"))

# Save scored data.tables for figures
saveRDS(list(wt = wt, tma = tma),
        file.path(OUT_DIR, "functional_survival_scored_cells.rds"))
message("Saved: functional_survival_scored_cells.rds")
