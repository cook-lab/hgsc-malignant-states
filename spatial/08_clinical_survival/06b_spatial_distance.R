# ============================================================================
# 06b_spatial_distance.R — macrophage vs lymphocyte distance to hypoxic anchors
# ----------------------------------------------------------------------------
# PURPOSE: Spatial-proximity test behind Fig 6 panel 29D. For each TME cell,
#   compute the distance to its nearest top-decile-hypoxia epithelial anchor,
#   then test whether macrophages sit systematically closer to those hypoxic
#   anchors than lymphocytes do. Per group (WT sample / TMA core->patient):
#     delta = median(D | macrophage) - median(D | lymphocyte)
#   Paired Wilcoxon (1-sided, H1: delta < 0) across samples / patients.
#   Complements the GLMM + paired-enrichment tests on a raw distance axis.
#
# INPUTS:
#   - 8 published whole-tissue SFEs + sfe_tma_filtered (load_sfe), each carrying
#     pathway_hypoxia (9b pathway scoring) and cell_label
#   - <output_root>/06f_reclassification_polarization/reclassified_xenium_scores.csv
#     (06f polarization override; written by spatial/03_annotation_polarization/04)
#
# OUTPUTS (<output_root>/29_macrophage_niche_survival/):
#   - spatial_distance_per_sample_wt.csv     (READ by 07_figures.R panel 29D)
#   - spatial_distance_per_patient_tma.csv   (READ by 07_figures.R panel 29D)
#   - spatial_distance_summary.csv, spatial_distance_per_cell.csv,
#     spatial_distance_results.rds
#
# MANUSCRIPT PANEL(S): Fig 6 panel 29D (WT + TMA). Must run BEFORE 07_figures.R.
# RUNTIME TIER: heavy
#
# Migrated from 2026_final_xenium_analysis/scripts/29d_spatial_distance.R.
# Analytical logic preserved; paths routed through central config; WT cohort
# pinned to the published 8 via CFG; epithelial label "Transitioning" ->
# "Intermediate" on read. Nearest-neighbour distances use FNN (already a
# pipeline dependency) instead of RANN — both are EXACT 1-NN, so the Euclidean
# distances are identical; this avoids introducing a new package.
# ============================================================================

# --- Config + shared setup (replaces hardcoded /Volumes/CookLab/Sarah paths) ---
here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
source(file.path(here, "..", "..", "config", "config.R"))   # CFG, cfg_obj, cfg_path
source(file.path(here, "..", "00_setup", "00_setup.R"))      # load_sfe, cohort, log_session
set.seed(CFG$seed)

suppressPackageStartupMessages({
  library(data.table); library(FNN)
})

OUT_DIR <- file.path(out_dir, "29_macrophage_niche_survival")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Published whole-tissue cohort (8) from config — excludes the FTE whole tissues.
WT_SAMPLES <- sfe_names_wt

EPI_LBLS   <- c("SecA epithelium", "Intermediate epithelium", "SecB epithelium")
LYMPH_LBLS <- c("T cell", "NK cell", "B cell", "Plasma cell")

# --- 06f override (same pattern as 03_compute_niche_metabolic_scores.R) -------
f06f <- file.path(out_dir, "06f_reclassification_polarization",
                  "reclassified_xenium_scores.csv")
stopifnot(file.exists(f06f))
recl <- fread(f06f, select = c("sample", "barcode_orig", "cell_label_06f"))

override_with_06f <- function(sfe, sample_key) {
  sub <- recl[sample == sample_key]
  if (nrow(sub) == 0) return(sfe)
  m <- match(colnames(sfe), sub$barcode_orig)
  hit <- !is.na(m)
  lab <- as.character(sfe$cell_label)
  lab[hit] <- sub$cell_label_06f[m[hit]]
  # Standardize legacy "Transitioning epithelium" -> "Intermediate epithelium"
  # so EPI_LBLS matching does not silently drop the Intermediate epitype.
  lab[lab == "Transitioning epithelium"] <- "Intermediate epithelium"
  sfe$cell_label <- lab
  sfe
}

# --- Nearest-anchor distance per query cell ----------------------------------
per_group <- function(coords, labels, hypo, group_name) {
  is_epi <- labels %in% EPI_LBLS
  if (sum(is_epi) == 0) return(NULL)
  epi_hypo <- hypo[is_epi]
  thr <- quantile(epi_hypo, 0.9, na.rm = TRUE)
  if (!is.finite(thr)) return(NULL)
  is_anchor <- is_epi & hypo >= thr
  if (sum(is_anchor) < 10) return(NULL)
  is_mac   <- labels == "Macrophage"
  is_lymph <- labels %in% LYMPH_LBLS
  if (sum(is_mac) < 10 || sum(is_lymph) < 10) return(NULL)

  anchor_coords <- coords[is_anchor, , drop = FALSE]
  mac_coords    <- coords[is_mac,    , drop = FALSE]
  lymph_coords  <- coords[is_lymph,  , drop = FALSE]

  # Nearest anchor distance per query cell (exact 1-NN; FNN == RANN here).
  d_mac   <- FNN::get.knnx(anchor_coords, mac_coords,   k = 1)$nn.dist[, 1]
  d_lymph <- FNN::get.knnx(anchor_coords, lymph_coords, k = 1)$nn.dist[, 1]

  data.table(
    group_id  = group_name,
    cell_type = c(rep("Macrophage", length(d_mac)),
                  rep("Lymphocyte", length(d_lymph))),
    dist      = c(d_mac, d_lymph)
  )
}

# ---------------------------------------------------------------------------
# Whole tissue (per sample)
# ---------------------------------------------------------------------------
wt_all <- list()
for (s in WT_SAMPLES) {
  samp <- sub("^sfe_", "", s)
  message("== WT distance ", samp, " ==")
  sfe <- load_sfe(s)
  sfe <- override_with_06f(sfe, s)
  coords <- SpatialFeatureExperiment::spatialCoords(sfe)
  labs   <- as.character(sfe$cell_label)
  hypo   <- as.numeric(sfe$pathway_hypoxia)
  res <- per_group(coords, labs, hypo, samp)
  if (!is.null(res)) {
    res[, cohort := "WT"]; res[, sample_id := samp]
    wt_all[[samp]] <- res
  }
  rm(sfe); gc(verbose = FALSE)
}
wt_all <- rbindlist(wt_all, fill = TRUE)
message("WT distance rows: ", format(nrow(wt_all), big.mark = ","))

# ---------------------------------------------------------------------------
# TMA (per core, attributed to patient)
# ---------------------------------------------------------------------------
message("\n== TMA distance (per-core) ==")
sfe_t <- load_sfe("sfe_tma_filtered")
sfe_t <- override_with_06f(sfe_t, "sfe_tma")
coords_t <- SpatialFeatureExperiment::spatialCoords(sfe_t)
labs_t   <- as.character(sfe_t$cell_label)
hypo_t   <- as.numeric(sfe_t$pathway_hypoxia)
core_ids <- as.character(sfe_t$core_id)
pat_ids  <- as.character(sfe_t$patient_id)

unique_cores <- unique(core_ids[!is.na(core_ids) & core_ids != ""])
tma_all <- list()
for (ci in unique_cores) {
  in_core <- core_ids == ci & !is.na(core_ids)
  if (sum(in_core) < 100) next
  res <- per_group(coords_t[in_core, , drop = FALSE],
                   labs_t[in_core], hypo_t[in_core], ci)
  if (!is.null(res)) {
    res[, cohort := "TMA"]
    res[, sample_id := "TMA"]
    # Attach patient_id (most common in this core)
    pat_in_core <- pat_ids[in_core]
    pat_in_core <- pat_in_core[!is.na(pat_in_core) & pat_in_core != ""]
    if (length(pat_in_core) > 0) {
      res[, patient_id := names(sort(table(pat_in_core), decreasing = TRUE))[1]]
    } else {
      res[, patient_id := NA_character_]
    }
    tma_all[[ci]] <- res
  }
}
tma_all <- rbindlist(tma_all, fill = TRUE)
rm(sfe_t); gc(verbose = FALSE)
message("TMA distance rows: ", format(nrow(tma_all), big.mark = ","))

# ---------------------------------------------------------------------------
# Summarise: median distance per cell type per group
# ---------------------------------------------------------------------------
wt_med  <- wt_all[, .(med_dist = median(dist, na.rm = TRUE), n_cells = .N),
                  by = .(sample_id, cell_type)]
tma_med <- tma_all[!is.na(patient_id) & patient_id != "",
                   .(med_dist = median(dist, na.rm = TRUE), n_cells = .N),
                   by = .(patient_id, cell_type)]

# Per-group delta: median(dist | Macrophage) - median(dist | Lymphocyte).
# Macrophages closer to hypoxia => delta < 0.
wt_delta <- dcast(wt_med, sample_id ~ cell_type, value.var = "med_dist")
wt_delta[, delta := Macrophage - Lymphocyte]
wt_delta[, cohort := "WT"]

tma_delta <- dcast(tma_med, patient_id ~ cell_type, value.var = "med_dist")
tma_delta[, delta := Macrophage - Lymphocyte]
tma_delta[, cohort := "TMA"]

cat("\n=== WT per-sample median distances (um) to top-decile-hypoxia anchors ===\n")
print(wt_delta)
cat("\n=== TMA per-patient median distances (um), most-negative deltas ===\n")
print(head(tma_delta[order(delta)], 10))

# ---------------------------------------------------------------------------
# Paired Wilcoxon on deltas (H0: delta == 0; expected: delta < 0)
# ---------------------------------------------------------------------------
test_wt <- wilcox.test(wt_delta$Macrophage, wt_delta$Lymphocyte,
                       paired = TRUE, alternative = "less")
test_tma <- wilcox.test(tma_delta[!is.na(delta), Macrophage],
                        tma_delta[!is.na(delta), Lymphocyte],
                        paired = TRUE, alternative = "less")

summary_dt <- data.table(
  cohort          = c("WT", "TMA"),
  n_groups        = c(sum(is.finite(wt_delta$delta)),
                      sum(is.finite(tma_delta$delta))),
  median_delta    = c(median(wt_delta$delta, na.rm = TRUE),
                      median(tma_delta$delta, na.rm = TRUE)),
  pct_mac_closer  = c(100 * mean(wt_delta$delta < 0, na.rm = TRUE),
                      100 * mean(tma_delta$delta < 0, na.rm = TRUE)),
  p_wilcox_1sided = c(test_wt$p.value, test_tma$p.value)
)
cat("\n=== SUMMARY — median distance of MAC vs LYMPH to top-hypoxia anchor ===\n")
print(summary_dt)
cat("\nInterpretation: negative median_delta + pct_mac_closer high means\n")
cat("macrophages sit closer to top-hypoxia anchors than lymphocytes.\n")

# --- Save --------------------------------------------------------------------
fwrite(summary_dt, file.path(OUT_DIR, "spatial_distance_summary.csv"))
fwrite(wt_delta,   file.path(OUT_DIR, "spatial_distance_per_sample_wt.csv"))
fwrite(tma_delta,  file.path(OUT_DIR, "spatial_distance_per_patient_tma.csv"))
fwrite(rbind(wt_all, tma_all, fill = TRUE),
       file.path(OUT_DIR, "spatial_distance_per_cell.csv"))
saveRDS(list(wt = wt_all, tma = tma_all,
             wt_delta = wt_delta, tma_delta = tma_delta,
             summary = summary_dt),
        file.path(OUT_DIR, "spatial_distance_results.rds"))
message("\nSaved spatial distance outputs to ", OUT_DIR)
log_session()
