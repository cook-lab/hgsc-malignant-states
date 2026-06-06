# ============================================================================
# 01_compute_cell_nucleus_morphometrics.R
# ----------------------------------------------------------------------------
# PURPOSE: Compute per-cell cell/nucleus morphometrics (area, N:C ratio) and eligibility (>=30 cells per epitype).
#
# INPUTS:
#   - SFEs (load_sfe) with cell_label + cellSeg/nucSeg geometries; 06f override
#
# OUTPUTS:
#   - output/33_morphometrics/ per-cell morphometrics + eligibility tables
#
# MANUSCRIPT PANEL(S): Fig 5D/5E (upstream)
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

suppressPackageStartupMessages({
  library(data.table); library(sf)
  library(SpatialFeatureExperiment); library(SummarizedExperiment)
})

OUT_DIR <- file.path(out_dir, "33_morphometrics")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE,
            showWarnings = FALSE)

WT_SAMPLES <- c("sfe_OTB_2384", "sfe_OTB_2417", "sfe_OTB_2432",
                "sfe_OTB_2454", "sfe_OTB_2457", "sfe_OTB_2461",
                "sfe_SP24_24824", "sfe_SP24_25573")
EPI_LBLS <- c("SecA epithelium", "Intermediate epithelium", "SecB epithelium")
CHUNK    <- 25000   # PCA-on-polygon chunk size

# ---------------------------------------------------------------------------
# 06f override
# ---------------------------------------------------------------------------
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
  sfe$cell_label <- lab
  sfe
}

# ---------------------------------------------------------------------------
# Eccentricity from polygon vertices via PCA
#  - Returns 1 - sqrt(lambda_min / lambda_max) where lambda are eigenvalues
#    of the covariance of vertex coords. 0 = circle, →1 = elongated.
# ---------------------------------------------------------------------------
ecc_from_polygon <- function(poly) {
  cc <- sf::st_coordinates(poly)
  if (nrow(cc) < 3) return(NA_real_)
  m <- cov(cc[, 1:2])
  if (any(!is.finite(m))) return(NA_real_)
  ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
  ev <- pmax(ev, 0)
  if (max(ev) <= 0) return(NA_real_)
  1 - sqrt(min(ev) / max(ev))
}

ecc_chunked <- function(sf_obj, chunk = CHUNK) {
  n <- nrow(sf_obj)
  ecc <- numeric(n)
  for (i in seq(1, n, by = chunk)) {
    j <- min(i + chunk - 1, n)
    ecc[i:j] <- vapply(seq_len(j - i + 1), function(k) {
      ecc_from_polygon(sf_obj$geometry[i + k - 1])
    }, numeric(1))
  }
  ecc
}

# ---------------------------------------------------------------------------
# Process one SFE → per-cell morphometric data.table
# ---------------------------------------------------------------------------
process_sfe <- function(sfe, cohort, sample_key, group_col = NULL) {
  message(sprintf("[%s] starting (n_cells=%d)", sample_key, ncol(sfe)))
  lbl <- as.character(sfe$cell_label)
  is_epi <- lbl %in% EPI_LBLS
  if (sum(is_epi) == 0) {
    message("  no epithelial cells; skipping")
    return(NULL)
  }
  message(sprintf("  epithelial cells: %d", sum(is_epi)))

  # Get cell + nucleus segmentation polygons
  cell_sf <- tryCatch(sf::st_as_sf(colGeometry(sfe, "cellSeg")),
                       error = function(e) NULL)
  nuc_sf  <- tryCatch(sf::st_as_sf(colGeometry(sfe, "nucSeg")),
                       error = function(e) NULL)
  if (is.null(cell_sf) || is.null(nuc_sf)) {
    message("  cellSeg / nucSeg unavailable; skipping")
    return(NULL)
  }

  # cellSeg / nucSeg should align row-wise with colData
  if (nrow(cell_sf) != ncol(sfe) || nrow(nuc_sf) != ncol(sfe)) {
    message(sprintf("  WARNING: polygon counts mismatch (cell=%d, nuc=%d, ncol=%d)",
                     nrow(cell_sf), nrow(nuc_sf), ncol(sfe)))
  }

  # Subset to epi cells (assume row alignment with colData)
  cell_e <- cell_sf[is_epi, ]
  nuc_e  <- nuc_sf[is_epi, ]
  lbl_e  <- lbl[is_epi]
  cells_x <- spatialCoords(sfe)[is_epi, 1]
  cells_y <- spatialCoords(sfe)[is_epi, 2]
  cell_id <- colnames(sfe)[is_epi]

  if (is.null(group_col)) {
    grp_id <- rep(sample_key, length(cell_id))
  } else {
    grp_id <- as.character(colData(sfe)[[group_col]][is_epi])
  }

  # Validate / repair polygons
  cell_e <- sf::st_make_valid(cell_e)
  nuc_e  <- sf::st_make_valid(nuc_e)
  v_cell <- sf::st_is_valid(cell_e)
  v_nuc  <- sf::st_is_valid(nuc_e)
  keep   <- v_cell & v_nuc

  message(sprintf("  valid polygons: cell=%.1f%% nuc=%.1f%% both=%.1f%%",
                   100 * mean(v_cell), 100 * mean(v_nuc), 100 * mean(keep)))
  if (sum(keep) == 0) return(NULL)

  cell_e  <- cell_e[keep, ]
  nuc_e   <- nuc_e[keep, ]
  lbl_e   <- lbl_e[keep]
  cells_x <- cells_x[keep]
  cells_y <- cells_y[keep]
  cell_id <- cell_id[keep]
  grp_id  <- grp_id[keep]

  # --- vector ops on full subset ---
  message("  computing area / perimeter / circularity / solidity ...")
  cell_area     <- as.numeric(sf::st_area(cell_e))
  nuc_area      <- as.numeric(sf::st_area(nuc_e))
  cell_perim    <- as.numeric(sf::st_length(sf::st_boundary(cell_e)))
  nuc_perim     <- as.numeric(sf::st_length(sf::st_boundary(nuc_e)))
  cell_circ     <- ifelse(cell_perim > 0,
                            (4 * pi * cell_area) / (cell_perim^2),
                            NA_real_)
  nuc_circ      <- ifelse(nuc_perim > 0,
                            (4 * pi * nuc_area) / (nuc_perim^2),
                            NA_real_)
  cell_hull     <- sf::st_convex_hull(cell_e)
  cell_hull_area <- as.numeric(sf::st_area(cell_hull))
  cell_solid    <- ifelse(cell_hull_area > 0,
                            cell_area / cell_hull_area, NA_real_)
  cell_cent <- sf::st_centroid(cell_e)
  nuc_cent  <- sf::st_centroid(nuc_e)
  nc_offset <- as.numeric(sf::st_distance(cell_cent, nuc_cent,
                                            by_element = TRUE))

  # PCA-based eccentricity (chunked)
  message("  computing eccentricity (PCA, chunked)...")
  cell_ecc <- ecc_chunked(cell_e)
  nuc_ecc  <- ecc_chunked(nuc_e)

  dt <- data.table(
    cohort     = cohort,
    sample_key = sample_key,
    group_id   = grp_id,
    cell_id    = cell_id,
    cell_label = lbl_e,
    x          = cells_x,
    y          = cells_y,
    cell_area, nuc_area,
    nc_ratio   = nuc_area / pmax(cell_area, 1e-6),
    cell_perimeter = cell_perim,
    nuc_perimeter  = nuc_perim,
    cell_circularity = cell_circ,
    nuc_circularity  = nuc_circ,
    cell_solidity    = cell_solid,
    cell_eccentricity = cell_ecc,
    nuc_eccentricity  = nuc_ecc,
    nc_centroid_offset = nc_offset
  )

  # Implausible-value filter
  before <- nrow(dt)
  dt <- dt[is.finite(cell_area) & is.finite(nuc_area) &
              cell_area >= 10 & cell_area <= 5000 &
              nuc_area  >= 5  & nuc_area  <= 1000 &
              nc_ratio  <= 1.0 &
              is.finite(cell_perimeter) & is.finite(nuc_perimeter)]
  after <- nrow(dt)
  message(sprintf("  filter: kept %d / %d (%.1f%%)",
                   after, before, 100 * after / max(1, before)))

  dt
}

# ---------------------------------------------------------------------------
# WT
# ---------------------------------------------------------------------------
wt_list <- list()
for (s in WT_SAMPLES) {
  sfe <- load_sfe(s)
  sfe <- override_with_06f(sfe, s)
  wt_list[[s]] <- process_sfe(sfe, "WT", s)
  rm(sfe); gc(verbose = FALSE)
}
wt_dt <- rbindlist(wt_list, fill = TRUE)
message(sprintf("\n[WT] total epithelial cells: %s",
                format(nrow(wt_dt), big.mark = ",")))

# ---------------------------------------------------------------------------
# TMA
# ---------------------------------------------------------------------------
sfe_t <- load_sfe("sfe_tma_filtered")
sfe_t <- override_with_06f(sfe_t, "sfe_tma")
tma_dt <- process_sfe(sfe_t, "TMA", "sfe_tma_filtered",
                       group_col = "patient_id")
rm(sfe_t); gc(verbose = FALSE)
message(sprintf("[TMA] total epithelial cells: %s",
                format(nrow(tma_dt), big.mark = ",")))

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
saveRDS(list(wt = wt_dt, tma = tma_dt),
        file.path(OUT_DIR, "per_cell_morphometrics.rds"))
message("Saved: ", file.path(OUT_DIR, "per_cell_morphometrics.rds"))

# ---------------------------------------------------------------------------
# Per-sample / per-patient summary
# ---------------------------------------------------------------------------
features <- c("cell_area","nuc_area","nc_ratio",
                "cell_perimeter","nuc_perimeter",
                "cell_circularity","nuc_circularity",
                "cell_solidity",
                "cell_eccentricity","nuc_eccentricity",
                "nc_centroid_offset")

per_group_summary <- function(dt) {
  dt[, c(list(n_cells = .N),
          lapply(.SD, function(x) round(median(x, na.rm = TRUE), 4))),
       by = .(cohort, sample_key, group_id, cell_label),
       .SDcols = features]
}

wt_sum  <- per_group_summary(wt_dt)
tma_sum <- per_group_summary(tma_dt)

# Wide form: cell counts per epitype per group
wt_counts <- dcast(wt_sum, sample_key ~ cell_label, value.var = "n_cells")
tma_counts <- dcast(tma_sum, group_id ~ cell_label, value.var = "n_cells")
setnames(wt_counts, EPI_LBLS, c("n_SecA","n_Int","n_SecB"),
          skip_absent = TRUE)
setnames(tma_counts, EPI_LBLS, c("n_SecA","n_Int","n_SecB"),
          skip_absent = TRUE)
wt_counts[, eligible_paired := !is.na(n_SecA)  & !is.na(n_Int) &
              !is.na(n_SecB)  & n_SecA  >= 30 & n_Int >= 30 & n_SecB >= 30]
tma_counts[, eligible_paired := !is.na(n_SecA) & !is.na(n_Int) &
              !is.na(n_SecB) & n_SecA  >= 30 & n_Int >= 30 & n_SecB >= 30]

fwrite(wt_sum,
       file.path(OUT_DIR, "per_sample_summary_wt.csv"))
fwrite(tma_sum,
       file.path(OUT_DIR, "per_patient_summary_tma.csv"))
fwrite(wt_counts,
       file.path(OUT_DIR, "wt_cell_counts.csv"))
fwrite(tma_counts,
       file.path(OUT_DIR, "tma_cell_counts.csv"))

cat("\n=== WT eligibility ===\n")
print(wt_counts)
cat(sprintf("\nWT samples eligible for paired tests: %d / %d\n",
             sum(wt_counts$eligible_paired, na.rm = TRUE),
             nrow(wt_counts)))

cat("\n=== TMA eligibility ===\n")
cat(sprintf("Total patients: %d\n", nrow(tma_counts)))
cat(sprintf("Eligible (>=30 cells in each of SecA/Trans/SecB): %d\n",
             sum(tma_counts$eligible_paired, na.rm = TRUE)))
cat("Top 10 by total cells:\n")
print(head(tma_counts[order(-rowSums(.SD, na.rm = TRUE)),
                        .(group_id, n_SecA, n_Int, n_SecB,
                          eligible_paired),
                        .SDcols = c("n_SecA","n_Int","n_SecB")], 10))

cat("\n=== Median feature values per epitype (WT pooled) ===\n")
wt_med <- wt_dt[, lapply(.SD, function(x) round(median(x, na.rm = TRUE), 3)),
                  by = cell_label, .SDcols = features]
print(wt_med)

message("\nStep 0 complete.")
