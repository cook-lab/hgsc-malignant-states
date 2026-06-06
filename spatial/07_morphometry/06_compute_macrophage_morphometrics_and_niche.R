# ============================================================================
# 06_compute_macrophage_morphometrics_and_niche.R
# ----------------------------------------------------------------------------
# PURPOSE: Compute per-macrophage morphometrics + 50um epithelial-neighbor niche class (SecA_dominant / SecB_dominant via RANN).
#
# INPUTS:
#   - SFEs (load_sfe) with cell_label, geometries; 06f override
#
# OUTPUTS:
#   - output/34_macrophage_morphometrics/per_cell_macrophage_morphometrics.rds
#
# MANUSCRIPT PANEL(S): Fig 6I/6J (upstream)
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
  library(data.table); library(sf); library(dbscan)
  library(SpatialFeatureExperiment); library(SummarizedExperiment)
})

OUT_DIR <- file.path(out_dir, "34_macrophage_morphometrics")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE,
            showWarnings = FALSE)

WT_SAMPLES <- c("sfe_OTB_2384", "sfe_OTB_2417", "sfe_OTB_2432",
                "sfe_OTB_2454", "sfe_OTB_2457", "sfe_OTB_2461",
                "sfe_SP24_24824", "sfe_SP24_25573")
EPI_LBLS <- c("SecA epithelium", "Intermediate epithelium", "SecB epithelium")
RADIUS_UM <- 50
DOM_RATIO <- 3   # 3:1 dominance threshold
DOM_MIN_N <- 5   # minimum count of dominant epitype to assign
CHUNK     <- 25000

# ---------------------------------------------------------------------------
# 06f override
# ---------------------------------------------------------------------------
f06f <- file.path(out_dir, "06f_reclassification_polarization",
                  "reclassified_xenium_scores.csv")
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
# Eccentricity helper (PCA on polygon coords; chunked for memory)
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
  if (n == 0) return(numeric(0))
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
# Per-SFE processor
# ---------------------------------------------------------------------------
process_sfe <- function(sfe, cohort, sample_key, group_col = NULL) {
  message(sprintf("[%s] starting (n_cells=%d)", sample_key, ncol(sfe)))
  lbl    <- as.character(sfe$cell_label)
  is_mac <- lbl == "Macrophage"
  is_epi <- lbl %in% EPI_LBLS
  if (sum(is_mac) == 0 || sum(is_epi) == 0) {
    message("  no macs or no epi; skipping"); return(NULL)
  }
  message(sprintf("  macrophages: %d  epithelial: %d",
                   sum(is_mac), sum(is_epi)))

  cell_sf <- tryCatch(sf::st_as_sf(colGeometry(sfe, "cellSeg")),
                       error = function(e) NULL)
  nuc_sf  <- tryCatch(sf::st_as_sf(colGeometry(sfe, "nucSeg")),
                       error = function(e) NULL)
  if (is.null(cell_sf) || is.null(nuc_sf)) {
    message("  segmentations unavailable"); return(NULL)
  }

  # Polygons may cover only the first N cells in colData (verified for TMA)
  n_with_poly <- nrow(cell_sf)
  if (n_with_poly < ncol(sfe)) {
    message(sprintf("  WARNING: %.1f%% of cells have polygons (%d / %d)",
                     100 * n_with_poly / ncol(sfe),
                     n_with_poly, ncol(sfe)))
  }
  # Restrict everything to the polygon-covered set
  is_mac_pl <- is_mac & seq_along(is_mac) <= n_with_poly
  is_epi_pl <- is_epi & seq_along(is_epi) <= n_with_poly
  message(sprintf("  with polygons: macs %d  epi %d",
                   sum(is_mac_pl), sum(is_epi_pl)))

  # Subset polygons + labels
  mac_sf <- cell_sf[is_mac_pl[seq_len(n_with_poly)], ]
  mac_nu <- nuc_sf[is_mac_pl[seq_len(n_with_poly)], ]
  mac_id <- colnames(sfe)[is_mac_pl]
  mac_x  <- spatialCoords(sfe)[is_mac_pl, 1]
  mac_y  <- spatialCoords(sfe)[is_mac_pl, 2]
  if (is.null(group_col)) {
    mac_grp <- rep(sample_key, length(mac_id))
  } else {
    mac_grp <- as.character(colData(sfe)[[group_col]][is_mac_pl])
  }

  # Validate / repair polygons
  mac_sf <- sf::st_make_valid(mac_sf)
  mac_nu <- sf::st_make_valid(mac_nu)
  v_cell <- sf::st_is_valid(mac_sf); v_nuc <- sf::st_is_valid(mac_nu)
  keep   <- v_cell & v_nuc
  message(sprintf("  valid polygons: cell=%.1f%% nuc=%.1f%% both=%.1f%%",
                   100 * mean(v_cell), 100 * mean(v_nuc), 100 * mean(keep)))
  if (sum(keep) == 0) return(NULL)
  mac_sf <- mac_sf[keep, ]; mac_nu <- mac_nu[keep, ]
  mac_id <- mac_id[keep]; mac_x <- mac_x[keep]; mac_y <- mac_y[keep]
  mac_grp <- mac_grp[keep]

  # --- 9 morphometric features ---
  message("  computing morphometric features...")
  cell_area  <- as.numeric(sf::st_area(mac_sf))
  nuc_area   <- as.numeric(sf::st_area(mac_nu))
  cell_per   <- as.numeric(sf::st_length(sf::st_boundary(mac_sf)))
  nuc_per    <- as.numeric(sf::st_length(sf::st_boundary(mac_nu)))
  cell_circ  <- ifelse(cell_per > 0, (4*pi*cell_area)/(cell_per^2), NA_real_)
  nuc_circ   <- ifelse(nuc_per  > 0, (4*pi*nuc_area )/(nuc_per^2 ), NA_real_)
  hull_area  <- as.numeric(sf::st_area(sf::st_convex_hull(mac_sf)))
  cell_solid <- ifelse(hull_area > 0, cell_area / hull_area, NA_real_)
  cc <- sf::st_centroid(mac_sf); nc <- sf::st_centroid(mac_nu)
  nc_off     <- as.numeric(sf::st_distance(cc, nc, by_element = TRUE))
  message("  computing eccentricity (chunked)...")
  cell_ecc   <- ecc_chunked(mac_sf)
  nuc_ecc    <- ecc_chunked(mac_nu)

  # --- Niche scoring ---
  message("  niche scoring (KNN within 50 µm)...")
  epi_coords <- spatialCoords(sfe)[is_epi_pl, , drop = FALSE]
  epi_label  <- lbl[is_epi_pl]
  epi_pol    <- as.numeric(sfe$polarization_UCell[is_epi_pl])
  mac_xy     <- cbind(mac_x, mac_y)

  nn <- dbscan::frNN(x = epi_coords, eps = RADIUS_UM, query = mac_xy,
                      sort = FALSE)
  n_macs <- nrow(mac_xy)
  n_secA <- integer(n_macs); n_trans <- integer(n_macs); n_secB <- integer(n_macs)
  pol_mean <- rep(NA_real_, n_macs)
  for (i in seq_len(n_macs)) {
    idx <- nn$id[[i]]
    if (length(idx) == 0) next
    nbr_lbl <- epi_label[idx]
    n_secA[i]  <- sum(nbr_lbl == "SecA epithelium")
    n_trans[i] <- sum(nbr_lbl == "Intermediate epithelium")
    n_secB[i]  <- sum(nbr_lbl == "SecB epithelium")
    pol_mean[i] <- mean(epi_pol[idx], na.rm = TRUE)
  }
  total_epi <- n_secA + n_trans + n_secB

  niche_class <- rep("mixed", n_macs)
  niche_class[n_secA >= DOM_MIN_N & n_secA > DOM_RATIO * n_secB] <- "SecA_dominant"
  niche_class[n_secB >= DOM_MIN_N & n_secB > DOM_RATIO * n_secA] <- "SecB_dominant"
  niche_class[total_epi == 0] <- "no_epi_neighbors"

  dt <- data.table(
    cohort     = cohort,
    sample_key = sample_key,
    group_id   = mac_grp,
    cell_id    = mac_id,
    x          = mac_x, y = mac_y,
    cell_area, nuc_area,
    nc_ratio   = nuc_area / pmax(cell_area, 1e-6),
    cell_perimeter   = cell_per,
    nuc_perimeter    = nuc_per,
    cell_circularity = cell_circ,
    nuc_circularity  = nuc_circ,
    cell_solidity    = cell_solid,
    cell_eccentricity = cell_ecc,
    nuc_eccentricity  = nuc_ecc,
    nc_centroid_offset = nc_off,
    n_SecA_neighbors = n_secA,
    n_Int_neighbors = n_trans,
    n_SecB_neighbors = n_secB,
    n_total_epi_neighbors = total_epi,
    niche_polarization_mean = pol_mean,
    niche_class = niche_class
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
message(sprintf("\n[WT] total macs after filter: %s",
                format(nrow(wt_dt), big.mark = ",")))

# ---------------------------------------------------------------------------
# TMA
# ---------------------------------------------------------------------------
sfe_t <- load_sfe("sfe_tma_filtered")
sfe_t <- override_with_06f(sfe_t, "sfe_tma")
tma_dt <- process_sfe(sfe_t, "TMA", "sfe_tma_filtered",
                       group_col = "patient_id")
rm(sfe_t); gc(verbose = FALSE)
message(sprintf("[TMA] total macs after filter: %s",
                format(nrow(tma_dt), big.mark = ",")))

# ---------------------------------------------------------------------------
# Save cache
# ---------------------------------------------------------------------------
saveRDS(list(wt = wt_dt, tma = tma_dt),
        file.path(OUT_DIR, "per_cell_macrophage_morphometrics.rds"))

# ---------------------------------------------------------------------------
# Eligibility tables
# ---------------------------------------------------------------------------
elig_table <- function(dt, group_col) {
  by_cols <- group_col
  dt[, .(n_total          = .N,
          n_SecA_dominant  = sum(niche_class == "SecA_dominant"),
          n_SecB_dominant  = sum(niche_class == "SecB_dominant"),
          n_mixed          = sum(niche_class == "mixed"),
          n_no_epi         = sum(niche_class == "no_epi_neighbors")),
       by = by_cols][, eligible_paired :=
                       n_SecA_dominant >= 30 & n_SecB_dominant >= 30][]
}
wt_elig  <- elig_table(wt_dt, "sample_key")
tma_elig <- elig_table(tma_dt, "group_id")

fwrite(wt_elig,  file.path(OUT_DIR, "wt_eligibility.csv"))
fwrite(tma_elig, file.path(OUT_DIR, "tma_eligibility.csv"))

cat("\n=== WT ELIGIBILITY ===\n")
print(wt_elig)
cat(sprintf("\nWT samples eligible (>=30 in each niche class): %d / %d\n",
             sum(wt_elig$eligible_paired, na.rm = TRUE), nrow(wt_elig)))

cat("\n=== TMA ELIGIBILITY ===\n")
cat(sprintf("Total patients: %d\n", nrow(tma_elig)))
cat(sprintf("Eligible (>=30 in each niche class): %d\n",
             sum(tma_elig$eligible_paired, na.rm = TRUE)))
cat("Top 10 by total macs:\n")
print(head(tma_elig[order(-n_total)], 10))

message("\nStep 0 complete.")
