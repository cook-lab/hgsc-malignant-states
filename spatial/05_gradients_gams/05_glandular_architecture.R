# ============================================================================
# 05_glandular_architecture.R
# ----------------------------------------------------------------------------
# PURPOSE: Glandular architecture metrics (local epithelial density, distance to non-epithelial boundary) and their association with SecA->SecB polarization (WT paired; TMA validation, patients with >=30 SecA & SecB).
#
# INPUTS:
#   - sfe_tma_filtered + 8 WT SFEs (load_sfe) with cell_label, 06f override
#
# OUTPUTS:
#   - output/28_glandular_architecture/per_cell_architecture_wt.rds
#   - per_sample_medians.csv, per_patient_medians_tma.csv
#
# MANUSCRIPT PANEL(S): SF12.
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
  library(SpatialFeatureExperiment); library(SummarizedExperiment)
  library(data.table); library(sf); library(terra)
  library(BiocNeighbors); library(mgcv); library(Matrix)
})

OUT_DIR <- cfg_path("output_root", "28_glandular_architecture")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
t0  <- Sys.time()
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n")

EPI_LABELS <- c("SecA epithelium", "Intermediate epithelium",
                "SecB epithelium")
WT_SAMPLES <- sfe_names_wt

NEIGHBOR_RADII  <- c(25, 50, 100)
RAST_RES_UM     <- 10
LUMEN_MIN_AREA  <- 500     # µm²
LUMEN_MAX_AREA  <- 50000   # µm²
LUMEN_EPI_FRAC  <- 0.70

MIN_CELLS_PER_PATIENT <- 30   # TMA filter

recl <- fread(cfg_path("output_root", "06f_reclassification_polarization", "reclassified_xenium_scores.csv"),
              select = c("sample", "barcode_orig", "cell_label_06f",
                          "polarization_UCell"))
override_with_06f <- function(sfe, sample_key) {
  sub <- recl[sample == sample_key]
  if (nrow(sub) == 0) return(sfe)
  m <- match(colnames(sfe), sub$barcode_orig)
  hit <- !is.na(m)
  lab <- as.character(sfe$cell_label)
  lab[hit] <- sub$cell_label_06f[m[hit]]
  # Idempotent legacy-label rename: the 06f reclassification cache (and the
  # deposited SFE) still carry the legacy "Transitioning epithelium";
  # downstream code keys on the canonical "Intermediate epithelium"
  # (EPI_LABELS, SecA/SecB matches). Rename here so every match captures the
  # Intermediate epitype. Harmless if already canonical.
  lab[lab == "Transitioning epithelium"] <- "Intermediate epithelium"
  sfe$cell_label <- lab
  # Polarization: also attach 06f-computed values
  sfe$polarization_UCell <- sub$polarization_UCell[m]
  sfe
}

# ---------------------------------------------------------------------------
# Feature extraction helpers
# ---------------------------------------------------------------------------
compute_density <- function(coords, radii) {
  # coords: matrix [N, 2]
  # returns named list radius -> integer vector of neighbors within radius
  # (exclusive of self)
  out <- vector("list", length(radii)); names(out) <- paste0("r", radii)
  # KmknnIndex + findDistance + findNeighbors... but simpler: use RANN
  # or BiocNeighbors findNeighbors with threshold. findNeighbors doesn't
  # support radius query on all versions; fall back to KD-tree via RANN.
  if (requireNamespace("RANN", quietly = TRUE)) {
    for (r in radii) {
      # k = 512 is a hard cap; we only need COUNT of cells within r.
      # Use RANN::nn2 with large k and then count those within r.
      k_use <- min(512, nrow(coords))
      nn <- RANN::nn2(coords, query = coords, k = k_use, searchtype = "radius",
                       radius = r)
      # nn$nn.idx: 0 where no hit; subtract 1 for self
      counts <- rowSums(nn$nn.idx > 0) - 1L
      out[[paste0("r", r)]] <- pmax(counts, 0L)
    }
  } else {
    stop("RANN is required for density computation")
  }
  out
}

detect_lumens <- function(cseg_sf_all, all_cell_labels_is_epi,
                          res_um, min_area, max_area, epi_frac_thresh) {
  # cseg_sf_all: sf of ALL cellSeg polygons (both epi and non-epi)
  # all_cell_labels_is_epi: logical vector (same length / order)
  bb <- sf::st_bbox(cseg_sf_all)
  tmpl <- terra::rast(xmin = bb["xmin"], xmax = bb["xmax"],
                      ymin = bb["ymin"], ymax = bb["ymax"],
                      resolution = res_um, crs = "")
  cell_rast <- terra::rasterize(terra::vect(cseg_sf_all), tmpl,
                                field = 1, background = 0, touches = TRUE)
  space <- cell_rast == 0
  cc <- terra::patches(space, directions = 8, zeroAsNA = TRUE,
                        allowGaps = FALSE)
  cc_freq <- as.data.table(terra::freq(cc, bylayer = FALSE))
  setnames(cc_freq, c("value", "count"))
  cc_freq <- cc_freq[!is.na(value)]
  cc_freq[, area_um2 := count * (res_um ^ 2)]
  keep_ids <- cc_freq[area_um2 >= min_area & area_um2 <= max_area, value]
  if (length(keep_ids) == 0) {
    return(list(lumens_sf = sf::st_sf(lumen_id = integer(0),
                                        area_um2 = numeric(0),
                                        geometry = sf::st_sfc()),
                 raster    = cell_rast,
                 lumen_rast = NULL))
  }
  cc_sub <- terra::ifel(cc %in% keep_ids, cc, NA)
  lumens_vect <- terra::as.polygons(cc_sub, dissolve = TRUE)
  lumens_sf   <- sf::st_as_sf(lumens_vect)
  id_col <- setdiff(names(lumens_sf), "geometry")[1]
  names(lumens_sf)[names(lumens_sf) == id_col] <- "lumen_id"
  lumens_sf <- suppressWarnings(sf::st_cast(sf::st_make_valid(lumens_sf),
                                              "POLYGON", warn = FALSE))
  lumens_sf$area_um2 <- as.numeric(sf::st_area(lumens_sf))
  lumens_sf <- lumens_sf[lumens_sf$area_um2 >= min_area, ]

  # Epithelial-boundary filter
  all_centroids <- sf::st_centroid(cseg_sf_all)
  lum_buf <- sf::st_buffer(lumens_sf, 15)
  touch <- sf::st_intersects(lum_buf, all_centroids)
  epi_frac <- sapply(touch, function(ix) {
    if (length(ix) == 0) return(0)
    mean(all_cell_labels_is_epi[ix])
  })
  lumens_sf$epi_boundary_frac <- epi_frac
  lumens_kept <- lumens_sf[epi_frac >= epi_frac_thresh, ]

  list(lumens_sf = lumens_kept, raster = cell_rast, lumen_rast = cc_sub)
}

dist_to_lumen <- function(cell_centroids_mat, lumens_sf, cell_rast,
                           res_um) {
  # Rasterize lumens into a binary mask on cell_rast's grid
  lumen_mask <- terra::rasterize(terra::vect(lumens_sf), cell_rast,
                                   field = 1, background = 0)
  # Distance transform in raster units (pixels); convert to µm
  lumen_bin <- lumen_mask == 1
  dist_rast <- terra::distance(lumen_bin, target = 0)  # µm already if res=10
  # terra::distance returns dist in the units of the CRS/resolution.
  # Here CRS is "", resolution is in µm, so the output is in µm units.
  # Sample at cell centroid coords
  vals <- terra::extract(dist_rast, cell_centroids_mat)
  as.numeric(vals[, 1])
}

dist_to_nearest_nonepi_cell <- function(epi_coords, nonepi_coords) {
  # For each epithelial cell, distance to the nearest NON-epithelial cell
  # (stromal / immune / vascular). Epi cells at the tumor-stromal boundary
  # get small distances; interior tumor cells get large distances. This
  # is a cell-centroid-based definition of "distance to epithelial edge"
  # that avoids the raster artefacts of the mask-based distance
  # transform (at 5 µm/px, per-cell rasterized polygons leave 1-px gaps
  # at cell-cell contacts and every centroid ends up on an "edge pixel").
  if (nrow(nonepi_coords) == 0)
    return(rep(NA_real_, nrow(epi_coords)))
  nn <- RANN::nn2(data = nonepi_coords, query = epi_coords, k = 1)
  as.numeric(nn$nn.dists[, 1])
}

# ---------------------------------------------------------------------------
# WT pipeline: per-sample feature extraction
# ---------------------------------------------------------------------------
wt_tables <- list()
lumens_per_sample <- list()

for (s in WT_SAMPLES) {
  msg("=== WT", s)
  sfe <- load_sfe(s)
  sfe <- override_with_06f(sfe, s)

  co <- spatialCoords(sfe)
  colnames(co) <- c("x", "y")
  all_sf   <- sf::st_as_sf(colGeometry(sfe, "cellSeg"))
  all_is_epi <- sfe$cell_label %in% EPI_LABELS

  # Lumen detection on WHOLE tissue (needs non-epi polygons for the
  # epithelial-boundary fraction test)
  msg("  detecting lumens...")
  lum <- detect_lumens(all_sf, all_is_epi,
                        res_um          = RAST_RES_UM,
                        min_area        = LUMEN_MIN_AREA,
                        max_area        = LUMEN_MAX_AREA,
                        epi_frac_thresh = LUMEN_EPI_FRAC)
  saveRDS(lum$lumens_sf, file.path(OUT_DIR, sprintf("lumens_%s.rds",
                                                     sub("^sfe_", "", s))))
  lumens_per_sample[[s]] <- lum$lumens_sf
  msg(sprintf("    lumens retained: %d", nrow(lum$lumens_sf)))

  # Subset to epithelial cells for density + distance metrics
  is_epi <- all_is_epi
  epi_coords <- co[is_epi, , drop = FALSE]
  epi_labels <- as.character(sfe$cell_label)[is_epi]
  epi_pol    <- as.numeric(sfe$polarization_UCell)[is_epi]

  msg(sprintf("  density at radii %s ...",
              paste(NEIGHBOR_RADII, collapse = ",")))
  dens <- compute_density(epi_coords, NEIGHBOR_RADII)

  msg("  dist_to_lumen...")
  if (nrow(lum$lumens_sf) > 0) {
    dtl <- dist_to_lumen(epi_coords, lum$lumens_sf, lum$raster, RAST_RES_UM)
  } else {
    dtl <- rep(NA_real_, nrow(epi_coords))
  }

  msg("  dist_to_nearest_nonepi...")
  nonepi_coords <- co[!is_epi, , drop = FALSE]
  dee <- dist_to_nearest_nonepi_cell(epi_coords, nonepi_coords)

  dt <- data.table(
    sample_id   = sub("^sfe_", "", s),
    cell_id     = colnames(sfe)[is_epi],
    cell_label  = epi_labels,
    polarization = epi_pol,
    x = epi_coords[, 1], y = epi_coords[, 2],
    dist_to_lumen_um       = dtl,
    dist_to_nonepi_cell_um = dee,
    epi_neighbors_25um  = dens[["r25"]],
    epi_neighbors_50um  = dens[["r50"]],
    epi_neighbors_100um = dens[["r100"]])
  wt_tables[[s]] <- dt

  rm(sfe, all_sf, lum, dens, dtl, dee); gc(verbose = FALSE)
  msg(sprintf("  done. cells retained: %d", nrow(dt)))
}
wt_cells <- rbindlist(wt_tables)
saveRDS(wt_cells, file.path(OUT_DIR, "per_cell_architecture_wt.rds"))
msg(sprintf("WT combined: %s cells", format(nrow(wt_cells), big.mark = ",")))

# ---------------------------------------------------------------------------
# TMA pipeline: density only
# ---------------------------------------------------------------------------
msg("=== TMA sfe_tma_filtered (density only)")
sfe_t <- load_sfe("sfe_tma_filtered")
sfe_t <- override_with_06f(sfe_t, "sfe_tma")
is_epi <- sfe_t$cell_label %in% EPI_LABELS &
          !is.na(sfe_t$patient_id) & sfe_t$patient_id != ""
co_t <- spatialCoords(sfe_t)[is_epi, , drop = FALSE]
colnames(co_t) <- c("x", "y")

# Density computed PER CORE to avoid cross-core false-neighbor inflation
core_id <- as.character(sfe_t$core_id)[is_epi]
dens_t  <- vector("list", length(NEIGHBOR_RADII))
names(dens_t) <- paste0("r", NEIGHBOR_RADII)
for (r in NEIGHBOR_RADII)
  dens_t[[paste0("r", r)]] <- integer(length(core_id))

msg("  per-core density...")
cores <- split(seq_along(core_id), core_id)
for (i in seq_along(cores)) {
  idx <- cores[[i]]
  if (length(idx) < 2) next
  d <- compute_density(co_t[idx, , drop = FALSE], NEIGHBOR_RADII)
  for (r in NEIGHBOR_RADII) dens_t[[paste0("r", r)]][idx] <- d[[paste0("r", r)]]
  if (i %% 25 == 0) msg(sprintf("    core %d / %d", i, length(cores)))
}

tma_cells <- data.table(
  patient_id  = as.character(sfe_t$patient_id)[is_epi],
  core_id     = core_id,
  cell_label  = as.character(sfe_t$cell_label)[is_epi],
  polarization = as.numeric(sfe_t$polarization_UCell)[is_epi],
  x = co_t[, 1], y = co_t[, 2],
  epi_neighbors_25um  = dens_t[["r25"]],
  epi_neighbors_50um  = dens_t[["r50"]],
  epi_neighbors_100um = dens_t[["r100"]])
saveRDS(tma_cells, file.path(OUT_DIR, "per_cell_architecture_tma.rds"))
msg(sprintf("TMA: %s cells across %d patients",
            format(nrow(tma_cells), big.mark = ","),
            length(unique(tma_cells$patient_id))))

rm(sfe_t, co_t); gc(verbose = FALSE)

# ---------------------------------------------------------------------------
# Per-sample / per-patient SecA vs SecB medians
# ---------------------------------------------------------------------------
msg("=== Per-sample / per-patient aggregation")
FEATS_WT  <- c("dist_to_lumen_um", "dist_to_nonepi_cell_um",
                "epi_neighbors_25um", "epi_neighbors_50um",
                "epi_neighbors_100um")
FEATS_TMA <- c("epi_neighbors_25um", "epi_neighbors_50um",
                "epi_neighbors_100um")

wt_med <- wt_cells[
  cell_label %in% c("SecA epithelium", "SecB epithelium"),
  c(list(n_cells = .N), lapply(.SD, median, na.rm = TRUE)),
  by = .(sample_id, cell_label), .SDcols = FEATS_WT]
fwrite(wt_med, file.path(OUT_DIR, "per_sample_medians.csv"))

tma_med <- tma_cells[
  cell_label %in% c("SecA epithelium", "SecB epithelium"),
  c(list(n_cells = .N), lapply(.SD, median, na.rm = TRUE)),
  by = .(patient_id, cell_label), .SDcols = FEATS_TMA]
tma_n <- dcast(tma_med, patient_id ~ cell_label, value.var = "n_cells",
                fill = 0)
setnames(tma_n, c("SecA epithelium", "SecB epithelium"),
         c("n_SecA", "n_SecB"), skip_absent = TRUE)
keep_pt <- tma_n[n_SecA >= MIN_CELLS_PER_PATIENT &
                  n_SecB >= MIN_CELLS_PER_PATIENT, patient_id]
tma_med <- tma_med[patient_id %in% keep_pt]
fwrite(tma_med,
       file.path(OUT_DIR, "per_patient_medians_tma.csv"))
msg(sprintf("  TMA patients passing >=30 SecA & >=30 SecB: %d",
            length(keep_pt)))

# ---------------------------------------------------------------------------
# Direction counts (paired Wilcoxon, two-sided; report n_fall / n_rise)
# ---------------------------------------------------------------------------
direction_row <- function(med_tbl, id_col, feat) {
  dw <- dcast(med_tbl[, c(id_col, "cell_label", feat), with = FALSE],
              as.formula(paste(id_col, "~ cell_label")),
              value.var = feat)
  setnames(dw, c("SecA epithelium", "SecB epithelium"), c("secA", "secB"),
           skip_absent = TRUE)
  dw <- dw[!is.na(secA) & !is.na(secB)]
  n <- nrow(dw)
  p <- tryCatch(wilcox.test(dw$secA, dw$secB, paired = TRUE)$p.value,
                 error = function(e) NA_real_)
  n_rise <- sum(dw$secB > dw$secA)
  n_fall <- sum(dw$secB < dw$secA)
  data.table(feature = feat, n = n, n_rise = n_rise, n_fall = n_fall,
             p_wilcox = p,
             median_secA = median(dw$secA, na.rm = TRUE),
             median_secB = median(dw$secB, na.rm = TRUE))
}
dir_wt  <- rbindlist(lapply(FEATS_WT,  direction_row,
                              med_tbl = wt_med,  id_col = "sample_id"))
dir_wt[, cohort := "WT"]
dir_tma <- rbindlist(lapply(FEATS_TMA, direction_row,
                              med_tbl = tma_med, id_col = "patient_id"))
dir_tma[, cohort := "TMA"]
dir_all <- rbind(dir_wt, dir_tma)
fwrite(dir_all, file.path(OUT_DIR, "per_sample_direction_counts.csv"))
msg("Direction counts:")
print(dir_all, digits = 3)

# ---------------------------------------------------------------------------
# Per-sample + pooled GAMs (polarization ~ s(feature))
# ---------------------------------------------------------------------------
msg("=== GAM fitting")
K_POOL <- 10; K_SAMP <- 5
N_CAP  <- 5000

fit_gams <- function(cells, feat, group_col) {
  groups <- unique(cells[[group_col]])
  per <- lapply(groups, function(g) {
    sub <- cells[get(group_col) == g &
                  !is.na(get(feat)) & !is.na(polarization)]
    if (nrow(sub) < 100 || sd(sub[[feat]]) < 1e-6) return(NULL)
    if (nrow(sub) > N_CAP) sub <- sub[sample(.N, N_CAP)]
    setnames(sub, feat, "x_feat")
    # Degenerate features (too few unique values) will crash the smoother
    if (length(unique(sub$x_feat)) < K_SAMP + 1) return(NULL)
    tryCatch(
      gam(polarization ~ s(x_feat, k = K_SAMP), data = sub, method = "REML"),
      error = function(e) NULL)
  })
  names(per) <- groups

  # Pooled fit (up to N_CAP per sample)
  pooled_sub <- cells[!is.na(get(feat)) & !is.na(polarization),
                       .SD[if (.N > N_CAP) sample(.N, N_CAP) else .I],
                       by = group_col]
  setnames(pooled_sub, feat, "x_feat")
  if (length(unique(pooled_sub$x_feat)) < K_POOL + 1) {
    pooled <- NULL
  } else {
    pooled <- tryCatch(
      gam(polarization ~ s(x_feat, k = K_POOL),
          data = pooled_sub, method = "REML"),
      error = function(e) NULL)
  }
  list(pooled = pooled, per = per)
}

gam_wt  <- lapply(FEATS_WT,  fit_gams, cells = wt_cells,  group_col = "sample_id")
names(gam_wt)  <- FEATS_WT
gam_tma <- lapply(FEATS_TMA, fit_gams, cells = tma_cells[
                   patient_id %in% keep_pt],
                   group_col = "patient_id")
names(gam_tma) <- FEATS_TMA
saveRDS(list(wt = gam_wt, tma = gam_tma),
        file.path(OUT_DIR, "gam_fits.rds"))
msg("  saved gam_fits.rds")

msg(sprintf("=== Phase 28 feature extraction done. Total elapsed: %s",
            format(round(difftime(Sys.time(), t0, units = "mins"), 1))))
