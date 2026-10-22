# ============================================================================
# 05_rebuild_existing_cohort_rerun.R — Build the whole-tissue and TMA A SFEs
#   from the raw Xenium runs.
# ============================================================================
# PURPOSE: Build one HDF5-backed SFE per sample from the raw Xenium runs into
#   RERUN_SFE_DIR (default OUTPUT_ROOT/sfe): the 8 whole tissues (6 single-tissue
#   runs; WT05 and WT01 share one slide and are split by DBSCAN on cell
#   coordinates) and TMA A (sfe_tma: slides CHUM_TMA_1 + CHUM_TMA_2, merged with
#   the second slide shifted in x). TMA A core / patient / sample-type labels are
#   copied by cell id from an earlier processed object (nearest matched cell for
#   cells without a match). Morphology images are not read (image-staging: the
#   reader cannot open the newer multi-channel morphology_focus/ folder).
#   TMA B is built by 03_build_sfe_new_coeur.R.
#
# INPUTS:
#   RAW_XENIUM_DIR/<sample>/ (workspace raw); TMA A prior object for metadata:
#   <data_dir>/previously processed/xenium_tma_final_sfe_hdf5/se.rds
# OUTPUTS: <RERUN_SFE_DIR>/{sfe_tma, sfe_<wt> x8}
# MANUSCRIPT PANEL(S): the whole-tissue and TMA A arms of Figs 4–7, Figures S11,
#   S12, S14–S18, Tables S8–S11.
# RUNTIME TIER: heavy
# ============================================================================

source("spatial/00_setup/00_setup.R")
suppressPackageStartupMessages({ library(dbscan); library(FNN) })

RAW <- Sys.getenv("RAW_XENIUM_DIR",
                  file.path(Sys.getenv("DATA_ROOT", "/path/to/deposit"), "202609_final_hgsc_spatial/data/HGSC xenium "))
sfe_out <- Sys.getenv("RERUN_SFE_DIR",
                      file.path(Sys.getenv("OUTPUT_ROOT", "./output"), "sfe"))
dir.create(sfe_out, recursive = TRUE, showWarnings = FALSE)
message("Re-run SFE output dir: ", sfe_out)

stage_without_images <- function(xdir) {
  st <- file.path(tempdir(), paste0("stage_", basename(xdir), "_", as.integer(runif(1,1e5,1e6))))
  if (dir.exists(st)) unlink(st, recursive = TRUE); dir.create(st, recursive = TRUE)
  for (f in list.files(xdir, full.names = TRUE)) {
    b <- basename(f)
    if (grepl("^morphology", b) || grepl("^\\._", b)) next
    file.symlink(normalizePath(f), file.path(st, b))
  }
  st
}
resolve_xenium <- function(sub) {
  d <- file.path(RAW, sub)
  if (file.exists(file.path(d, "cells.parquet"))) return(d)
  n <- list.dirs(d, recursive = FALSE); n <- n[file.exists(file.path(n, "cells.parquet"))]
  if (length(n) >= 1) return(n[[1]]); stop("No cells.parquet under: ", d)
}
read_staged <- function(sub, sid, segs = c("cell","nucleus")) {
  st <- stage_without_images(resolve_xenium(sub))
  s <- readXenium(st, row.names = "symbol", sample_id = sid, segmentations = segs)
  unlink(st, recursive = TRUE); s
}

# ── 1. Whole tissues (8): 6 single-tissue + WT05_WT01 dual split ──
single_wt <- c("WT02","WT03","WT04","WT06","WT07","WT08")
for (sid in single_wt) {
  out <- file.path(sfe_out, paste0("sfe_", sid))
  if (dir.exists(out)) { message(sid, ": exists, skip"); next }
  message("\n=== WT ", sid, " ===")
  sfe <- read_staged(sid, sid)
  sfe$patient_id <- NA_character_; sfe$sample_type <- NA_character_; sfe$cohort <- "whole_tissue"
  saveHDF5SummarizedExperiment(sfe, dir = out, replace = TRUE); rm(sfe); gc(FALSE)
}

# Dual tissue WT05_WT01 -> WT05 (upper) + WT01 (lower), DBSCAN split
if (!all(dir.exists(file.path(sfe_out, c("sfe_WT01","sfe_WT05"))))) {
  message("\n=== WT WT05_WT01 (DBSCAN split) ===")
  sfe <- read_staged("WT05_WT01", "WT05_WT01")
  coords <- sf::st_coordinates(colGeometry(sfe, "centroids"))
  db <- dbscan(coords, eps = 200, minPts = 50)
  # two largest clusters = the two tissues; assign by median y (upper=WT05, lower=WT01)
  tab <- sort(table(db$cluster[db$cluster != 0]), decreasing = TRUE)
  two <- as.integer(names(tab)[1:2])
  ymed <- sapply(two, function(k) median(coords[db$cluster == k, 2]))
  upper <- two[which.min(ymed)]; lower <- two[which.max(ymed)]
  assign_map <- ifelse(db$cluster == upper, "WT05",
                ifelse(db$cluster == lower, "WT01", NA))
  for (nm in c("WT01","WT05")) {
    out <- file.path(sfe_out, paste0("sfe_", nm))
    if (dir.exists(out)) next
    sub <- sfe[, which(assign_map == nm)]
    sub$sample_id <- nm; sub$patient_id <- NA_character_; sub$sample_type <- NA_character_
    sub$cohort <- "whole_tissue"
    message("  ", nm, ": ", ncol(sub), " cells")
    saveHDF5SummarizedExperiment(sub, dir = out, replace = TRUE); rm(sub); gc(FALSE)
  }
  rm(sfe); gc(FALSE)
}

# ── 2. TMA A (sfe_tma): load TMA_1 + TMA_2, transfer core/patient from the
#       prior processed object (cell-ID match + spatial-NN fallback), shift, merge ─
out_tma <- file.path(sfe_out, "sfe_tma")
if (!dir.exists(out_tma)) {
  message("\n=== CHUM/STP TMA -> sfe_tma ===")
  tma1 <- read_staged("CHUM_TMA_1", "TMA_1")
  tma2 <- read_staged("CHUM_TMA_2", "TMA_2")
  message("  TMA_1: ", ncol(tma1), "; TMA_2: ", ncol(tma2))

  old <- readRDS(file.path(data_dir, "previously processed",
                           "xenium_tma_final_sfe_hdf5", "se.rds"))
  om <- data.frame(cell_id = colnames(old), core_id = old$core_id,
                   patient_id = old$patient_id, sample_type = old$sample_type,
                   core_flag = old$core_flag, stringsAsFactors = FALSE)
  rm(old); gc(FALSE)

  transfer <- function(sfe, om, prefix = NULL) {
    ids <- if (is.null(prefix)) colnames(sfe) else paste0(prefix, colnames(sfe))
    idx <- match(ids, om$cell_id); m <- !is.na(idx)
    for (col in c("core_id","patient_id","sample_type","core_flag"))
      sfe[[col]] <- ifelse(m, om[[col]][idx], NA_character_)
    sfe$meta_source <- ifelse(m, "id_match", NA_character_)
    if (any(!m)) {
      cen <- sf::st_coordinates(colGeometry(sfe, "centroids"))
      nn <- FNN::get.knnx(cen[m,,drop=FALSE], cen[!m,,drop=FALSE], k = 1)$nn.index[,1]
      donor <- which(m)[nn]
      for (col in c("core_id","patient_id","sample_type","core_flag"))
        sfe[[col]][!m] <- sfe[[col]][donor]
      sfe$meta_source[!m] <- "spatial_nn"
    }
    sfe
  }
  tma1 <- transfer(tma1, om)
  tma2 <- transfer(tma2, om, prefix = "tma2_")
  rm(om); gc(FALSE)

  # shift TMA_2 to avoid overlap
  bb <- sf::st_bbox(colGeometry(tma1, "centroids")); xoff <- ceiling(bb["xmax"]) + 5000
  for (g in names(colGeometries(tma2))) {
    gg <- colGeometry(tma2, g); gg$geometry <- gg$geometry + c(xoff, 0); colGeometry(tma2, g) <- gg
  }
  shared <- intersect(names(colData(tma1)), names(colData(tma2)))
  colData(tma1) <- colData(tma1)[, shared]; colData(tma2) <- colData(tma2)[, shared]
  dup <- intersect(colnames(tma1), colnames(tma2))
  if (length(dup)) { ids <- colnames(tma2); ids[ids %in% dup] <- paste0(ids[ids %in% dup], "_TMA2")
    colnames(tma2) <- ids }
  tma1$cohort <- "CHUM"; tma2$cohort <- "CHUM"
  tma <- cbind(tma1, tma2); rm(tma1, tma2); gc(FALSE)
  message("  Merged sfe_tma: ", ncol(tma), " cells; cores: ",
          length(unique(na.omit(tma$core_id))), "; patients: ",
          length(unique(na.omit(tma$patient_id))))
  saveHDF5SummarizedExperiment(tma, dir = out_tma, replace = TRUE); rm(tma); gc(FALSE)
}

message("\nDone. Existing cohort rebuilt into the re-run dir.")
log_session()
