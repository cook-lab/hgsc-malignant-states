# ============================================================================
# 01_build_sfe_tma.R — Build merged TMA SFE from raw Xenium data
# ============================================================================
# PURPOSE: Load fresh TMA_1 + TMA_2 Xenium slides, transfer core/patient
#   metadata from the previously processed SFE (cell-ID match with a spatial
#   nearest-neighbor fallback), shift TMA_2 coordinates to avoid overlap, and
#   merge into a single TMA SFE. All cells kept; only raw counts stored
#   (normalization happens in 02_qc/04_qc_filter_normalize.R).
#
# INPUTS:
#   - <data_root>/2026_final_xenium_analysis/data/xenium/tma/HGSC_STP_TMA_{1,2}
#   - <data_root>/2026_final_xenium_analysis/data/previously processed/
#       xenium_tma_final_sfe_hdf5/se.rds   (core/patient metadata source)
#
# OUTPUTS:
#   - <sfe_dir>/sfe_tma   (HDF5-backed merged TMA SFE, raw counts)
#
# MANUSCRIPT PANEL(S): backend for all TMA panels (Fig 4–7, SF10–SF14);
#   upstream of sfe_tma_filtered (the canonical TMA entry-point).
#
# RUNTIME TIER: heavy (loads two full Xenium slides; ~671K cells)
# ============================================================================

source("spatial/00_setup/00_setup.R")
library(FNN)

# ── 1. Load fresh TMA SFEs ──────────────────────────────────────────────────

message("Loading fresh TMA_1...")
tma1 <- readXenium(file.path(data_dir, "xenium/tma/HGSC_STP_TMA_1"),
                   row.names = "symbol", sample_id = "TMA_1")
message("  ", ncol(tma1), " cells, ", nrow(tma1), " genes")

message("Loading fresh TMA_2...")
tma2 <- readXenium(file.path(data_dir, "xenium/tma/HGSC_STP_TMA_2"),
                   row.names = "symbol", sample_id = "TMA_2")
message("  ", ncol(tma2), " cells, ", nrow(tma2), " genes")

# ── 2. Load old processed SFE for metadata transfer ─────────────────────────

message("Loading old processed TMA SFE for metadata...")
old_sfe <- readRDS(file.path(data_dir, "previously processed",
                              "xenium_tma_final_sfe_hdf5", "se.rds"))
message("  ", ncol(old_sfe), " cells in old SFE")

old_meta <- data.frame(
  cell_id      = colnames(old_sfe),
  core_id      = old_sfe$core_id,
  patient_id   = old_sfe$patient_id,
  sample_type  = old_sfe$sample_type,
  core_flag    = old_sfe$core_flag,
  stringsAsFactors = FALSE
)

rm(old_sfe); gc(verbose = FALSE)
message("  Extracted metadata for ", nrow(old_meta), " cells")

# ── 3. Helper: transfer metadata with spatial NN fallback ────────────────────

#' For each fresh SFE, first match cell IDs directly to old metadata.
#' Unmatched cells (likely QC-filtered previously but still physically on a
#' core) get assigned the nearest matched neighbor's core spatially. A
#' `meta_source` column tracks how each cell was assigned:
#'   "id_match"   = direct cell ID match to old SFE
#'   "spatial_nn" = assigned via nearest spatial neighbor

transfer_metadata <- function(sfe, old_meta, label) {

  message("Matching ", label, " cell IDs to old SFE...")
  idx <- match(colnames(sfe), old_meta$cell_id)
  matched  <- sum(!is.na(idx))
  unmatched <- sum(is.na(idx))
  message("  Direct match: ", matched, " / ", ncol(sfe),
          " (", round(matched / ncol(sfe) * 100, 1), "%)")
  message("  Unmatched: ", unmatched)

  sfe$core_id     <- ifelse(!is.na(idx), old_meta$core_id[idx],     NA_character_)
  sfe$patient_id  <- ifelse(!is.na(idx), old_meta$patient_id[idx],  NA_character_)
  sfe$sample_type <- ifelse(!is.na(idx), old_meta$sample_type[idx], NA_character_)
  sfe$core_flag   <- ifelse(!is.na(idx), old_meta$core_flag[idx],   NA_character_)
  sfe$meta_source <- ifelse(!is.na(idx), "id_match", NA_character_)

  if (unmatched > 0) {
    message("  Assigning ", unmatched, " unmatched cells via spatial nearest neighbor...")

    centroids <- sf::st_coordinates(colGeometry(sfe, "centroids"))

    matched_mask   <- !is.na(idx)
    unmatched_mask <- is.na(idx)

    coords_matched   <- centroids[matched_mask, , drop = FALSE]
    coords_unmatched <- centroids[unmatched_mask, , drop = FALSE]

    nn <- FNN::get.knnx(data = coords_matched, query = coords_unmatched, k = 1)
    nn_idx <- nn$nn.index[, 1]
    nn_dist <- nn$nn.dist[, 1]

    matched_positions <- which(matched_mask)
    donor_positions <- matched_positions[nn_idx]

    sfe$core_id[unmatched_mask]     <- sfe$core_id[donor_positions]
    sfe$patient_id[unmatched_mask]  <- sfe$patient_id[donor_positions]
    sfe$sample_type[unmatched_mask] <- sfe$sample_type[donor_positions]
    sfe$core_flag[unmatched_mask]   <- sfe$core_flag[donor_positions]
    sfe$meta_source[unmatched_mask] <- "spatial_nn"

    message("  Spatial NN distances — median: ", round(median(nn_dist), 1),
            " um, max: ", round(max(nn_dist), 1), " um")
  }

  sfe
}

# ── 4. Transfer metadata to fresh TMA_1 and TMA_2 ───────────────────────────

# TMA_1: cell IDs match directly
tma1 <- transfer_metadata(tma1, old_meta, "TMA_1")

# TMA_2: old SFE used "tma2_" prefix on TMA_2 cell IDs
tma2_orig_ids <- colnames(tma2)
colnames(tma2) <- paste0("tma2_", colnames(tma2))
tma2 <- transfer_metadata(tma2, old_meta, "TMA_2")
colnames(tma2) <- tma2_orig_ids
rm(tma2_orig_ids)

rm(old_meta); gc(verbose = FALSE)

# ── 5. Shift TMA_2 coordinates to avoid overlap ─────────────────────────────

tma1_bbox <- sf::st_bbox(colGeometry(tma1, "centroids"))
x_offset <- ceiling(tma1_bbox["xmax"]) + 5000  # 5000 um gap
rm(tma1_bbox)

message("Shifting TMA_2 x-coordinates by +", x_offset, " um")

for (geom_name in names(colGeometries(tma2))) {
  g <- colGeometry(tma2, geom_name)
  g$geometry <- g$geometry + c(x_offset, 0)
  colGeometry(tma2, geom_name) <- g
}

# ── 6. Ensure colData columns match before merge ────────────────────────────

shared_cols <- intersect(names(colData(tma1)), names(colData(tma2)))
colData(tma1) <- colData(tma1)[, shared_cols]
colData(tma2) <- colData(tma2)[, shared_cols]

message("Shared colData columns: ", paste(shared_cols, collapse = ", "))

# ── 7. Make column names unique before merge ─────────────────────────────────

# Cell IDs are unique within each slide but 33 overlap between slides.
dupes <- intersect(colnames(tma1), colnames(tma2))
message("Overlapping cell IDs between slides: ", length(dupes))
if (length(dupes) > 0) {
  tma2_ids <- colnames(tma2)
  dupe_mask <- tma2_ids %in% dupes
  tma2_ids[dupe_mask] <- paste0(tma2_ids[dupe_mask], "_TMA2")
  colnames(tma2) <- tma2_ids
  for (geom_name in names(colGeometries(tma2))) {
    g <- colGeometry(tma2, geom_name)
    rn <- rownames(g)
    rn_mask <- rn %in% dupes
    rn[rn_mask] <- paste0(rn[rn_mask], "_TMA2")
    rownames(g) <- rn
    colGeometry(tma2, geom_name) <- g
  }
  message("  Suffixed ", sum(dupe_mask), " TMA_2 cell IDs with '_TMA2'")
}

# ── 8. Merge ─────────────────────────────────────────────────────────────────

message("Merging TMA_1 + TMA_2...")

shared_geoms <- intersect(names(colGeometries(tma1)), names(colGeometries(tma2)))
message("Shared geometries: ", paste(shared_geoms, collapse = ", "))

tma <- cbind(tma1, tma2)
message("Merged TMA SFE: ", ncol(tma), " cells, ", nrow(tma), " genes")

rm(tma1, tma2); gc(verbose = FALSE)

# ── 9. Summary diagnostics ──────────────────────────────────────────────────

message("\n=== TMA SFE Summary ===")
message("Total cells: ", ncol(tma))
message("Genes: ", nrow(tma))
message("Assays: ", paste(assayNames(tma), collapse = ", "))
message("colGeometries: ", paste(names(colGeometries(tma)), collapse = ", "))

message("\nCells by sample_id (slide):")
print(table(tma$sample_id))

message("\nCells by sample_type:")
print(table(tma$sample_type, useNA = "ifany"))

message("\nMetadata assignment method:")
print(table(tma$meta_source))

message("\ncore_flag breakdown:")
print(table(tma$core_flag, useNA = "ifany"))

message("\nUnique cores: ", length(unique(tma$core_id[tma$core_id != "Off core"])))
message("Unique patients: ", length(unique(tma$patient_id[!is.na(tma$patient_id)])))

# ── 10. Save ─────────────────────────────────────────────────────────────────

out_path <- file.path(sfe_dir, "sfe_tma")
message("\nSaving to: ", out_path)
saveHDF5SummarizedExperiment(tma, dir = out_path, replace = TRUE)
message("Saved HDF5-backed SFE.")

message("\nDone. TMA SFE build complete.")
log_session()
