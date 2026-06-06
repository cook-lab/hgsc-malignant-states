# ============================================================================
# 02_build_sfe_whole_tissue.R — Build individual whole-tissue SFEs
# ============================================================================
# PURPOSE: Load each published whole-tissue Xenium sample and attach sample_id.
#   OTB_2457_2384 contains two tissues on one slide — split via DBSCAN into
#   OTB_2457 (upper) + OTB_2384 (lower). All cells kept; only raw counts stored
#   (normalization happens in 02_qc/04_qc_filter_normalize.R).
#
# COHORT PIN: the published whole-tissue arm is exactly the 8 samples in
#   CFG$cohort$whole_tissue. The two FTE whole-tissue samples
#   (CFG$cohort$fte_exclude_wt = FT1-1, EAOC-1-FTE) are EXCLUDED here — they
#   were added post-preprint and caused cohort drift. FTE TMA cores stay in the
#   TMA (built by 01_build_sfe_tma.R). See docs/REPRODUCIBILITY.md.
#
# INPUTS:
#   - <data_root>/2026_final_xenium_analysis/data/xenium/whole_tissue/<sample>
#       (raw single-tissue samples + the dual-tissue OTB_2457_2384)
#
# OUTPUTS:
#   - <sfe_dir>/sfe_<sample>   for each of the 8 published whole tissues
#
# MANUSCRIPT PANEL(S): backend for whole-tissue panels (Fig 4C/E/F/G/H/I/K,
#   Fig 5, Fig 6, SF10A, SF12–SF14).
#
# RUNTIME TIER: heavy (loads multiple full whole-tissue Xenium runs)
# ============================================================================

source("spatial/00_setup/00_setup.R")
library(dbscan)

# ── 1. Define samples (published whole-tissue cohort only) ────────────────────

# Published whole-tissue arm (the 8 in config). OTB_2457 + OTB_2384 are derived
# from the dual-tissue slide below; the remaining 6 are single-tissue raw runs.
wt_cohort <- CFG$cohort$whole_tissue

# Dual-tissue sample that splits into OTB_2457 + OTB_2384.
dual_tissue       <- "OTB_2457_2384"
dual_components   <- c("OTB_2457", "OTB_2384")

# Single-tissue raw samples = published cohort minus the dual-derived pair.
single_tissue <- setdiff(wt_cohort, dual_components)

message("Published whole-tissue cohort (", length(wt_cohort), "): ",
        paste(wt_cohort, collapse = ", "))
message("Excluded from whole-tissue arm (cohort PIN): ",
        paste(CFG$cohort$fte_exclude_wt, collapse = ", "))

# ── 2. Process single-tissue samples ─────────────────────────────────────────

for (sid in single_tissue) {
  message("\n=== ", sid, " ===")

  out_path <- file.path(sfe_dir, paste0("sfe_", sid))
  if (dir.exists(out_path)) {
    message("  Already built: ", out_path, " — skipping.")
    next
  }

  sfe <- readXenium(file.path(data_dir, "xenium/whole_tissue", sid),
                    row.names = "symbol", sample_id = sid)
  message("  Loaded: ", ncol(sfe), " cells, ", nrow(sfe), " genes")

  # Whole tissue: no patient_id; sample_type NA (no clinical metadata)
  sfe$patient_id  <- NA_character_
  sfe$sample_type <- NA_character_

  saveHDF5SummarizedExperiment(sfe, dir = out_path, replace = TRUE)
  message("  Saved: ", out_path)

  rm(sfe); gc(verbose = FALSE)
}

# ── 3. Process dual-tissue sample (OTB_2457_2384) ────────────────────────────

dual_already_built <- dir.exists(file.path(sfe_dir, "sfe_OTB_2457")) &&
                      dir.exists(file.path(sfe_dir, "sfe_OTB_2384"))

if (!dual_already_built) {

message("\n=== ", dual_tissue, " (DBSCAN split) ===")

sfe <- readXenium(file.path(data_dir, "xenium/whole_tissue", dual_tissue),
                  row.names = "symbol", sample_id = dual_tissue)
message("  Loaded: ", ncol(sfe), " cells, ", nrow(sfe), " genes")

coords <- sf::st_coordinates(colGeometry(sfe, "centroids"))

# DBSCAN: eps=200 um bridges within-tissue gaps, separates between-tissue gap
db <- dbscan(coords, eps = 200, minPts = 50)
message("  DBSCAN clusters: ", paste(names(table(db$cluster)),
                                      table(db$cluster), sep = "=", collapse = ", "))

main_clusters <- sort(setdiff(unique(db$cluster), 0))
if (length(main_clusters) != 2) {
  stop("Expected 2 tissue clusters, found ", length(main_clusters),
       ". Adjust DBSCAN parameters.")
}

# Higher y-centroid median = upper tissue (OTB_2457); lower = OTB_2384
median_y <- sapply(main_clusters, function(cl) median(coords[db$cluster == cl, 2]))
upper_cluster <- main_clusters[which.max(median_y)]
lower_cluster <- main_clusters[which.min(median_y)]

message("  Upper tissue (OTB_2457): cluster ", upper_cluster, " (",
        sum(db$cluster == upper_cluster), " cells)")
message("  Lower tissue (OTB_2384): cluster ", lower_cluster, " (",
        sum(db$cluster == lower_cluster), " cells)")
message("  Noise (cluster 0): ", sum(db$cluster == 0), " cells — discarded")

sfe$dbscan_cluster <- db$cluster

mask_2457 <- db$cluster == upper_cluster
sfe_2457 <- sfe[, mask_2457]
sfe_2457$sample_id   <- "OTB_2457"
sfe_2457$patient_id  <- NA_character_
sfe_2457$sample_type <- NA_character_

mask_2384 <- db$cluster == lower_cluster
sfe_2384 <- sfe[, mask_2384]
sfe_2384$sample_id   <- "OTB_2384"
sfe_2384$patient_id  <- NA_character_
sfe_2384$sample_type <- NA_character_

rm(sfe, coords, db); gc(verbose = FALSE)

out_2457 <- file.path(sfe_dir, "sfe_OTB_2457")
saveHDF5SummarizedExperiment(sfe_2457, dir = out_2457, replace = TRUE)
message("  Saved OTB_2457: ", ncol(sfe_2457), " cells -> ", out_2457)

out_2384 <- file.path(sfe_dir, "sfe_OTB_2384")
saveHDF5SummarizedExperiment(sfe_2384, dir = out_2384, replace = TRUE)
message("  Saved OTB_2384: ", ncol(sfe_2384), " cells -> ", out_2384)

rm(sfe_2457, sfe_2384); gc(verbose = FALSE)

} else {
  message("\n=== ", dual_tissue, " — DBSCAN-split outputs already exist, skipping ===")
}

# ── 4. Summary ───────────────────────────────────────────────────────────────

message("\n=== Whole-tissue build summary ===")
expected <- paste0("sfe_", wt_cohort)
built    <- expected[dir.exists(file.path(sfe_dir, expected))]
message("Whole-tissue SFEs present (", length(built), "/", length(expected), "):")
for (f in sort(built)) message("  ", f)
missing <- setdiff(expected, built)
if (length(missing) > 0) message("MISSING: ", paste(missing, collapse = ", "))

message("\nDone. Whole-tissue build complete.")
log_session()
