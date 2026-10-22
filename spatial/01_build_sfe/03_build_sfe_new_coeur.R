# ============================================================================
# 03_build_sfe_new_coeur.R — Build the TMA B SFEs (COEUR blocks G and H)
# ============================================================================
# PURPOSE: Load the Xenium runs of TMA B (COEUR_TMA_G, COEUR_TMA_H) and attach
#   the de-array core/patient metadata (exported by the de-array editor to
#   <DEARRAY_DIR>/<TAG>_cell_core_assignment.parquet) by cell-ID join. All cells
#   kept; only raw counts stored (normalization happens in
#   spatial/02_qc/04_qc_filter_normalize.R).
#
# COHORT: TMA B block G (149 patients) + block H (168 patients), patient ids from
#   the COEUR-TFRI column of the de-array tables. Off-core cells keep core_id NA.
#
# INPUTS (paths configurable via env; defaults = the final workspace):
#   - RAW_XENIUM_DIR/<sample>/[output-.../]cells.parquet + full Xenium bundle
#   - DEARRAY_DIR/<TAG>_cell_core_assignment.parquet   (cell_id -> core/patient)
#
# OUTPUTS:
#   - <RERUN_SFE_DIR>/sfe_<sample>   (HDF5-backed SFE, raw counts + de-array colData)
#
# MANUSCRIPT PANEL(S): the TMA B arm of Figs 4–7, Figures S13–S15, Tables S8 and S11.
# RUNTIME TIER: heavy (loads two full Xenium TMA runs)
# ============================================================================

source("spatial/00_setup/00_setup.R")
suppressPackageStartupMessages({ library(arrow) })

# SFEs are written to RERUN_SFE_DIR (default OUTPUT_ROOT/sfe); a sample whose
# folder already exists is skipped.
sfe_out <- Sys.getenv("RERUN_SFE_DIR",
                      file.path(Sys.getenv("OUTPUT_ROOT", "./output"), "sfe"))
dir.create(sfe_out, recursive = TRUE, showWarnings = FALSE)
message("Re-run SFE output dir: ", sfe_out)

# ── Paths (override via env for the deposit; default = final workspace) ───────
RAW  <- Sys.getenv("RAW_XENIUM_DIR",
                   file.path(Sys.getenv("DATA_ROOT", "/path/to/deposit"), "202609_final_hgsc_spatial/data/HGSC xenium "))
DEAR <- Sys.getenv("DEARRAY_DIR",
                   file.path(Sys.getenv("DATA_ROOT", "/path/to/deposit"), "202609_final_hgsc_spatial/data/dearray"))

# sample -> (raw subdir, de-array tag, has patient map)
samples <- list(
  COEUR_TMA_G   = list(sub = "COEUR_TMA_G",   tag = "G",    patient = TRUE,  segs = c("cell","nucleus")),
  COEUR_TMA_H   = list(sub = "COEUR_TMA_H",   tag = "H",    patient = TRUE,  segs = c("cell","nucleus"))
)

# readXenium can't read the newer multi-channel morphology_focus/ format and
# match.arg() rejects an empty image list, so we read from a STAGING dir that
# symlinks every file EXCEPT the morphology images -> readXenium then skips
# images with a warning and builds counts + segmentation geometries from parquet.
stage_without_images <- function(xdir) {
  st <- file.path(tempdir(), paste0("stage_", basename(xdir), "_", as.integer(runif(1,1e5,1e6))))
  if (dir.exists(st)) unlink(st, recursive = TRUE)
  dir.create(st, recursive = TRUE)
  for (f in list.files(xdir, full.names = TRUE)) {
    b <- basename(f)
    if (grepl("^morphology", b) || grepl("^\\._", b)) next
    file.symlink(normalizePath(f), file.path(st, b))
  }
  st
}

# Resolve the Xenium run dir (some samples nest under an output-XETG* folder)
resolve_xenium <- function(sub) {
  d <- file.path(RAW, sub)
  if (file.exists(file.path(d, "cells.parquet"))) return(d)
  nested <- list.dirs(d, recursive = FALSE)
  nested <- nested[file.exists(file.path(nested, "cells.parquet"))]
  if (length(nested) >= 1) return(nested[[1]])
  stop("No cells.parquet under: ", d)
}

for (sid in names(samples)) {
  spec <- samples[[sid]]
  message("\n=== ", sid, " (Bloc ", spec$tag, ") ===")

  out_path <- file.path(sfe_out, paste0("sfe_", sid))
  if (dir.exists(out_path)) { message("  Already built: ", out_path, " — skipping."); next }

  xdir  <- resolve_xenium(spec$sub)
  stage <- stage_without_images(xdir)
  message("  Reading Xenium (staged, no images): ", xdir)
  sfe <- readXenium(stage, row.names = "symbol", sample_id = sid,
                    segmentations = spec$segs)
  unlink(stage, recursive = TRUE)
  message("  Loaded: ", ncol(sfe), " cells, ", nrow(sfe), " genes")

  # ── attach de-array metadata by cell-ID join ──────────────────────────────
  da <- as.data.frame(arrow::read_parquet(
    file.path(DEAR, paste0(spec$tag, "_cell_core_assignment.parquet"))))
  idx <- match(colnames(sfe), da$cell_id)
  matched <- sum(!is.na(idx))
  message("  De-array cell-ID match: ", matched, " / ", ncol(sfe),
          " (", round(100 * matched / ncol(sfe), 1), "%)")

  sfe$core_id     <- ifelse(!is.na(idx), da$core_id[idx],   NA_character_)
  sfe$core_rank   <- ifelse(!is.na(idx), da$core_rank[idx], NA)
  sfe$core_flag   <- ifelse(!is.na(idx) & !is.na(da$core_id[idx]) & da$core_id[idx] != "",
                            "On core", "Off core")
  if (isTRUE(spec$patient)) {
    sfe$patient_id  <- ifelse(!is.na(idx), as.character(da$patient_coeur_tfri[idx]), NA_character_)
    sfe$sample_type <- ifelse(sfe$core_flag == "On core", "tumour", NA_character_)  # COEUR TMA = HGSC tumour cores
    sfe$cohort      <- "COEUR"
  }

  message("  On-core cells: ", sum(sfe$core_flag == "On core"),
          "; unique cores: ", length(unique(na.omit(sfe$core_id[sfe$core_id != ""]))),
          "; unique patients: ", length(unique(na.omit(sfe$patient_id))))

  saveHDF5SummarizedExperiment(sfe, dir = out_path, replace = TRUE)
  message("  Saved: ", out_path)
  rm(sfe, da); gc(verbose = FALSE)
}

message("\nDone. TMA B (COEUR) SFEs built.")
log_session()
