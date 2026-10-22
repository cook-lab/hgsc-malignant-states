#!/usr/bin/env Rscript
# ============================================================================
# 04c_pathway_scoring_allcells.R
# ----------------------------------------------------------------------------
# PURPOSE: Per-cell pathway activity scoring (UCell, the 37 modules of
#   9b_scoring/pathway_gene_sets_v2.csv) written into each SFE as pathway_*
#   colData; any pathway_* columns already in the SFE are dropped first.
#
# SUBSTRATE (author decision 2026-09-21):
#   counts(sfe): the FULL 477-gene panel, over ALL CELLS. UCell ranks genes
#   WITHIN each cell and maxRank is fixed at 477 for every sample, so scores are
#   per-cell independent and comparable across samples.
#   The 28 probe-QC exclusions (genes_exclude_singler.txt) apply to annotation
#   only: flagged probes should not drive cell-type calls, but they are retained
#   for expression and pathway scoring.
#
# PROBE-QC TRANSPARENCY: pathway_coverage_allcells.csv records, per set, which genes
#   carry a probe-QC FAIL or WARN, so any caption resting on such a set can say so.
#   Sets affected by a FAIL gene: matrix_remodeling (MMP7), adc_targets
#   (ERBB2, VTCN1), apoptosis_pro (BAX), ecm_adhesion (CLDN3).
#
# INPUTS : SFEs via load_sfe() (the 8 whole tissues, TMA A, TMA B blocks G and H);
#          9b_scoring/pathway_gene_sets_v2.csv (04_pathway_scoring.R or 04b);
#          probe-QC flag table <data_root>/2026_final_xenium_analysis/output/05_probe_qc/flag_qc.csv
#          (optional: without it the FAIL / WARN columns stay empty)
# OUTPUTS: SFEs updated with pathway_* colData (all cells scored);
#          9b_scoring/pathway_coverage_allcells.csv, scoring_summary_allcells.csv
#          (both read by tables/09_tableS9_ucell_gene_sets.py)
# MANUSCRIPT PANEL(S): the pathway scores of Figs 4f, 4j and 5 and Figures
#   S16–S18; Table S9.
# RUNTIME TIER: heavy (~1.5 h over the 11 SFEs)
# ============================================================================
source("spatial/00_setup/00_setup.R")
suppressPackageStartupMessages({
  library(SpatialFeatureExperiment); library(SummarizedExperiment); library(HDF5Array)
  library(data.table); library(Seurat); library(UCell)
})
R      <- cfg_path("output_root")
OUT9B  <- file.path(R, "9b_scoring"); dir.create(OUT9B, showWarnings = FALSE, recursive = TRUE)
# the 8 whole tissues, TMA A and TMA B blocks G and H (every SFE used in the manuscript)
SAMPLES <- c("sfe_WT01","sfe_WT02","sfe_WT03","sfe_WT04",
             "sfe_WT05","sfe_WT06","sfe_WT07","sfe_WT08",
             "sfe_tma","sfe_COEUR_TMA_G","sfe_COEUR_TMA_H")

gs   <- fread(file.path(OUT9B, "pathway_gene_sets_v2.csv"))
sets <- split(gs$gene, gs$pathway)
fq   <- tryCatch(fread(cfg_path("data_root","2026_final_xenium_analysis","output",
                                "05_probe_qc","flag_qc.csv")), error = function(e) NULL)
qcmap <- if (!is.null(fq)) setNames(fq$qc_flag, fq$gene) else character(0)
message(sprintf("Pathway sets: %d (%d unique genes)", length(sets), uniqueN(gs$gene)))

save_sfe_safely <- function(sfe, name) {
  base <- file.path(cfg_path("output_root"), "sfe", name)
  tmp  <- paste0(base, "__new")
  if (dir.exists(tmp)) unlink(tmp, recursive = TRUE)
  saveHDF5SummarizedExperiment(sfe, dir = tmp, replace = TRUE)
  chk <- tryCatch({ z <- loadHDF5SummarizedExperiment(tmp); ncol(z) }, error = function(e) -1L)
  if (chk != ncol(sfe)) stop("verification failed for ", name, " (reloaded ncol=", chk, ")")
  old <- paste0(base, "__old")
  if (dir.exists(old)) unlink(old, recursive = TRUE)
  file.rename(base, old); file.rename(tmp, base); unlink(old, recursive = TRUE)
  message("    saved + verified: ", base)
}

cov_rows <- list(); sum_rows <- list()
for (nm in SAMPLES) {
  t0 <- Sys.time(); message("\n[", format(t0), "] === ", nm, " ===")
  sfe <- load_sfe(nm)
  M <- counts(sfe)
  M <- as(M, "CsparseMatrix")
  ng <- nrow(M)
  message(sprintf("  full panel: %d genes x %d cells (ALL cells; no SPLIT filter)",
                  ng, ncol(M)))

  present <- lapply(sets, function(g) intersect(g, rownames(M)))
  for (p in names(sets)) {
    miss <- setdiff(sets[[p]], rownames(M))
    fl   <- intersect(sets[[p]], names(qcmap)[qcmap == "FAIL"])
    wn   <- intersect(sets[[p]], names(qcmap)[qcmap == "WARN"])
    cov_rows[[length(cov_rows)+1]] <- data.table(
      sample = nm, pathway = p, n_requested = length(sets[[p]]),
      n_found = length(present[[p]]), missing_genes = paste(miss, collapse=";"),
      probeqc_FAIL_genes = paste(fl, collapse=";"),
      probeqc_WARN_genes = paste(wn, collapse=";"))
  }
  partial <- names(present)[lengths(present) < lengths(sets)]
  if (length(partial)) message("  partial sets: ", paste(sprintf("%s(%d/%d)",
      partial, lengths(present)[partial], lengths(sets)[partial]), collapse=", "))
  else message("  all 37 sets complete on the full panel")
  scoreable <- present[lengths(present) >= 2]

  so <- CreateSeuratObject(counts = M)
  so <- AddModuleScore_UCell(so, features = scoreable, maxRank = ng, name = NULL)
  sc <- so@meta.data[, names(scoreable), drop = FALSE]
  rownames(sc) <- colnames(M)
  rm(so, M); gc(verbose = FALSE)

  idx <- match(colnames(sfe), rownames(sc))
  for (c0 in grep("^pathway_", colnames(colData(sfe)), value = TRUE)) colData(sfe)[[c0]] <- NULL
  for (p in names(scoreable)) colData(sfe)[[paste0("pathway_", p)]] <- sc[[p]][idx]
  message(sprintf("  scored %d sets on %d cells (0 NA expected)",
                  length(scoreable), sum(!is.na(idx))))
  save_sfe_safely(sfe, nm)
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  sum_rows[[nm]] <- data.table(sample=nm, n_genes=ng, n_cells_scored=sum(!is.na(idx)),
                               n_cells_sfe=ncol(sfe), n_sets=length(scoreable),
                               elapsed_min=round(el,1))
  message(sprintf("  done in %.1f min", el))
  rm(sfe, sc); gc(verbose = FALSE)
}
fwrite(rbindlist(cov_rows), file.path(OUT9B, "pathway_coverage_allcells.csv"))
fwrite(rbindlist(sum_rows), file.path(OUT9B, "scoring_summary_allcells.csv"))
message("\nALL DONE")
