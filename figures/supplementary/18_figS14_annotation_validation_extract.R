#!/usr/bin/env Rscript
# =============================================================================
# Figure S14 (step 1/2) — Xenium annotation validation: per-cell-type expression
#                         aggregates for the canonical-marker dot plot, read from the SFEs
# -----------------------------------------------------------------------------
# PURPOSE
#   Per-cell-type sums of logcounts (the exact SingleR test input: scuttle
#   library-size logNormCounts over the 477-gene panel) and the number of
#   expressing cells (counts > 0), over ALL cells, grouped by
#   sample x patient x core x sample_type x cell_label_singler x pruned
#   (pruned = SingleR::pruneScores() low-confidence flag, colData singler_pruned,
#   written by spatial/03_annotation_polarization/01_annotation_singler.R; no SingleR
#   re-run). The Python renderer (19_figS14_annotation_validation_plot.py) applies
#   the cohort (final TMA core set) and pools.
#   Only READS the SFEs (no writes to sfe/).
#
# INPUTS:   output_root/sfe/sfe_<sample>          (tma, COEUR_TMA_G/H, core-8 WT)
#           output_root/singler_recompute/<sample>_singler_cells.csv.gz
#               (cell_label_singler = SingleR label, Secretory split SecA/Int/SecB
#                by the frozen 06f polarization thresholds; patient_id, sample_type)
# OUTPUTS:  output_root/_prep_caches/sf12/<sample>_agg_groups.csv
#           output_root/_prep_caches/sf12/<sample>_agg_sumlog.csv.gz   (groups x genes)
#           output_root/_prep_caches/sf12/<sample>_agg_nexpr.csv.gz    (groups x genes)
#           (cache folder and the SF12_SAMPLES variable carry the earlier number SF12 = Figure S14)
# MANUSCRIPT PANEL(S): Figure S14 (data).
# RUNTIME:  moderate (~10-20 min; HDF5 logcounts read in 50k-cell blocks).
#
# Run from repo root:
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OUTPUT_ROOT=<outputs>
#   Rscript figures/supplementary/18_figS14_annotation_validation_extract.R
# =============================================================================

source("spatial/00_setup/00_setup.R")
suppressPackageStartupMessages({library(data.table); library(Matrix)})

OUT <- file.path(out_dir, "_prep_caches", "sf12")
SGR <- file.path(out_dir, "singler_recompute")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

samples <- c("tma", "COEUR_TMA_G", "COEUR_TMA_H",
             "WT01", "WT02", "WT03", "WT04",
             "WT05", "WT06", "WT07", "WT08")
if (nzchar(Sys.getenv("SF12_SAMPLES"))) samples <- strsplit(Sys.getenv("SF12_SAMPLES"), ",")[[1]]
BLOCK <- 50000L

for (s in samples) {
  t0 <- Sys.time()
  sfe <- load_sfe(paste0("sfe_", s))
  cd <- colData(sfe)
  sc <- data.table(barcode        = colnames(sfe),
                   singler_label  = as.character(cd$singler_label),
                   singler_pruned = as.character(cd$singler_pruned),
                   singler_score  = as.numeric(cd$singler_score),
                   size_factor    = as.numeric(cd$sizeFactor))

  # labels / patient from the canonical SingleR per-cell table
  lab <- fread(file.path(SGR, paste0(s, "_singler_cells.csv.gz")))
  stopifnot(nrow(lab) == ncol(sfe))
  m <- match(colnames(sfe), lab$barcode)
  stopifnot(!anyNA(m))
  lab <- lab[m]
  stopifnot(identical(lab$singler_label, sc$singler_label))
  pid <- if ("patient_id" %in% names(lab)) as.character(lab$patient_id) else rep(NA_character_, nrow(lab))
  pid <- sub("\\.0$", "", pid)
  st  <- if ("sample_type" %in% names(lab)) as.character(lab$sample_type) else rep(NA_character_, nrow(lab))
  core <- if ("core_id" %in% names(lab)) as.character(lab$core_id) else rep(NA_character_, nrow(lab))
  key <- data.table(sample = s, patient_id = pid, core_id = core, sample_type = st,
                    cell_label = lab$cell_label_singler,
                    pruned = is.na(sc$singler_pruned))
  key[, gkey := paste(sample, patient_id, core_id, sample_type, cell_label, pruned, sep = "|")]
  groups <- unique(key[, .(gkey, sample, patient_id, core_id, sample_type, cell_label, pruned)])
  groups[, gi := .I]
  g_idx <- groups$gi[match(key$gkey, groups$gkey)]
  groups[, n_cells := tabulate(g_idx, nbins = nrow(groups))]

  genes <- rownames(sfe)
  SUM <- matrix(0, nrow(groups), length(genes), dimnames = list(NULL, genes))
  NEX <- matrix(0, nrow(groups), length(genes), dimnames = list(NULL, genes))
  lc <- SummarizedExperiment::assay(sfe, "logcounts")
  n <- ncol(sfe)
  for (b0 in seq(1L, n, by = BLOCK)) {
    idx <- b0:min(n, b0 + BLOCK - 1L)
    X <- as(as.matrix(lc[, idx]), "dgCMatrix")                 # genes x cells
    G <- sparseMatrix(i = g_idx[idx], j = seq_along(idx), x = 1,
                      dims = c(nrow(groups), length(idx)))      # groups x cells
    SUM <- SUM + as.matrix(G %*% t(X))
    Xb <- X; Xb@x <- rep(1, length(Xb@x))                       # logcounts > 0 <=> counts > 0
    NEX <- NEX + as.matrix(G %*% t(Xb))
  }
  fwrite(groups[, .(gi, sample, patient_id, core_id, sample_type, cell_label, pruned, n_cells)],
         file.path(OUT, paste0(s, "_agg_groups.csv")))
  fwrite(as.data.table(SUM), file.path(OUT, paste0(s, "_agg_sumlog.csv.gz")))
  fwrite(as.data.table(NEX), file.path(OUT, paste0(s, "_agg_nexpr.csv.gz")))
  message(sprintf("  %-12s %8d cells  groups=%d  pruned=%d  %.1f min", s, n, nrow(groups),
                  sum(key$pruned), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(sfe, cd, sc, lab, key, SUM, NEX, lc); gc(verbose = FALSE)
}
message("EXTRACT_DONE")
