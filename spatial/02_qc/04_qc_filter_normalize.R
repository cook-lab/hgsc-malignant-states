# ============================================================================
# 04_qc_filter_normalize.R — Filter cells, remove controls, normalize
# ============================================================================
# PURPOSE: Apply QC filters to each SFE independently, then:
#   - remove the 64 control-probe features (keep 477 biological genes)
#   - normalize: logNormCounts (library-size) + an area-adjusted variant
#   Thresholds (uniform across samples): counts >= 10, genes >= 5,
#   neg_pct <= 5%, cell area 10–500 um^2. Lenient — removes only clear
#   technical artifacts (neg controls ~0%, no bimodality).
#
# COHORT PIN: operates over sfe_tma + the 8 published whole tissues.
#
# INPUTS:
#   - <sfe_dir>/sfe_tma, <sfe_dir>/sfe_<wt>  (raw SFEs)
#
# OUTPUTS:
#   - <sfe_dir>/sfe_*   updated in place (controls removed; logcounts +
#       logcounts_area assays added; HDF5-backed via _v2 swap)
#   - <output_root>/03_04_qc/filtering_summary.csv
#
# MANUSCRIPT PANEL(S): QC backend; produces the normalized SFEs consumed by
#   annotation (06_annotation.R) and all downstream spatial analyses.
#
# RUNTIME TIER: heavy
# ============================================================================

source("spatial/00_setup/00_setup.R")

# ── 1. Define filtering thresholds ─────────────────────────────────────────

THRESH <- list(
  min_counts  = 10,
  min_genes   = 5,
  max_neg_pct = 5,
  min_area    = 10,
  max_area    = 500
)

message("=== Filter, Remove Controls, Normalize ===")
message("Thresholds:")
for (nm in names(THRESH)) message("  ", nm, " = ", THRESH[[nm]])

# ── 2. Enumerate SFE directories ──────────────────────────────────────────

sfe_names <- sfe_names_all   # sfe_tma + published 8 whole tissues (cohort PIN)

# ── 3. Process each SFE ───────────────────────────────────────────────────

results <- list()

for (sfe_name in sfe_names) {

  message("\n", "=" |> rep(60) |> paste(collapse = ""))
  message("Processing: ", sfe_name)
  message("=" |> rep(60) |> paste(collapse = ""))

  sfe_path <- file.path(sfe_dir, sfe_name)

  # Skip if already processed (has logcounts assay)
  sfe <- loadHDF5SummarizedExperiment(dir = sfe_path)
  if ("logcounts" %in% assayNames(sfe)) {
    message("  Already processed (logcounts found). Skipping.")
    results[[sfe_name]] <- data.table(
      sfe_name = sfe_name, cells_before = ncol(sfe), cells_after = ncol(sfe),
      cells_removed = 0, pct_kept = 100, fail_counts = NA_integer_,
      fail_genes = NA_integer_, fail_neg = NA_integer_,
      fail_min_area = NA_integer_, fail_max_area = NA_integer_,
      genes_before = nrow(sfe), genes_after = nrow(sfe),
      controls_removed = 0, assays = paste(assayNames(sfe), collapse = ", ")
    )
    rm(sfe); gc(verbose = FALSE)
    next
  }

  n_before <- ncol(sfe)
  g_before <- nrow(sfe)
  message("  Loaded: ", format(n_before, big.mark = ","), " cells x ", g_before, " genes")

  # Fix spatial graph mismatch (DBSCAN-split samples have parent sample_id)
  if (!is.null(int_metadata(sfe)$spatialGraphs)) {
    graph_samples <- tryCatch(
      names(int_metadata(sfe)$spatialGraphs),
      error = function(e) NULL
    )
    cd_samples <- unique(colData(sfe)$sample_id)
    if (!is.null(graph_samples) && !all(graph_samples %in% cd_samples)) {
      int_metadata(sfe)$spatialGraphs <- NULL
      message("  Cleared mismatched spatial graphs")
    }
  }

  # ── 3a. Compute QC metrics for filtering ─────────────────────────────────

  rn <- rownames(sfe)
  is_neg_probe    <- grepl("^NegControlProbe_", rn)
  is_neg_codeword <- grepl("^NegControlCodeword_", rn)
  is_unassigned   <- grepl("^UnassignedCodeword_", rn)
  is_any_control  <- is_neg_probe | is_neg_codeword | is_unassigned
  is_biological   <- !is_any_control

  message("  Features: ", sum(is_biological), " biological, ",
          sum(is_any_control), " controls")

  qc_metrics <- perCellQCMetrics(sfe, subsets = list(
    any_neg = which(is_any_control)
  ))

  sfe$qc_sum      <- qc_metrics$sum
  sfe$qc_detected <- qc_metrics$detected
  sfe$qc_neg_pct  <- qc_metrics$subsets_any_neg_percent

  # ── 3b. Apply cell filters ───────────────────────────────────────────────

  f_counts  <- sfe$qc_sum >= THRESH$min_counts
  f_genes   <- sfe$qc_detected >= THRESH$min_genes
  f_neg     <- sfe$qc_neg_pct <= THRESH$max_neg_pct
  f_min_area <- sfe$cell_area >= THRESH$min_area
  f_max_area <- sfe$cell_area <= THRESH$max_area

  n_fail_counts  <- sum(!f_counts)
  n_fail_genes   <- sum(!f_genes)
  n_fail_neg     <- sum(!f_neg)
  n_fail_minarea <- sum(!f_min_area)
  n_fail_maxarea <- sum(!f_max_area)

  pass_all <- f_counts & f_genes & f_neg & f_min_area & f_max_area
  n_after <- sum(pass_all)
  n_removed <- n_before - n_after

  message(sprintf("  Filtering: %s -> %s cells (removed %s, %.2f%%)",
                  format(n_before, big.mark = ","),
                  format(n_after, big.mark = ","),
                  format(n_removed, big.mark = ","),
                  100 * n_removed / n_before))
  message(sprintf("    counts < %d:  %d", THRESH$min_counts, n_fail_counts))
  message(sprintf("    genes < %d:   %d", THRESH$min_genes, n_fail_genes))
  message(sprintf("    neg > %.0f%%:    %d", THRESH$max_neg_pct, n_fail_neg))
  message(sprintf("    area < %d:    %d", THRESH$min_area, n_fail_minarea))
  message(sprintf("    area > %d:   %d", THRESH$max_area, n_fail_maxarea))

  sfe <- sfe[, pass_all]

  if (!is.null(int_metadata(sfe)$spatialGraphs)) {
    int_metadata(sfe)$spatialGraphs <- NULL
    message("  Cleared spatial graphs (invalidated by cell subsetting)")
  }

  sfe$qc_sum      <- NULL
  sfe$qc_detected <- NULL
  sfe$qc_neg_pct  <- NULL

  # ── 3c. Remove control probes ─────────────────────────────────────────────

  rn2 <- rownames(sfe)
  is_bio2 <- !grepl("^(NegControlProbe_|NegControlCodeword_|UnassignedCodeword_)", rn2)
  n_controls_removed <- sum(!is_bio2)

  sfe <- sfe[is_bio2, ]
  g_after <- nrow(sfe)
  message(sprintf("  Control removal: %d -> %d genes (removed %d controls)",
                  g_before, g_after, n_controls_removed))

  # ── 3d. Normalize ──────────────────────────────────────────────────────────

  # Library-size normalization (standard for SingleR)
  sfe <- logNormCounts(sfe, assay.type = "counts")
  message("  Added assay: logcounts (library-size normalized)")

  # Area-adjusted normalization: size factors proportional to cell area
  area_sf <- sfe$cell_area / mean(sfe$cell_area)
  sizeFactors(sfe) <- area_sf
  sfe <- logNormCounts(sfe, assay.type = "counts", name = "logcounts_area")
  message("  Added assay: logcounts_area (area-adjusted normalized)")

  # Reset size factors to library-size (default = standard normalization)
  sfe <- computeLibraryFactors(sfe)

  message("  Assays: ", paste(assayNames(sfe), collapse = ", "))

  # ── 3e. Save ───────────────────────────────────────────────────────────────
  # Realize the HDF5-backed counts assay in memory before the _v2 swap.
  assay(sfe, "counts") <- as(assay(sfe, "counts"), "dgCMatrix")
  message("  Realized counts matrix in memory")

  sfe_path_new <- paste0(sfe_path, "_v2")
  if (dir.exists(sfe_path_new)) unlink(sfe_path_new, recursive = TRUE)
  saveHDF5SummarizedExperiment(sfe, dir = sfe_path_new)

  unlink(sfe_path, recursive = TRUE)
  file.rename(sfe_path_new, sfe_path)
  message("  Saved: ", sfe_path)

  results[[sfe_name]] <- data.table(
    sfe_name        = sfe_name,
    cells_before    = n_before,
    cells_after     = n_after,
    cells_removed   = n_removed,
    pct_kept        = round(100 * n_after / n_before, 2),
    fail_counts     = n_fail_counts,
    fail_genes      = n_fail_genes,
    fail_neg        = n_fail_neg,
    fail_min_area   = n_fail_minarea,
    fail_max_area   = n_fail_maxarea,
    genes_before    = g_before,
    genes_after     = g_after,
    controls_removed = n_controls_removed,
    assays          = paste(assayNames(sfe), collapse = ", ")
  )

  rm(sfe, qc_metrics)
  gc(verbose = FALSE)
}

# ── 4. Summary table ──────────────────────────────────────────────────────

summary_dt <- rbindlist(results)

message("\n", "=" |> rep(60) |> paste(collapse = ""))
message("FILTERING SUMMARY")
message("=" |> rep(60) |> paste(collapse = ""))
print(summary_dt[, .(sfe_name, cells_before, cells_after, pct_kept,
                      genes_after, assays)])

message(sprintf("\nTotal: %s -> %s cells (%.1f%% retained)",
                format(sum(summary_dt$cells_before), big.mark = ","),
                format(sum(summary_dt$cells_after), big.mark = ","),
                100 * sum(summary_dt$cells_after) / sum(summary_dt$cells_before)))
message(sprintf("Genes: %d -> %d (removed %d controls)",
                summary_dt$genes_before[1],
                summary_dt$genes_after[1],
                summary_dt$controls_removed[1]))

out_path <- file.path(out_dir, "03_04_qc", "filtering_summary.csv")
fwrite(summary_dt, out_path)
message("\nSaved: ", out_path)

message("\nDone. All SFEs filtered, control probes removed, and normalized.")
log_session()
