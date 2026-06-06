# ============================================================================
# 03_qc_metrics.R — Compute QC metrics across all Xenium samples
# ============================================================================
# PURPOSE: Compute per-cell QC metrics on each SFE independently, then collect
#   into a single summary table for cross-sample comparison. No filtering — this
#   script only computes and saves metrics.
#   Metrics: total counts (sum), genes detected, negative-control subset sums/%
#   (negProbe, negCodeword, unassigned), cell/nucleus area, prop_nuc, density.
#
# COHORT PIN: operates over sfe_tma + the 8 published whole tissues
#   (CFG$cohort$whole_tissue). The FTE whole-tissue samples are excluded.
#
# INPUTS:
#   - <sfe_dir>/sfe_tma, <sfe_dir>/sfe_<wt> (raw SFEs from 01_build_sfe/)
#
# OUTPUTS:
#   - <output_root>/03_04_qc/qc_metrics_all_samples.csv
#   - <output_root>/03_04_qc/qc_summary_per_sample.csv
#
# MANUSCRIPT PANEL(S): QC backend; feeds 04_qc_filter_normalize.R.
#
# RUNTIME TIER: moderate
# ============================================================================

source("spatial/00_setup/00_setup.R")

# ── 1. Define samples and control probe sets ─────────────────────────────────

sfe_names <- sfe_names_all   # sfe_tma + published 8 whole tissues (cohort PIN)

# Append-only guard: existing samples were already filtered (controls removed),
# so recomputing QC on their post-filter SFEs would change values already in the
# CSV. Load the existing CSV and only process samples missing from it.
out_csv     <- file.path(out_dir, "03_04_qc", "qc_metrics_all_samples.csv")
out_summary <- file.path(out_dir, "03_04_qc", "qc_summary_per_sample.csv")

existing_qc <- if (file.exists(out_csv)) {
  fread(out_csv)
} else {
  data.table()
}
existing_samples <- if (nrow(existing_qc) > 0) unique(existing_qc$sample_id) else character(0)
message("Existing samples in qc_metrics_all_samples.csv: ",
        if (length(existing_samples) > 0) paste(existing_samples, collapse = ", ") else "(none)")

# ── 2. Compute QC metrics per sample ─────────────────────────────────────────

all_metrics <- list()

for (sfe_name in sfe_names) {
  message("\n=== ", sfe_name, " ===")

  sfe_path_probe <- file.path(sfe_dir, sfe_name)
  if (length(existing_samples) > 0) {
    sfe_probe <- loadHDF5SummarizedExperiment(dir = sfe_path_probe)
    probe_sample_ids <- unique(as.character(colData(sfe_probe)$sample_id))
    rm(sfe_probe); gc(verbose = FALSE)
    if (all(probe_sample_ids %in% existing_samples)) {
      message("  Already in CSV (sample_id=", paste(probe_sample_ids, collapse = ","),
              ") — skipping.")
      all_metrics[[sfe_name]] <- existing_qc[sample_id %in% probe_sample_ids]
      next
    }
  }

  sfe <- loadHDF5SummarizedExperiment(dir = sfe_path_probe)

  # Fix sample_id mismatch in spatial graphs for DBSCAN-split samples
  int_meta <- int_metadata(sfe)
  if (!is.null(int_meta$spatialGraphs)) {
    graph_samples <- names(int_meta$spatialGraphs)
    cd_samples <- unique(sfe$sample_id)
    if (!all(graph_samples %in% cd_samples)) {
      int_metadata(sfe)$spatialGraphs <- NULL
    }
  }

  message("  Cells: ", ncol(sfe), " | Genes: ", nrow(sfe))

  rn <- rownames(sfe)
  is_neg_probe    <- grepl("^NegControlProbe_", rn)
  is_neg_codeword <- grepl("^NegControlCodeword_", rn)
  is_unassigned   <- grepl("^UnassignedCodeword_", rn)
  is_any_control  <- is_neg_probe | is_neg_codeword | is_unassigned
  is_biological   <- !is_any_control

  message("  Features: ", sum(is_biological), " biological, ",
          sum(is_neg_probe), " negProbe, ",
          sum(is_neg_codeword), " negCodeword, ",
          sum(is_unassigned), " unassigned")

  subsets <- list(
    negProbe    = which(is_neg_probe),
    negCodeword = which(is_neg_codeword),
    unassigned  = which(is_unassigned),
    any_neg     = which(is_any_control)
  )

  qc <- perCellQCMetrics(sfe, subsets = subsets)

  cd <- as.data.frame(colData(sfe))

  metrics <- data.frame(
    cell_id       = colnames(sfe),
    sample_id     = cd$sample_id,
    sum           = qc$sum,
    detected      = qc$detected,
    subsets_negProbe_sum     = qc$subsets_negProbe_sum,
    subsets_negProbe_percent = qc$subsets_negProbe_percent,
    subsets_negCodeword_sum     = qc$subsets_negCodeword_sum,
    subsets_negCodeword_percent = qc$subsets_negCodeword_percent,
    subsets_unassigned_sum     = qc$subsets_unassigned_sum,
    subsets_unassigned_percent = qc$subsets_unassigned_percent,
    subsets_any_neg_sum     = qc$subsets_any_neg_sum,
    subsets_any_neg_percent = qc$subsets_any_neg_percent,
    cell_area     = cd$cell_area,
    nucleus_area  = cd$nucleus_area,
    stringsAsFactors = FALSE
  )

  metrics$prop_nuc       <- metrics$nucleus_area / metrics$cell_area
  metrics$counts_per_um2 <- metrics$sum / metrics$cell_area

  if ("core_id" %in% names(cd)) {
    metrics$core_id     <- cd$core_id
    metrics$patient_id  <- cd$patient_id
    metrics$sample_type <- cd$sample_type
  }

  message("  sum:      median=", round(median(metrics$sum)),
          ", range=[", min(metrics$sum), ", ", max(metrics$sum), "]")
  message("  detected: median=", round(median(metrics$detected)),
          ", range=[", min(metrics$detected), ", ", max(metrics$detected), "]")
  message("  cell_area: median=", round(median(metrics$cell_area, na.rm = TRUE), 1),
          " um2")
  message("  neg_probe_%: median=", round(median(metrics$subsets_negProbe_percent, na.rm = TRUE), 2),
          "%, 95th=", round(quantile(metrics$subsets_negProbe_percent, 0.95, na.rm = TRUE), 2), "%")

  all_metrics[[sfe_name]] <- metrics
  rm(sfe, qc, cd, metrics); gc(verbose = FALSE)
}

# ── 3. Combine and save ──────────────────────────────────────────────────────

qc_all <- rbindlist(all_metrics, fill = TRUE)
qc_all <- as.data.frame(qc_all)

message("\n=== Combined QC Summary ===")
message("Total cells: ", nrow(qc_all))
message("\nCells per sample_id:")
print(table(qc_all$sample_id))

fwrite(qc_all, out_csv)
message("\nSaved: ", out_csv)
message("File size: ", format(file.size(out_csv), big.mark = ","), " bytes")

# ── 4. Per-sample summary table ──────────────────────────────────────────────

summary_dt <- qc_all |>
  as.data.table() |>
  _[, .(
    n_cells           = .N,
    median_counts     = as.double(median(sum, na.rm = TRUE)),
    median_genes      = as.double(median(detected, na.rm = TRUE)),
    median_cell_area  = round(median(cell_area, na.rm = TRUE), 1),
    median_nuc_area   = round(median(nucleus_area, na.rm = TRUE), 1),
    median_prop_nuc   = round(median(prop_nuc, na.rm = TRUE), 3),
    median_neg_pct    = round(median(subsets_any_neg_percent, na.rm = TRUE), 2),
    p95_neg_pct       = round(quantile(subsets_any_neg_percent, 0.95, na.rm = TRUE), 2),
    median_density    = round(median(counts_per_um2, na.rm = TRUE), 2)
  ), by = sample_id]

message("\nPer-sample summary:")
print(summary_dt, nrows = 20)

fwrite(summary_dt, out_summary)
message("\nSaved: ", out_summary)

message("\nDone. QC metrics computed for all samples.")
log_session()
