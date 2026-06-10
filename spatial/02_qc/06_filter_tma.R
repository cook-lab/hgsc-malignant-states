# ============================================================================
# 06_filter_tma.R — Remove QC-flagged cores and save filtered TMA SFE
# ============================================================================
# PURPOSE: Apply the core-level QC exclusions (from the core_qc analysis) to the
#   merged TMA SFE. Removes tumour cores with >= 2 QC flags (non-representative
#   fragments dominated by stroma/immune rather than the expected HGSC
#   epithelial composition). Writes the canonical TMA entry-point SFE plus
#   exclusion documentation. FTE cores and off-core cells are retained.
#
# INPUTS:
#   - <sfe_dir>/sfe_tma                                   (unfiltered TMA)
#   - <data_root>/.../07_core_qc/core_qc_summary.csv      (FROZEN deposited per-core
#       decisions = the PUBLISHED 97-patient cohort, 16 cores excluded; CANONICAL.
#       Re-running 05b_core_qc on the current (drifted) sfe_tma gives 14 — the
#       deposited frozen table is authoritative for the paper.)
#   - <data_root>/.../07_core_qc/patient_core_status.csv  (frozen; output_root fallback)
#   - <data_root>/2026_final_xenium_analysis/data/clinical_data_clean.csv
#
# OUTPUTS:
#   - <sfe_dir>/sfe_tma_filtered                          (canonical TMA entry-point)
#   - <output_root>/07_core_qc/excluded_cores_documentation.csv
#
# MANUSCRIPT PANEL(S): produces sfe_tma_filtered, the canonical TMA object
#   feeding Fig 4–7 TMA panels and SF10B/SF12.
#
# RUNTIME TIER: moderate
# ============================================================================

source("spatial/00_setup/00_setup.R")

# ── 1. Load core QC results ──────────────────────────────────────────────────

# The PUBLISHED 97-patient TMA cohort comes from the FROZEN deposited core-QC table
# (16 cores excluded). Re-running 05b_core_qc on the CURRENT (drifted) sfe_tma yields
# 14 excludes (cores 107/126 sit at the composition-outlier threshold and the object's
# cell_labels have drifted since the table was frozen).
# Prefer the deposited frozen table (canonical = paper); fall back to a locally
# regenerated copy under output_root only if the deposited one is absent.
.qc_deposited <- cfg_path("data_root", "2026_final_xenium_analysis", "output",
                          "07_core_qc", "core_qc_summary.csv")
.qc_path <- if (file.exists(.qc_deposited)) .qc_deposited else
            file.path(out_dir, "07_core_qc", "core_qc_summary.csv")
if (!identical(.qc_path, .qc_deposited))
  warning("core_qc_summary.csv: using REGENERATED copy (", .qc_path, "); the published ",
          "97-patient cohort uses the deposited frozen table (16 excludes).")
core_qc <- read.csv(.qc_path, stringsAsFactors = FALSE)

excluded_cores <- core_qc$core_id[core_qc$recommendation == "exclude"]
cat("Cores to exclude:", length(excluded_cores), "\n")
cat("Core IDs:", paste(excluded_cores, collapse = ", "), "\n\n")

# ── 2. Document exclusion reasons ────────────────────────────────────────────

exc <- core_qc[core_qc$recommendation == "exclude", ]

exc$exclusion_reason <- vapply(seq_len(nrow(exc)), function(i) {
  r <- exc[i, ]
  reasons <- c()
  if (r$flag_small)                reasons <- c(reasons, sprintf("small core (n=%d cells)", r$n_cells))
  if (r$flag_low_epi)              reasons <- c(reasons, sprintf("low epithelial (%.1f%%)", r$pct_epithelial))
  if (r$flag_low_confidence)       reasons <- c(reasons, sprintf("low annotation confidence (median score=%.3f, %.1f%% pruned)", r$median_singler_score, r$pct_pruned))
  if (r$flag_composition_outlier)  reasons <- c(reasons, sprintf("composition outlier (Mahalanobis=%.1f)", r$maha_dist))
  if (r$flag_discordant_pair)      reasons <- c(reasons, "discordant with replicate core")
  paste(reasons, collapse = "; ")
}, character(1))

exc$dominant_lineage <- vapply(seq_len(nrow(exc)), function(i) {
  r <- exc[i, ]
  comps <- c(epithelial  = r$pct_epithelial,
             immune      = r$pct_immune,
             stromal     = r$pct_stromal,
             mesothelial = r$pct_mesothelial,
             vascular    = r$pct_vascular)
  paste0(names(which.max(comps)), "-dominated (",
         sprintf("%.0f%%", max(comps)), ")")
}, character(1))

doc <- exc[, c("core_id", "patient_id", "n_cells",
               "pct_epithelial", "pct_immune", "pct_stromal", "pct_mesothelial",
               "dominant_lineage", "n_flags", "exclusion_reason")]
doc <- doc[order(doc$core_id), ]

cat("=== Excluded cores documentation ===\n\n")
for (i in seq_len(nrow(doc))) {
  r <- doc[i, ]
  cat(sprintf("Core %s (patient %s, n=%d): %s\n  %s\n\n",
              r$core_id, r$patient_id, r$n_cells,
              r$dominant_lineage, r$exclusion_reason))
}

write.csv(doc, file.path(out_dir, "07_core_qc", "excluded_cores_documentation.csv"),
          row.names = FALSE)
cat("Saved: 07_core_qc/excluded_cores_documentation.csv\n\n")

# ── 3. Document lost patients ────────────────────────────────────────────────

.ps_deposited <- cfg_path("data_root", "2026_final_xenium_analysis", "output",
                          "07_core_qc", "patient_core_status.csv")
patient_status <- read.csv(
  if (file.exists(.ps_deposited)) .ps_deposited else
    file.path(out_dir, "07_core_qc", "patient_core_status.csv"),
  stringsAsFactors = FALSE)
lost_patients <- patient_status$patient_id[patient_status$passing_cores == 0 &
                                            patient_status$review_cores == 0]

cat("=== Patients losing ALL cores (n=", length(lost_patients), ") ===\n")
cat("Patient IDs:", paste(lost_patients, collapse = ", "), "\n\n")

clinical <- read.csv(file.path(data_dir, "clinical_data_clean.csv"),
                     stringsAsFactors = FALSE)
lost_clin <- clinical[clinical$patient_id %in% lost_patients,
                      c("patient_id", "age", "stage_figo", "treatment_status",
                        "survival_months", "survival_outcome", "residual_binary",
                        "chemo_status_6months")]

cat("Clinical context for lost patients:\n")
for (i in seq_len(nrow(lost_clin))) {
  r <- lost_clin[i, ]
  cat(sprintf("  Patient %s: age %d, stage %s, %s, %s resection, %s months OS (%s)\n",
              r$patient_id, r$age, r$stage_figo, r$treatment_status,
              ifelse(is.na(r$residual_binary), "unknown", r$residual_binary),
              r$survival_months,
              ifelse(r$survival_outcome == 1, "deceased", "alive")))
}

cat("\nNote: lost patients span stages I-IV and include both alive and deceased")
cat("\n      outcomes — no systematic clinical bias from their loss.\n\n")

# ── 4. Load and filter TMA SFE ──────────────────────────────────────────────

message("Loading TMA SFE...")
sfe_tma <- load_sfe("sfe_tma")
cat("Original TMA:", ncol(sfe_tma), "cells\n")

cells_in_excluded <- sfe_tma$core_id %in% as.character(excluded_cores)
cat("Cells in excluded cores:", sum(cells_in_excluded), "\n")

# Keep everything else: passing/review tumour cores, all FT cores, off-core cells
sfe_filtered <- sfe_tma[, !cells_in_excluded]
cat("Filtered TMA:", ncol(sfe_filtered), "cells\n")
cat("Cells removed:", sum(cells_in_excluded), "\n\n")

rm(sfe_tma); gc(verbose = FALSE)

# ── 5. Verify ────────────────────────────────────────────────────────────────

cat("=== Filtered TMA summary ===\n")
cat("Total cells:", ncol(sfe_filtered), "\n")
cat("Genes:", nrow(sfe_filtered), "\n")

cat("\nCells by sample_type:\n")
print(table(sfe_filtered$sample_type, useNA = "ifany"))

tumour_mask <- sfe_filtered$sample_type == "tumour" &
               !is.na(sfe_filtered$core_id) &
               sfe_filtered$core_id != "Off core"
remaining_cores <- length(unique(sfe_filtered$core_id[tumour_mask]))
remaining_patients <- length(unique(sfe_filtered$patient_id[tumour_mask]))

cat("\nRemaining tumour cores:", remaining_cores, "\n")
cat("Remaining tumour patients:", remaining_patients, "\n")

# ── 6. Save filtered SFE ────────────────────────────────────────────────────

message("Saving filtered TMA SFE...")
save_sfe(sfe_filtered, "sfe_tma_filtered")

rm(sfe_filtered); gc(verbose = FALSE)

cat("\nDone. Filtered TMA build complete.\n")
cat("Downstream analyses should use load_sfe('sfe_tma_filtered').\n")
log_session()
