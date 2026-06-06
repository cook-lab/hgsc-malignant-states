#!/usr/bin/env Rscript
# ============================================================================
# Supp Data 5 — Xenium gene panel + cross-platform probe QC
# ============================================================================
# PURPOSE
#   Assemble the supplemental table describing the Xenium gene panel after probe
#   QC: which genes pass into the analysis panel vs are excluded from annotation,
#   with cross-platform (scRNA-seq vs Xenium) concordance metrics. Logic preserved
#   verbatim from the canonical generate_xenium_gene_panel_table.R.
#
# INPUTS  (under output_root)
#   - 05_probe_qc/gene_exclusion_decisions.csv   per-gene QC flags + include/exclude tier
#   - 05_probe_qc/probe_qc_full.csv              per-gene cross-platform QC metrics
#
# OUTPUTS
#   - supplemental/Supplemental_Table_5_Xenium_Gene_Panel.csv
#
# MANUSCRIPT PANEL(S)
#   Supp Data 5 (Xenium gene panel + probe QC).
#
# RUNTIME TIER
#   fast (two CSV reads + a merge).
# ============================================================================

# --- central config (tables/ is 1 level below repo root) ---
# Resolve this script's directory so config can be sourced via a relative path
# regardless of the working directory.
`%||%` <- function(a, b) if (is.null(a)) b else a
.this_file <- tryCatch(
  normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  error = function(e) NA_character_)
.script_dir <- if (!is.na(.this_file)) dirname(.this_file) else "tables"
source(file.path(.script_dir, "..", "config", "config.R"))

library(data.table)

# ---- Inputs (config-resolved) ----------------------------------------------
qc_decisions <- fread(cfg_path("output_root", "05_probe_qc", "gene_exclusion_decisions.csv"))
qc_full      <- fread(cfg_path("output_root", "05_probe_qc", "probe_qc_full.csv"))

# ---- Build panel table ------------------------------------------------------
# Drop overlapping QC columns from decisions (keep from qc_full for metrics)
decisions_slim <- qc_decisions[, .(gene, qc_flag, n_flags, marker_for, tier,
                                    decision, reason)]

# Merge: decisions (all genes) + QC metrics (NA for xenium-only)
panel <- merge(decisions_slim, qc_full, by = "gene", all.x = TRUE)

# Classify: genes in final analysis panel vs excluded
panel[, in_analysis_panel := !grepl("EXCLUDE|HARD_EXCLUDE", tier)]

# Clean up decision labels
panel[, status := fcase(
  tier == "INCLUDE",         "Pass",
  tier == "INCLUDE_MONITOR", "Pass (monitored)",
  tier == "KEEP_CAUTION",    "Pass (caution)",
  tier == "OK",              "Pass",
  grepl("TRUE", tier),       "Pass",
  tier == "EXCLUDE",         "Excluded from annotation",
  tier == "HARD_EXCLUDE",    "Excluded from annotation",
  default = "Pass"
)]

# Add xenium-only flag
panel[, xenium_only := (qc_flag == "XENIUM_ONLY")]

# Select and rename columns for supplemental table
supp_table <- panel[, .(
  gene,
  in_analysis_panel,
  qc_status         = status,
  qc_flag,
  n_qc_flags        = n_flags,
  xenium_only,
  marker_for,
  corr_of_corr      = round(rho_corr_of_corr, 3),
  jaccard_top20     = round(top20_jaccard, 3),
  sc_detect_rate    = round(sc_detect_rate, 3),
  xe_detect_rate    = round(xe_detect_rate, 3),
  exclusion_reason  = reason
)]

# Sort: included genes first (alphabetical), then excluded (alphabetical)
setorder(supp_table, -in_analysis_panel, gene)

# ---- Summary ----------------------------------------------------------------
message("Xenium Gene Panel Summary")
message("  Total genes (post control removal):  ", nrow(supp_table))
message("  Included in analysis panel:           ", sum(supp_table$in_analysis_panel))
message("  Excluded from annotation:             ", sum(!supp_table$in_analysis_panel))
message("  Xenium-only (no scRNA-seq match):     ", sum(supp_table$xenium_only))

# ---- Write ------------------------------------------------------------------
out_path <- cfg_path("output_root", "supplemental", "Supplemental_Table_5_Xenium_Gene_Panel.csv")
fwrite(supp_table, out_path)
message("\nWritten: ", out_path)
