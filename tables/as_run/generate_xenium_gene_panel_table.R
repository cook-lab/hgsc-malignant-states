# ============================================================================
# generate_xenium_gene_panel_table.R
# ============================================================================
# Generates Supplemental Table: Xenium gene panel after probe QC filtering
#
# Source data:  output/05_probe_qc/gene_exclusion_decisions.csv
#               output/05_probe_qc/probe_qc_full.csv
#
# Output:       supplemental_material/Supplemental_Table_Xenium_Gene_Panel.csv
# ============================================================================

library(data.table)

# ── Paths ──────────────────────────────────────────────────────────────────

proj_root    <- "/Volumes/CookLab/Sarah/2026_final_xenium_analysis"
ms_root      <- "/Volumes/CookLab/Sarah/202605_epitype_manuscript"
out_dir      <- file.path(ms_root, "supplemental_material")

qc_decisions <- fread(file.path(proj_root, "output/05_probe_qc/gene_exclusion_decisions.csv"))
qc_full      <- fread(file.path(proj_root, "output/05_probe_qc/probe_qc_full.csv"))

# ── Build panel table ──────────────────────────────────────────────────────

# Drop overlapping QC columns from decisions (keep from qc_full for metrics)
decisions_slim <- qc_decisions[, .(gene, qc_flag, n_flags, marker_for, tier,
                                    decision, reason)]

# Merge: decisions (all 477 genes) + QC metrics (457 genes; NA for xenium-only)
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

# ── Summary ────────────────────────────────────────────────────────────────

n_total    <- nrow(supp_table)
n_included <- sum(supp_table$in_analysis_panel)
n_excluded <- sum(!supp_table$in_analysis_panel)
n_xenium   <- sum(supp_table$xenium_only)

message("Xenium Gene Panel Summary")
message("  Total genes (post control removal):  ", n_total)
message("  Included in analysis panel:           ", n_included)
message("  Excluded from annotation:             ", n_excluded)
message("  Xenium-only (no scRNA-seq match):     ", n_xenium)

# ── Write ──────────────────────────────────────────────────────────────────

out_path <- file.path(out_dir, "Supplemental_Table_Xenium_Gene_Panel.csv")
fwrite(supp_table, out_path)
message("\nWritten: ", out_path)
