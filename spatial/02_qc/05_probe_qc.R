# ============================================================================
# 05_probe_qc.R — Probe QC (cross-platform validation)
# ============================================================================
# PURPOSE: Compare each of the 477 Xenium panel genes against a matched
#   scRNA-seq reference (16K cells, 16 cell types) to flag probes with aberrant
#   cross-platform behavior. Core analysis (gene-gene Spearman
#   correlation-of-correlations, top-20 partner Jaccard, per-gene summary stats,
#   multi-criteria PASS/WARN/FAIL flagging, cell-type marker coverage) is
#   implemented in the probe-QC Rmd, which this staged script renders and then
#   verifies/summarizes.
#
# INPUTS:
#   - <sfe_dir>/sfe_*  (normalized SFEs)
#   - <data_root>/.../xenium_celltype reference (read inside the Rmd)
#   - probe-QC analysis notebook (see NOTE below)
#
# OUTPUTS:
#   - <output_root>/05_probe_qc/probe_qc_full.csv          (Supp Data 5)
#   - <output_root>/05_probe_qc/flag_qc.csv,
#       gene_exclusion_decisions.csv, celltype_coverage.csv
#   - <output_root>/05_probe_qc/explore_probe_qc.html
#   - <output_root>/05_probe_qc/genes_exclude_singler.txt,
#       genes_monitor_singler.txt   (consumed by 06_annotation.R / 06g)
#
# MANUSCRIPT PANEL(S): Supp Data 5 (Xenium gene-panel QC table).
#
# RUNTIME TIER: moderate
#
# NOTE: the heavy interactive computation lives in the probe-QC analysis
# notebook. Point PROBE_QC_RMD at it (the original lived in the exploratory
# sandbox, which is not part of the canonical migrate-set). The notebook itself
# uses the central config for paths.
# ============================================================================

source("spatial/00_setup/00_setup.R")

message("=== Probe QC ===")

# ── 1. Render the HTML report ────────────────────────────────────────────────

rmd_path <- Sys.getenv("PROBE_QC_RMD", unset = "spatial/02_qc/explore_probe_qc.Rmd")
out_path <- file.path(out_dir, "05_probe_qc")

if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)

if (!file.exists(rmd_path)) {
  stop("Probe-QC notebook not found: ", rmd_path,
       "\n  Set the PROBE_QC_RMD env var to its path before running.")
}

message("Rendering: ", rmd_path)
rmarkdown::render(
  input       = rmd_path,
  output_dir  = out_path,
  output_file = "explore_probe_qc.html",
  quiet       = FALSE
)
message("Report saved: ", file.path(out_path, "explore_probe_qc.html"))

# ── 2. Verify outputs ───────────────────────────────────────────────────────

expected_files <- c("probe_qc_full.csv", "flag_qc.csv",
                    "celltype_coverage.csv", "explore_probe_qc.html")
for (f in expected_files) {
  fp <- file.path(out_path, f)
  if (file.exists(fp)) {
    message("  OK: ", f, " (", format(file.size(fp), big.mark = ","), " bytes)")
  } else {
    warning("  MISSING: ", f)
  }
}

# ── 3. Print summary ────────────────────────────────────────────────────────

flag_qc <- data.table::fread(file.path(out_path, "flag_qc.csv"))
message("\nProbe QC Summary:")
message("  Genes assessed: ", nrow(flag_qc))
message("  PASS: ", sum(flag_qc$qc_flag == "PASS"))
message("  WARN: ", sum(flag_qc$qc_flag == "WARN"))
message("  FAIL: ", sum(flag_qc$qc_flag == "FAIL"))

coverage <- data.table::fread(file.path(out_path, "celltype_coverage.csv"))
message("\nCell Type Coverage Risk:")
for (r in c("HIGH RISK", "MODERATE", "OK")) {
  cts <- coverage$cell_type[coverage$risk == r]
  if (length(cts) > 0) {
    message("  ", r, " (", length(cts), "): ", paste(cts, collapse = ", "))
  }
}

message("\nDone. Review report: ", file.path(out_path, "explore_probe_qc.html"))
log_session()
