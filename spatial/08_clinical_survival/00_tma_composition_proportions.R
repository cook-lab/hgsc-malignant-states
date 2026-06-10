# ============================================================================
# 00_tma_composition_proportions.R
# ----------------------------------------------------------------------------
# PURPOSE: Per-TMA-core cell-type composition and wide proportion table. For
#   each core, count cells by cell_label and express each cell type as a percent
#   of that core's cells, then reshape to one row per core (wide). This is the
#   per_core_proportions_wide.csv that the Fig 7C/7D protein-correlation script
#   consumes for per-core SecB% / SecA% and total cell counts.
#
# INPUTS:
#   - cfg_obj("sfe_tma_filtered")  (TMA SFE: colData core_id, patient_id,
#     sample_type, cell_label; loaded via load_sfe)
#
# OUTPUTS:
#   - output/38_FTE_baseline_TMA/composition_per_core.csv   (long: per core x cell_label)
#   - output/38_FTE_baseline_TMA/per_core_proportions_wide.csv
#       columns: core_id, patient_id, sample_type, total, <cell_label %> ...
#
# PIPELINE ROLE: producer for 08_xenium_protein_correlation.R (Fig 7C/7D).
#   Number 00 so it runs before 08.
#
# RUNTIME TIER: moderate (SFE load + reshape).
#
# Ported from 2026_final_xenium_analysis/scripts/38b_TMA_composition_and_proportions.R
#   (composition + per_core_proportions_wide.csv build only). Analytical logic
#   preserved; paths routed through central config, seed from CFG$seed. The
#   "Intermediate epithelium" proportion column is also emitted under its legacy
#   name "Transitioning epithelium" so the Fig 7 consumer resolves unchanged.
# ============================================================================

# --- Config + shared setup (replaces hardcoded /Volumes/CookLab/Sarah paths) ---
here <- local({                       # robust: works via `Rscript <path>`, source(), or from-dir
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa)) dirname(sub("^--file=", "", fa[1]))
  else tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
})
source(file.path(here, "..", "..", "config", "config.R"))   # CFG, cfg_obj, cfg_path
source(file.path(here, "..", "00_setup", "00_setup.R"))      # load_sfe, theme_lab, palettes
set.seed(CFG$seed)

suppressPackageStartupMessages({
  library(data.table)
  library(SpatialFeatureExperiment); library(SummarizedExperiment); library(HDF5Array)
})

message("\n=== 38: TMA per-core composition + wide proportions ===")

OUT_DIR <- cfg_path("output_root", "38_FTE_baseline_TMA")

# --- Load TMA SFE ------------------------------------------------------------
sfe <- load_sfe("sfe_tma_filtered")
cd  <- as.data.table(as.data.frame(
  colData(sfe))[, c("core_id", "patient_id", "sample_type", "cell_label")])
rm(sfe); gc(verbose = FALSE)

# Restrict to the depth-matched FTE / HGSC cohorts (drop NA / unknown).
cd <- cd[sample_type %in% c("fallopian", "tumour")]
cd[, cell_label := as.character(cell_label)]

# --- Per-core composition (long) --------------------------------------------
comp <- cd[, .(N = .N), by = .(core_id, patient_id, sample_type, cell_label)]
comp[, total := sum(N), by = core_id]
comp[, pct := 100 * N / total]
fwrite(comp, file.path(OUT_DIR, "composition_per_core.csv"))
message("[saved] composition_per_core.csv")

# --- Wide per-core proportions ----------------------------------------------
core_pcts <- dcast(comp, core_id + patient_id + sample_type + total ~ cell_label,
                   value.var = "pct", fill = 0)

# Backward-compat: the Fig 7 consumer reads a "Transitioning epithelium" column;
# the repo standardizes this label to "Intermediate epithelium". Mirror it so
# 08_xenium_protein_correlation.R resolves unchanged.
if ("Intermediate epithelium" %in% names(core_pcts) &&
    !("Transitioning epithelium" %in% names(core_pcts))) {
  core_pcts[, `Transitioning epithelium` := `Intermediate epithelium`]
}

fwrite(core_pcts, file.path(OUT_DIR, "per_core_proportions_wide.csv"))
message("[saved] per_core_proportions_wide.csv  (",
        nrow(core_pcts), " cores x ", ncol(core_pcts), " cols)")

message("\n[done]")
log_session()
