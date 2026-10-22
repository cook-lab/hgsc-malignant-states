# Trimmed copy (2026-10-21) of 2026_final_xenium_analysis/scripts/00_setup.R, the script that ran (sourced by the
# probe-QC and clinical notebooks). Kept: libraries, project paths, output folders, load_sfe() and log_session().
# Removed: the sample sheet, plot theme, palettes, plot helpers, save_sfe() and the startup message.
# The untrimmed script is in the lab project folder.
# ============================================================================
# 00_setup.R — Shared setup for 2026 Final Xenium Analysis
# ============================================================================
# Source this at the top of every analysis script:
#   source("scripts/00_setup.R")
# ============================================================================

# --- Libraries ---------------------------------------------------------------

# Core spatial / single-cell
library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scuttle)
library(scran)
library(HDF5Array)

# Data handling
library(arrow)
library(data.table)

# Plotting
library(ggplot2)
library(ggrastr)
library(patchwork)
library(scico)
library(viridis)
library(ComplexHeatmap)

# --- Project paths -----------------------------------------------------------

proj_dir  <- "/Volumes/CookLab/Sarah/2026_final_xenium_analysis"
data_dir  <- file.path(proj_dir, "data")
meta_dir  <- file.path(proj_dir, "metadata")
out_dir   <- file.path(proj_dir, "output")
fig_dir   <- file.path(proj_dir, "output", "figures")
sfe_dir   <- file.path(proj_dir, "output", "sfe")

# Create output directories if needed
for (d in c(out_dir, fig_dir, sfe_dir,
            file.path(out_dir, "03_04_qc"),
            file.path(out_dir, "05_probe_qc"),
            file.path(out_dir, "06_annotation"),
            file.path(out_dir, "06b_adaptive_secretory"),
            file.path(out_dir, "06c_annotation_evaluation"),
            file.path(out_dir, "07_core_qc"),
            file.path(out_dir, "08_colocalization"),
            file.path(out_dir, "09_neighborhood"),
            file.path(out_dir, "10_clinical"),
            file.path(out_dir, "11_visium"),
            file.path(out_dir, "9b_scoring"),
            file.path(out_dir, "9c_neighborhood_network"),
            file.path(out_dir, "06d_annotation_noBCAM"),
            file.path(out_dir, "13_macrophage_niche"))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# --- Utility functions -------------------------------------------------------

# Load an HDF5-backed SFE saved with saveHDF5SummarizedExperiment
load_sfe <- function(name) {
  sfe_path <- file.path(sfe_dir, name)
  if (!dir.exists(sfe_path)) stop("SFE directory not found: ", sfe_path)
  loadHDF5SummarizedExperiment(dir = sfe_path)
}

# --- Session info helper -----------------------------------------------------

log_session <- function() {
  si <- sessionInfo()
  cat("R version:", si$R.version$version.string, "\n")
  pkgs <- c(si$otherPkgs, si$loadedOnly)
  key_pkgs <- c("SpatialFeatureExperiment", "SingleCellExperiment", "scuttle",
                "scran", "ggplot2", "ComplexHeatmap", "SingleR")
  for (p in key_pkgs) {
    if (p %in% names(pkgs)) cat(" ", p, ":", pkgs[[p]]$Version, "\n")
  }
}
