# ============================================================================
# 00_setup.R — Shared setup for the HGSC Xenium spatial pipeline
# ============================================================================
# PURPOSE: Single shared setup sourced at the top of every spatial/ script.
#   Loads libraries, the central config (paths/seed/cohort), the canonical
#   color palettes / ggplot theme, and SFE load/save helpers.
#
# INPUTS:
#   - config/config.{R,yml}            (paths, seed, cohort definition)
#   - <data_root>/2026_final_xenium_analysis/metadata/samples.csv  (sample table)
#
# OUTPUTS:
#   - none (defines objects + helpers in the calling environment; creates the
#     output_root subdirectories used by downstream stages)
#
# MANUSCRIPT PANEL(S): none directly — supports the entire Xenium backend
#   (Fig 4–7, SF10–SF14).
#
# RUNTIME TIER: fast
#
# USAGE (from a stage script two levels deep, e.g. spatial/02_qc/):
#   source(file.path(dirname(dirname(sys.frame(1)$ofile)), "00_setup", "00_setup.R"))
#   # or, when run with Rscript from the repo root:
#   source("spatial/00_setup/00_setup.R")
# ============================================================================

# --- Central config ----------------------------------------------------------
# This file lives at spatial/00_setup/00_setup.R (2 levels under the repo root),
# so config/config.R is two directories up.

.setup_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(.setup_dir) || !nzchar(.setup_dir)) .setup_dir <- "spatial/00_setup"
.config_path <- normalizePath(file.path(.setup_dir, "..", "..", "config", "config.R"),
                              mustWork = FALSE)
if (!file.exists(.config_path)) .config_path <- "config/config.R"
source(.config_path)

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

# --- Project paths (resolved from central config) ----------------------------
# data_root holds the deposited input objects; output_root receives all
# pipeline outputs. Both are overridable via DATA_ROOT / OUTPUT_ROOT env vars.

xen_root  <- file.path(path.expand(CFG$paths$data_root), "2026_final_xenium_analysis")
data_dir  <- file.path(xen_root, "data")
meta_dir  <- file.path(xen_root, "metadata")

out_dir   <- path.expand(CFG$paths$output_root)
fig_dir   <- file.path(out_dir, "figures")
# Canonical SFE directory (config entry-point objects live here).
sfe_dir   <- file.path(path.expand(CFG$paths$data_root), CFG$objects$sfe_dir)

# Create output directories if needed
for (d in c(out_dir, fig_dir, sfe_dir,
            file.path(out_dir, "03_04_qc"),
            file.path(out_dir, "05_probe_qc"),
            file.path(out_dir, "06_annotation"),
            file.path(out_dir, "06b_adaptive_secretory_noBCAM"),
            file.path(out_dir, "06d_annotation_noBCAM"),
            file.path(out_dir, "06f_reclassification_polarization"),
            file.path(out_dir, "06g_clean_split"),
            file.path(out_dir, "07_core_qc"))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# --- Cohort (published whole-tissue arm) -------------------------------------
# Pinned to the published 8 whole tissues. The 2 FTE whole-tissue samples
# (FT1-1, EAOC-1-FTE) are EXCLUDED from the whole-tissue arm (cohort drift fix —
# see docs/REPRODUCIBILITY.md). FTE TMA cores (n=15) remain in the TMA.
whole_tissue_samples <- CFG$cohort$whole_tissue
sfe_names_wt         <- paste0("sfe_", whole_tissue_samples)
# Full SFE set used downstream: merged TMA + the 8 published whole tissues.
sfe_names_all        <- c("sfe_tma", sfe_names_wt)

# --- Sample metadata ---------------------------------------------------------

samples_csv <- file.path(meta_dir, "samples.csv")
if (file.exists(samples_csv)) {
  samples <- read.csv(samples_csv, stringsAsFactors = FALSE)
  xenium_samples <- samples[samples$platform == "xenium", ]
} else {
  samples <- NULL
  xenium_samples <- NULL
}

# --- theme_lab() — Cook Lab style guide v0.6 ---------------------------------

theme_lab <- function(base_size = 8, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      text             = element_text(color = "black"),
      axis.text        = element_text(size = rel(0.9), color = "black"),
      axis.title       = element_text(size = rel(1.05)),
      axis.title.x     = element_text(margin = margin(t = 6)),
      axis.title.y     = element_text(margin = margin(r = 6), angle = 90),
      axis.line        = element_line(color = "black", linewidth = 0.5),
      legend.background = element_blank(),
      legend.key       = element_blank(),
      legend.text      = element_text(size = rel(0.85)),
      legend.title     = element_text(size = rel(0.9)),
      legend.position  = "right",
      legend.key.size  = unit(0.8, "lines"),
      panel.background = element_blank(),
      panel.border     = element_blank(),
      panel.grid       = element_blank(),
      strip.background = element_blank(),
      strip.text       = element_text(size = rel(1), margin = margin(b = 4)),
      plot.margin      = margin(8, 8, 8, 8),
      plot.title       = element_text(size = rel(1.15), margin = margin(b = 8)),
      plot.subtitle    = element_text(size = rel(0.9), color = "gray40",
                                      margin = margin(b = 8))
    )
}

# --- Color palettes ----------------------------------------------------------

# Cell type palette v1.2 (lab standard)
celltype_palette <- c(
  "Epi_Secretory"       = "#E6A141",
  "Epi_Secretory_Polar" = "#9A7D55",
  "Epi_Ciliated"        = "#E07850",
  "Mesothelial"         = "#D4A574",
  "Fibroblast"          = "#DDD5CA",
  "Smooth_Muscle"       = "#D14E6C",
  "Pericyte"            = "#B87A7A",
  "Endothelial"         = "#7D4E4E",
  "T_cell"              = "#87CEFA",
  "NK_cell"             = "#56AFC4",
  "B_cell"              = "#5665B6",
  "Plasma_cell"         = "#8A5DAF",
  "Macrophage"          = "#8FBC8F",
  "DC"                  = "#2E8B57",
  "Neutrophil"          = "#6B8E23",
  "Mast"                = "#8B9B6B",
  "Erythrocyte"         = "#CD5C5C",
  "Other"               = "#A0A0A0"
)

# Reference palette: 18 cell_label types (canonical annotation palette).
# Epithelial polarization label standardized: "Intermediate" (was "Transitioning").
ref_palette <- c(
  "SecA epithelium"                      = "#E6A141",
  "Intermediate epithelium"              = "#C08E48",
  "SecB epithelium"                      = "#9A7D55",
  "Ciliated epithelium"                  = "#E07850",
  "Mesothelial"                          = "#D4A574",
  "Fibroblast"                           = "#DDD5CA",
  "Smooth muscle"                        = "#D14E6C",
  "Pericyte"                             = "#B87A7A",
  "Endothelial"                          = "#7D4E4E",
  "T cell"                               = "#87CEFA",
  "NK cell"                              = "#56AFC4",
  "B cell"                               = "#5665B6",
  "Plasma cell"                          = "#8A5DAF",
  "Macrophage"                           = "#8FBC8F",
  "Conventional dendritic cell"          = "#2E8B57",
  "Neutrophil"                           = "#6B8E23",
  "Mast cell"                            = "#8B9B6B",
  "Plasmacytoid dendritic cell"          = "#A0A0A0"
)

# Neighborhood names and palette (10 neighborhoods from k=10 with 06f
# polarization-based reclassification; identities from
# output/09_neighborhood/neighborhood_composition.csv).
nb_names <- c(
  nb_1  = "SecA-dominant epithelium",
  nb_2  = "SecA-stroma border",
  nb_3  = "Fibroblast-rich stroma",
  nb_4  = "Immune-rich niche",
  nb_5  = "Intermediate-SecB epithelium",
  nb_6  = "Early intermediate epithelium",
  nb_7  = "SecA-mixed epithelium",
  nb_8  = "Mesothelial niche",
  nb_9  = "Stromal-vascular niche",
  nb_10 = "Ciliated niche"
)

nb_palette <- c(
  "SecA-dominant epithelium"      = "#E6A141",
  "SecA-stroma border"            = "#F2C98C",
  "Fibroblast-rich stroma"        = "#DDD5CA",
  "Immune-rich niche"             = "#87CEFA",
  "Intermediate-SecB epithelium"  = "#9A7D55",
  "Early intermediate epithelium" = "#C08E48",
  "SecA-mixed epithelium"         = "#F5D78E",
  "Mesothelial niche"             = "#D4A574",
  "Stromal-vascular niche"        = "#B87A7A",
  "Ciliated niche"                = "#E07850"
)

# Okabe-Ito for general discrete (<=8 categories)
okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
)

# Expression color scales
expr_umap    <- c("lightgrey", rev(viridis::magma(100)))
expr_spatial <- scico(30, palette = "lapaz", direction = -1)

# --- SecA/SecB signatures (shared, noBCAM 7-gene set) ------------------------
# Single source of truth: shared/signatures.yml. Loaded here so every spatial
# script uses identical gene lists (no inlined divergent copies).

.sig_path <- normalizePath(file.path(.setup_dir, "..", "..", "shared", "signatures.yml"),
                           mustWork = FALSE)
if (!file.exists(.sig_path)) .sig_path <- "shared/signatures.yml"
SIGNATURES   <- yaml::read_yaml(.sig_path)
SECA_GENES   <- SIGNATURES$SecA   # MECOM FBXO21 LGR5 LPAR3 PBX1 SOX17 RCN2 (noBCAM)
SECB_GENES   <- SIGNATURES$SecB   # KRT17 KRT19 KRT7 LCN2 PRSS22 SLPI TACSTD2

# --- Spatial plot defaults ---------------------------------------------------

# Standard rasterized point layer for spatial / UMAP plots
spatial_point_layer <- function(mapping = aes(), ...) {
  ggrastr::geom_point_rast(
    mapping  = mapping,
    shape    = 16,
    size     = 0.1,
    alpha    = 0.8,
    raster.dpi = 300,
    ...
  )
}

# --- Utility functions -------------------------------------------------------

# Load an HDF5-backed SFE saved with saveHDF5SummarizedExperiment
load_sfe <- function(name) {
  sfe_path <- file.path(sfe_dir, name)
  if (!dir.exists(sfe_path)) stop("SFE directory not found: ", sfe_path)
  loadHDF5SummarizedExperiment(dir = sfe_path)
}

# Save SFE as HDF5-backed, using a _v2 swap for crash safety
save_sfe <- function(sfe, name) {
  sfe_path     <- file.path(sfe_dir, name)
  sfe_path_new <- paste0(sfe_path, "_v2")
  if (dir.exists(sfe_path_new)) unlink(sfe_path_new, recursive = TRUE)
  saveHDF5SummarizedExperiment(sfe, dir = sfe_path_new)
  if (dir.exists(sfe_path)) unlink(sfe_path, recursive = TRUE)
  file.rename(sfe_path_new, sfe_path)
  message("Saved: ", sfe_path)
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

# --- Startup message ---------------------------------------------------------

message("=== HGSC Xenium spatial pipeline ===")
message("data_root:   ", path.expand(CFG$paths$data_root))
message("output_root: ", out_dir)
message("seed:        ", CFG$seed)
message("Whole-tissue cohort (published 8): ", paste(whole_tissue_samples, collapse = ", "))
message("Setup loaded successfully.")
