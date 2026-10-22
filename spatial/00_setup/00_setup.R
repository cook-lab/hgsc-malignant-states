# ============================================================================
# 00_setup.R — Shared setup for the HGSC Xenium spatial pipeline
# ============================================================================
# PURPOSE: Single shared setup sourced at the top of every spatial/ script and
#   of the Xenium figure R scripts. Loads libraries, the central config
#   (paths/seed/cohort), the colour palettes / ggplot theme, the SecA/SecB
#   genes, and SFE load/save helpers.
#
# INPUTS:
#   - config/config.{R,yml}            (paths, seed, cohort definition)
#   - shared/signatures.yml            (SecA / SecB genes)
#
# OUTPUTS:
#   - none (defines objects + helpers in the calling environment; creates the
#     output_root subdirectories used by downstream stages)
#
# MANUSCRIPT PANEL(S): none directly — supports the entire Xenium backend
#   (Figs 4–7, Figures S11–S18, Tables S7–S11).
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
library(viridis)
library(ComplexHeatmap)

# --- Project paths (resolved from central config) ----------------------------
# data_root holds the deposited input objects; output_root receives all
# pipeline outputs. Both are overridable via DATA_ROOT / OUTPUT_ROOT env vars.

xen_root  <- file.path(path.expand(CFG$paths$data_root), "2026_final_xenium_analysis")
data_dir  <- file.path(xen_root, "data")

out_dir   <- path.expand(CFG$paths$output_root)
fig_dir   <- file.path(out_dir, "figures")
# SFE input directory: RERUN_SFE_DIR if set (the final run: OUTPUT_ROOT/sfe);
# otherwise config objects.sfe_dir under data_root, which holds the SFEs from
# before the final re-run.
sfe_dir   <- Sys.getenv("RERUN_SFE_DIR",
                        file.path(path.expand(CFG$paths$data_root), CFG$objects$sfe_dir))
# Regenerated SFEs (downstream re-runs that add colData) are written HERE, under
# output_root — never into the read-only deposited data_root.
sfe_out_dir <- file.path(out_dir, "sfe")

# Create output directories if needed
for (d in c(out_dir, fig_dir, sfe_out_dir,
            file.path(out_dir, "03_04_qc"),
            file.path(out_dir, "05_probe_qc"),
            file.path(out_dir, "06d_annotation_noBCAM"),
            file.path(out_dir, "06f_reclassification_polarization"),
            file.path(out_dir, "06g_clean_split"),
            file.path(out_dir, "07_core_qc"))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# --- Cohort (whole-tissue arm and arrays) ------------------------------------
# The 8 whole tissues of the paper (config cohort.whole_tissue, WT01–WT08). FTE
# whole-tissue sections are not part of the whole-tissue arm; the 15 FTE cores
# lie on the TMA A slides.
whole_tissue_samples <- CFG$cohort$whole_tissue
sfe_names_wt         <- paste0("sfe_", whole_tissue_samples)
# Arrays in config cohort.tma_new: TMA B blocks G and H (separate objects,
# de-arrayed with the interactive editor).
sfe_names_tma_new    <- if (!is.null(CFG$cohort$tma_new)) paste0("sfe_", CFG$cohort$tma_new) else character(0)
# Full SFE set looped over by the backend steps: TMA A + cohort.tma_new + whole tissues.
sfe_names_all        <- c("sfe_tma", sfe_names_tma_new, sfe_names_wt)

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

# --- svg_inline_fill(): make svglite SVGs editable in Affinity / Illustrator -
# svglite declares element defaults (notably `fill: none` for stroked geoms) in
# a <style type='text/css'> block scoped to the `.svglite` class, NOT as inline
# presentation attributes. Affinity Designer and Adobe Illustrator ignore SVG
# <style> blocks, so stroked geoms (polyline/line/path: KM step curves, forest
# error bars, gridlines) fall back to the SVG default fill:black and render as
# SOLID BLOCKS, with the inline stroke barely visible. This rewrites the CSS
# default as an INLINE `fill: none;` prepended to each such element's style; any
# genuine inline fill later in the same style string still wins (last-declared
# property wins within one inline style), so points / filled shapes are
# unchanged. Call once on each *.svg written via ggsave()/svglite.
svg_inline_fill <- function(path) {
  if (!file.exists(path)) return(invisible(FALSE))
  x <- paste(readLines(path, warn = FALSE), collapse = "\n")
  # svglite emits the element defaults `fill: none; stroke: #000000;` ONLY inside
  # a `.svglite {...}` CSS block, which Affinity Designer / Adobe Illustrator
  # ignore. Without those defaults inlined:
  #   * black axis lines + ticks (no inline stroke colour) -> stroke:none -> AXES
  #     DISAPPEAR;
  #   * coloured step curves / gridlines (a <polyline>, no inline fill) -> SVG
  #     default fill:black -> the area under the curve fills as a SOLID BLOCK.
  # We inline each missing default ONCE, and ONLY when the element does not
  # already declare it. (The previous version prepended stroke:#000000 to every
  # stroked geom unconditionally, so coloured lines carried TWO stroke: declns —
  # Affinity drew both, painting a faint #000/near-white line THROUGH every
  # coloured line. Adding only-when-absent removes that artifact.)
  # `fixed = TRUE` on the tests distinguishes the `stroke:` colour property from
  # `stroke-width:` / `stroke-linecap:` etc. `line/polyline/path/polygon` only;
  # `rect`/`circle` (filled area shapes) always specify their own fill/stroke.
  # Line-like geoms: inject BOTH missing defaults (a stroked line with no fill
  # would fill as a block; a black line with no inline stroke would vanish).
  for (tag in c("polyline", "line", "path", "polygon")) {
    patt <- sprintf("<%s\\b[^>]*?style='[^']*'", tag)
    m <- gregexpr(patt, x, perl = TRUE)
    old <- regmatches(x, m)[[1]]
    if (!length(old)) next
    new <- vapply(old, function(el) {
      sty <- sub(".*style='([^']*)'.*", "\\1", el)
      add <- ""
      if (!grepl("fill:",   sty, fixed = TRUE)) add <- paste0(add, "fill: none; ")
      if (!grepl("stroke:", sty, fixed = TRUE)) add <- paste0(add, "stroke: #000000; ")
      if (nzchar(add)) el <- sub("style='", paste0("style='", add), el, fixed = TRUE)
      el
    }, character(1), USE.NAMES = FALSE)
    regmatches(x, m) <- list(new)
  }
  # <circle> (geom_point): inject ONLY `fill: none` when absent — this is the
  # open-marker case (shape 1/hollow rings: stroke but no fill), which otherwise
  # takes the SVG default fill:black in Affinity and renders as SOLID BLACK dots.
  # NEVER add a stroke to a circle: a filled dot legitimately has no stroke, and
  # a black stroke would ring every point.
  {
    m <- gregexpr("<circle\\b[^>]*?style='[^']*'", x, perl = TRUE)
    old <- regmatches(x, m)[[1]]
    if (length(old)) {
      new <- vapply(old, function(el) {
        sty <- sub(".*style='([^']*)'.*", "\\1", el)
        if (!grepl("fill:", sty, fixed = TRUE))
          el <- sub("style='", "style='fill: none; ", el, fixed = TRUE)
        el
      }, character(1), USE.NAMES = FALSE)
      regmatches(x, m) <- list(new)
    }
  }
  writeLines(x, path)
  invisible(TRUE)
}

# --- Color palettes ----------------------------------------------------------

# Reference palette: the 18 Xenium cell types (SingleR labels, with Secretory
# epithelium split into SecA / Intermediate / SecB by polarization).
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

# Load an HDF5-backed SFE. Prefer the copy under output_root/sfe; fall back to
# sfe_dir (RERUN_SFE_DIR, else config objects.sfe_dir under data_root). SFEs are
# only ever written under output_root/sfe, so data_root inputs are never mutated.
load_sfe <- function(name) {
  for (base in c(sfe_out_dir, sfe_dir)) {
    sfe_path <- file.path(base, name)
    if (dir.exists(sfe_path)) return(loadHDF5SummarizedExperiment(dir = sfe_path))
  }
  stop("SFE directory not found in output_root/sfe or data_root: ", name)
}

# Save SFE as HDF5-backed under output_root/sfe, using a _v2 swap for crash safety.
# Writes/unlinks ONLY under output_root — the deposited data_root inputs are never
# mutated.
save_sfe <- function(sfe, name) {
  if (!dir.exists(sfe_out_dir)) dir.create(sfe_out_dir, recursive = TRUE)
  sfe_path     <- file.path(sfe_out_dir, name)
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
message("Whole-tissue cohort (core ", length(whole_tissue_samples), "): ",
        paste(whole_tissue_samples, collapse = ", "))
message("Setup loaded successfully.")
