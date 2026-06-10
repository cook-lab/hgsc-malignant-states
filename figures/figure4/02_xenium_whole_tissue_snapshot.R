#!/usr/bin/env Rscript
# ============================================================================
# Figure 4C / SF10A — Whole-tissue cell-type spatial maps (06f labels)
# ----------------------------------------------------------------------------
# PURPOSE: Spatial cell-type maps for whole-tissue samples (and the side-by-side
#   TMA), using the 06f polarization-based epithelial reclassification
#   (SecA / Intermediate / SecB). Default = OTB_2384 (Fig 4C); `all` renders all
#   8 published whole tissues + TMA (SF10A).
#
# INPUTS:
#   - SFE objects : <data_root>/<sfe_dir>/sfe_<sample>  (coords + non-epi labels)
#   - output_root/06f_reclassification_polarization/reclassified_xenium_scores.csv
#   - ref_palette / load_sfe from spatial/00_setup/00_setup.R
#
# OUTPUTS:
#   - figures_dir/figure4/xenium_whole_tissue_snapshot_<sample>.{png,svg}
#   - figures_dir/figure4/xenium_whole_tissue_snapshot_TMA.{png,svg}
#
# MANUSCRIPT PANEL(S): Fig 4C (OTB_2384), SF10A (all whole tissues).
#
# RUNTIME TIER: heavy (loads each whole-tissue SFE; rasterized rendering).
#
# CLI: no args / "trial" = OTB_2384 only; "all" = all WT + TMA; or sample names.
# ============================================================================

# --- Shared spatial setup (provides config, load_sfe, ref_palette, seed) -----
.fig_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(.fig_dir) || !nzchar(.fig_dir)) .fig_dir <- "figures/figure4"
.setup_path <- normalizePath(file.path(.fig_dir, "..", "..", "spatial", "00_setup", "00_setup.R"),
                             mustWork = FALSE)
if (!file.exists(.setup_path)) .setup_path <- "spatial/00_setup/00_setup.R"
source(.setup_path)

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

fig_dir <- cfg_path("figures_dir", "figure4")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
f06f <- cfg_path("data_root", "2026_final_xenium_analysis", "output", "06f_reclassification_polarization",
                 "reclassified_xenium_scores.csv")
stopifnot(file.exists(f06f))

recl <- fread(f06f, select = c("sample", "barcode_orig", "cell_label_06f"))
setnames(recl, c("sample", "barcode_orig", "cell_label_06f"))

args <- commandArgs(trailingOnly = TRUE)
wt_all <- sfe_names_wt   # published 8 whole tissues, from 00_setup.R cohort

if (length(args) == 0 || identical(args[1], "trial")) {
  wt_names <- "sfe_OTB_2384"; run_tma <- FALSE
  cat("Mode: trial (OTB_2384 only, no TMA)\n")
} else if (identical(args[1], "all")) {
  wt_names <- wt_all; run_tma <- TRUE
  cat("Mode: all WT + TMA\n")
} else {
  wt_names <- paste0("sfe_", args)
  run_tma  <- any(args == "TMA")
  wt_names <- setdiff(wt_names, "sfe_TMA")
  cat(sprintf("Mode: custom (%d WT, tma=%s)\n", length(wt_names), run_tma))
}

spatial_point_layer <- function(mapping = aes(), ...) {
  ggrastr::geom_point_rast(mapping = mapping, shape = 16, size = 0.1,
                           alpha = 0.8, raster.dpi = 300, dev = "ragg", ...)
}

plot_tissue <- function(cd, title) {
  p <- ggplot(cd, aes(x = x, y = y)) +
    spatial_point_layer(aes(color = cell_label)) +
    scale_color_manual(values = ref_palette, drop = TRUE) +
    coord_fixed() +
    theme_void(base_size = 8) +
    theme(legend.position = "none",
          plot.title = element_text(size = 9, hjust = 0.5)) +
    labs(title = title)
  xr <- diff(range(cd$x)); yr <- diff(range(cd$y))
  h <- 3; w <- h * (xr / yr)
  list(plot = p, w = w, h = h)
}

override_with_06f <- function(cd, sample_name) {
  sub <- recl[sample == sample_name]
  if (nrow(sub) == 0) {
    warning("No 06f rows for ", sample_name, "; cell_label unchanged")
    return(cd)
  }
  m <- match(cd$cell_id, sub$barcode_orig)
  n_hit <- sum(!is.na(m))
  cat(sprintf("  06f override: %d / %d cells relabeled (%s)\n", n_hit, nrow(cd), sample_name))
  cd[!is.na(m), cell_label := sub$cell_label_06f[m[!is.na(m)]]]
  # Deposited 06f cache still carries the legacy literal "Transitioning
  # epithelium"; rename to the standardized "Intermediate epithelium" so the
  # value matches ref_palette (otherwise these cells render with no colour).
  cd[cell_label == "Transitioning epithelium", cell_label := "Intermediate epithelium"]
  cd
}

# -------------------- Whole-tissue samples --------------------
for (nm in wt_names) {
  cat(sprintf("\n=== %s ===\n", nm))
  sfe <- load_sfe(nm)
  cd  <- as.data.table(as.data.frame(colData(sfe)))
  if (!"cell_id" %in% names(cd)) cd[, cell_id := colnames(sfe)]
  coords <- spatialCoords(sfe)
  cd[, x := coords[, 1]]; cd[, y := coords[, 2]]
  cd <- override_with_06f(cd, nm)
  set.seed(CFG$seed)
  cd <- cd[sample(.N)]
  label <- sub("^sfe_", "", nm)
  pk <- plot_tissue(cd, label)
  out_png <- file.path(fig_dir, sprintf("xenium_whole_tissue_snapshot_%s.png", label))
  out_svg <- sub("\\.png$", ".svg", out_png)
  ggsave(out_png, pk$plot, width = pk$w, height = pk$h, dpi = 300)
  ggsave(out_svg, pk$plot, width = pk$w, height = pk$h, dpi = 300)
  cat("  Saved:", basename(out_png), "\n")
  rm(sfe, cd); gc(verbose = FALSE)
}

# -------------------- TMA (side-by-side) --------------------
if (run_tma) {
  cat("\n=== TMA ===\n")
  sfe_tma <- load_sfe("sfe_tma_filtered")
  cd <- as.data.table(as.data.frame(colData(sfe_tma)))
  if (!"cell_id" %in% names(cd)) cd[, cell_id := colnames(sfe_tma)]
  coords <- spatialCoords(sfe_tma)
  cd[, x := coords[, 1]]; cd[, y := coords[, 2]]
  cd <- override_with_06f(cd, "sfe_tma")
  gap <- 1500
  tma1_xmax <- cd[sample_id == "TMA_1", max(x)]
  tma2_xmin <- cd[sample_id == "TMA_2", min(x)]
  cd[sample_id == "TMA_2", x := x + (tma1_xmax - tma2_xmin + gap)]
  set.seed(CFG$seed); cd <- cd[sample(.N)]
  pk <- plot_tissue(cd, "TMA")
  out_png <- file.path(fig_dir, "xenium_whole_tissue_snapshot_TMA.png")
  out_svg <- sub("\\.png$", ".svg", out_png)
  ggsave(out_png, pk$plot, width = pk$w, height = pk$h, dpi = 300)
  ggsave(out_svg, pk$plot, width = pk$w, height = pk$h, dpi = 300)
  cat("  Saved:", basename(out_png), "\n")
}

cat("\nDone.\n")
