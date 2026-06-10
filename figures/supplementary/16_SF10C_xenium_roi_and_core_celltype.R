#!/usr/bin/env Rscript
# =============================================================================
# SF10C — Xenium ROI / per-core cell-type maps (06f reclassification)
# -----------------------------------------------------------------------------
# Purpose:  Voyager cellSeg-polygon cell-type maps using the current 06f
#           polarization-based reclassification. Styling matches Fig 4D
#           (plotSpatialFeature, aes_use="fill", linewidth=0.15, color="grey50",
#           ref_palette, theme_void).
#
#           SF10C uses the per-core TMA gallery (mode "cores"). The "trial"/"rois"
#           modes render OTB_2384 ROIs (shared utility; not the SF10C panel).
#
# INPUTS:
#   cfg_obj("sfe_dir")/sfe_OTB_2384         (ROIs A/B/C)        [modes trial/rois]
#   cfg_obj("sfe_tma_filtered")             (per-core gallery)  [mode cores]
#   output_root/06f_reclassification_polarization/reclassified_xenium_scores.csv
#       (06f epithelial polarization override of cell_label)
# OUTPUTS:
#   output_root/figures/supplementary/SF10C_xenium_by_core_celltype/core_<id>_celltype.png
#   output_root/figures/supplementary/SF10C_xenium_roi_<sample>_<roi>_celltype.{png,svg}
# MANUSCRIPT PANEL(S):  SF10C (per-core FTE/TMA cell-type gallery).
# RUNTIME TIER:  heavy (per-core SFE subsetting + polygon rendering across cores).
# =============================================================================

# Config + shared spatial setup (2 levels deep -> ../../).
source(file.path("..", "..", "config", "config.R"))
source(file.path("..", "..", "spatial", "00_setup", "00_setup.R"))  # load_sfe, ref_palette

suppressPackageStartupMessages({
  library(Voyager); library(data.table); library(ggplot2)
})

set.seed(CFG$seed)

OUT_DIR  <- cfg_path("output_root", "figures", "supplementary")
core_dir <- file.path(OUT_DIR, "SF10C_xenium_by_core_celltype")

f06f <- cfg_path("data_root", "2026_final_xenium_analysis", "output",
                 "06f_reclassification_polarization",
                 "reclassified_xenium_scores.csv")
stopifnot(file.exists(f06f))
recl <- fread(f06f, select = c("sample", "barcode_orig", "cell_label_06f"))
# Deposited 06f cache still carries the legacy literal "Transitioning
# epithelium"; rename to the standardized "Intermediate epithelium" so the
# value (later written onto sfe$cell_label) matches ref_palette; otherwise
# these cells render with no fill colour in plotSpatialFeature.
recl[cell_label_06f == "Transitioning epithelium", cell_label_06f := "Intermediate epithelium"]

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) == 0) "cores" else args[1]
cat("Mode:", mode, "\n")

# ROI bounding boxes (OTB_2384; per Fig 4 ROI definitions)
ROI_BBOX <- list(
  ROI_A = c(xmin = 3200, xmax = 4000, ymin = -7200, ymax = -6400),
  ROI_B = c(xmin = 3400, xmax = 4200, ymin = -6800, ymax = -6000),
  ROI_C = c(xmin = 2800, xmax = 4000, ymin = -7800, ymax = -6600)
)

# Override cell_label in-place on an SFE using 06f
override_with_06f <- function(sfe, sample_name) {
  sub <- recl[sample == sample_name]
  if (nrow(sub) == 0) { warning("no 06f rows for ", sample_name); return(sfe) }
  m <- match(colnames(sfe), sub$barcode_orig)
  hit <- !is.na(m)
  cat(sprintf("  06f override: %d / %d cells relabeled\n", sum(hit), length(hit)))
  new_lab <- as.character(sfe$cell_label)
  new_lab[hit] <- sub$cell_label_06f[m[hit]]
  sfe$cell_label <- new_lab
  sfe
}

# Cell-type polygon plot matching Fig 4D
plot_ct_polygons <- function(sfe, show_legend = TRUE, base_size = 10) {
  p <- plotSpatialFeature(sfe, "cell_label", colGeometryName = "cellSeg",
                          aes_use = "fill", linewidth = 0.15, color = "grey50") +
    scale_fill_manual(values = ref_palette, name = "Cell type", drop = TRUE) +
    theme_void(base_size = base_size)
  if (show_legend) {
    p <- p + theme(legend.position = "right",
                   legend.text = element_text(size = 8),
                   legend.title = element_text(size = 9),
                   legend.key.size = unit(0.5, "lines")) +
      guides(fill = guide_legend(override.aes = list(linewidth = 0)))
  } else {
    p <- p + theme(legend.position = "none")
  }
  p
}

render_roi <- function(sfe_full, roi_name, sample_name) {
  bb <- ROI_BBOX[[roi_name]]
  co <- spatialCoords(sfe_full)
  in_roi <- co[, 1] >= bb["xmin"] & co[, 1] <= bb["xmax"] &
            co[, 2] >= bb["ymin"] & co[, 2] <= bb["ymax"]
  sfe_roi <- sfe_full[, in_roi]
  cat(sprintf("  %s: %d cells in ROI\n", roi_name, ncol(sfe_roi)))
  p <- plot_ct_polygons(sfe_roi, show_legend = TRUE)
  base <- sprintf("SF10C_xenium_roi_%s_%s_celltype", sample_name, roi_name)
  ggsave(file.path(OUT_DIR, paste0(base, ".png")), p, width = 8, height = 7, dpi = 300)
  ggsave(file.path(OUT_DIR, paste0(base, ".svg")), p, width = 8, height = 7, dpi = 300)
  cat("  Saved:", paste0(base, ".png"), "\n")
}

# ---------------------------------------------------------------------------
# Mode: trial | rois  (OTB_2384 ROIs — shared utility, not the SF10C panel)
# ---------------------------------------------------------------------------
if (mode %in% c("trial", "rois")) {
  cat("\n=== OTB_2384 ROI render ===\n")
  sfe <- load_sfe("sfe_OTB_2384")
  sfe <- override_with_06f(sfe, "sfe_OTB_2384")
  rois <- if (mode == "trial") "ROI_C" else c("ROI_A", "ROI_B", "ROI_C")
  for (r in rois) render_roi(sfe, r, "OTB_2384")
  rm(sfe); gc(verbose = FALSE)
}

# ---------------------------------------------------------------------------
# Mode: cores  (per-core gallery — the SF10C panel)
# ---------------------------------------------------------------------------
if (mode == "cores") {
  dir.create(core_dir, showWarnings = FALSE, recursive = TRUE)
  cat("\n=== Per-core gallery ->", core_dir, "===\n")
  sfe_tma <- load_sfe("sfe_tma_filtered")
  sfe_tma <- override_with_06f(sfe_tma, "sfe_tma")

  cd <- as.data.table(as.data.frame(colData(sfe_tma)))
  if (!"core_id" %in% names(cd)) stop("core_id column missing from TMA colData")

  core_ids <- sort(unique(cd$core_id))
  core_ids <- core_ids[!is.na(core_ids) & core_ids != "Off core"]
  cat(sprintf("  %d cores to render\n", length(core_ids)))

  W <- 3; H <- 3; DPI <- 150
  for (cid in core_ids) {
    idx <- which(cd$core_id == cid)
    if (length(idx) < 50) { cat("  skip", cid, "(", length(idx), "cells)\n"); next }
    sfe_c <- sfe_tma[, idx]
    p <- plot_ct_polygons(sfe_c, show_legend = FALSE, base_size = 8) +
      labs(title = as.character(cid)) +
      theme(plot.title = element_text(size = 7, hjust = 0.5))
    out <- file.path(core_dir, sprintf("core_%s_celltype.png", cid))
    ggsave(out, p, width = W, height = H, dpi = DPI)
  }
  leg <- plot_ct_polygons(sfe_tma[, 1:min(1000, ncol(sfe_tma))], show_legend = TRUE)
  ggsave(file.path(core_dir, "_legend.png"), leg, width = 4, height = 6, dpi = 200)
  cat("  Legend saved to _legend.png\n")
  rm(sfe_tma); gc(verbose = FALSE)
}

cat("\nDone.\n")
