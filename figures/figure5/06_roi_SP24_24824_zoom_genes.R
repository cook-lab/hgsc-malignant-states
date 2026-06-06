#!/usr/bin/env Rscript
# ============================================================================
# Figure 5H,5I — SP24_24824 ROI zoom: cell type + gene expression panels
# ----------------------------------------------------------------------------
# PURPOSE
#   Zoomed spatial panels for SP24_24824 ROI 06g slit region. Gene-expression
#   maps (CDH2, CTNNB1, ITGB5, TGM2, MMP7, ICAM1) over cell-segmentation
#   polygons, colorbar outside the panel.
#   Full ROI 06g: x=[7800,9000] y=[-7100,-5900]; zoom slit: x=[7980,8700]
#   y=[-6860,-6140].
#
# INPUTS
#   data_root/2026_final_xenium_analysis/output/sfe/sfe_SP24_24824 (via load_sfe)
#   data_root/2026_final_xenium_analysis/output/06f_reclassification_polarization/
#     reclassified_xenium_scores.csv  (06f cell-label override)
#   Shared helpers: config/config.R, spatial/00_setup/00_setup.R (load_sfe).
#
# OUTPUTS
#   figures_dir/figure5/ROI_figure_5/SP24_24824_roi_<gene>_zoom.{svg,png} (5I)
#
# MANUSCRIPT PANEL(S): Fig 5I (gene maps). LINEAGE attributes the Fig 5H
#   cell-type map to "the same zoom_genes script", but the canonical source
#   file only renders the gene-expression panels below; the 5H cell-type
#   panel is not present in this generator (see manifest note).
# RUNTIME TIER: moderate (SFE load + crop)
# ============================================================================

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

.here     <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))
source(file.path(.here, "..", "..", "spatial", "00_setup", "00_setup.R"))

suppressPackageStartupMessages({
  library(SpatialFeatureExperiment); library(SummarizedExperiment)
  library(data.table); library(ggplot2); library(sf); library(scico)
})

set.seed(CFG$seed)

FIG_DIR <- cfg_path("figures_dir", "figure5", "ROI_figure_5")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
SAMPLE  <- "sfe_SP24_24824"
GENES   <- c("CDH2", "CTNNB1", "ITGB5", "TGM2", "MMP7", "ICAM1")

BB_FULL <- c(xmin = 7800, xmax = 9000, ymin = -7100, ymax = -5900)
BB_ZOOM <- c(xmin = 7980, xmax = 8700, ymin = -6860, ymax = -6140)

GRAD_COLS  <- c("#F6EFE5", "#ECDDD0", "#D89E97", "#A03A4A", "#3A111A")
GRAD_STOPS <- c(0.00,      0.28,      0.55,      0.80,      1.00)

# ── Load SFE + 06f override ─────────────────────────────────────────────
cat("=== Loading", SAMPLE, "===\n")
sfe <- load_sfe(SAMPLE)

f06f <- cfg_path("data_root", "2026_final_xenium_analysis", "output",
                 "06f_reclassification_polarization",
                 "reclassified_xenium_scores.csv")
recl <- fread(f06f, select = c("sample", "barcode_orig", "cell_label_06f"))
m <- match(colnames(sfe), recl[sample == SAMPLE, barcode_orig])
hit <- !is.na(m)
lab <- as.character(sfe$cell_label)
lab[hit] <- recl[sample == SAMPLE, cell_label_06f][m[hit]]
sfe$cell_label <- lab
cat(sprintf("  06f override: %d / %d cells\n", sum(hit), length(hit)))

# ── Subset to full ROI (superset of zoom) ────────────────────────────────
co <- spatialCoords(sfe)
in_roi <- co[, 1] >= BB_FULL["xmin"] & co[, 1] <= BB_FULL["xmax"] &
          co[, 2] >= BB_FULL["ymin"] & co[, 2] <= BB_FULL["ymax"]
sfe_roi <- sfe[, in_roi]
cat(sprintf("  Full ROI cells: %d\n", ncol(sfe_roi)))

cseg <- colGeometry(sfe_roi, "cellSeg")
cseg$cell_label <- as.character(sfe_roi$cell_label)
for (g in GENES) {
  cseg[[g]] <- as.numeric(assay(sfe_roi, "counts")[g, ])
}
cseg_sf <- sf::st_as_sf(cseg)

# ── Gene panels (Fig 5I) ─────────────────────────────────────────────────
for (g in GENES) {
  cat(sprintf("\n--- %s zoom ---\n", g))

  vals <- cseg_sf[[g]]
  qlim <- quantile(vals, c(0.02, 0.98), na.rm = TRUE)
  vmin <- as.numeric(qlim[1]); vmax <- as.numeric(qlim[2])
  if (vmin == vmax) { vmin <- 0; vmax <- max(vals, na.rm = TRUE) }
  if (vmax == 0) vmax <- 1

  cseg_sf$val <- pmin(pmax(vals, vmin), vmax)
  plot_data <- cseg_sf[order(cseg_sf$val), ]

  p <- ggplot(plot_data) +
    geom_sf(aes(fill = val), colour = "grey55", linewidth = 0.10) +
    scale_fill_gradientn(
      colours = GRAD_COLS,
      values  = GRAD_STOPS,
      name    = g,
      limits  = c(vmin, vmax),
      breaks  = pretty(c(vmin, vmax), n = 3),
      labels  = function(x) round(x, 2),
      oob     = scales::squish,
      guide   = guide_colorbar(
        barwidth       = unit(0.4, "cm"),
        barheight      = unit(2.2, "cm"),
        ticks          = TRUE,
        ticks.linewidth = unit(0.6, "mm"),
        frame.colour   = "grey30",
        frame.linewidth = unit(0.5, "mm"),
        ticks.colour   = "grey30",
        title.position = "top",
        title.hjust    = 0.5,
        title.theme    = element_text(size = 11, face = "bold.italic",
                                       colour = "grey10"),
        label.theme    = element_text(size = 9, colour = "grey10",
                                       face = "bold"))) +
    coord_sf(xlim = c(BB_ZOOM["xmin"], BB_ZOOM["xmax"]),
             ylim = c(BB_ZOOM["ymin"], BB_ZOOM["ymax"]),
             expand = FALSE) +
    theme_void(base_size = 8) +
    theme(
      legend.position  = "right",
      legend.justification = c(0, 1),
      legend.margin    = margin(0, 0, 0, 4),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin      = margin(4, 4, 4, 4)
    )

  svg_file <- file.path(FIG_DIR, sprintf("SP24_24824_roi_%s_zoom.svg", g))
  ggsave(svg_file, p, width = 3.8, height = 3, bg = "transparent")
  cat(sprintf("  SVG: %s\n", basename(svg_file)))

  png_file <- file.path(FIG_DIR, sprintf("SP24_24824_roi_%s_zoom.png", g))
  ggsave(png_file, p, width = 3.8, height = 3, dpi = 450,
         bg = "transparent", device = ragg::agg_png)
  cat(sprintf("  PNG: %s\n", basename(png_file)))
}

cat("\n=== ALL ZOOM PANELS DONE ===\n")
