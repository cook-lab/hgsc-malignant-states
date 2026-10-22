#!/usr/bin/env Rscript
# ============================================================================
# Figure 4h — WT01 ROI_C distance-to-vasculature spatial map
# ----------------------------------------------------------------------------
# PURPOSE: Per-cell distance to nearest vascular cell (Pericyte/Endothelial),
#   plotted as a continuous gradient on cellSeg polygons within WT01 ROI_C
#   (the ROI of 4c and 4j). Distance via RANN::nn2 on spatialCoords centroids
#   (microns), over the whole section; cells without a SingleR label are dropped.
#
# INPUTS:
#   - SFE object : sfe_WT01 via load_sfe() (spatialCoords, cellSeg)
#   - output_root/singler_recompute/WT01_singler_cells.csv.gz (barcode, cell_label_singler)
#   - load_sfe from spatial/00_setup/00_setup.R; figures/figure4/scalebar.R
#
# OUTPUTS:
#   - figures_dir/figure4/xenium_roi_WT01_vascular_distance.{png,svg}
#
# MANUSCRIPT PANEL(S): Fig 4h.
#
# RUNTIME TIER: moderate.
# ============================================================================

# --- Shared spatial setup (provides config, load_sfe) ------------------------
.fig_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(.fig_dir) || !nzchar(.fig_dir)) .fig_dir <- "figures/figure4"
.setup_path <- normalizePath(file.path(.fig_dir, "..", "..", "spatial", "00_setup", "00_setup.R"),
                             mustWork = FALSE)
if (!file.exists(.setup_path)) .setup_path <- "spatial/00_setup/00_setup.R"
source(.setup_path)

suppressPackageStartupMessages({
  library(SpatialFeatureExperiment); library(SummarizedExperiment)
  library(data.table); library(ggplot2); library(sf); library(scales)
  library(grid); library(RANN)
})

fig_dir <- cfg_path("figures_dir", "figure4")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
SAMPLE  <- "sfe_WT01"
# ROI_C SELECTION bbox, in the spatialCoords frame. The SFE stores spatialCoords
# y POSITIVE (~+5600..+11900) whereas the cellSeg polygons are stored NEGATIVE
# (~-11946..-5592) -- the two frames are mirrored in y. So cells are SELECTED
# here in the spatialCoords frame (the same box as 4c in
# 03_fig4bc_wt01_celltype_maps.R), while the coord_sf DISPLAY limits are taken
# from the subset polygons' own bbox further down.
BBOX    <- c(xmin = 2800, xmax = 4000, ymin = 6600, ymax = 7800)
VASCULAR_TYPES <- c("Pericyte", "Endothelial")

DIST_COLS  <- c("#F6EFE5", "#ECDDD0", "#D89E97", "#A03A4A", "#3A111A")
DIST_STOPS <- c(0.00,      0.28,      0.55,      0.80,      1.00)

# --- Load SFE + SingleR labels (cell_label_singler) -----------------------------

source(file.path("figures", "figure4", "scalebar.R"))   # 4b/c/h/j scale bars
cat("\n=== Loading", SAMPLE, "===\n")
# SingleR label table of WT01
.sgl <- fread(file.path(cfg_path("output_root"),"singler_recompute","WT01_singler_cells.csv.gz"), select=c("barcode","cell_label_singler"))
sm <- data.table(barcode=.sgl$barcode, lab=.sgl$cell_label_singler)
sfe <- load_sfe(SAMPLE)
lab <- sm$lab[match(colnames(sfe), sm$barcode)]
sfe$cell_label <- lab
cat(sprintf("  SPLIT label on %d / %d cells (%d rejected -> dropped)\n",
            sum(!is.na(lab)), length(lab), sum(is.na(lab))))
sfe <- sfe[, !is.na(lab)]

# --- Distance to nearest vascular cell ---------------------------------------
co <- spatialCoords(sfe)
vasc_idx <- which(sm$lab[match(colnames(sfe),sm$barcode)] %in% VASCULAR_TYPES)
vasc_coords <- co[vasc_idx, , drop = FALSE]
nn <- nn2(vasc_coords, co, k = 1)
sfe$dist_to_vascular <- as.numeric(nn$nn.dists)

# --- Subset to ROI -----------------------------------------------------------
in_roi <- co[, 1] >= BBOX["xmin"] & co[, 1] <= BBOX["xmax"] &
          co[, 2] >= BBOX["ymin"] & co[, 2] <= BBOX["ymax"]
sfe_roi <- sfe[, in_roi]
cat(sprintf("  ROI cells: %d\n", ncol(sfe_roi)))

cseg <- colGeometry(sfe_roi, "cellSeg")
cseg$dist_to_vascular <- as.numeric(sfe_roi$dist_to_vascular)
cseg$cell_label <- as.character(sfe_roi$cell_label)
cseg_sf <- sf::st_as_sf(cseg)
# display limits from the polygons themselves (cellSeg frame), not from BBOX
GBB <- sf::st_bbox(cseg_sf)
cat(sprintf("  cellSeg display bbox: x [%.0f, %.0f]  y [%.0f, %.0f]\n",
            GBB[["xmin"]], GBB[["xmax"]], GBB[["ymin"]], GBB[["ymax"]]))

# --- Plot --------------------------------------------------------------------
qlim <- quantile(cseg_sf$dist_to_vascular, c(0.02, 0.98), na.rm = TRUE)
vmin <- as.numeric(qlim[1]); vmax <- as.numeric(qlim[2])
cseg_sf$val <- pmin(pmax(cseg_sf$dist_to_vascular, vmin), vmax)
cseg_sf <- cseg_sf[order(cseg_sf$val), ]

p <- ggplot(cseg_sf) +
  geom_sf(aes(fill = val), colour = "grey55", linewidth = 0.10) +
  scale_fill_gradientn(
    colours = DIST_COLS, values = DIST_STOPS,
    name = expression(paste("Distance (", mu, "m)")),
    limits = c(vmin, vmax), breaks = pretty(c(vmin, vmax), n = 3),
    labels = function(x) round(x), oob = scales::squish,
    guide = guide_colorbar(barwidth = unit(0.22, "cm"), barheight = unit(1.6, "cm"),
                           ticks = TRUE, frame.colour = "grey40", ticks.colour = "grey40",
                           title.theme = element_text(size = 5, colour = "grey20"),
                           label.theme = element_text(size = 5, colour = "grey20"))) +
  coord_sf(xlim = c(GBB[["xmin"]], GBB[["xmax"]]), ylim = c(GBB[["ymin"]], GBB[["ymax"]]),
           expand = FALSE) +
  theme_void(base_size = 6) +
  theme(legend.position = "right", legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 1), plot.margin = margin(1, 2, 1, 1))

aspect <- as.numeric(GBB[["xmax"]] - GBB[["xmin"]]) / as.numeric(GBB[["ymax"]] - GBB[["ymin"]])
panel_h_mm <- 60
p <- p + scalebar_layers(GBB[["xmin"]], GBB[["xmax"]], GBB[["ymin"]], GBB[["ymax"]], 200, panel_h_mm / 25.4)
panel_w_mm <- panel_h_mm * aspect + 9
w_in <- panel_w_mm / 25.4; h_in <- panel_h_mm / 25.4

png_out <- file.path(fig_dir, "xenium_roi_WT01_vascular_distance.png")
svg_out <- file.path(fig_dir, "xenium_roi_WT01_vascular_distance.svg")
ggsave(png_out, p, width = w_in, height = h_in, dpi = 450, bg = "white",
       device = ragg::agg_png, limitsize = FALSE)
# SVG: the cell polygons go in as one 600-dpi image (2026-10-21, author); colour bar,
# scale bar and text stay vector
p_svg <- p; p_svg$layers[[1]] <- ggrastr::rasterise(p_svg$layers[[1]], dpi = 600, dev = "ragg")
ggsave(svg_out, p_svg, width = w_in, height = h_in, bg = "white", limitsize = FALSE)
cat("\nSaved:\n  ", png_out, "\n  ", svg_out, "\n\nDONE\n")
