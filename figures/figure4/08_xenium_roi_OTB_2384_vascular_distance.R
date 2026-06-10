#!/usr/bin/env Rscript
# ============================================================================
# Figure 4I — OTB_2384 ROI_C distance-to-vasculature spatial map
# ----------------------------------------------------------------------------
# PURPOSE: Per-cell distance to nearest vascular cell (Pericyte/Endothelial),
#   plotted as a continuous gradient on cellSeg polygons within OTB_2384 ROI_C.
#   Distance via RANN::nn2 on spatialCoords centroids (microns).
#
# INPUTS:
#   - SFE object : <data_root>/<sfe_dir>/sfe_OTB_2384  (cellSeg, cell_label)
#   - output_root/06f_reclassification_polarization/reclassified_xenium_scores.csv
#   - load_sfe from spatial/00_setup/00_setup.R
#
# OUTPUTS:
#   - figures_dir/figure4/xenium_roi_OTB_2384_vascular_distance.{png,svg}
#
# MANUSCRIPT PANEL(S): Fig 4I.
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
SAMPLE  <- "sfe_OTB_2384"
BBOX    <- c(xmin = 2800, xmax = 4000, ymin = -7800, ymax = -6600)
VASCULAR_TYPES <- c("Pericyte", "Endothelial")

DIST_COLS  <- c("#F6EFE5", "#ECDDD0", "#D89E97", "#A03A4A", "#3A111A")
DIST_STOPS <- c(0.00,      0.28,      0.55,      0.80,      1.00)

# --- Load + 06f override -----------------------------------------------------
f06f <- cfg_path("data_root", "2026_final_xenium_analysis", "output", "06f_reclassification_polarization",
                 "reclassified_xenium_scores.csv")
recl <- fread(f06f, select = c("sample", "barcode_orig", "cell_label_06f"))

cat("\n=== Loading", SAMPLE, "===\n")
sfe <- load_sfe(SAMPLE)
m <- match(colnames(sfe), recl[sample == SAMPLE, barcode_orig])
hit <- !is.na(m)
lab <- as.character(sfe$cell_label)
lab[hit] <- recl[sample == SAMPLE, cell_label_06f][m[hit]]
sfe$cell_label <- lab
cat(sprintf("  06f override applied to %d / %d cells\n", sum(hit), length(hit)))

# --- Distance to nearest vascular cell ---------------------------------------
co <- spatialCoords(sfe)
vasc_idx <- which(sfe$cell_label %in% VASCULAR_TYPES)
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
  coord_sf(xlim = c(BBOX["xmin"], BBOX["xmax"]), ylim = c(BBOX["ymin"], BBOX["ymax"]),
           expand = FALSE) +
  theme_void(base_size = 6) +
  theme(legend.position = "right", legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 1), plot.margin = margin(1, 2, 1, 1))

aspect <- as.numeric(BBOX["xmax"] - BBOX["xmin"]) / as.numeric(BBOX["ymax"] - BBOX["ymin"])
panel_h_mm <- 60
panel_w_mm <- panel_h_mm * aspect + 9
w_in <- panel_w_mm / 25.4; h_in <- panel_h_mm / 25.4

png_out <- file.path(fig_dir, "xenium_roi_OTB_2384_vascular_distance.png")
svg_out <- file.path(fig_dir, "xenium_roi_OTB_2384_vascular_distance.svg")
ggsave(png_out, p, width = w_in, height = h_in, dpi = 450, bg = "white",
       device = ragg::agg_png, limitsize = FALSE)
ggsave(svg_out, p, width = w_in, height = h_in, bg = "white", limitsize = FALSE)
cat("\nSaved:\n  ", png_out, "\n  ", svg_out, "\n\nDONE\n")
