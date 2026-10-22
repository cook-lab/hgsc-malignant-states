#!/usr/bin/env Rscript
# ============================================================================
# Figure 4b / 4c — WT01 cell-type maps (flattened high-resolution rasters)
# ----------------------------------------------------------------------------
# PURPOSE
#   4b: WT01 whole-tissue cell-type map (1 mm scale bar). 4c: cell-type map of
#   ROI_C in WT01 (cell polygons, 200 um scale bar; the same ROI as 4h and 4j).
#   Cells are coloured by the SingleR label cell_label_singler (secretory cells
#   split into SecA / Intermediate / SecB); cells without a label are not drawn.
#
#   OUTPUT FORM: PNG @ 600 dpi ONLY. A PNG is a single flattened raster, so
#   Affinity embeds one compact image instead of a million-point SVG — small
#   file, publication-grade resolution. (No SVG is written for these panels.)
#
# INPUTS
#   - SFE object : sfe_WT01 via load_sfe() (coords, cellSeg, cell_id)
#   - SingleR    : output_root/singler_recompute/WT01_singler_cells.csv.gz
#                    (barcode, cell_label_singler)
#   - ref_palette / load_sfe from spatial/00_setup/00_setup.R; scalebar.R
#
# OUTPUTS (figures_dir/figure4/; "4C"/"4D"/"split" in the names are earlier labels)
#   - fig4C_WT01_split_celltype.png            (Fig 4b, whole tissue, 600 dpi)
#   - fig4D_ROI_C_split_celltype.png               (Fig 4c, ROI_C cell type, 600 dpi)
#
# MANUSCRIPT PANEL(S): Fig 4b, 4c. RUNTIME TIER: heavy.
# ============================================================================

# --- Shared spatial setup (config, load_sfe, ref_palette, seed) --------------
.fig_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(.fig_dir) || !nzchar(.fig_dir)) .fig_dir <- "figures/figure4"
.setup_path <- normalizePath(file.path(.fig_dir, "..", "..", "spatial", "00_setup", "00_setup.R"),
                             mustWork = FALSE)
if (!file.exists(.setup_path)) .setup_path <- "spatial/00_setup/00_setup.R"
source(.setup_path)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(Voyager)
})

fig_dir <- cfg_path("figures_dir", "figure4")
source(file.path("figures", "figure4", "scalebar.R"))   # 4b/c/h/j scale bars
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

SAMPLE   <- "WT01"
SFE_NAME <- paste0("sfe_", SAMPLE)
DPI      <- 600

# ---- palette and SingleR labels -------------------------------------------
pal <- c(ref_palette, "Secretory epithelium" = "#E6D2A8")

# SingleR label table of WT01; the column split_label holds cell_label_singler
.sgl <- fread(file.path(cfg_path("output_root"),"singler_recompute","WT01_singler_cells.csv.gz"), select=c("barcode","cell_label_singler"))
sm <- .sgl[, .(barcode, split_label=cell_label_singler)]

# ---- load SFE, attach the SingleR label (split_label) ------------------------
cat("\n=== load", SFE_NAME, "===\n")
sfe <- load_sfe(SFE_NAME)
cid <- colnames(sfe)
lab <- sm$split_label[match(cid, sm$barcode)]
sfe$split_label <- lab
cat(sprintf("  matched %d / %d SFE cells to SPLIT label\n", sum(!is.na(lab)), length(cid)))

spatial_point_layer <- function(mapping = aes(), ...) {
  ggrastr::geom_point_rast(mapping = mapping, shape = 16, size = 0.1,
                           alpha = 0.8, raster.dpi = DPI, dev = "ragg", ...)
}

# ============================ Panel 4b ======================================
cd <- as.data.table(as.data.frame(colData(sfe)))
cd[, cell_id := colnames(sfe)]
coords <- spatialCoords(sfe)
cd[, x := coords[, 1]]; cd[, y := coords[, 2]]
cd <- cd[!is.na(split_label)]                 # cells with a SingleR label
set.seed(CFG$seed); cd <- cd[sample(.N)]

xr <- diff(range(cd$x)); yr <- diff(range(cd$y))
h_in <- 3; w_in <- h_in * (xr / yr)
p_c <- ggplot(cd, aes(x = x, y = y)) +
  spatial_point_layer(aes(color = split_label)) +
  scale_color_manual(values = pal, drop = TRUE) +
  coord_fixed() + theme_void(base_size = 8) +
  theme(legend.position = "none") +
  scalebar_layers(min(cd$x), max(cd$x), min(cd$y), max(cd$y), 1000, h_in)   # 1 mm
out_c <- file.path(fig_dir, "fig4C_WT01_split_celltype.png")
ggsave(out_c, p_c, width = w_in, height = h_in, dpi = DPI)
cat("  Saved:", basename(out_c), sprintf("(%.0f x %.0f px)\n", w_in*DPI, h_in*DPI))

# ============================ Panel 4c ======================================
# ROI_C, the ROI shared with 4h and 4j. Cells are selected in the spatialCoords
# frame (y positive); the cellSeg polygons are stored with y negative, so the
# scale bar takes its limits from the polygons' own bounding box.
bb <- c(xmin = 2800, xmax = 4000, ymin = 6600, ymax = 7800)     # ROI_C (spatialCoords frame)
in_roi <- coords[, 1] >= bb["xmin"] & coords[, 1] <= bb["xmax"] &
          coords[, 2] >= bb["ymin"] & coords[, 2] <= bb["ymax"]
sfe_roi <- sfe[, in_roi & !is.na(sfe$split_label)]
cat("\n=== ROI_C cells:", ncol(sfe_roi), "===\n")
RB <- sf::st_bbox(colGeometry(sfe_roi, "cellSeg"))   # cellSeg (display) frame
sb_roi <- scalebar_layers(RB[["xmin"]], RB[["xmax"]], RB[["ymin"]], RB[["ymax"]], 200, 7)   # 200 um

p_ct <- plotSpatialFeature(sfe_roi, "split_label", colGeometryName = "cellSeg",
                           aes_use = "fill", linewidth = 0.15, color = "grey50") +
  scale_fill_manual(values = pal, name = "Cell type", drop = TRUE) +
  theme_void(base_size = 10) +
  theme(legend.position = "none") + sb_roi
out_d1 <- file.path(fig_dir, "fig4D_ROI_C_split_celltype.png")
ggsave(out_d1, p_ct, width = 8, height = 7, dpi = DPI)
cat("  Saved:", basename(out_d1), "\n")

cat("\nDone. 4C/4D recolored by clean-SPLIT label; flattened PNG @", DPI, "dpi (no SVG).\n")
