#!/usr/bin/env Rscript
# ============================================================================
# Figure 4K — OTB_2384 ROI_C pathway UCell spatial maps (3 panels)
# ----------------------------------------------------------------------------
# PURPOSE: Per-pathway ROI spotlight on OTB_2384 ROI_C: every cell filled by its
#   per-cell UCell pathway score on a beige->burgundy ramp, one panel per
#   pathway (matrix remodeling, hypoxia, oxidative stress). Fig 4K uses the
#   hypoxia panel; the others support adjacent ROI figures.
#
# INPUTS:
#   - SFE object : <data_root>/<sfe_dir>/sfe_OTB_2384  (cellSeg + UCell pathway cols)
#   - output_root/06f_reclassification_polarization/reclassified_xenium_scores.csv
#     (06f override used for the cell census only; fill is by pathway score)
#   - load_sfe from spatial/00_setup/00_setup.R
#
# OUTPUTS:
#   - figures_dir/xenium_roi_OTB_2384_{matrix_remodeling,hypoxia,oxidative_stress}.{png,svg}
#
# MANUSCRIPT PANEL(S): Fig 4K (hypoxia).
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
  library(data.table); library(ggplot2); library(sf); library(scales); library(grid)
})

fig_dir <- path.expand(CFG$paths$figures_dir)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
SAMPLE  <- "sfe_OTB_2384"
BBOX    <- c(xmin = 2800, xmax = 4000, ymin = -7800, ymax = -6600)

PATHWAYS <- list(
  matrix_remodeling = list(col = "pathway_matrix_remodeling", label = "Matrix remodeling"),
  hypoxia           = list(col = "pathway_hypoxia",           label = "Hypoxia"),
  oxidative_stress  = list(col = "pathway_oxidative_stress",  label = "Oxidative stress")
)

EXPR_COLS  <- c("#F6EFE5", "#ECDDD0", "#D89E97", "#A03A4A", "#3A111A")
EXPR_STOPS <- c(0.00,      0.28,      0.55,      0.80,      1.00)

# --- Load + 06f override (cell census only) ----------------------------------
f06f <- cfg_path("output_root", "06f_reclassification_polarization",
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

co <- spatialCoords(sfe)
in_roi <- co[, 1] >= BBOX["xmin"] & co[, 1] <= BBOX["xmax"] &
          co[, 2] >= BBOX["ymin"] & co[, 2] <= BBOX["ymax"]
sfe_roi <- sfe[, in_roi]
cat(sprintf("  ROI cells: %d\n", ncol(sfe_roi)))

cseg <- colGeometry(sfe_roi, "cellSeg")
cseg$cell_label <- as.character(sfe_roi$cell_label)
for (k in names(PATHWAYS)) {
  cseg[[PATHWAYS[[k]]$col]] <- as.numeric(colData(sfe_roi)[[PATHWAYS[[k]]$col]])
}
cseg_sf <- sf::st_as_sf(cseg)

render_panel <- function(sf_layer, pw_col, label) {
  vals <- sf_layer[[pw_col]]
  qlim <- quantile(vals, c(0.02, 0.98), na.rm = TRUE)
  vmin <- as.numeric(qlim[1]); vmax <- as.numeric(qlim[2])
  sf_layer$val <- pmin(pmax(vals, vmin), vmax)
  sf_layer <- sf_layer[order(sf_layer$val), ]
  ggplot(sf_layer) +
    geom_sf(aes(fill = val), colour = "grey55", linewidth = 0.10) +
    scale_fill_gradientn(
      colours = EXPR_COLS, values = EXPR_STOPS, name = "UCell",
      limits = c(vmin, vmax), breaks = c(vmin, (vmin + vmax) / 2, vmax),
      labels = function(x) format(round(x, 2), nsmall = 2), oob = scales::squish,
      guide = guide_colorbar(barwidth = unit(0.22, "cm"), barheight = unit(1.6, "cm"),
                             ticks = TRUE, frame.colour = "grey40", ticks.colour = "grey40",
                             title.theme = element_text(size = 5, colour = "grey20"),
                             label.theme = element_text(size = 5, colour = "grey20"))) +
    coord_sf(xlim = c(BBOX["xmin"], BBOX["xmax"]), ylim = c(BBOX["ymin"], BBOX["ymax"]),
             expand = FALSE) +
    theme_void(base_size = 6) +
    theme(legend.position = "right", legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 1), plot.margin = margin(1, 2, 1, 1))
}

aspect <- as.numeric(BBOX["xmax"] - BBOX["xmin"]) / as.numeric(BBOX["ymax"] - BBOX["ymin"])
panel_h_mm <- 60
panel_w_mm <- panel_h_mm * aspect + 9
w_in <- panel_w_mm / 25.4; h_in <- panel_h_mm / 25.4

for (k in names(PATHWAYS)) {
  p <- render_panel(cseg_sf, PATHWAYS[[k]]$col, PATHWAYS[[k]]$label)
  out_stem <- sprintf("xenium_roi_OTB_2384_%s", k)
  png_out <- file.path(fig_dir, paste0(out_stem, ".png"))
  svg_out <- file.path(fig_dir, paste0(out_stem, ".svg"))
  ggsave(png_out, p, width = w_in, height = h_in, dpi = 450, bg = "white",
         device = ragg::agg_png, limitsize = FALSE)
  ggsave(svg_out, p, width = w_in, height = h_in, bg = "white", limitsize = FALSE)
  cat("  saved: ", basename(png_out), " | ", basename(svg_out), "\n", sep = "")
}
