#!/usr/bin/env Rscript
# ============================================================================
# Fig 4j -- WT01 ROI_C hypoxia UCell map
# ----------------------------------------------------------------------------
# PURPOSE: the ROI of 4c / 4h coloured by the per-cell pathway_hypoxia UCell
#   score (12 genes; 2nd-98th percentile colour limits), on cellSeg polygons.
#   Reads pathway_* straight from the SFE colData (written by
#   spatial/04_pathway_scoring/04c_pathway_scoring_allcells.R: UCell, full
#   477-gene panel, all cells, maxRank = 477); cells without a score are not drawn.
#   ROI is SELECTED in the spatialCoords frame and DISPLAYED in the cellSeg frame
#   (the two are mirrored in y); limits come from the subset geometry's own bbox.
# INPUTS: sfe_WT01 (pathway_*, spatialCoords, cellSeg); output_root/9b_scoring/
#   pathway_gene_sets_v2.csv (gene list printed for the log); figures/figure4/scalebar.R
# OUTPUTS: figures_dir/figure4/xenium_roi_WT01_hypoxia.{png,svg}
# MANUSCRIPT PANEL(S): Fig 4j. RUNTIME TIER: moderate.
# ============================================================================
source("spatial/00_setup/00_setup.R")
suppressPackageStartupMessages({
  library(SpatialFeatureExperiment); library(SummarizedExperiment)
  library(data.table); library(ggplot2); library(sf); library(scales); library(grid)
})
fig_dir <- cfg_path("figures_dir","figure4"); R <- cfg_path("output_root")
SAMPLE <- "sfe_WT01"
SELBOX <- c(xmin=2800, xmax=4000, ymin=6600, ymax=7800)
WANT <- c("hypoxia")
LABS <- c(hypoxia="Hypoxia")
EXPR_COLS  <- c("#F6EFE5","#ECDDD0","#D89E97","#A03A4A","#3A111A")
EXPR_STOPS <- c(0.00,0.28,0.55,0.80,1.00)

gs <- fread(file.path(R,"9b_scoring","pathway_gene_sets_v2.csv"))
for (p in WANT) cat(sprintf("  %-18s set: %s\n", p, paste(gs[pathway==p, gene], collapse=", ")))

sfe <- load_sfe(SAMPLE)
cols <- paste0("pathway_", WANT)
stopifnot(all(cols %in% colnames(colData(sfe))))
cat("\npathway columns found in SFE:", paste(cols, collapse=", "), "\n")
co <- spatialCoords(sfe)
scored <- !is.na(colData(sfe)[["pathway_hypoxia"]])
cat("scored (SPLIT-kept) cells:", sum(scored), "/", ncol(sfe), "\n")

in_roi <- co[,1] >= SELBOX["xmin"] & co[,1] <= SELBOX["xmax"] &
          co[,2] >= SELBOX["ymin"] & co[,2] <= SELBOX["ymax"] & scored
sfe_roi <- sfe[, in_roi]
cat("ROI cells:", ncol(sfe_roi), "\n")

cseg <- colGeometry(sfe_roi, "cellSeg")
for (c0 in cols) cseg[[c0]] <- as.numeric(sfe_roi[[c0]])
cseg_sf <- sf::st_as_sf(cseg); GBB <- sf::st_bbox(cseg_sf)

source(file.path("figures", "figure4", "scalebar.R"))   # 4b/c/h/j scale bars
for (p in WANT) {
  v <- cseg_sf[[paste0("pathway_", p)]]
  ql <- quantile(v, c(0.02, 0.98), na.rm = TRUE)
  d <- cseg_sf; d$val <- pmin(pmax(v, ql[1]), ql[2]); d <- d[order(d$val), ]
  pl <- ggplot(d) + geom_sf(aes(fill = val), colour = "grey55", linewidth = 0.10) +
    scale_fill_gradientn(colours = EXPR_COLS, values = EXPR_STOPS, name = LABS[[p]],
      limits = as.numeric(ql), breaks = pretty(as.numeric(ql), n = 3),
      labels = function(x) sprintf("%.2f", x), oob = scales::squish,
      guide = guide_colorbar(barwidth = unit(0.22,"cm"), barheight = unit(1.6,"cm"),
        ticks = TRUE, frame.colour = "grey40", ticks.colour = "grey40",
        title.theme = element_text(size = 5, colour = "grey20"),
        label.theme = element_text(size = 5, colour = "grey20"))) +
    coord_sf(xlim = c(GBB[["xmin"]], GBB[["xmax"]]),
             ylim = c(GBB[["ymin"]], GBB[["ymax"]]), expand = FALSE) +
    theme_void(base_size = 6) +
    theme(legend.position = "right", legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(0,0,0,1), plot.margin = margin(1,2,1,1))
  asp <- as.numeric(GBB[["xmax"]]-GBB[["xmin"]])/as.numeric(GBB[["ymax"]]-GBB[["ymin"]])
  h_mm <- 60; w_in <- (h_mm*asp+9)/25.4; h_in <- h_mm/25.4
  pl <- pl + scalebar_layers(GBB[["xmin"]], GBB[["xmax"]], GBB[["ymin"]], GBB[["ymax"]], 200, h_in)
  f <- file.path(fig_dir, paste0("xenium_roi_WT01_", p))
  ggsave(paste0(f,".png"), pl, width=w_in, height=h_in, dpi=450, bg="white",
         device=ragg::agg_png, limitsize=FALSE)
  # SVG: cell polygons as one 600-dpi image (2026-10-21, author); colour bar, scale bar, text vector
  pl_svg <- pl; pl_svg$layers[[1]] <- ggrastr::rasterise(pl_svg$layers[[1]], dpi = 600, dev = "ragg")
  ggsave(paste0(f,".svg"), pl_svg, width=w_in, height=h_in, bg="white", limitsize=FALSE)
  cat(sprintf("  saved %-34s range %.3f-%.3f\n", basename(f), ql[1], ql[2]))
}
cat("\nDONE\n")
