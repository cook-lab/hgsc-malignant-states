#!/usr/bin/env Rscript
# ============================================================================
# Figure 5H,5I — SP24_24824 ROI zoom: cell type + gene expression panels
# ----------------------------------------------------------------------------
# PURPOSE
#   SP24_24824 ROI 06g spatial panels: cell-type map (5H) + gene-expression
#   maps (5I) over cell-segmentation polygons, colorbar outside the panel.
#   CANONICAL = the PUBLISHED Fig 5I: the FULL ROI (x=[7800,9000] y=[-7100,-5900])
#   for 4 genes (CTNNB1, ITGB5, MMP7, ICAM1). CDH2/TGM2 appear only in the
#   Fig 5G GAM curves, not as 5I maps; the published 5H/5I overlay two ROI
#   rectangles in Illustrator (an assembly-time annotation, not rendered here).
#
# INPUTS
#   data_root/2026_final_xenium_analysis/output/sfe/sfe_SP24_24824 (via load_sfe)
#   data_root/2026_final_xenium_analysis/output/06f_reclassification_polarization/
#     reclassified_xenium_scores.csv  (06f cell-label override)
#   Shared helpers: config/config.R, spatial/00_setup/00_setup.R (load_sfe).
#
# OUTPUTS
#   figures_dir/figure5/ROI_figure_5/SP24_24824_roi_celltype_full.{svg,png} (5H)
#   figures_dir/figure5/ROI_figure_5/SP24_24824_roi_<gene>.{svg,png}        (5I; full ROI)
#
# MANUSCRIPT PANEL(S): Fig 5H (cell-type map, full ROI) + Fig 5I (gene maps).
#   LINEAGE attributes both panels to this generator. The canonical 5H render
#   is the full-ROI cellSeg polygons filled by cell_label with the canonical
#   ref_palette (matches 20260429_figures/ROI figure 5/SP24_24824_roi_celltype.svg;
#   the published panel adds two ROI rectangles downstream in Illustrator).
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
# Published Fig 5I = 4 genes at the full ROI (CDH2/TGM2 are Fig 5G GAM-only).
GENES   <- c("CTNNB1", "ITGB5", "MMP7", "ICAM1")

BB_FULL <- c(xmin = 7800, xmax = 9000, ymin = -7100, ymax = -5900)  # published 5H/5I extent
# BB_ZOOM marks one of the two ROI rectangles the published panels overlay in
# Illustrator; rendering is at BB_FULL (the boxes are an assembly-time annotation).
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

# ── Cell-type panel (Fig 5H) ─────────────────────────────────────────────
# Same full-ROI crop and cellSeg polygons as the canonical 5H render
# (202605_epitype_manuscript/20260429_figures/ROI figure 5/
#  SP24_24824_roi_celltype.svg). Fill polygons by cell_label using the
# canonical ref_palette (sourced from spatial/00_setup/00_setup.R). The
# published panel adds two ROI rectangles in Illustrator; here we emit the
# clean cell-type map (boxes are a downstream annotation, not data).
cat("\n--- cell type (Fig 5H) ---\n")

ct_data <- cseg_sf
# Reconcile naming drift: deposited SFE / 06f labels say "Transitioning
# epithelium"; ref_palette uses the standardized "Intermediate epithelium".
# Rename before the fill so those cells are coloured (not dropped to grey).
ct_data$cell_label[ct_data$cell_label == "Transitioning epithelium"] <-
  "Intermediate epithelium"
# A small number of ROI cells carry no annotation (empty cell_label) in the
# deposited data; draw them in neutral grey but keep them out of the legend.
ct_data$cell_label[!nzchar(ct_data$cell_label)] <- NA
ct_data$cell_label <- factor(ct_data$cell_label, levels = names(ref_palette))

p_ct <- ggplot(ct_data) +
  geom_sf(aes(fill = cell_label), colour = "grey55", linewidth = 0.10) +
  scale_fill_manual(values = ref_palette, name = "Cell type", drop = TRUE,
                    na.value = "grey85", na.translate = FALSE) +
  coord_sf(xlim = c(BB_FULL["xmin"], BB_FULL["xmax"]),
           ylim = c(BB_FULL["ymin"], BB_FULL["ymax"]),
           expand = FALSE) +
  theme_void(base_size = 8) +
  theme(
    legend.position  = "right",
    legend.justification = c(0, 1),
    legend.title     = element_text(size = 11, face = "bold", colour = "grey10"),
    legend.text      = element_text(size = 9, colour = "grey10"),
    legend.key.size  = unit(0.4, "cm"),
    legend.margin    = margin(0, 0, 0, 4),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin      = margin(4, 4, 4, 4)
  ) +
  guides(fill = guide_legend(override.aes = list(linewidth = 0)))

ct_svg <- file.path(FIG_DIR, "SP24_24824_roi_celltype_full.svg")
ggsave(ct_svg, p_ct, width = 5.0, height = 3, bg = "transparent")
cat(sprintf("  SVG: %s\n", basename(ct_svg)))

ct_png <- file.path(FIG_DIR, "SP24_24824_roi_celltype_full.png")
ggsave(ct_png, p_ct, width = 5.0, height = 3, dpi = 450,
       bg = "transparent", device = ragg::agg_png)
cat(sprintf("  PNG: %s\n", basename(ct_png)))

# ── Gene panels (Fig 5I — full ROI, published 4-gene set) ────────────────
for (g in GENES) {
  cat(sprintf("\n--- %s (full ROI) ---\n", g))

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
    coord_sf(xlim = c(BB_FULL["xmin"], BB_FULL["xmax"]),
             ylim = c(BB_FULL["ymin"], BB_FULL["ymax"]),
             expand = FALSE) +
    theme_void(base_size = 8) +
    theme(
      legend.position  = "right",
      legend.justification = c(0, 1),
      legend.margin    = margin(0, 0, 0, 4),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin      = margin(4, 4, 4, 4)
    )

  svg_file <- file.path(FIG_DIR, sprintf("SP24_24824_roi_%s.svg", g))
  ggsave(svg_file, p, width = 3.8, height = 3, bg = "transparent")
  cat(sprintf("  SVG: %s\n", basename(svg_file)))

  png_file <- file.path(FIG_DIR, sprintf("SP24_24824_roi_%s.png", g))
  ggsave(png_file, p, width = 3.8, height = 3, dpi = 450,
         bg = "transparent", device = ragg::agg_png)
  cat(sprintf("  PNG: %s\n", basename(png_file)))
}

cat("\n=== ALL PANELS DONE (5H cell-type + 5I full-ROI gene maps) ===\n")
