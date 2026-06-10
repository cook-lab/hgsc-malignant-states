#!/usr/bin/env Rscript
# ============================================================================
# Figure 4D — OTB_2384 ROI_C cell-type + macrophage-overlay spatial maps
# ----------------------------------------------------------------------------
# PURPOSE: Magnified ROI_C (OTB_2384) cell-type map and a macrophage+epitype
#   overlay, via Voyager plotSpatialFeature on cellSeg polygons. Shows the
#   SecA->SecB spatial gradient.
#
# INPUTS:
#   - SFE object : <data_root>/<sfe_dir>/sfe_OTB_2384  (cell_label, cellSeg)
#   - ref_palette / load_sfe from spatial/00_setup/00_setup.R
#
# OUTPUTS:
#   - figures_dir/figure4/fig2d_ROI_C_wide_celltype.{svg,png}
#   - figures_dir/figure4/fig2d_ROI_C_wide_macrophage.{svg,png}
#
# MANUSCRIPT PANEL(S): Fig 4D.
#
# RUNTIME TIER: moderate.
#
# FLAG (LINEAGE.md): this generator colours by the in-object `cell_label`
#   (SingleR / pre-06f), NOT the 06f polarization reclassification used by the
#   rest of Fig 4. If the published 4D used 06f labels, apply the 06f override
#   as in 02_xenium_whole_tissue_snapshot.R before plotting. Migrated faithfully;
#   resolution is the authors' call (see docs/REPRODUCIBILITY.md).
#
# NOTE: epithelial label standardized "Transitioning" -> "Intermediate".
# ============================================================================

# --- Shared spatial setup (provides config, load_sfe, ref_palette) -----------
.fig_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(.fig_dir) || !nzchar(.fig_dir)) .fig_dir <- "figures/figure4"
.setup_path <- normalizePath(file.path(.fig_dir, "..", "..", "spatial", "00_setup", "00_setup.R"),
                             mustWork = FALSE)
if (!file.exists(.setup_path)) .setup_path <- "spatial/00_setup/00_setup.R"
source(.setup_path)
library(Voyager)

fig_dir <- cfg_path("figures_dir", "figure4")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

epi_types <- c("SecA epithelium", "Intermediate epithelium", "SecB epithelium")

spatial_point_layer <- function(mapping = aes(), ...) {
  ggrastr::geom_point_rast(mapping = mapping, shape = 16, size = 0.1,
                           alpha = 0.8, raster.dpi = 300, dev = "ragg", ...)
}

cat("=== Fig 4D ===\n")
sfe <- load_sfe("sfe_OTB_2384")
# Rename-mismatch fix (idempotent): deposited SFE still carries the legacy
# "Transitioning epithelium"; standardize to "Intermediate epithelium" on the
# SAME cell_label vector the colour/filter below keys on (lines ~59, ~73).
.lab <- as.character(sfe$cell_label)
.lab[.lab == "Transitioning epithelium"] <- "Intermediate epithelium"
sfe$cell_label <- .lab
bb_2d <- c(xmin = 2800, xmax = 4000, ymin = -7800, ymax = -6600)
coords <- spatialCoords(sfe)
in_roi <- coords[, 1] >= bb_2d["xmin"] & coords[, 1] <= bb_2d["xmax"] &
          coords[, 2] >= bb_2d["ymin"] & coords[, 2] <= bb_2d["ymax"]
sfe_roi <- sfe[, in_roi]
cat("  ROI cells:", ncol(sfe_roi), "\n")

# Cell type
p_ct <- plotSpatialFeature(sfe_roi, "cell_label", colGeometryName = "cellSeg",
                           aes_use = "fill", linewidth = 0.15, color = "grey50") +
  scale_fill_manual(values = ref_palette, name = "Cell type", drop = TRUE) +
  theme_void(base_size = 10) +
  theme(legend.position = "right", legend.text = element_text(size = 8),
        legend.title = element_text(size = 9), legend.key.size = unit(0.5, "lines")) +
  guides(fill = guide_legend(override.aes = list(linewidth = 0)))

ggsave(file.path(fig_dir, "fig2d_ROI_C_wide_celltype.png"), p_ct, width = 8, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "fig2d_ROI_C_wide_celltype.svg"), p_ct, width = 8, height = 7, dpi = 300)

# Epitype + Macrophage overlay
sfe_roi$mac_group <- ifelse(
  sfe_roi$cell_label == "Macrophage", "Macrophage",
  ifelse(sfe_roi$cell_label %in% epi_types, as.character(sfe_roi$cell_label), "Other"))

mac_pal <- c("SecA epithelium" = "#E6A141", "Intermediate epithelium" = "#C08E48",
             "SecB epithelium" = "#6B5530", "Macrophage" = "#8FBC8F", "Other" = "#E0E0E0")

p_mac <- plotSpatialFeature(sfe_roi, "mac_group", colGeometryName = "cellSeg",
                            aes_use = "fill", linewidth = 0.15, color = "grey50") +
  scale_fill_manual(values = mac_pal, name = "Cell type", drop = TRUE) +
  theme_void(base_size = 10) +
  theme(legend.position = "right", legend.text = element_text(size = 8),
        legend.title = element_text(size = 9), legend.key.size = unit(0.5, "lines")) +
  guides(fill = guide_legend(override.aes = list(linewidth = 0)))

ggsave(file.path(fig_dir, "fig2d_ROI_C_wide_macrophage.png"), p_mac, width = 8, height = 7, dpi = 300)
ggsave(file.path(fig_dir, "fig2d_ROI_C_wide_macrophage.svg"), p_mac, width = 8, height = 7, dpi = 300)
cat("Saved fig2d\n")
