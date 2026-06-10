#!/usr/bin/env Rscript
# =============================================================================
# SF10B — TMA spatial composition map
# -----------------------------------------------------------------------------
# Purpose:  All TMA cells plotted at spatial coordinates, coloured by cell_label
#           (ref_palette). TMA_1 and TMA_2 shown side by side with a gap; legend
#           included. High-resolution raster for publication (half letter page).
#
# INPUTS:   cfg_obj("sfe_tma_filtered")   (canonical filtered TMA SFE; cell_label
#                                           includes the 06f epithelial polarization)
# OUTPUTS:  output_root/figures/supplementary/SF10B_tma_composition_spatial.{pdf,png,svg}
# MANUSCRIPT PANEL(S):  SF10B.
# RUNTIME TIER:  moderate (600-dpi rasterized full-TMA point cloud).
# =============================================================================

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

# Config + shared spatial setup (2 levels deep -> ../../).
source(file.path("..", "..", "config", "config.R"))
source(file.path("..", "..", "spatial", "00_setup", "00_setup.R"))  # load_sfe, ref_palette

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
  library(ggrastr); library(ragg); library(svglite)
})

set.seed(CFG$seed)

OUT_STEM <- cfg_path("output_root", "figures", "supplementary", "SF10B_tma_composition_spatial")

# -- Load TMA -----------------------------------------------------------------
message("[1] Loading TMA ...")
sfe <- load_sfe("sfe_tma_filtered")
cd <- as.data.table(as.data.frame(colData(sfe)))
cd[, cell_id := colnames(sfe)]
coords <- spatialCoords(sfe)
cd[, x := coords[, 1]]
cd[, y := coords[, 2]]

message(sprintf("  %s cells, %d cell types",
                format(nrow(cd), big.mark = ","),
                length(unique(cd$cell_label))))

# -- Arrange TMA_1 and TMA_2 side by side -------------------------------------
gap <- 1500
tma1_xmax <- cd[sample_id == "TMA_1", max(x)]
tma2_xmin <- cd[sample_id == "TMA_2", min(x)]
cd[sample_id == "TMA_2", x := x + (tma1_xmax - tma2_xmin + gap)]

# Shuffle draw order
cd <- cd[sample(.N)]

rm(sfe); gc(verbose = FALSE)

# -- Order cell types for legend: epithelial first, then stromal, immune ------
# "Intermediate epithelium" (was "Transitioning epithelium").
epi_types <- c("SecA epithelium", "Intermediate epithelium",
               "SecB epithelium", "Ciliated epithelium")
stromal_types <- c("Fibroblast", "Smooth muscle", "Pericyte",
                   "Endothelial", "Mesothelial")
immune_types <- c("Macrophage", "T cell", "NK cell", "B cell",
                  "Plasma cell", "Conventional dendritic cell",
                  "Plasmacytoid dendritic cell", "Neutrophil", "Mast cell")

all_types <- c(epi_types, stromal_types, immune_types)
present <- intersect(all_types, unique(cd$cell_label))
# SFE colData still carries the legacy literal "Transitioning epithelium";
# rename to the standardized "Intermediate epithelium" before factor() so these
# cells are retained (otherwise they fall outside `present` levels and drop).
cd[cell_label == "Transitioning epithelium", cell_label := "Intermediate epithelium"]
cd[, cell_label := factor(cell_label, levels = present)]
pal <- ref_palette[present]

# -- Compute figure dimensions ------------------------------------------------
xr <- diff(range(cd$x, na.rm = TRUE))
yr <- diff(range(cd$y, na.rm = TRUE))
plot_w <- 8.5
plot_h <- plot_w * (yr / xr)
if (plot_h > 5.5) {
  plot_h <- 5.5
  plot_w <- plot_h * (xr / yr)
}
total_w <- plot_w + 1.5
message(sprintf("  Figure: %.1f x %.1f in (plot) + legend", plot_w, plot_h))

# -- Plot ---------------------------------------------------------------------
message("\n[2] Plotting ...")
p <- ggplot(cd, aes(x = x, y = y)) +
  geom_point_rast(aes(color = cell_label),
                  shape = 16, size = 0.05, alpha = 0.8,
                  raster.dpi = 600, dev = "ragg") +
  scale_color_manual(values = pal, name = "Cell type", drop = TRUE) +
  coord_fixed() +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 1)) +
  theme_void(base_size = 7) +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 7, face = "bold"),
    legend.text      = element_text(size = 5.5),
    legend.key.size  = unit(0.3, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    plot.margin      = margin(4, 4, 4, 4)
  )

# -- Save ---------------------------------------------------------------------
message("[3] Saving ...")
ggsave(paste0(OUT_STEM, ".png"), p, width = total_w, height = plot_h, dpi = 600,
       bg = "white", device = ragg::agg_png)
ggsave(paste0(OUT_STEM, ".svg"), p, width = total_w, height = plot_h, bg = "white")
ggsave(paste0(OUT_STEM, ".pdf"), p, width = total_w, height = plot_h, bg = "white")

message("\n  ", OUT_STEM, ".{png,svg,pdf}")
message("\nDONE")
