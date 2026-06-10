#!/usr/bin/env Rscript
# ============================================================================
# Figure 5B,5C — OTB_2384 ROI 06 spatial panels
# ----------------------------------------------------------------------------
# PURPOSE
#   Individual SVG + PNG panels for OTB_2384 ROI 06 (800 x 800 µm). Cell-type
#   map (5B) plus epithelial-only pathway / gene expression maps (5C):
#   glycolysis (UCell scored fresh on-ROI), proliferation, NF-κB, LDHA,
#   nuclear area. Non-epithelial cells drawn grey.
#   ROI 06 bbox: x=[1499,2299] y=[-7933,-7133].
#
# INPUTS
#   data_root/2026_final_xenium_analysis/output/sfe/sfe_OTB_2384  (via load_sfe)
#   data_root/2026_final_xenium_analysis/output/06f_reclassification_polarization/
#     reclassified_xenium_scores.csv  (06f cell-label override)
#   data_root/2026_final_xenium_analysis/output/9b_scoring/pathway_gene_sets_v2.csv
#   Shared helpers: config/config.R, spatial/00_setup/00_setup.R (load_sfe, ref_palette).
#
# OUTPUTS
#   figures_dir/figure5/roi_exploration_OTB_2384_secb/OTB_2384_roi06_*.{svg,png}
#
# MANUSCRIPT PANEL(S): Fig 5B (cell type), Fig 5C (4 expression/UCell maps)
# RUNTIME TIER: moderate (SFE load + on-ROI UCell)
# ============================================================================

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

.here     <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))
source(file.path(.here, "..", "..", "spatial", "00_setup", "00_setup.R"))

suppressPackageStartupMessages({
  library(SpatialFeatureExperiment); library(SummarizedExperiment)
  library(data.table); library(ggplot2); library(sf)
  library(UCell); library(Matrix)
})

set.seed(CFG$seed)

FIG_DIR <- cfg_path("figures_dir", "figure5", "roi_exploration_OTB_2384_secb")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
SAMPLE  <- "sfe_OTB_2384"

EPI_TYPES <- c("SecA epithelium", "Intermediate epithelium",
               "SecB epithelium", "Ciliated epithelium")

BB <- c(xmin = 1499, xmax = 2299, ymin = -7933, ymax = -7133)

# Beige-burgundy gradient
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
# Rename-mismatch fix (idempotent): the 06f cache (and SFE) still carry the
# legacy "Transitioning epithelium"; standardize on the SAME 'lab' vector that
# EPI_TYPES filtering (~line 100) and the cell_label colour scale key on.
lab[lab == "Transitioning epithelium"] <- "Intermediate epithelium"
sfe$cell_label <- lab

# ── Subset to ROI ────────────────────────────────────────────────────────
co <- spatialCoords(sfe)
in_roi <- co[, 1] >= BB["xmin"] & co[, 1] <= BB["xmax"] &
          co[, 2] >= BB["ymin"] & co[, 2] <= BB["ymax"]
sfe_roi <- sfe[, in_roi]
cat(sprintf("  ROI cells: %d\n", ncol(sfe_roi)))

# ── UCell glycolysis scoring on ROI ──────────────────────────────────────
pw_dt <- fread(cfg_path("data_root", "2026_final_xenium_analysis", "output",
                        "9b_scoring", "pathway_gene_sets_v2.csv"))
glyc_genes <- intersect(pw_dt[pathway == "glycolysis"]$gene, rownames(sfe_roi))
pw_sets <- list(glycolysis = glyc_genes)

mat <- as(assay(sfe_roi, "counts"), "dgCMatrix")
cat(sprintf("[ucell] scoring %d cells ...\n", ncol(mat)))
ucell_scores <- ScoreSignatures_UCell(matrix = mat,
                                       features = pw_sets,
                                       maxRank  = nrow(mat),
                                       name     = "",
                                       chunk.size = 500)
ucell_dt <- as.data.table(ucell_scores, keep.rownames = "cell_id")
m_uc <- match(colnames(sfe_roi), ucell_dt$cell_id)
sfe_roi$pw_glycolysis <- ucell_dt$glycolysis[m_uc]

# ── Extract polygons ─────────────────────────────────────────────────────
cseg <- colGeometry(sfe_roi, "cellSeg")
cseg$cell_label            <- as.character(sfe_roi$cell_label)
cseg$pw_glycolysis         <- as.numeric(sfe_roi$pw_glycolysis)
cseg$pathway_proliferation <- as.numeric(colData(sfe_roi)[["pathway_proliferation"]])
cseg$pathway_nfkb          <- as.numeric(colData(sfe_roi)[["pathway_nfkb"]])
cseg$LDHA                  <- as.numeric(assay(sfe_roi, "counts")["LDHA", ])
cseg$nucleus_area          <- as.numeric(sfe_roi$nucleus_area)
cseg$is_epi                <- cseg$cell_label %in% EPI_TYPES
cseg_sf <- sf::st_as_sf(cseg)

cat(sprintf("  Epithelial: %d / %d\n", sum(cseg_sf$is_epi), nrow(cseg_sf)))

# ── Panel dimensions ─────────────────────────────────────────────────────
w_in <- 3.8
h_in <- 3.0

# ── 1. Cell type panel (Fig 5B) ──────────────────────────────────────────
cat("\n--- Cell type ---\n")
p_ct <- ggplot(cseg_sf) +
  geom_sf(aes(fill = cell_label), colour = "grey50", linewidth = 0.10) +
  scale_fill_manual(values = ref_palette, name = "Cell type", drop = TRUE) +
  coord_sf(xlim = c(BB["xmin"], BB["xmax"]),
           ylim = c(BB["ymin"], BB["ymax"]), expand = FALSE) +
  theme_void(base_size = 6) +
  theme(legend.position = "right",
        legend.text  = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.4, "lines"),
        plot.margin = margin(2, 2, 2, 2)) +
  guides(fill = guide_legend(override.aes = list(linewidth = 0)))

ggsave(file.path(FIG_DIR, "OTB_2384_roi06_celltype.svg"),
       p_ct, width = w_in + 1.0, height = h_in, bg = "transparent")
ggsave(file.path(FIG_DIR, "OTB_2384_roi06_celltype.png"),
       p_ct, width = w_in + 1.0, height = h_in, dpi = 450,
       bg = "transparent", device = ragg::agg_png)
cat("  OK\n")

# ── Epithelial-only feature panels (Fig 5C) ──────────────────────────────
epi_sf   <- cseg_sf[cseg_sf$is_epi, ]
other_sf <- cseg_sf[!cseg_sf$is_epi, ]

nfkb_label <- "NF-κB"

PANEL_SPECS <- list(
  list(col = "pw_glycolysis",       label = "Glycolysis",       filename = "glycolysis",
       fmt = function(x) round(x, 2), italic = FALSE),
  list(col = "pathway_proliferation", label = "Proliferation",  filename = "proliferation",
       fmt = function(x) round(x, 2), italic = FALSE),
  list(col = "pathway_nfkb",        label = nfkb_label,         filename = "nfkb",
       fmt = function(x) round(x, 2), italic = FALSE),
  list(col = "LDHA",                label = "LDHA",             filename = "LDHA",
       fmt = function(x) round(x, 1), italic = TRUE),
  list(col = "nucleus_area",        label = "Nuclear area",     filename = "nuclear_area",
       fmt = function(x) round(x),    italic = FALSE)
)

for (spec in PANEL_SPECS) {
  cat(sprintf("\n--- %s ---\n", spec$label))

  vals <- epi_sf[[spec$col]]
  qlim <- quantile(vals, c(0.02, 0.98), na.rm = TRUE)
  vmin <- as.numeric(qlim[1]); vmax <- as.numeric(qlim[2])
  if (vmin == vmax) { vmin <- min(vals, na.rm = TRUE); vmax <- max(vals, na.rm = TRUE) }
  if (vmax == vmin) { vmin <- 0; vmax <- 1 }

  epi_sf$val <- pmin(pmax(vals, vmin), vmax)
  plot_epi <- epi_sf[order(epi_sf$val), ]

  title_face <- if (spec$italic) "bold.italic" else "bold"

  p <- ggplot() +
    geom_sf(data = other_sf, fill = "#E8E8E8", colour = "grey70",
            linewidth = 0.06) +
    geom_sf(data = plot_epi, aes(fill = val), colour = "grey55",
            linewidth = 0.08) +
    scale_fill_gradientn(
      colours = GRAD_COLS,
      values  = GRAD_STOPS,
      name    = spec$label,
      limits  = c(vmin, vmax),
      breaks  = pretty(c(vmin, vmax), n = 3),
      labels  = spec$fmt,
      oob     = scales::squish,
      guide   = guide_colorbar(
        barwidth        = unit(0.4, "cm"),
        barheight       = unit(2.2, "cm"),
        ticks           = TRUE,
        ticks.linewidth = unit(0.6, "mm"),
        frame.colour    = "grey30",
        frame.linewidth = unit(0.5, "mm"),
        ticks.colour    = "grey30",
        title.position  = "top",
        title.hjust     = 0.5,
        title.theme     = element_text(size = 11, face = title_face,
                                        colour = "grey10"),
        label.theme     = element_text(size = 9, colour = "grey10",
                                        face = "bold"))) +
    coord_sf(xlim = c(BB["xmin"], BB["xmax"]),
             ylim = c(BB["ymin"], BB["ymax"]), expand = FALSE) +
    theme_void(base_size = 8) +
    theme(
      legend.position      = "right",
      legend.justification = c(0, 1),
      legend.margin        = margin(0, 0, 0, 4),
      legend.box.margin    = margin(0, 0, 0, 0),
      plot.margin          = margin(4, 4, 4, 4)
    )

  svg_out <- file.path(FIG_DIR, sprintf("OTB_2384_roi06_%s.svg", spec$filename))
  png_out <- file.path(FIG_DIR, sprintf("OTB_2384_roi06_%s.png", spec$filename))
  ggsave(svg_out, p, width = w_in, height = h_in, bg = "transparent")
  ggsave(png_out, p, width = w_in, height = h_in, dpi = 450,
         bg = "transparent", device = ragg::agg_png)
  cat(sprintf("  SVG: %s\n  PNG: %s\n", basename(svg_out), basename(png_out)))
}

cat("\n=== ALL PANELS DONE ===\n")
