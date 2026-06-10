#!/usr/bin/env Rscript
# ============================================================================
# Figure 6K,6L — SP24_24824 ROI cell-type + gene maps (whole-tissue investigation)
# ----------------------------------------------------------------------------
# PURPOSE
#   Per-WT-sample spatial gene panels (light-grey-floor cellSeg polygons,
#   Voyager fill). One sample per invocation, optionally cropped to a
#   rectangular ROI and/or restricted to a gene subset. For the published
#   Fig 6K/6L the relevant ROIs of SP24_24824 are rendered for CTSL / MMP7 /
#   ICAM1 (plus cell-type context); see Usage.
#
#   *** REQUIRES A <sample> CLI ARGUMENT *** (and optional ROI bbox / label /
#   gene subset). The script renders nothing without a sample. See Usage.
#
# INPUTS
#   data_root/2026_final_xenium_analysis/output/sfe/sfe_<sample> (load_sfe)
#     samples = CFG$cohort$whole_tissue (published 8 whole-tissue).
#   06f polarization reclassification is baked into the SFE cell_label.
#   Assay: logcounts.
#   Shared helpers: config/config.R, spatial/00_setup/00_setup.R (load_sfe).
#
# OUTPUTS
#   figures_dir/figure6/roi_whole_tissue_investigation/
#     {sample}_{label}_{gene}.{png,svg}      (gene-expression maps; Fig 6L)
#     {sample}_{label}_celltype.{png,svg}    (cell-type map;        Fig 6K)
#
# MANUSCRIPT PANEL(S): Fig 6K (SP24_24824 ROI cell-type, 2 boxes),
#                      Fig 6L (SP24_24824 ROI CTSL/MMP7/ICAM1, 2 boxes)
# RUNTIME TIER: heavy (full 173-gene loop unless a gene subset is supplied)
#
# USAGE
#   Rscript 09_xenium_roi_whole_tissue_investigation.R <sample>
#     -> whole-tissue render of all panel genes for <sample>
#   Rscript ... <sample> <xmin,xmax,ymin,ymax> <label>
#     -> ROI render, stored as {sample}_{label}_{gene}.png
#   Rscript ... <sample> whole "" CTSL,MMP7,ICAM1
#     -> whole-tissue render of a gene subset (arg 4 comma-separated)
#   Rscript ... SP24_24824 <xmin,xmax,ymin,ymax> upper_right CTSL,MMP7,ICAM1
#     -> the Fig 6L production call (per-ROI gene subset)
#   Rscript ... SP24_24824 <xmin,xmax,ymin,ymax> upper_right celltype
#     -> the Fig 6K production call: cell-type map coloured by `cell_label`
#        with ref_palette (instead of gene expression). The reserved gene arg
#        "celltype" switches the renderer; the gene-expression path is unchanged.
# ============================================================================

.here     <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))
source(file.path(.here, "..", "..", "spatial", "00_setup", "00_setup.R"))

suppressPackageStartupMessages({
  library(Voyager); library(data.table); library(ggplot2); library(scales)
  library(SummarizedExperiment)
})

OUT_DIR <- cfg_path("figures_dir", "figure6", "roi_whole_tissue_investigation")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

WT_SAMPLES <- CFG$cohort$whole_tissue

# Full 173-gene panel (union of worked-on + doc-mentioned genes, intersected
# with rownames(sfe) at run time).
GENES <- c(
  "A2M", "ACTA2", "AKT1", "ALDOA", "APC", "AXL", "BAX", "BCAM",
  "BCL2", "BIRC3", "BTLA", "C1QC", "C1S", "C3", "C7", "CALR",
  "CASP3", "CASP8", "CCNB1", "CCR5", "CCR7", "CD14", "CD226", "CD24",
  "CD274", "CD276", "CD4", "CD40", "CD44", "CD55", "CD68", "CD69",
  "CD79A", "CD86", "CDC20", "CDH2", "CDK1", "CDKN2A", "CIITA", "CLDN3",
  "COTL1", "CPT1A", "CTHRC1", "CTNNB1", "CTSL", "CTSS", "CXCL10", "CXCL11",
  "CXCL3", "CXCL9", "CXCR4", "CXCR6", "DCN", "DKK3", "DLL4", "EGFR",
  "ENG", "ENO1", "ENTPD1", "EPAS1", "EPCAM", "F3", "FAS", "FASLG",
  "FBXO21", "FCGR1A", "FCGR3A", "FGFR1", "GPX3", "GZMA", "HAVCR2", "HIF1A",
  "HSPG2", "ICAM1", "IDO1", "IFI44L", "IFI6", "IFIT1", "IFIT3", "IFNG",
  "IFNGR1", "IL10", "IL2RG", "IL6", "IL6R", "IL7R", "INHBA", "IRF1",
  "ISG20", "ITGAM", "ITGAV", "ITGB2", "ITGB5", "JAG2", "KRT17", "KRT19",
  "KRT7", "KRT8", "LAG3", "LCN2", "LDHA", "LGALS3", "LGALS9", "LGR5",
  "LPAR3", "LTB", "MDM2", "MDM4", "MECOM", "MET", "MKI67", "MMP11",
  "MMP7", "MRC1", "MS4A1", "MUC1", "MX1", "MYC", "NCR1", "NDRG1",
  "NECTIN4", "NFKB1", "NFKB2", "NFKBIA", "NOTCH1", "NOTCH2", "NOTCH3", "NRP2",
  "NT5E", "PBX1", "PCNA", "PDCD1", "PDGFB", "PDGFRA", "PDK1", "PECAM1",
  "PGK1", "PIK3R1", "POSTN", "PRF1", "PRKCB", "PRSS22", "PTEN", "RARRES2",
  "RCN2", "RUNX1", "SERPINE1", "SFRP4", "SLC16A3", "SLC2A1", "SLPI", "SOX17",
  "STAT1", "STMN1", "TACSTD2", "TAGLN", "TAP1", "TAPBP", "TCF7", "TGFB1",
  "TGFBR2", "TGM2", "TIGIT", "TNF", "TOP2A", "TP53", "TREM2", "VCAM1",
  "VCAN", "VEGFA", "VWF", "WNT7A", "XBP1"
)

# Light-gray-floor -> dusty rose -> burgundy
expr_cols <- c("#EDEDED", "#F0D8D8", "#C68B8B", "#7A2B35")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Pass a sample, e.g. SP24_24824.  See script header (USAGE) for full CLI.")
}
SAMPLE_ARG <- args[1]
ROI_ARG    <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "whole"
LABEL_ARG  <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "whole"
GENE_ARG   <- if (length(args) >= 4 && nzchar(args[4])) args[4] else NA_character_

# Reserved gene arg "celltype" switches to the Fig 6K cell-type renderer.
CELLTYPE_MODE <- !is.na(GENE_ARG) && tolower(GENE_ARG) == "celltype"

sfe_name <- if (startsWith(SAMPLE_ARG, "sfe_")) SAMPLE_ARG else paste0("sfe_", SAMPLE_ARG)
base_sample <- sub("^sfe_", "", sfe_name)
if (!base_sample %in% WT_SAMPLES) {
  warning("Sample not in canonical WT list; proceeding: ", base_sample)
}

if (!is.na(GENE_ARG) && !CELLTYPE_MODE) {
  sel <- strsplit(GENE_ARG, ",")[[1]]
  not_in_union <- setdiff(sel, GENES)
  if (length(not_in_union)) {
    cat("!! requested genes not in union list (will still try if on panel):\n   ",
        paste(not_in_union, collapse = ", "), "\n")
  }
  GENES <- sel
}

if (CELLTYPE_MODE) {
  cat(sprintf("Sample : %s\nROI    : %s\nLabel  : %s\nMode   : CELLTYPE (Fig 6K)\n",
              sfe_name, ROI_ARG, LABEL_ARG))
} else {
  cat(sprintf("Sample : %s\nROI    : %s\nLabel  : %s\nGenes  : %d (%s%s)\n",
              sfe_name, ROI_ARG, LABEL_ARG, length(GENES),
              paste(head(GENES, 5), collapse = ", "),
              ifelse(length(GENES) > 5, ", ...", "")))
}

# --- Load SFE ---------------------------------------------------------------
cat(sprintf("Loading %s ...\n", sfe_name))
sfe <- load_sfe(sfe_name)
co  <- spatialCoords(sfe)
cat(sprintf("  %d cells  x=[%.0f,%.0f]  y=[%.0f,%.0f]\n",
            ncol(sfe), min(co[,1]), max(co[,1]), min(co[,2]), max(co[,2])))

# --- Optional ROI crop ------------------------------------------------------
if (ROI_ARG != "whole") {
  bb <- as.numeric(strsplit(ROI_ARG, ",")[[1]])
  if (length(bb) != 4) stop("ROI must be xmin,xmax,ymin,ymax")
  names(bb) <- c("xmin","xmax","ymin","ymax")
  in_roi <- co[,1] >= bb["xmin"] & co[,1] <= bb["xmax"] &
            co[,2] >= bb["ymin"] & co[,2] <= bb["ymax"]
  sfe <- sfe[, in_roi]
  co  <- spatialCoords(sfe)
  cat(sprintf("  ROI (%s) bbox=%s -> %d cells\n",
              LABEL_ARG, paste(bb, collapse=","), ncol(sfe)))
}

if (!CELLTYPE_MODE) {
  missing_genes <- setdiff(GENES, rownames(sfe))
  if (length(missing_genes) > 0) {
    cat("!! genes not in panel (skipped):", paste(missing_genes, collapse=", "), "\n")
    GENES <- setdiff(GENES, missing_genes)
  }
}

aspect <- (max(co[,1]) - min(co[,1])) / (max(co[,2]) - min(co[,2]))
LONG_IN <- 10
if (aspect >= 1) { w_in <- LONG_IN; h_in <- LONG_IN / aspect } else {
                   h_in <- LONG_IN; w_in <- LONG_IN * aspect }
w_in <- w_in + 1.5  # room for legend
cat(sprintf("  figure size: %.1f x %.1f in   aspect %.2f\n", w_in, h_in, aspect))

render_gene <- function(sfe_x, gene) {
  feat <- paste0(gene, "_log")
  colData(sfe_x)[[feat]] <- as.numeric(assay(sfe_x, "logcounts")[gene, ])
  v <- colData(sfe_x)[[feat]]
  nz <- sum(v > 0, na.rm = TRUE)
  lims <- as.numeric(quantile(v, c(0.01, 0.99), na.rm = TRUE))
  if (lims[1] == lims[2]) lims[2] <- lims[1] + 1e-6
  cat(sprintf("  [%s] nonzero %d / %d (%.1f%%)  lim %.2f..%.2f\n",
              gene, nz, length(v), 100*nz/length(v), lims[1], lims[2]))

  p <- plotSpatialFeature(sfe_x, feat, colGeometryName = "cellSeg",
                          aes_use = "fill", linewidth = 0.1,
                          color = "grey50") +
    scale_fill_gradientn(colors = expr_cols,
                          name = sprintf("%s\n(log counts)", gene),
                          limits = lims, oob = squish,
                          na.value = "#BBBBBB") +
    theme_void(base_size = 10) +
    theme(legend.position   = "right",
          legend.text       = element_text(size = 8, color = "black"),
          legend.title      = element_text(size = 9, color = "black"),
          legend.key.height = unit(1, "lines"),
          plot.background   = element_rect(fill = "white", color = NA),
          panel.background  = element_rect(fill = "white", color = NA),
          legend.background = element_rect(fill = "white", color = NA),
          legend.key        = element_rect(fill = "white", color = NA))

  stem <- sprintf("%s_%s_%s", base_sample, LABEL_ARG, gene)
  ggsave(file.path(OUT_DIR, paste0(stem, ".png")), p,
         width = w_in, height = h_in, dpi = 450, bg = "white",
         limitsize = FALSE)
  ggsave(file.path(OUT_DIR, paste0(stem, ".svg")), p,
         width = w_in, height = h_in, bg = "white", limitsize = FALSE)
}

# --- Cell-type render (Fig 6K) ----------------------------------------------
# Colour the cropped cells by `cell_label` with ref_palette, matching the
# Fig 4D / SF10C cellSeg-polygon styling. The deposited SFE carries the legacy
# literal "Transitioning epithelium"; rename to the standardized "Intermediate
# epithelium" so values match ref_palette (otherwise those cells would render
# with no fill). The gene-expression path above is untouched.
if (CELLTYPE_MODE) {
  lab <- as.character(sfe$cell_label)
  n_trans <- sum(lab == "Transitioning epithelium", na.rm = TRUE)
  lab[lab == "Transitioning epithelium"] <- "Intermediate epithelium"
  sfe$cell_label <- lab
  cat(sprintf("  renamed Transitioning -> Intermediate: %d cells\n", n_trans))
  cat("  cell_label tally:\n")
  print(sort(table(lab), decreasing = TRUE))

  p <- plotSpatialFeature(sfe, "cell_label", colGeometryName = "cellSeg",
                          aes_use = "fill", linewidth = 0.15, color = "grey50") +
    scale_fill_manual(values = ref_palette, name = "Cell type", drop = TRUE) +
    theme_void(base_size = 10) +
    theme(legend.position  = "right",
          legend.text      = element_text(size = 8, color = "black"),
          legend.title     = element_text(size = 9, color = "black"),
          legend.key.size  = unit(0.5, "lines"),
          plot.background  = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA)) +
    guides(fill = guide_legend(override.aes = list(linewidth = 0)))

  stem <- sprintf("%s_%s_celltype", base_sample, LABEL_ARG)
  ggsave(file.path(OUT_DIR, paste0(stem, ".png")), p,
         width = w_in, height = h_in, dpi = 450, bg = "white", limitsize = FALSE)
  ggsave(file.path(OUT_DIR, paste0(stem, ".svg")), p,
         width = w_in, height = h_in, bg = "white", limitsize = FALSE)
  cat(sprintf("\n=== Saved cell-type map: %s.{png,svg} ===\n", stem))
  cat("\nDone. Output in:", OUT_DIR, "\n")
  quit(save = "no", status = 0)
}

cat(sprintf("\n=== Rendering %d genes to %s ===\n", length(GENES), OUT_DIR))
for (g in GENES) render_gene(sfe, g)

cat("\nDone. Output in:", OUT_DIR, "\n")
