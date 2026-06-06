#!/usr/bin/env Rscript
# ============================================================================
# Figure 6D (SecB) — TMA gradient-core SecB-epithelium highlight panels
# ----------------------------------------------------------------------------
# PURPOSE
#   SecB-epithelium highlight SVG/PNG panels for TMA cores 2 and 161 (background
#   cells grey, SecB epithelium highlighted). Sibling of the immune-highlight
#   script; together they make the highlighted-population row of Fig 6D.
#   Requires gradient_metrics.csv from 01_tma_hypoxia_gradient_cores.R.
#
# INPUTS
#   data_root/2026_final_xenium_analysis/output/sfe/sfe_tma_filtered (load_sfe)
#   data_root/2026_final_xenium_analysis/output/06f_reclassification_polarization/
#     reclassified_xenium_scores.csv
#   figures_dir/figure6/tma_hypoxia_gradient_cores/gradient_metrics.csv
#   Shared helpers: config/config.R, spatial/00_setup/00_setup.R (load_sfe, ref_palette).
#
# OUTPUTS
#   figures_dir/figure6/tma_hypoxia_gradient_cores/core{ID}_{patient}_secb.{svg,png}
#
# MANUSCRIPT PANEL(S): Fig 6D (SecB highlight, core 2)
# RUNTIME TIER: moderate (SFE load)
# ============================================================================

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

.here     <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))
source(file.path(.here, "..", "..", "spatial", "00_setup", "00_setup.R"))

suppressPackageStartupMessages({
  library(sf); library(data.table); library(ggplot2)
  library(ragg); library(svglite)
})

OUT_DIR <- cfg_path("figures_dir", "figure6", "tma_hypoxia_gradient_cores")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
BG_COL  <- "#E8E8E8"

CORE_IDS <- c("2", "161")

# --- Load SFE + 06f --------------------------------------------------------
message("[1/2] Loading sfe_tma_filtered ...")
sfe <- load_sfe("sfe_tma_filtered")

f06f <- cfg_path("data_root", "2026_final_xenium_analysis", "output",
                 "06f_reclassification_polarization",
                 "reclassified_xenium_scores.csv")
recl <- fread(f06f, select = c("sample", "barcode_orig", "cell_label_06f"))
sample_recl <- recl[sample == "sfe_tma_filtered"]
m <- match(colnames(sfe), sample_recl$barcode_orig)
hit <- !is.na(m)
lab <- as.character(sfe$cell_label)
lab[hit] <- as.character(sample_recl$cell_label_06f[m[hit]])
if ("singler_label" %in% colnames(colData(sfe))) {
  needs_fallback <- is.na(lab) | lab == "" | lab == "NA"
  singler <- as.character(sfe$singler_label)
  lab[needs_fallback & !is.na(singler)] <- singler[needs_fallback & !is.na(singler)]
}
sfe$cell_label <- lab
rm(recl, sample_recl)

grad_met <- fread(file.path(OUT_DIR, "gradient_metrics.csv"))

# --- Render -----------------------------------------------------------------
message("[2/2] Rendering SecB highlight panels ...")

for (core in CORE_IDS) {
  pid <- grad_met[core_id == core, patient_id][1]
  message(sprintf("\n  Core %s (patient %s)", core, pid))

  core_mask <- sfe$core_id == core & !is.na(sfe$core_id)
  sfe_core  <- sfe[, core_mask]
  cseg_sf   <- colGeometry(sfe_core, "cellSeg")
  cseg_sf$cell_label <- as.character(sfe_core$cell_label)

  bbox_data <- sf::st_bbox(cseg_sf)
  mg <- 20
  BBOX <- c(xmin = as.numeric(bbox_data["xmin"]) - mg,
            xmax = as.numeric(bbox_data["xmax"]) + mg,
            ymin = as.numeric(bbox_data["ymin"]) - mg,
            ymax = as.numeric(bbox_data["ymax"]) + mg)

  aspect  <- (BBOX["xmax"] - BBOX["xmin"]) / (BBOX["ymax"] - BBOX["ymin"])
  panel_h <- 70
  panel_w <- panel_h * as.numeric(aspect) + 12
  w_in    <- panel_w / 25.4
  h_in    <- panel_h / 25.4

  is_secb <- cseg_sf$cell_label == "SecB epithelium"
  n_secb  <- sum(is_secb, na.rm = TRUE)
  message(sprintf("    %d SecB epithelium cells", n_secb))

  bg_sf <- cseg_sf[!is_secb | is.na(is_secb), ]
  fg_sf <- cseg_sf[which(is_secb), ]

  secb_col <- unname(ref_palette["SecB epithelium"])

  p <- ggplot() +
    geom_sf(data = bg_sf, fill = BG_COL, colour = "grey70", linewidth = 0.06) +
    geom_sf(data = fg_sf, fill = secb_col, colour = "grey55", linewidth = 0.10) +
    coord_sf(xlim = c(BBOX["xmin"], BBOX["xmax"]),
             ylim = c(BBOX["ymin"], BBOX["ymax"]),
             expand = FALSE, default_crs = NULL) +
    theme_void(base_size = 6) +
    theme(plot.margin = margin(4, 4, 4, 4))

  ggsave(file.path(OUT_DIR, sprintf("core%s_%s_secb.svg", core, pid)),
         p, width = w_in, height = h_in, bg = "transparent")
  ggsave(file.path(OUT_DIR, sprintf("core%s_%s_secb.png", core, pid)),
         p, width = w_in, height = h_in, dpi = 450,
         bg = "transparent", device = ragg::agg_png, limitsize = FALSE)
  message("    saved secb SVG + PNG")
}

message("\nDone.")
