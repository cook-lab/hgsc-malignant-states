#!/usr/bin/env Rscript
# ============================================================================
# Figure 6D (immune) — TMA gradient-core highlighted-population panels
# ----------------------------------------------------------------------------
# PURPOSE
#   SVG + PNG panels for TMA cores 2 and 161: background cells grey, a single
#   immune cell type highlighted (macrophage, T cell, plasma, NK, B cell).
#   Also emits the matching hypoxia and cell-type panels as SVG/PNG.
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
#   figures_dir/figure6/tma_hypoxia_gradient_cores/
#     core{ID}_{patient}_{macrophage,tcell,plasmacell,nkcell,bcell,hypoxia,celltype}.{svg,png}
#
# MANUSCRIPT PANEL(S): Fig 6D (highlighted immune populations, core 2)
# RUNTIME TIER: moderate (SFE load + per-core rendering)
# ============================================================================

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

.here     <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))
source(file.path(.here, "..", "..", "spatial", "00_setup", "00_setup.R"))

suppressPackageStartupMessages({
  library(sf)
  library(data.table)
  library(ggplot2)
  library(ragg)
  library(svglite)
  library(scales)
})

OUT_DIR <- cfg_path("figures_dir", "figure6", "tma_hypoxia_gradient_cores")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CORE_IDS <- c("2", "161")

HIGHLIGHTS <- list(
  list(ct = "Macrophage",  col = unname(ref_palette["Macrophage"]),
       suffix = "macrophage",  subtitle = "Macrophages highlighted"),
  list(ct = "T cell",      col = unname(ref_palette["T cell"]),
       suffix = "tcell",       subtitle = "T cells highlighted"),
  list(ct = "Plasma cell", col = unname(ref_palette["Plasma cell"]),
       suffix = "plasmacell",  subtitle = "Plasma cells highlighted"),
  list(ct = "NK cell",     col = unname(ref_palette["NK cell"]),
       suffix = "nkcell",      subtitle = "NK cells highlighted"),
  list(ct = "B cell",      col = unname(ref_palette["B cell"]),
       suffix = "bcell",       subtitle = "B cells highlighted")
)

BG_COL <- "#E8E8E8"

GRAD_COLS  <- c("#F6EFE5", "#ECDDD0", "#D89E97", "#A03A4A", "#3A111A")
GRAD_STOPS <- c(0.00,      0.28,      0.55,      0.80,      1.00)

# --- Load SFE ---------------------------------------------------------------
message("[1/3] Loading sfe_tma_filtered ...")
sfe <- load_sfe("sfe_tma_filtered")

# --- Apply 06f reclassification ---------------------------------------------
message("[2/3] Applying 06f reclassification ...")

f06f <- cfg_path("data_root", "2026_final_xenium_analysis", "output",
                 "06f_reclassification_polarization",
                 "reclassified_xenium_scores.csv")
if (file.exists(f06f)) {
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
} else {
  message("  WARNING: 06f file not found, using cell_label as-is")
}

# Standardize the epithelial polarization label: the deposited SFE (and the
# 06f override) carry the legacy literal "Transitioning epithelium", which is
# absent from ref_palette — so those cells would render with no fill / drop
# from the legend. Rename to the canonical "Intermediate epithelium" before
# colouring (consistent with sibling scripts 06/09 and 00_setup ref_palette).
lab <- as.character(sfe$cell_label)
n_trans <- sum(lab == "Transitioning epithelium", na.rm = TRUE)
lab[lab == "Transitioning epithelium"] <- "Intermediate epithelium"
sfe$cell_label <- lab
message(sprintf("  renamed Transitioning -> Intermediate: %d cells", n_trans))

# --- Gradient metrics --------------------------------------------------------
grad_met <- fread(file.path(OUT_DIR, "gradient_metrics.csv"))

# --- Render panels -----------------------------------------------------------
message("[3/3] Rendering SVG + PNG panels ...")

for (core in CORE_IDS) {
  pid  <- grad_met[core_id == core, patient_id][1]
  grad <- grad_met[core_id == core, grad_max][1]

  message(sprintf("\n  === Core %s (patient %s, grad=%.4f) ===", core, pid, grad))

  core_mask <- sfe$core_id == core & !is.na(sfe$core_id)
  if (sum(core_mask) < 50) {
    message("    Skipping: too few cells"); next
  }

  sfe_core <- sfe[, core_mask]
  cseg_sf  <- colGeometry(sfe_core, "cellSeg")
  cseg_sf$cell_label      <- as.character(sfe_core$cell_label)
  cseg_sf$pathway_hypoxia <- as.numeric(sfe_core$pathway_hypoxia)

  bbox_data <- sf::st_bbox(cseg_sf)
  mg <- 20
  BBOX <- c(
    xmin = as.numeric(bbox_data["xmin"]) - mg,
    xmax = as.numeric(bbox_data["xmax"]) + mg,
    ymin = as.numeric(bbox_data["ymin"]) - mg,
    ymax = as.numeric(bbox_data["ymax"]) + mg
  )

  aspect  <- (BBOX["xmax"] - BBOX["xmin"]) / (BBOX["ymax"] - BBOX["ymin"])
  panel_h <- 70
  panel_w <- panel_h * as.numeric(aspect) + 12
  w_in    <- panel_w / 25.4
  h_in    <- panel_h / 25.4

  # ── Hypoxia panel ──────────────────────────────────────────────────────
  message("    Hypoxia ...")
  vals <- cseg_sf$pathway_hypoxia
  vals[is.na(vals)] <- 0
  qlim <- quantile(vals, c(0.02, 0.98), na.rm = TRUE)
  vmin <- as.numeric(qlim[1]); vmax <- as.numeric(qlim[2])
  if (vmin == vmax) { vmin <- vmin - 0.01; vmax <- vmax + 0.01 }

  cseg_sf$val <- pmin(pmax(vals, vmin), vmax)
  cseg_ordered <- cseg_sf[order(cseg_sf$val), ]

  p_hyp <- ggplot(cseg_ordered) +
    geom_sf(aes(fill = val), colour = "grey55", linewidth = 0.10) +
    scale_fill_gradientn(
      colours = GRAD_COLS, values = GRAD_STOPS,
      name    = "Hypoxia\nscore",
      limits  = c(vmin, vmax),
      breaks  = c(vmin, (vmin + vmax) / 2, vmax),
      labels  = function(x) format(round(x, 2), nsmall = 2),
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
        title.theme     = element_text(size = 11, face = "bold",
                                        colour = "grey10"),
        label.theme     = element_text(size = 9, colour = "grey10",
                                        face = "bold"))
    ) +
    coord_sf(xlim = c(BBOX["xmin"], BBOX["xmax"]),
             ylim = c(BBOX["ymin"], BBOX["ymax"]),
             expand = FALSE, default_crs = NULL) +
    theme_void(base_size = 6) +
    theme(legend.position      = "right",
          legend.justification = c(0, 1),
          legend.margin        = margin(0, 0, 0, 4),
          legend.box.margin    = margin(0, 0, 0, 0),
          plot.margin          = margin(4, 4, 4, 4))

  ggsave(file.path(OUT_DIR, sprintf("core%s_%s_hypoxia.svg", core, pid)),
         p_hyp, width = w_in, height = h_in, bg = "transparent")
  ggsave(file.path(OUT_DIR, sprintf("core%s_%s_hypoxia.png", core, pid)),
         p_hyp, width = w_in, height = h_in, dpi = 450,
         bg = "transparent", device = ragg::agg_png, limitsize = FALSE)
  message("      saved hypoxia SVG + PNG")

  # ── Cell type panel ────────────────────────────────────────────────────
  message("    Cell type ...")
  palette <- c(ref_palette, "Secretory epithelium" = "#E6A141")

  p_ct <- ggplot(cseg_sf) +
    geom_sf(aes(fill = cell_label), colour = "grey55", linewidth = 0.10) +
    scale_fill_manual(values = palette, name = "Cell type",
                      guide = guide_legend(
                        override.aes = list(linewidth = 0),
                        ncol = 1, keyheight = unit(0.4, "cm"))) +
    coord_sf(xlim = c(BBOX["xmin"], BBOX["xmax"]),
             ylim = c(BBOX["ymin"], BBOX["ymax"]),
             expand = FALSE, default_crs = NULL) +
    theme_void(base_size = 6) +
    theme(legend.position   = "right",
          legend.text       = element_text(size = 6),
          legend.title      = element_text(size = 7),
          legend.key.size   = unit(0.4, "lines"),
          legend.margin     = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 1),
          plot.margin       = margin(4, 4, 4, 4))

  w_ct <- (panel_w + 25) / 25.4
  ggsave(file.path(OUT_DIR, sprintf("core%s_%s_celltype.svg", core, pid)),
         p_ct, width = w_ct, height = h_in, bg = "transparent")
  ggsave(file.path(OUT_DIR, sprintf("core%s_%s_celltype.png", core, pid)),
         p_ct, width = w_ct, height = h_in, dpi = 450,
         bg = "transparent", device = ragg::agg_png, limitsize = FALSE)
  message("      saved celltype SVG + PNG")

  # ── Immune highlight panels ────────────────────────────────────────────
  for (hl in HIGHLIGHTS) {
    message(sprintf("    %s ...", hl$ct))

    is_highlight <- cseg_sf$cell_label == hl$ct
    n_hl <- sum(is_highlight, na.rm = TRUE)
    message(sprintf("      %d %s cells", n_hl, hl$ct))

    bg_sf <- cseg_sf[!is_highlight | is.na(is_highlight), ]
    fg_sf <- cseg_sf[which(is_highlight), ]

    p <- ggplot() +
      geom_sf(data = bg_sf, fill = BG_COL, colour = "grey70",
              linewidth = 0.06) +
      geom_sf(data = fg_sf, fill = hl$col, colour = "grey55",
              linewidth = 0.10) +
      coord_sf(xlim = c(BBOX["xmin"], BBOX["xmax"]),
               ylim = c(BBOX["ymin"], BBOX["ymax"]),
               expand = FALSE, default_crs = NULL) +
      theme_void(base_size = 6) +
      theme(plot.margin = margin(4, 4, 4, 4))

    ggsave(file.path(OUT_DIR, sprintf("core%s_%s_%s.svg", core, pid, hl$suffix)),
           p, width = w_in, height = h_in, bg = "transparent")
    ggsave(file.path(OUT_DIR, sprintf("core%s_%s_%s.png", core, pid, hl$suffix)),
           p, width = w_in, height = h_in, dpi = 450,
           bg = "transparent", device = ragg::agg_png, limitsize = FALSE)
    message(sprintf("      saved %s SVG + PNG", hl$suffix))
  }
}

message("\nDone.")
