#!/usr/bin/env Rscript
# ============================================================================
# Figure 6A,6C — TMA hypoxia-gradient cores (cell type + hypoxia maps)
# ----------------------------------------------------------------------------
# PURPOSE
#   Identify TMA cores with strong spatial hypoxia gradients and render
#   per-core panels: pathway_hypoxia (beige-burgundy gradient) and cell-type
#   labels (ref_palette), both over cellSeg polygons. Gradient detection
#   splits each core along 4 axes (L/R, T/B, two diagonals); the axis with
#   the largest mean-hypoxia difference defines the gradient strength. Top 6
#   cores are rendered (Fig 6A uses core 2 / patient 2054).
#
# INPUTS
#   data_root/2026_final_xenium_analysis/output/sfe/sfe_tma_filtered (load_sfe)
#   data_root/2026_final_xenium_analysis/output/06f_reclassification_polarization/
#     reclassified_xenium_scores.csv  (06f cell-label override)
#   Shared helpers: config/config.R, spatial/00_setup/00_setup.R
#     (load_sfe, ref_palette, log_session).
#
# OUTPUTS
#   figures_dir/figure6/tma_hypoxia_gradient_cores/
#     core{ID}_{patient}_hypoxia.png, core{ID}_{patient}_celltype.png,
#     gradient_metrics.csv
#
# MANUSCRIPT PANEL(S): Fig 6A (core 2 cell-type), Fig 6C (core 2 hypoxia)
# RUNTIME TIER: moderate (SFE load + per-core rendering)
# ============================================================================

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

.here     <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))
source(file.path(.here, "..", "..", "spatial", "00_setup", "00_setup.R"))

suppressPackageStartupMessages({
  library(sf)
  library(scales)
  library(grid)
  library(data.table)
})

message("\n", strrep("=", 70))
message("TMA Hypoxia Gradient Cores — ", Sys.Date())
message(strrep("=", 70))

# --- Paths ------------------------------------------------------------------
OUT_DIR <- cfg_path("figures_dir", "figure6", "tma_hypoxia_gradient_cores")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

GRAD_COLS  <- c("#F6EFE5", "#ECDDD0", "#D89E97", "#A03A4A", "#3A111A")
GRAD_STOPS <- c(0.00,      0.28,      0.55,      0.80,      1.00)

N_TOP <- 6  # number of gradient cores to select

# --- Load SFE ---------------------------------------------------------------
message("[1/6] Loading sfe_tma_filtered ...")
sfe <- load_sfe("sfe_tma_filtered")
message(sprintf("  %d cells, %d cores", ncol(sfe), length(unique(sfe$core_id))))

# --- Apply 06f reclassification ---------------------------------------------
message("[2/6] Applying 06f reclassification ...")

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
    message(sprintf("  06f override: %d cells; singler fallback: %d cells",
                    sum(hit), sum(needs_fallback & !is.na(singler))))
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

# --- Compute gradient metrics -----------------------------------------------
message("[3/6] Computing per-core hypoxia gradient metrics ...")

centroids <- colGeometry(sfe, "centroids")
coords <- sf::st_coordinates(sf::st_as_sf(centroids))

dt <- data.table(
  core_id    = as.character(sfe$core_id),
  patient_id = as.character(sfe$patient_id),
  x          = coords[, "X"],
  y          = coords[, "Y"],
  hypoxia    = as.numeric(sfe$pathway_hypoxia)
)

dt <- dt[core_id != "Off core" & !is.na(core_id)]

gradient_results <- dt[, {
  n <- .N
  if (n < 200) {
    list(n_cells = n, grad_x = NA_real_, grad_y = NA_real_,
         grad_d1 = NA_real_, grad_d2 = NA_real_,
         grad_max = NA_real_, mean_hyp = NA_real_, sd_hyp = NA_real_,
         mean_hi = NA_real_, mean_lo = NA_real_, split_dir = NA_character_)
  } else {
    cx <- x - median(x)
    cy <- y - median(y)

    left  <- hypoxia[cx <= 0]; right  <- hypoxia[cx > 0]
    top   <- hypoxia[cy <= 0]; bottom <- hypoxia[cy > 0]
    d1p   <- hypoxia[(cx + cy) > 0]; d1n <- hypoxia[(cx + cy) <= 0]
    d2p   <- hypoxia[(cx - cy) > 0]; d2n <- hypoxia[(cx - cy) <= 0]

    grad_x  <- abs(mean(right, na.rm = TRUE) - mean(left, na.rm = TRUE))
    grad_y  <- abs(mean(bottom, na.rm = TRUE) - mean(top, na.rm = TRUE))
    grad_d1 <- abs(mean(d1p, na.rm = TRUE) - mean(d1n, na.rm = TRUE))
    grad_d2 <- abs(mean(d2p, na.rm = TRUE) - mean(d2n, na.rm = TRUE))

    all_grads <- c(x = grad_x, y = grad_y, d1 = grad_d1, d2 = grad_d2)
    best_dir  <- names(which.max(all_grads))
    grad_max  <- max(all_grads)

    if (best_dir == "x") {
      mean_hi <- max(mean(left, na.rm = TRUE), mean(right, na.rm = TRUE))
      mean_lo <- min(mean(left, na.rm = TRUE), mean(right, na.rm = TRUE))
    } else if (best_dir == "y") {
      mean_hi <- max(mean(top, na.rm = TRUE), mean(bottom, na.rm = TRUE))
      mean_lo <- min(mean(top, na.rm = TRUE), mean(bottom, na.rm = TRUE))
    } else if (best_dir == "d1") {
      mean_hi <- max(mean(d1p, na.rm = TRUE), mean(d1n, na.rm = TRUE))
      mean_lo <- min(mean(d1p, na.rm = TRUE), mean(d1n, na.rm = TRUE))
    } else {
      mean_hi <- max(mean(d2p, na.rm = TRUE), mean(d2n, na.rm = TRUE))
      mean_lo <- min(mean(d2p, na.rm = TRUE), mean(d2n, na.rm = TRUE))
    }

    list(n_cells = n, grad_x = grad_x, grad_y = grad_y,
         grad_d1 = grad_d1, grad_d2 = grad_d2,
         grad_max = grad_max, mean_hyp = mean(hypoxia, na.rm = TRUE),
         sd_hyp = sd(hypoxia, na.rm = TRUE),
         mean_hi = mean_hi, mean_lo = mean_lo, split_dir = best_dir)
  }
}, by = .(core_id, patient_id)]

gradient_results <- gradient_results[!is.na(grad_max)]
gradient_results <- gradient_results[order(-grad_max)]

fwrite(gradient_results, file.path(OUT_DIR, "gradient_metrics.csv"))
message(sprintf("  %d cores scored; saved gradient_metrics.csv", nrow(gradient_results)))

top_cores <- gradient_results[1:N_TOP]

message("\n  Selected cores:")
for (i in 1:N_TOP) {
  r <- top_cores[i]
  message(sprintf("    Core %4s (patient %5s): grad=%.4f  split=%2s  hi_half=%.3f  lo_half=%.3f  sd=%.3f  n=%d",
                  r$core_id, r$patient_id, r$grad_max, r$split_dir,
                  r$mean_hi, r$mean_lo, r$sd_hyp, r$n_cells))
}

# --- Render panels ----------------------------------------------------------
message("\n[4/6] Rendering panels ...")

for (i in 1:N_TOP) {
  core  <- top_cores$core_id[i]
  pid   <- top_cores$patient_id[i]
  grad  <- top_cores$grad_max[i]

  message(sprintf("\n  --- Core %s (patient %s, grad=%.4f) ---", core, pid, grad))

  core_mask <- sfe$core_id == core & !is.na(sfe$core_id)

  if (sum(core_mask) < 50) {
    message("    Skipping: too few cells")
    next
  }

  sfe_core <- sfe[, core_mask]

  cseg_sf <- colGeometry(sfe_core, "cellSeg")
  cseg_sf$cell_label      <- as.character(sfe_core$cell_label)
  cseg_sf$pathway_hypoxia <- as.numeric(sfe_core$pathway_hypoxia)

  bbox_data <- sf::st_bbox(cseg_sf)
  margin <- 20
  BBOX <- c(
    xmin = as.numeric(bbox_data["xmin"]) - margin,
    xmax = as.numeric(bbox_data["xmax"]) + margin,
    ymin = as.numeric(bbox_data["ymin"]) - margin,
    ymax = as.numeric(bbox_data["ymax"]) + margin
  )

  aspect   <- (BBOX["xmax"] - BBOX["xmin"]) / (BBOX["ymax"] - BBOX["ymin"])
  panel_h  <- 70  # mm
  panel_w  <- panel_h * as.numeric(aspect) + 12
  w_in     <- panel_w / 25.4
  h_in     <- panel_h / 25.4

  # --- Hypoxia panel ---
  vals <- cseg_sf$pathway_hypoxia
  vals[is.na(vals)] <- 0
  qlim <- quantile(vals, c(0.02, 0.98), na.rm = TRUE)
  vmin <- as.numeric(qlim[1]); vmax <- as.numeric(qlim[2])
  if (vmin == vmax) { vmin <- vmin - 0.01; vmax <- vmax + 0.01 }

  cseg_sf$val <- pmin(pmax(vals, vmin), vmax)
  cseg_sf <- cseg_sf[order(cseg_sf$val), ]

  p_hyp <- ggplot(cseg_sf) +
    geom_sf(aes(fill = val), colour = "grey55", linewidth = 0.10) +
    scale_fill_gradientn(
      colours = GRAD_COLS, values = GRAD_STOPS,
      name    = "Hypoxia\nscore",
      limits  = c(vmin, vmax),
      breaks  = c(vmin, (vmin + vmax) / 2, vmax),
      labels  = function(x) format(round(x, 2), nsmall = 2),
      oob     = scales::squish,
      guide   = guide_colorbar(
        barwidth     = unit(0.22, "cm"),
        barheight    = unit(1.6, "cm"),
        ticks        = TRUE,
        frame.colour = "grey40",
        ticks.colour = "grey40",
        title.theme  = element_text(size = 5, colour = "grey20"),
        label.theme  = element_text(size = 5, colour = "grey20"))
    ) +
    coord_sf(xlim = c(BBOX["xmin"], BBOX["xmax"]),
             ylim = c(BBOX["ymin"], BBOX["ymax"]),
             expand = FALSE, default_crs = NULL) +
    labs(title = sprintf("Core %s (patient %s)", core, pid),
         subtitle = sprintf("Hypoxia gradient = %.4f", grad)) +
    theme_void(base_size = 6) +
    theme(legend.position   = "right",
          legend.margin     = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 1),
          plot.title        = element_text(size = 7, colour = "black"),
          plot.subtitle     = element_text(size = 5, colour = "grey40"),
          plot.margin       = margin(2, 2, 2, 2))

  hyp_file <- file.path(OUT_DIR,
                         sprintf("core%s_%s_hypoxia.png", core, pid))
  ggsave(hyp_file, p_hyp, width = w_in, height = h_in, dpi = 450,
         bg = "white", device = ragg::agg_png, limitsize = FALSE)
  message(sprintf("    saved: %s", basename(hyp_file)))

  # --- Cell type panel ---
  palette <- c(ref_palette, "Secretory epithelium" = "#E6A141")

  p_ct <- ggplot(cseg_sf) +
    geom_sf(aes(fill = cell_label), colour = "grey55", linewidth = 0.10) +
    scale_fill_manual(values = palette, name = "Cell type",
                      guide = guide_legend(override.aes = list(linewidth = 0.3),
                                           ncol = 1, keyheight = unit(0.4, "cm"))) +
    coord_sf(xlim = c(BBOX["xmin"], BBOX["xmax"]),
             ylim = c(BBOX["ymin"], BBOX["ymax"]),
             expand = FALSE, default_crs = NULL) +
    labs(title = sprintf("Core %s (patient %s)", core, pid),
         subtitle = "Cell type labels") +
    theme_void(base_size = 6) +
    theme(legend.position   = "right",
          legend.text       = element_text(size = 4),
          legend.title      = element_text(size = 5),
          legend.margin     = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 1),
          plot.title        = element_text(size = 7, colour = "black"),
          plot.subtitle     = element_text(size = 5, colour = "grey40"),
          plot.margin       = margin(2, 2, 2, 2))

  w_ct <- (panel_w + 25) / 25.4
  ct_file <- file.path(OUT_DIR,
                        sprintf("core%s_%s_celltype.png", core, pid))
  ggsave(ct_file, p_ct, width = w_ct, height = h_in, dpi = 450,
         bg = "white", device = ragg::agg_png, limitsize = FALSE)
  message(sprintf("    saved: %s", basename(ct_file)))
}

# --- Summary ----------------------------------------------------------------
message("\n[5/6] Summary of selected cores ...")
for (i in 1:N_TOP) {
  r <- top_cores[i]
  message(sprintf("  Core %4s | patient %5s | grad %.4f | split %2s | hi=%.3f lo=%.3f | sd=%.3f | n=%d",
                  r$core_id, r$patient_id, r$grad_max, r$split_dir,
                  r$mean_hi, r$mean_lo, r$sd_hyp, r$n_cells))
}

message("\n[6/6] Done.")
message("Output: ", OUT_DIR)
message(strrep("=", 70))

log_session()
