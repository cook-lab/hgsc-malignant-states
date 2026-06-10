#!/usr/bin/env Rscript
# ============================================================================
# Figure 6I,6J — Macrophage morphology / phenotype in hypoxic niches
# ----------------------------------------------------------------------------
# PURPOSE
#   Four paired-slope SVGs supporting macrophage phenotype across niches:
#     Apoptosis (lower in hypoxic niche, stress decile 1 vs 10)
#     Proliferation (higher in hypoxic niche, stress decile 1 vs 10)
#     Cell area (increased in SecB-rich neighborhoods, polarization tertile) — Fig 6I
#     Circularity (increased in SecB-rich neighborhoods, polarization tertile) — Fig 6J
#   Per-sample paired medians + cell-level violins; Wilcoxon signed-rank.
#
# INPUTS
#   data_root/2026_final_xenium_analysis/output/29_macrophage_niche_survival/
#     functional_survival_scored_cells.rds, functional_survival_per_patient.csv
#   data_root/2026_final_xenium_analysis/output/34_macrophage_morphometrics/
#     per_cell_macrophage_morphometrics.rds
#   Shared helpers: config/config.R, spatial/00_setup/00_setup.R (theme_lab).
#
# OUTPUTS
#   figures_dir/figure6/fig_macrophage_apoptosis_wt.{pdf,png,svg}
#   figures_dir/figure6/fig_macrophage_proliferation_wt.{pdf,png,svg}
#   figures_dir/figure6/fig_macrophage_cell_area_wt.{pdf,png,svg}    (Fig 6I)
#   figures_dir/figure6/fig_macrophage_circularity_wt.{pdf,png,svg}  (Fig 6J)
#
# MANUSCRIPT PANEL(S): Fig 6I (cell area), Fig 6J (circularity); apoptosis &
#   proliferation are supporting panels.
# RUNTIME TIER: moderate (reads caches)
#
# NOTE: Fig 6I/6J header p-values differ between header text and published PDF
#   (6I ~p=0.014 vs PDF 0.008; 6J 6/8 p=0.042 vs PDF 8/8 p=0.008). Code is
#   migrated faithfully; reproduced p should be
#   verified against the live-filtered data, not silently changed.
# ============================================================================

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

.here     <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))
source(file.path(.here, "..", "..", "spatial", "00_setup", "00_setup.R"))

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
  library(ragg); library(svglite)
})

set.seed(CFG$seed)

FIG_DIR <- cfg_path("figures_dir", "figure6")
XEN_OUT <- function(...) cfg_path("data_root", "2026_final_xenium_analysis", "output", ...)

COL_MAC_LO <- "#B8D8B8"
COL_MAC_HI <- "#4A8F4A"

MIN_CELLS <- 15

# =========================================================================
# Apoptosis & Proliferation (stress decile)
# =========================================================================
message("\n[load] functional_survival_scored_cells.rds ...")
scored <- readRDS(XEN_OUT("29_macrophage_niche_survival",
                          "functional_survival_scored_cells.rds"))
wt_cells <- scored$wt
setDT(wt_cells)

per <- fread(XEN_OUT("29_macrophage_niche_survival",
                     "functional_survival_per_patient.csv"))

build_slope_violin <- function(cells_dt, score_col, per_dt,
                                y_lab, fill_lo, fill_hi) {

  plot_cells <- cells_dt[stress_decile %in% c(1, 10)]
  plot_cells[, x_num := fifelse(stress_decile == 1, 1, 2)]
  plot_cells[, score := get(score_col)]

  paired <- rbind(
    per_dt[, .(sample = sample_key, x = 1, y = med_bot)],
    per_dt[, .(sample = sample_key, x = 2, y = med_top)]
  )
  paired[, direction := fifelse(
    y[x == 2] > y[x == 1], "up", "down"), by = sample]

  w <- wilcox.test(per_dt$delta, mu = 0)
  if (median(per_dt$delta, na.rm = TRUE) < 0) {
    n_conc <- sum(per_dt$delta < 0)
    conc_dir <- "decrease"
    line_col <- fill_lo
  } else {
    n_conc <- sum(per_dt$delta > 0)
    conc_dir <- "increase"
    line_col <- fill_hi
  }

  p <- ggplot() +
    geom_violin(data = plot_cells,
                aes(x = x_num, y = score,
                    fill = factor(x_num), group = x_num),
                colour = NA, alpha = 0.15,
                scale = "width", width = 0.6) +
    geom_line(data = paired,
              aes(x = x, y = y, group = sample, colour = direction),
              linewidth = 0.6, alpha = 0.8) +
    geom_point(data = paired,
               aes(x = x, y = y, fill = factor(x)),
               shape = 21, size = 2.5, colour = "grey30", stroke = 0.3) +
    scale_fill_manual(values = c("1" = fill_lo, "2" = fill_hi),
                      guide = "none") +
    scale_colour_manual(
      values = c("up"   = if (conc_dir == "increase") fill_hi else "grey70",
                 "down" = if (conc_dir == "decrease") fill_hi else "grey70"),
      guide = "none") +
    scale_x_continuous(breaks = c(1, 2),
                       labels = c("Low stress\n(bottom decile)",
                                  "High stress\n(top decile)"),
                       expand = expansion(mult = 0.3)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    labs(x = NULL, y = y_lab) +
    annotate("text", x = 1.5, y = Inf, vjust = 1.3,
             label = sprintf("n=%d samples\n%d/%d %s\np = %.3f",
                             nrow(per_dt), n_conc, nrow(per_dt),
                             conc_dir, w$p.value),
             size = 1.8, colour = "grey30") +
    theme_lab() +
    theme(
      axis.title      = element_text(size = 6),
      axis.text.x     = element_text(size = 5.5),
      axis.text.y     = element_text(size = 5.5),
      plot.margin     = margin(4, 4, 4, 4)
    )

  p
}

# ── Macrophage apoptosis ─────────────────────────────────────────────────
message("[panel] Macrophage apoptosis ...")

mac_cells <- wt_cells[cell_label == "Macrophage"]

per_mac_apop <- per[cell_type == "Macrophage" &
                    score == "apoptosis_score" &
                    cohort == "WT" &
                    n_top >= MIN_CELLS & n_bot >= MIN_CELLS]

pK <- build_slope_violin(
  cells_dt  = mac_cells,
  score_col = "apoptosis_score",
  per_dt    = per_mac_apop,
  y_lab     = "Macrophage apoptosis score\n(mean logcounts)",
  fill_lo   = COL_MAC_LO,
  fill_hi   = COL_MAC_HI
)

stem_K <- file.path(FIG_DIR, "fig_macrophage_apoptosis_wt")
ggsave(paste0(stem_K, ".pdf"), pK, width = 2.2, height = 2.5, bg = "white")
ggsave(paste0(stem_K, ".png"), pK, width = 2.2, height = 2.5, dpi = 450,
       bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_K, ".svg"), pK, width = 2.2, height = 2.5, bg = "white")
message("  Saved: ", stem_K)

# ── Macrophage proliferation ─────────────────────────────────────────────
message("[panel] Macrophage proliferation ...")

per_mac_prolif <- per[cell_type == "Macrophage" &
                      score == "proliferation_score" &
                      cohort == "WT" &
                      n_top >= MIN_CELLS & n_bot >= MIN_CELLS]

pL <- build_slope_violin(
  cells_dt  = mac_cells,
  score_col = "proliferation_score",
  per_dt    = per_mac_prolif,
  y_lab     = "Macrophage proliferation score\n(mean logcounts)",
  fill_lo   = COL_MAC_LO,
  fill_hi   = COL_MAC_HI
)

stem_L <- file.path(FIG_DIR, "fig_macrophage_proliferation_wt")
ggsave(paste0(stem_L, ".pdf"), pL, width = 2.2, height = 2.5, bg = "white")
ggsave(paste0(stem_L, ".png"), pL, width = 2.2, height = 2.5, dpi = 450,
       bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_L, ".svg"), pL, width = 2.2, height = 2.5, bg = "white")
message("  Saved: ", stem_L)

# =========================================================================
# Cell area & Circularity (polarization tertile) — Fig 6I / 6J
# =========================================================================
message("\n[load] per_cell_macrophage_morphometrics.rds ...")
morph <- readRDS(XEN_OUT("34_macrophage_morphometrics",
                         "per_cell_macrophage_morphometrics.rds"))
wt_morph <- morph$wt
setDT(wt_morph)

wt_morph <- wt_morph[n_total_epi_neighbors >= 5 & !is.na(niche_polarization_mean)]
wt_morph[, pol_tertile := cut(niche_polarization_mean,
                               breaks = quantile(niche_polarization_mean,
                                                  probs = c(0, 1/3, 2/3, 1),
                                                  na.rm = TRUE),
                               labels = c("low_secA", "mid", "high_secB"),
                               include.lowest = TRUE),
          by = sample_key]

wt_morph_plot <- wt_morph[pol_tertile %in% c("low_secA", "high_secB")]
wt_morph_plot[, pol_bin := factor(
  fifelse(pol_tertile == "low_secA",
          "SecA-rich\n(low polarization)",
          "SecB-rich\n(high polarization)"),
  levels = c("SecA-rich\n(low polarization)",
             "SecB-rich\n(high polarization)")
)]

get_morph_paired <- function(morph_dt, metric_col) {
  morph_dt[, metric := get(metric_col)]
  sample_meds <- morph_dt[, .(med = median(metric, na.rm = TRUE)),
                           by = .(sample_key, pol_tertile)]
  sample_meds_wide <- dcast(sample_meds, sample_key ~ pol_tertile,
                             value.var = "med")
  sample_meds_wide <- sample_meds_wide[!is.na(low_secA) & !is.na(high_secB)]
  sample_meds_wide[, delta := high_secB - low_secA]

  paired <- rbind(
    sample_meds_wide[, .(sample = sample_key, x = 1, y = low_secA)],
    sample_meds_wide[, .(sample = sample_key, x = 2, y = high_secB)]
  )

  w <- wilcox.test(sample_meds_wide$delta, mu = 0)
  pct_pos <- round(100 * mean(sample_meds_wide$delta > 0), 1)

  list(paired = paired, w = w, pct_pos = pct_pos,
       n_samples = nrow(sample_meds_wide))
}

build_slope_plot <- function(morph_dt, metric_col, y_lab,
                              col_lo, col_hi) {

  morph_dt[, metric := get(metric_col)]

  sample_meds <- morph_dt[, .(med = median(metric, na.rm = TRUE)),
                           by = .(sample_key, pol_tertile)]
  sample_meds_wide <- dcast(sample_meds, sample_key ~ pol_tertile,
                             value.var = "med")
  sample_meds_wide <- sample_meds_wide[!is.na(low_secA) & !is.na(high_secB)]
  sample_meds_wide[, delta := high_secB - low_secA]
  sample_meds_wide[, direction := fifelse(delta > 0, "up", "down")]

  paired <- rbind(
    sample_meds_wide[, .(sample = sample_key, x = 1, y = low_secA, direction)],
    sample_meds_wide[, .(sample = sample_key, x = 2, y = high_secB, direction)]
  )

  w <- wilcox.test(sample_meds_wide$delta, mu = 0)
  pct_pos <- round(100 * mean(sample_meds_wide$delta > 0), 1)
  n_up   <- sum(sample_meds_wide$delta > 0)
  n_down <- sum(sample_meds_wide$delta <= 0)

  morph_dt[, x_num := fifelse(pol_tertile == "low_secA", 1, 2)]

  p <- ggplot() +
    geom_violin(data = morph_dt,
                aes(x = x_num, y = metric,
                    fill = factor(x_num), group = x_num),
                colour = NA, alpha = 0.15,
                scale = "width", width = 0.6) +
    geom_line(data = paired,
              aes(x = x, y = y, group = sample, colour = direction),
              linewidth = 0.6, alpha = 0.8) +
    geom_point(data = paired,
               aes(x = x, y = y, fill = factor(x)),
               shape = 21, size = 2.5, colour = "grey30", stroke = 0.3) +
    scale_colour_manual(values = c("up" = col_hi, "down" = "grey70"),
                        guide = "none") +
    scale_fill_manual(values = c("1" = col_lo, "2" = col_hi),
                      guide = "none") +
    scale_x_continuous(breaks = c(1, 2),
                       labels = c("SecA-rich\n(low pol.)",
                                  "SecB-rich\n(high pol.)"),
                       expand = expansion(mult = 0.3)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    labs(x = NULL, y = y_lab) +
    annotate("text", x = 1.5, y = Inf, vjust = 1.3,
             label = sprintf("n=%d samples\n%d/%d increase\np = %.3f",
                             nrow(sample_meds_wide),
                             n_up, nrow(sample_meds_wide),
                             w$p.value),
             size = 1.8, colour = "grey30") +
    theme_lab() +
    theme(
      axis.title      = element_text(size = 6),
      axis.text.x     = element_text(size = 5.5),
      axis.text.y     = element_text(size = 5.5),
      plot.margin     = margin(4, 4, 4, 4)
    )

  p
}

# ── Fig 6I: Cell area (paired slope) ─────────────────────────────────────
message("[panel I] Macrophage cell area (paired slope) ...")

pM <- build_slope_plot(
  morph_dt   = wt_morph_plot,
  metric_col = "cell_area",
  y_lab      = expression(paste("Median cell area (", mu, "m"^2, ")")),
  col_lo     = COL_MAC_LO,
  col_hi     = COL_MAC_HI
) + coord_cartesian(ylim = c(50, 100))

stem_M <- file.path(FIG_DIR, "fig_macrophage_cell_area_wt")
ggsave(paste0(stem_M, ".pdf"), pM, width = 2.2, height = 2.5, bg = "white")
ggsave(paste0(stem_M, ".png"), pM, width = 2.2, height = 2.5, dpi = 450,
       bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_M, ".svg"), pM, width = 2.2, height = 2.5, bg = "white")
message("  Saved: ", stem_M)

# ── Fig 6J: Circularity (paired slope) ───────────────────────────────────
message("[panel J] Macrophage circularity (paired slope) ...")

pN <- build_slope_plot(
  morph_dt   = wt_morph_plot,
  metric_col = "cell_circularity",
  y_lab      = "Median circularity",
  col_lo     = COL_MAC_LO,
  col_hi     = COL_MAC_HI
) + coord_cartesian(ylim = c(0.7, 0.9))

stem_N <- file.path(FIG_DIR, "fig_macrophage_circularity_wt")
ggsave(paste0(stem_N, ".pdf"), pN, width = 2.2, height = 2.5, bg = "white")
ggsave(paste0(stem_N, ".png"), pN, width = 2.2, height = 2.5, dpi = 450,
       bg = "white", device = ragg::agg_png)
ggsave(paste0(stem_N, ".svg"), pN, width = 2.2, height = 2.5, bg = "white")
message("  Saved: ", stem_N)

message("\n\nAll saved.")
message("DONE")
