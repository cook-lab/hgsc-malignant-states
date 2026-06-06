# ============================================================================
# 10_figures.R
# ----------------------------------------------------------------------------
# PURPOSE: Render macrophage morphology figures (cell area, circularity) by niche.
#
# INPUTS:
#   - output/34_macrophage_morphometrics/ summaries + cross-cohort stats
#
# OUTPUTS:
#   - output/figures/34_macrophage_morphometrics/ figures
#
# MANUSCRIPT PANEL(S): Fig 6I, Fig 6J
# RUNTIME TIER: moderate
#
# Migrated from 2026_final_xenium_analysis/scripts/. Analytical logic preserved;
# paths routed through central config, seed from CFG$seed, epithelial label
# "Transitioning" -> "Intermediate", SecA/SecB from shared/signatures.yml.
# ============================================================================

# --- Config + shared setup (replaces hardcoded /Volumes/CookLab/Sarah paths) ---
here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
source(file.path(here, "..", "..", "config", "config.R"))   # CFG, cfg_obj, cfg_path
source(file.path(here, "..", "00_setup", "00_setup.R"))      # load_sfe, save_sfe, theme_lab, nb_names, palettes
set.seed(CFG$seed)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scico)
  library(ggnewscale)
})

OUT_DIR <- file.path(out_dir, "34_macrophage_morphometrics")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cache <- readRDS(file.path(OUT_DIR, "per_cell_macrophage_morphometrics.rds"))
wt <- cache$wt; tma <- cache$tma
xc  <- fread(file.path(OUT_DIR, "cross_cohort_summary.csv"))

NICHE_COLS <- c(SecA_dominant = "#E6A141", SecB_dominant = "#6B5530")

INCREASE <- "#D65146"; DECREASE <- "#4575B4"

fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) sprintf("p = %.1e", p) else sprintf("p = %.3f", p)
}

# Pick top 3 features by absolute WT Cliff's δ + valid (or top 3 by absolute
# magnitude for headline)
xc_sub <- xc[!is.na(wt_cliffs_delta)]
xc_sub[, abs_cd := abs(wt_cliffs_delta)]
setorder(xc_sub, -abs_cd)
HEADLINE <- head(xc_sub$feature, 3)
HEADLINE_LABELS <- setNames(c("Cell area (µm²)", "Nuclear area (µm²)",
                                "Nuclear : Cytoplasmic ratio",
                                "Cell perimeter (µm)", "Nuc perimeter (µm)",
                                "Cell circularity", "Nuc circularity",
                                "Cell solidity", "Cell eccentricity",
                                "Nuc eccentricity", "NC centroid offset (µm)"),
                              c("cell_area","nuc_area","nc_ratio",
                                "cell_perimeter","nuc_perimeter",
                                "cell_circularity","nuc_circularity",
                                "cell_solidity","cell_eccentricity",
                                "nuc_eccentricity","nc_centroid_offset"))

build_wt_panel <- function(feat) {
  # Tertile design: split per sample on niche_polarization_mean
  wt_sub <- wt[!is.na(niche_polarization_mean) & n_total_epi_neighbors >= 5]
  wt_sub[, pol_tertile := cut(niche_polarization_mean,
                                  breaks = quantile(niche_polarization_mean,
                                                      c(0, 1/3, 2/3, 1),
                                                      na.rm = TRUE),
                                  labels = c("low_secA","mid","high_secB"),
                                  include.lowest = TRUE),
            by = sample_key]
  d <- wt_sub[pol_tertile %in% c("low_secA","high_secB"),
                lapply(.SD, median, na.rm = TRUE),
                by = .(sample_key, pol_tertile),
                .SDcols = feat]
  setnames(d, feat, "y")
  per <- dcast(d, sample_key ~ pol_tertile, value.var = "y")
  per <- per[!is.na(low_secA) & !is.na(high_secB)]
  if (nrow(per) == 0) return(ggplot() + theme_void())

  segs <- per[, .(id = sample_key, x = 1, xend = 2,
                    y = low_secA, yend = high_secB,
                    col = ifelse(high_secB - low_secA > 0,
                                  "up", "down"))]
  d[, x := match(pol_tertile, c("low_secA","high_secB"))]

  p_two <- xc[feature == feat, wt_p_bonf]
  p_text <- sprintf("WT n=%d • Bonferroni p = %s", nrow(per), fmt_p(p_two))

  TERTILE_COLS <- c(low_secA = "#E6A141", high_secB = "#6B5530")
  ggplot() +
    geom_segment(data = segs,
                 aes(x = x, xend = xend, y = y, yend = yend, colour = col),
                 linewidth = 0.7, alpha = 0.85, show.legend = FALSE) +
    scale_colour_manual(values = c(up = INCREASE, down = DECREASE)) +
    ggnewscale::new_scale_colour() +
    geom_point(data = d, aes(x = x, y = y, colour = pol_tertile),
                size = 3, alpha = 0.95) +
    scale_colour_manual(values = TERTILE_COLS, guide = "none") +
    scale_x_continuous(breaks = 1:2,
                       labels = c("Low niche-pol\n(SecA-rich)",
                                    "High niche-pol\n(SecB-rich)")) +
    labs(x = NULL, y = HEADLINE_LABELS[[feat]],
          title    = sprintf("WT — %s", HEADLINE_LABELS[[feat]]),
          subtitle = p_text) +
    theme_lab(base_size = 9.5) +
    theme(panel.grid.major.x = element_blank(),
           plot.title = element_text(size = rel(1.05)),
           plot.subtitle = element_text(size = rel(0.85), colour = "grey35"))
}

build_tma_panel <- function(feat) {
  tma_sub <- tma[!is.na(niche_polarization_mean) & n_total_epi_neighbors >= 5]
  patient_n <- tma_sub[, .N, by = group_id]
  big_p <- patient_n[N >= 30, group_id]
  tma_sub <- tma_sub[group_id %in% big_p]
  tma_sub[, pol_rank := frank(niche_polarization_mean, ties.method = "average"),
            by = group_id]
  tma_sub[, pol_pct := pol_rank / max(pol_rank), by = group_id]
  tma_sub[, pol_tertile := fcase(
    pol_pct <= 1/3, "low_secA",
    pol_pct >= 2/3, "high_secB",
    default = "mid")]

  d <- tma_sub[pol_tertile %in% c("low_secA","high_secB"),
                  lapply(.SD, median, na.rm = TRUE),
                  by = .(group_id, pol_tertile),
                  .SDcols = feat]
  setnames(d, feat, "y")
  per <- dcast(d, group_id ~ pol_tertile, value.var = "y")
  per <- per[!is.na(low_secA) & !is.na(high_secB)]
  if (nrow(per) == 0) return(ggplot() + theme_void())

  segs <- per[, .(id = group_id, x = 1, xend = 2,
                    y = low_secA, yend = high_secB,
                    col = ifelse(high_secB - low_secA > 0,
                                  "up", "down"))]
  d_sub <- d[group_id %in% per$group_id]
  d_sub[, x := match(pol_tertile, c("low_secA","high_secB"))]

  p_bonf <- xc[feature == feat, tma_p_bonf]
  pct <- xc[feature == feat, tma_pct_concordant]
  p_text <- sprintf("TMA n=%d • %.0f%% concordant • Bonferroni p = %s",
                     nrow(per), pct, fmt_p(p_bonf))

  TERTILE_COLS <- c(low_secA = "#E6A141", high_secB = "#6B5530")
  ggplot() +
    geom_violin(data = d_sub, aes(x = x, y = y, group = x,
                                       fill = pol_tertile),
                colour = NA, alpha = 0.18, width = 0.6,
                trim = FALSE, scale = "width") +
    scale_fill_manual(values = TERTILE_COLS, guide = "none") +
    geom_segment(data = segs,
                 aes(x = x, xend = xend, y = y, yend = yend, colour = col),
                 linewidth = 0.25, alpha = 0.5, show.legend = FALSE) +
    scale_colour_manual(values = c(up = INCREASE, down = DECREASE)) +
    ggnewscale::new_scale_colour() +
    geom_boxplot(data = d_sub,
                 aes(x = x, y = y, group = x, colour = pol_tertile),
                 width = 0.18, fill = NA, outlier.shape = NA,
                 linewidth = 0.45, show.legend = FALSE) +
    scale_colour_manual(values = TERTILE_COLS, guide = "none") +
    scale_x_continuous(breaks = 1:2,
                       labels = c("Low niche-pol\n(SecA-rich)",
                                    "High niche-pol\n(SecB-rich)")) +
    labs(x = NULL, y = HEADLINE_LABELS[[feat]],
          title    = sprintf("TMA — %s", HEADLINE_LABELS[[feat]]),
          subtitle = p_text) +
    theme_lab(base_size = 9.5) +
    theme(panel.grid.major.x = element_blank(),
           plot.title = element_text(size = rel(1.05)),
           plot.subtitle = element_text(size = rel(0.85), colour = "grey35"))
}

panels <- list()
for (f in HEADLINE) {
  panels[[paste0(f, "_wt")]]  <- build_wt_panel(f)
  panels[[paste0(f, "_tma")]] <- build_tma_panel(f)
}

fig_A <- (panels[[paste0(HEADLINE[1], "_wt")]]  | panels[[paste0(HEADLINE[1], "_tma")]]) /
          (panels[[paste0(HEADLINE[2], "_wt")]]  | panels[[paste0(HEADLINE[2], "_tma")]]) /
          (panels[[paste0(HEADLINE[3], "_wt")]]  | panels[[paste0(HEADLINE[3], "_tma")]]) +
  plot_annotation(
    title    = "Phase 34 — Macrophage morphology by epithelial-niche dominance",
    subtitle = "Top 3 features (by |Cliff's δ| in WT) for SecA-dominant vs SecB-dominant niche macrophages.\nWT (n=8 samples paired by sample) + TMA (eligible patients paired by patient).",
    theme = theme(plot.title    = element_text(size = 12, face = "bold"),
                   plot.subtitle = element_text(size = 9, colour = "grey30"))
  )

ggsave(file.path(FIG_DIR, "fig34A_headline.png"),
       fig_A, width = 12, height = 11, dpi = 400, bg = "white",
       device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "fig34A_headline.svg"),
       fig_A, width = 12, height = 11, bg = "white")

# ---------------------------------------------------------------------------
# Heatmap
# ---------------------------------------------------------------------------
hm <- xc[, .(feature,
               WT  = wt_cliffs_delta,
               TMA = tma_cliffs_delta,
               WT_sig  = wt_p_bonf < 0.05,
               TMA_sig = tma_p_bonf < 0.05)]
hm_long <- melt(hm, id.vars = "feature",
                  measure.vars = patterns(value = "^(WT|TMA)$",
                                            sig   = "_sig$"),
                  variable.factor = FALSE)
hm_long[, cohort := ifelse(variable == "1", "WT", "TMA")]
hm_long[, cohort := factor(cohort, levels = c("WT", "TMA"))]
ord <- xc[order(abs(wt_cliffs_delta), decreasing = TRUE), feature]
hm_long[, feature := factor(feature, levels = rev(ord))]
hm_long[, label := ifelse(sig, "*", "")]

heatmap_p <- ggplot(hm_long, aes(x = cohort, y = feature, fill = value)) +
  geom_tile(colour = "grey95", linewidth = 0.5) +
  geom_text(aes(label = label), size = 5, fontface = "bold",
             colour = "grey15") +
  scale_fill_scico(palette = "vik", midpoint = 0,
                    name = "Cliff's δ\n(neg = SecB-niche < SecA-niche)",
                    limits = c(-0.2, 0.2), oob = scales::squish) +
  labs(x = NULL, y = NULL,
       title = "Phase 34 — Macrophage morphometric effect sizes by niche dominance",
       subtitle = "SecA-dominant niche macrophages vs SecB-dominant niche macrophages.\nCliff's δ signed; * = Bonferroni-corrected paired Wilcoxon p < 0.05.") +
  theme_lab(base_size = 10) +
  theme(axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank(),
        plot.title = element_text(size = rel(1.15), face = "bold"),
        plot.subtitle = element_text(size = rel(0.85), colour = "grey35"))

ggsave(file.path(FIG_DIR, "fig34B_summary_heatmap.png"),
       heatmap_p, width = 7, height = 7, dpi = 400, bg = "white",
       device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "fig34B_summary_heatmap.svg"),
       heatmap_p, width = 7, height = 7, bg = "white")

message("\nSaved figures.")
