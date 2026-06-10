# ============================================================================
# 05_morphometric_figures.R
# ----------------------------------------------------------------------------
# PURPOSE: Render paired-slope morphometric figures (nuclear area, N:C ratio) across epitypes.
#
# INPUTS:
#   - output/33_morphometrics/ per-sample/per-patient summaries + cross-cohort stats
#
# OUTPUTS:
#   - output/figures/33_morphometrics/ paired-slope figures
#
# MANUSCRIPT PANEL(S): Fig 5D, Fig 5E
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

OUT_DIR <- file.path(out_dir, "33_morphometrics")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Load per-cell cache + summaries
cache <- readRDS(file.path(OUT_DIR, "per_cell_morphometrics.rds"))
wt <- cache$wt
tma <- cache$tma
xc  <- fread(file.path(OUT_DIR, "cross_cohort_summary.csv"))

# Idempotent rename: the per-cell morphometrics cache may have been written from
# SFE/06f sources still carrying the legacy epithelial label "Transitioning
# epithelium"; standardize to "Intermediate epithelium" before the EPI_LBLS
# match/colour so the Intermediate epitype is not silently dropped/greyed.
# Harmless no-op if the cache was already written with the standardized label.
if ("cell_label" %in% names(wt))
  wt[cell_label == "Transitioning epithelium", cell_label := "Intermediate epithelium"]
if ("cell_label" %in% names(tma))
  tma[cell_label == "Transitioning epithelium", cell_label := "Intermediate epithelium"]

EPI_LBLS <- c("SecA epithelium", "Intermediate epithelium", "SecB epithelium")
EPI_COLS <- ref_palette[EPI_LBLS]

# ---------------------------------------------------------------------------
# Headline panels: paired WT + violin/box TMA for nc_ratio, nuc_area,
# cell_area (the three features with strongest gradient signal)
# ---------------------------------------------------------------------------
HEADLINE_FEATURES <- c(
  nc_ratio  = "Nuclear : Cytoplasmic ratio",
  nuc_area  = "Nuclear area (µm²)",
  cell_area = "Cell area (µm²)"
)

INCREASE <- "#D65146"; DECREASE <- "#4575B4"

fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) sprintf("p = %.1e", p) else sprintf("p = %.3f", p)
}

build_wt_panel <- function(feat) {
  feat_label <- HEADLINE_FEATURES[[feat]]
  d <- wt[cell_label %in% EPI_LBLS,
            lapply(.SD, median, na.rm = TRUE),
            by = .(sample_key, cell_label),
            .SDcols = feat]
  setnames(d, feat, "y")
  d <- d[!is.na(y)]
  d[, x := match(cell_label, EPI_LBLS)]
  per <- dcast(d, sample_key ~ cell_label, value.var = "y")
  ok <- !is.na(per[["SecA epithelium"]]) &
        !is.na(per[["Intermediate epithelium"]]) &
        !is.na(per[["SecB epithelium"]])
  per <- per[ok]
  delta <- per[["SecB epithelium"]] - per[["SecA epithelium"]]
  dir_col <- ifelse(delta > 0, INCREASE, DECREASE)
  seg <- per[, .(id = sample_key,
                   x = 1, xend = 3,
                   y = get("SecA epithelium"),
                   yend = get("SecB epithelium"),
                   col = ifelse(get("SecB epithelium") - get("SecA epithelium") > 0,
                                  "up", "down"))]
  # Three-segment line per sample (1→2→3)
  segs <- rbindlist(list(
    per[, .(id = sample_key, x = 1, xend = 2,
              y = get("SecA epithelium"),
              yend = get("Intermediate epithelium"))],
    per[, .(id = sample_key, x = 2, xend = 3,
              y = get("Intermediate epithelium"),
              yend = get("SecB epithelium"))]
  ))
  segs[, col := ifelse(yend - y > 0, "up", "down")]

  # Pull p-value for SecA_to_SecB
  p_secA_secB <- xc[feature == feat & comparison == "SecA_to_SecB",
                     wt_p_bonf_two]
  p_text <- sprintf("WT n=8 • Bonferroni p = %s", fmt_p(p_secA_secB))

  ggplot() +
    geom_segment(data = segs,
                 aes(x = x, xend = xend, y = y, yend = yend, colour = col),
                 linewidth = 0.7, alpha = 0.78,
                 show.legend = FALSE) +
    scale_colour_manual(values = c(up = INCREASE, down = DECREASE)) +
    ggnewscale::new_scale_colour() +
    geom_point(data = d, aes(x = x, y = y, colour = cell_label),
                size = 2.8, alpha = 0.95) +
    scale_colour_manual(values = EPI_COLS, guide = "none") +
    scale_x_continuous(breaks = 1:3,
                       labels = c("SecA", "Intermediate", "SecB")) +
    labs(x = NULL, y = feat_label,
          title    = sprintf("WT — %s", feat_label),
          subtitle = p_text) +
    theme_lab(base_size = 9.5) +
    theme(panel.grid.major.x = element_blank(),
           plot.title = element_text(size = rel(1.05)),
           plot.subtitle = element_text(size = rel(0.85), colour = "grey35"))
}

build_tma_panel <- function(feat) {
  feat_label <- HEADLINE_FEATURES[[feat]]
  # Per-patient median per epitype
  d <- tma[cell_label %in% EPI_LBLS,
              lapply(.SD, median, na.rm = TRUE),
              by = .(group_id, cell_label),
              .SDcols = feat]
  setnames(d, feat, "y")
  per <- dcast(d, group_id ~ cell_label, value.var = "y")
  ok <- !is.na(per[["SecA epithelium"]]) &
        !is.na(per[["Intermediate epithelium"]]) &
        !is.na(per[["SecB epithelium"]])
  per <- per[ok]
  d_sub <- d[group_id %in% per$group_id]
  d_sub[, x := match(cell_label, EPI_LBLS)]

  # Build segments
  segs <- rbindlist(list(
    per[, .(id = group_id, x = 1, xend = 2,
              y = get("SecA epithelium"),
              yend = get("Intermediate epithelium"))],
    per[, .(id = group_id, x = 2, xend = 3,
              y = get("Intermediate epithelium"),
              yend = get("SecB epithelium"))]
  ))
  segs[, col := ifelse(yend - y > 0, "up", "down")]

  p_secA_secB <- xc[feature == feat & comparison == "SecA_to_SecB",
                      tma_p_bonf_two]
  pct_concordant <- xc[feature == feat & comparison == "SecA_to_SecB",
                         tma_pct_concordant]
  p_text <- sprintf("TMA n = %d patients • %.0f%% concordant • Bonferroni p = %s",
                     nrow(per), pct_concordant, fmt_p(p_secA_secB))

  ggplot() +
    geom_violin(data = d_sub, aes(x = x, y = y, group = x,
                                       fill = cell_label),
                colour = NA, alpha = 0.18, width = 0.7,
                trim = FALSE, scale = "width") +
    scale_fill_manual(values = EPI_COLS, guide = "none") +
    geom_segment(data = segs,
                 aes(x = x, xend = xend, y = y, yend = yend, colour = col),
                 linewidth = 0.25, alpha = 0.4,
                 show.legend = FALSE) +
    scale_colour_manual(values = c(up = INCREASE, down = DECREASE)) +
    ggnewscale::new_scale_colour() +
    geom_boxplot(data = d_sub,
                 aes(x = x, y = y, group = x, colour = cell_label),
                 width = 0.18, fill = NA, outlier.shape = NA,
                 linewidth = 0.45, show.legend = FALSE) +
    scale_colour_manual(values = EPI_COLS, guide = "none") +
    scale_x_continuous(breaks = 1:3,
                       labels = c("SecA", "Intermediate", "SecB")) +
    labs(x = NULL, y = feat_label,
          title    = sprintf("TMA — %s", feat_label),
          subtitle = p_text) +
    theme_lab(base_size = 9.5) +
    theme(panel.grid.major.x = element_blank(),
           plot.title = element_text(size = rel(1.05)),
           plot.subtitle = element_text(size = rel(0.85), colour = "grey35"))
}

panels <- list()
for (f in names(HEADLINE_FEATURES)) {
  panels[[paste0(f, "_wt")]]  <- build_wt_panel(f)
  panels[[paste0(f, "_tma")]] <- build_tma_panel(f)
}

fig_A <- (panels$nc_ratio_wt | panels$nc_ratio_tma) /
          (panels$nuc_area_wt | panels$nuc_area_tma) /
          (panels$cell_area_wt | panels$cell_area_tma) +
  plot_annotation(
    title    = "Phase 33 — Morphometric reversal: SecB cells are SMALLER (cell + nucleus) and have LOWER N/C ratio than SecA",
    subtitle = "WT (n=8 samples, paired by sample) and TMA (n=63 patients, paired by patient) both confirm the gradient is opposite to the histopath textbook prediction.\nProliferating SecA cells have LARGER nuclei + HIGHER N/C ratio than quiescent SecB cells. Standard H&E grading may systematically under-score SecB-rich regions.",
    theme = theme(plot.title    = element_text(size = 12, face = "bold"),
                   plot.subtitle = element_text(size = 9, colour = "grey30"))
  )

ggsave(file.path(FIG_DIR, "fig33A_headline.png"),
       fig_A, width = 12, height = 11, dpi = 400, bg = "white",
       device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "fig33A_headline.svg"),
       fig_A, width = 12, height = 11, bg = "white")

# ---------------------------------------------------------------------------
# Summary heatmap: 9 features × 3 pairwise × 2 cohorts (= 6 cols)
# Cell value: signed Cliff's δ (negative = b < a, positive = b > a)
# Asterisk: Bonferroni p < 0.05 (two-sided)
# ---------------------------------------------------------------------------
heat_dt <- xc[, .(feature, comparison,
                    wt_cd = wt_cliffs_delta,
                    tma_cd = tma_cliffs_delta,
                    wt_sig = wt_p_bonf_two < 0.05,
                    tma_sig = tma_p_bonf_two < 0.05)]
hm <- melt(heat_dt, id.vars = c("feature","comparison"),
            measure.vars = patterns(value = "_cd$", sig = "_sig$"),
            variable.factor = FALSE)
hm[, cohort := ifelse(variable == "1", "WT", "TMA")]
hm[, cohort := factor(cohort, levels = c("WT", "TMA"))]
hm[, comparison := factor(comparison,
                            levels = c("SecA_to_Trans", "Trans_to_SecB",
                                        "SecA_to_SecB"),
                            labels = c("SecA→Trans", "Trans→SecB",
                                        "SecA→SecB"))]
# Order features by absolute value of WT SecA→SecB Cliff's δ
ord <- xc[comparison == "SecA_to_SecB"][order(abs(wt_cliffs_delta), decreasing = TRUE), feature]
hm[, feature := factor(feature, levels = rev(ord))]
hm[, label := ifelse(sig, "*", "")]
hm[, x_grp := factor(paste(cohort, comparison, sep = " | "),
                       levels = c("WT | SecA→Trans", "WT | Trans→SecB",
                                    "WT | SecA→SecB",
                                    "TMA | SecA→Trans", "TMA | Trans→SecB",
                                    "TMA | SecA→SecB"))]

heatmap_p <- ggplot(hm, aes(x = x_grp, y = feature, fill = value)) +
  geom_tile(colour = "grey95", linewidth = 0.5) +
  geom_text(aes(label = label), size = 5, colour = "grey15",
             fontface = "bold") +
  scale_fill_scico(palette = "vik", midpoint = 0, name = "Cliff's δ\n(neg = SecB < SecA)",
                    limits = c(-0.2, 0.2), oob = scales::squish) +
  labs(x = NULL, y = NULL,
       title = "Phase 33 — Morphometric effect sizes across cohorts",
       subtitle = "Cliff's δ (signed); negative = feature value lower in B than A. * = Bonferroni-corrected paired Wilcoxon p < 0.05.") +
  theme_lab(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = rel(0.9)),
        axis.text.y = element_text(face = "bold", size = rel(0.95)),
        panel.grid = element_blank(),
        plot.title    = element_text(size = rel(1.15), face = "bold"),
        plot.subtitle = element_text(size = rel(0.85), colour = "grey35"))

ggsave(file.path(FIG_DIR, "fig33B_summary_heatmap.png"),
       heatmap_p, width = 10, height = 7, dpi = 400, bg = "white",
       device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "fig33B_summary_heatmap.svg"),
       heatmap_p, width = 10, height = 7, bg = "white")

message("\nSaved:")
message("  ", file.path(FIG_DIR, "fig33A_headline.png"))
message("  ", file.path(FIG_DIR, "fig33B_summary_heatmap.png"))
