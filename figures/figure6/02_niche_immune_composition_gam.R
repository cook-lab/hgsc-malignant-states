#!/usr/bin/env Rscript
# ============================================================================
# Figure 6B — Immune composition vs niche hypoxia / glycolysis (GAM)
# ----------------------------------------------------------------------------
# PURPOSE
#   Immune cell composition (% of TME cells) as a function of niche hypoxia
#   and niche glycolysis. All immune cell types overlaid as GAM lines
#   (binomial REML) on a single plot; 2 columns (hypoxia, glycolysis) ×
#   2 rows (WT, TMA).
#
# INPUTS
#   data_root/2026_final_xenium_analysis/output/29_macrophage_niche_survival/
#     per_cell_niche_scores.rds
#   Shared helpers: config/config.R, spatial/00_setup/00_setup.R
#     (ref_palette, theme_lab).
#
# OUTPUTS
#   figures_dir/figure6/niche_immune_composition_gam.{pdf,png,svg}
#
# MANUSCRIPT PANEL(S): Fig 6B
# RUNTIME TIER: moderate (binomial GAM fits over many cells)
# ============================================================================

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

.here     <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))
source(file.path(.here, "..", "..", "spatial", "00_setup", "00_setup.R"))

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(mgcv)
  library(patchwork); library(ragg)
})
set.seed(CFG$seed)

OUT_STEM <- cfg_path("figures_dir", "figure6", "niche_immune_composition_gam")

CT_ORDER <- c("Macrophage", "T cell", "NK cell", "B cell", "Plasma cell")
CT_PAL   <- ref_palette[CT_ORDER]

N_GRID <- 200
K_POOL <- 10

# ── Load data ────────────────────────────────────────────────────────────
message("[load] per_cell_niche_scores.rds")
cache <- readRDS(cfg_path("data_root", "2026_final_xenium_analysis", "output",
                          "29_macrophage_niche_survival",
                          "per_cell_niche_scores.rds"))
wt  <- cache$wt
tma <- cache$tma[!is.na(patient_id) & patient_id != ""]

for (ct in CT_ORDER) {
  flag <- paste0("is_", gsub(" ", "_", ct))
  wt[,  (flag) := as.integer(cell_label == ct)]
  tma[, (flag) := as.integer(cell_label == ct)]
}
message(sprintf("  WT: %s cells, TMA: %s cells",
                format(nrow(wt), big.mark = ","),
                format(nrow(tma), big.mark = ",")))

n_wt  <- length(unique(wt$sample_key))
n_tma <- length(unique(tma$patient_id))

# ── GAM fitting helper ───────────────────────────────────────────────────
fit_immune_gam <- function(dt, x_col, ct) {
  flag <- paste0("is_", gsub(" ", "_", ct))
  d <- data.frame(y = dt[[flag]], x = dt[[x_col]])
  d <- d[is.finite(d$x), ]
  if (nrow(d) < 500 || sum(d$y) < 20) return(NULL)

  x_range <- quantile(d$x, c(0.02, 0.98))
  x_grid  <- seq(x_range[1], x_range[2], length.out = N_GRID)

  fit <- tryCatch(
    gam(y ~ s(x, k = K_POOL), data = d, family = binomial, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  pr <- predict(fit, data.frame(x = x_grid), type = "link", se.fit = TRUE)
  data.table(
    x       = x_grid,
    fitted  = plogis(pr$fit) * 100,
    lower   = plogis(pr$fit - 1.96 * pr$se.fit) * 100,
    upper   = plogis(pr$fit + 1.96 * pr$se.fit) * 100,
    feature = ct
  )
}

# ── Build one plot panel ─────────────────────────────────────────────────
build_niche_panel <- function(dt, x_col, x_label, cohort_label, n_label) {

  pooled_all <- list()
  for (ct in CT_ORDER) {
    res <- fit_immune_gam(dt, x_col, ct)
    if (!is.null(res)) {
      pooled_all[[ct]] <- res
      message(sprintf("  %s / %s: OK", cohort_label, ct))
    }
  }
  pooled <- rbindlist(pooled_all)
  pooled[, feature := factor(feature, levels = CT_ORDER)]

  label_df <- pooled[, .SD[which.max(x)], by = feature]
  label_df[, color := CT_PAL[as.character(feature)]]
  label_df <- label_df[order(fitted)]
  min_gap <- 1.5
  label_df[, fitted_orig := fitted]
  if (nrow(label_df) > 1) {
    for (i in 2:nrow(label_df)) {
      if (label_df$fitted[i] - label_df$fitted[i - 1] < min_gap) {
        label_df$fitted[i] <- label_df$fitted[i - 1] + min_gap
      }
    }
  }

  x_range <- range(pooled$x)
  x_nudge <- diff(x_range) * 0.012

  p <- ggplot() +
    geom_ribbon(data = pooled,
                aes(x = x, ymin = lower, ymax = upper, fill = feature),
                alpha = 0.12, show.legend = FALSE) +
    geom_line(data = pooled,
              aes(x = x, y = fitted, color = feature),
              linewidth = 0.8, show.legend = FALSE) +
    geom_segment(data = label_df[abs(fitted - fitted_orig) > 0.3],
                 aes(x = x, xend = x + x_nudge * 0.6,
                     y = fitted_orig, yend = fitted, color = feature),
                 linewidth = 0.3, show.legend = FALSE) +
    geom_text(data = label_df,
              aes(x = x + x_nudge, y = fitted,
                  label = feature, color = feature),
              hjust = 0, size = 2.0, fontface = "bold",
              show.legend = FALSE) +
    scale_color_manual(values = CT_PAL) +
    scale_fill_manual(values = CT_PAL) +
    labs(x = x_label,
         y = "% of TME cells") +
    coord_cartesian(clip = "off") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.22))) +
    theme_lab() +
    theme(
      plot.margin  = margin(4, 55, 4, 4),
      axis.title.y = element_text(size = 6, angle = 90, margin = margin(r = 4)),
      axis.title.x = element_text(size = 6, margin = margin(t = 4)),
      axis.text    = element_text(size = 5.5)
    )

  p <- p + annotate("text", x = x_range[2], y = -Inf,
                     hjust = 1, vjust = -0.5,
                     label = sprintf("%s (n=%d)", cohort_label, n_label),
                     size = 1.8, colour = "grey40")
  p
}

# ── Build 4 panels ───────────────────────────────────────────────────────
message("\n[panels] fitting GAMs ...")

p_wt_hyp  <- build_niche_panel(wt,  "niche_hypoxia_z",   "Niche hypoxia (z)",    "WT",  n_wt)
p_wt_glyc <- build_niche_panel(wt,  "niche_glycolysis_z", "Niche glycolysis (z)", "WT",  n_wt)
p_tma_hyp <- build_niche_panel(tma, "niche_hypoxia_z",   "Niche hypoxia (z)",    "TMA", n_tma)
p_tma_glyc <- build_niche_panel(tma, "niche_glycolysis_z", "Niche glycolysis (z)", "TMA", n_tma)

# ── Compose 2×2 ──────────────────────────────────────────────────────────
fig <- (p_wt_hyp + p_wt_glyc) /
       (p_tma_hyp + p_tma_glyc)

# ── Save ─────────────────────────────────────────────────────────────────
w <- 7.5
h <- 5.0

ggsave(paste0(OUT_STEM, ".pdf"), fig, width = w, height = h, bg = "white")
ggsave(paste0(OUT_STEM, ".png"), fig, width = w, height = h, dpi = 450,
       bg = "white", device = ragg::agg_png)
ggsave(paste0(OUT_STEM, ".svg"), fig, width = w, height = h, bg = "white")

message("\nDONE")
