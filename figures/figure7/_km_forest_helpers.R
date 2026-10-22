# ============================================================================
# Figure 7 — shared survival-row helpers (KM curves + forest, one SVG per row)
# ----------------------------------------------------------------------------
# Native-ggplot KM panels drawn with geom_step under theme_lab(), so line
# thickness / fonts / axes match the manuscript branding used in Figures 1-3
# (8pt labels, 6pt ticks, 0.4pt axes); survminer/ggsurvplot is not used, as its
# defaults don't honour theme base_size.
#
# A "survival row" = two KM panels on the left + one forest on the right,
# assembled with patchwork and written as ONE combined *.svg (+ *.png). The same
# functions build the TMA row (7a-7d), the TCGA row (7g-7j) and the ICON7 row
# (7k-7n) so the styling is identical. No sub-panel letters (added in Affinity).
#
# Requires (sourced by the caller before this file): config/config.R and
# spatial/00_setup/00_setup.R (theme_lab, svg_inline_fill, ref_palette).
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# House typography for the whole figure-7 survival block.
F7_FA <- 8    # axis titles / annotations
F7_FK <- 6    # tick labels
F7_FL <- 7    # in-panel KM legend
F7_CURVE_LW <- 0.7   # KM step-curve weight (manuscript branding)
F7_AXIS_LW  <- 0.4   # axis + tick weight

# ---------------------------------------------------------------------------
# km_panel(): one Kaplan-Meier panel from a survfit object.
#   fit      survfit(Surv(time, event) ~ grp)  (strata in factor-level order)
#   cols     stroke colours, in that stratum order
#   labs     legend labels, in that stratum order
#   p_text   log-rank annotation string already formatted (e.g. "p < 0.001")
#   subtitle short endpoint tag shown top-left (e.g. "5-yr OS")
#   show_y   FALSE drops the y axis text/title (2nd KM shares the 1st's scale)
# ---------------------------------------------------------------------------
km_panel <- function(fit, cols, labs, p_text, subtitle,
                     xlim = c(0, 60), xbreak = 12, show_y = TRUE) {
  if (is.null(fit$strata)) {
    strata <- factor(rep(1L, length(fit$time)))
  } else {
    strata <- factor(rep(seq_along(fit$strata), fit$strata))
  }
  d  <- data.frame(time = fit$time, surv = fit$surv,
                   ncens = fit$n.censor, g = strata)
  lv <- levels(d$g)
  d0 <- do.call(rbind, lapply(lv, function(g)
          data.frame(time = 0, surv = 1, ncens = 0, g = g)))
  d0$g <- factor(d0$g, levels = lv)
  d  <- rbind(d0, d)
  d  <- d[order(d$g, d$time), ]

  col_v <- setNames(cols, lv)
  lab_v <- setNames(labs, lv)
  cens  <- d[d$ncens > 0, ]

  p <- ggplot(d, aes(time, surv, colour = g)) +
    geom_step(linewidth = F7_CURVE_LW)
  if (nrow(cens)) {
    p <- p + geom_point(data = cens, aes(colour = g),
                        shape = 124, size = 1.3, show.legend = FALSE)
  }
  p <- p +
    scale_colour_manual(values = col_v, labels = lab_v, name = NULL) +
    scale_x_continuous(limits = xlim,
                       breaks = seq(xlim[1], xlim[2], xbreak),
                       expand = expansion(mult = c(0.01, 0.03))) +
    scale_y_continuous(limits = c(0, 1.02),
                       breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                       expand = expansion(mult = c(0.01, 0))) +
    annotate("text", x = xlim[1] + diff(xlim) * 0.03, y = 0.05,
             label = p_text, hjust = 0, size = F7_FA / .pt, colour = "grey30") +
    labs(x = "Months",
         y = if (show_y) "Survival probability" else NULL,
         subtitle = subtitle) +
    theme_lab(base_size = F7_FA) +
    theme(
      legend.position      = c(0.99, 1.00),
      legend.justification = c(1, 1),
      legend.text          = element_text(size = F7_FL),
      legend.key.height    = unit(0.30, "cm"),
      legend.key.width     = unit(0.44, "cm"),
      legend.background     = element_blank(),
      axis.line            = element_line(colour = "black", linewidth = F7_AXIS_LW),
      axis.ticks           = element_line(colour = "black", linewidth = F7_AXIS_LW),
      axis.text            = element_text(size = F7_FK, colour = "black"),
      axis.title           = element_text(size = F7_FA),
      plot.subtitle        = element_text(size = F7_FA, colour = "gray30", hjust = 0),
      plot.margin          = margin(10, 6, 4, 4)
    ) +
    guides(colour = guide_legend(ncol = 1, byrow = TRUE))

  if (!show_y) {
    p <- p + theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.title.y = element_blank())
  }
  p
}

# ---------------------------------------------------------------------------
# forest_theme(): shared forest look (thin axes, no CI text column — HR, p and q text only).
# ---------------------------------------------------------------------------
forest_theme <- function() {
  theme_lab(base_size = F7_FA) +
    theme(
      plot.title         = element_blank(),
      axis.line          = element_line(colour = "black", linewidth = F7_AXIS_LW),
      axis.ticks         = element_line(colour = "black", linewidth = F7_AXIS_LW),
      axis.text.y        = element_text(size = F7_FA, colour = "black"),
      axis.text.x        = element_text(size = F7_FK, colour = "black"),
      axis.title.x       = element_text(size = F7_FA),
      panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.15),
      plot.subtitle      = element_text(size = F7_FA, colour = "gray30", hjust = 0),
      plot.margin        = margin(10, 14, 4, 4)
    )
}

# ---------------------------------------------------------------------------
# Shared forest geometry, so the three survival rows (TMA a-d, TCGA g-j, ICON7 k-n) have the same layout:
#   * forest_ylab_pad(labels): extra left plot.margin (pt) that pads the widest y label to
#     F7_FOREST_YLAB_PT, so every forest panel starts at the same x and has the same width;
#   * forest_geometry(lo, data_hi): log-axis limits that put the data (lo .. data_hi) in the left
#     F7_FOREST_DATA_FRAC of the axis and the HR / p / q text columns at the same fractions in every forest.
# ---------------------------------------------------------------------------
F7_FOREST_YLAB_PT   <- 66                                # widest label: "+ Epithelial fraction" (65 pt at 8 pt)
F7_FOREST_DATA_FRAC <- 0.50
F7_FOREST_COL_FRAC  <- c(HR = 0.54, p = 0.71, q = 0.88)
forest_geometry <- function(lo, data_hi) {
  span <- log10(data_hi / lo) / F7_FOREST_DATA_FRAC
  list(limits = c(lo, lo * 10^span),
       x = setNames(lo * 10^(span * F7_FOREST_COL_FRAC), names(F7_FOREST_COL_FRAC)))
}
forest_ylab_pad <- function(labels, size_pt = F7_FA) {
  w <- systemfonts::string_width(as.character(labels), family = "sans", size = size_pt, res = 72)
  max(0, F7_FOREST_YLAB_PT - max(w))
}

# ---------------------------------------------------------------------------
# save_survival_row(): assemble [km1 | km2 | forest] into ONE svg (+ png),
# Affinity-safe. widths default = KM, KM, wider forest.
# ---------------------------------------------------------------------------
save_survival_row <- function(km1, km2, forest, stem, fig_dir,
                              width_in = 7.0, height_in = 2.5,
                              widths = c(1, 1, 1.75)) {
  row <- km1 + km2 + forest + plot_layout(widths = widths, nrow = 1)
  ggsave(file.path(fig_dir, paste0(stem, ".svg")), row,
         width = width_in, height = height_in, bg = "white")
  ggsave(file.path(fig_dir, paste0(stem, ".png")), row,
         width = width_in, height = height_in, dpi = 450, bg = "white")
  svg_inline_fill(file.path(fig_dir, paste0(stem, ".svg")))
  message("  Saved: ", stem, ".{svg,png}")
  invisible(row)
}
