# ============================================================================
# 36_helpers.R
# ----------------------------------------------------------------------------
# Phase 36 figure-series helpers. Standardised polarization-axis GAM overlay
# style used by every figure in the series:
#   * x-axis = epithelial polarization_UCell
#   * y-axis = z-scored pooled GAM fit (per-curve normalised)
#   * pooled GAM (mgcv, k=10, REML) + 95% CI ribbon (alpha 0.15)
#   * horizontal reference at z = 0
#   * end labels via ggrepel
#   * theme_lab base
# ============================================================================

suppressPackageStartupMessages({
  library(mgcv); library(data.table); library(ggplot2)
  library(ggrepel); library(patchwork)
})

# Fit pooled GAM and return tidy result on a 200-point x grid.
fit_pooled_gam <- function(x, y, k = 10, n_grid = 200, x_range = NULL) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(y) < 1000 || sd(y) == 0) return(NULL)
  d <- data.table(x = x, y = y)
  f <- mgcv::gam(y ~ s(x, k = k), data = d, method = "REML")
  if (is.null(x_range)) x_range <- range(d$x)
  xg <- seq(x_range[1], x_range[2], length.out = n_grid)
  pr <- predict(f, data.frame(x = xg), se.fit = TRUE)
  data.table(
    x   = xg,
    fit = as.numeric(pr$fit),
    lo  = as.numeric(pr$fit - 1.96 * pr$se.fit),
    hi  = as.numeric(pr$fit + 1.96 * pr$se.fit),
    dev_expl = summary(f)$dev.expl,
    n_cells  = nrow(d)
  )
}

# z-score a pooled fit using its own mean/sd; ribbon is z-scored consistently
zscore_fit <- function(d) {
  mu <- mean(d$fit); sd_ <- stats::sd(d$fit)
  if (sd_ == 0) sd_ <- 1
  d[, fit_z := (fit - mu) / sd_]
  d[, lo_z  := (lo  - mu) / sd_]
  d[, hi_z  := (hi  - mu) / sd_]
  d
}

# UCell-like score (rank-based mean signature) for a small gene set.
# Uses simple rank-mean (not the bidirectional UCell exact form, but matches
# the Phase 9b approximation pattern used elsewhere).
compute_signature_score <- function(logc_mat, gene_set) {
  hit <- intersect(gene_set, rownames(logc_mat))
  if (length(hit) < 2) {
    warning("Fewer than 2 genes available for signature; returning NA")
    return(rep(NA_real_, ncol(logc_mat)))
  }
  sub <- as.matrix(logc_mat[hit, , drop = FALSE])
  # rank within cell, rescale to [0,1], mean across genes
  per_cell_score <- function(j) {
    v <- sub[, j]
    if (all(v == v[1])) return(0.5)
    r <- rank(v) / length(v)
    mean(r)
  }
  vapply(seq_len(ncol(sub)), per_cell_score, numeric(1))
}

# Build the standardised overlay figure from a list of curves.
# `curves` is a list of data.tables with columns:
#   x, fit_z, lo_z, hi_z, group, name, color
# Optionally each curve can have `direction` ("up","down","flat") for an arrow.
phase36_overlay_plot <- function(
  curves,
  x_label = "polarization_UCell  (SecA → SecB)",
  y_label = "z-scored pooled GAM  (per-curve)",
  title    = NULL,
  subtitle = NULL,
  x_extra_right = 0.10,    # extra x space to fit ggrepel end labels
  ribbon_alpha  = 0.15,
  line_size     = 1.05,
  ref_line      = 0,
  legend_position = "bottom",
  base_size = 10
) {
  d <- rbindlist(curves, fill = TRUE)
  # End-label table: rightmost x for each curve
  label_dt <- d[d[, .I[which.max(x)], by = name]$V1]
  x_range  <- range(d$x)
  x_max    <- x_range[2] + x_extra_right * diff(x_range)

  # Color vector, group → color (constant within group)
  group_cols <- unique(d[, .(group, color)])
  group_pal  <- setNames(group_cols$color, group_cols$group)

  # Per-curve color for end-label text
  curve_pal <- setNames(d[!duplicated(name)]$color, d[!duplicated(name)]$name)

  p <- ggplot() +
    geom_hline(yintercept = ref_line, colour = "grey85", linewidth = 0.3) +
    # 95% CI ribbon per curve, low alpha
    geom_ribbon(data = d,
                aes(x = x, ymin = lo_z, ymax = hi_z, fill = group,
                    group = name),
                alpha = ribbon_alpha, colour = NA) +
    geom_line(data = d,
              aes(x = x, y = fit_z, colour = group, group = name),
              linewidth = line_size, alpha = 0.95) +
    scale_fill_manual(values = group_pal, name = NULL) +
    scale_colour_manual(values = group_pal, name = NULL) +
    ggrepel::geom_text_repel(data = label_dt,
                              aes(x = x, y = fit_z, label = name,
                                  colour = group),
                              size = 3.0, fontface = "bold",
                              direction = "y",
                              nudge_x = 0.04 * diff(x_range),
                              hjust = 0,
                              segment.size = 0.2,
                              segment.colour = "grey60",
                              min.segment.length = 0,
                              max.overlaps = Inf, seed = 1,
                              show.legend = FALSE) +
    coord_cartesian(xlim = c(x_range[1], x_max), clip = "off") +
    labs(x = x_label, y = y_label, title = title, subtitle = subtitle) +
    theme_lab(base_size = base_size) +
    theme(
      legend.position = legend_position,
      legend.box      = "horizontal",
      legend.key.height = unit(0.4, "cm"),
      legend.text     = element_text(size = rel(0.85)),
      plot.title      = element_text(size = rel(1.1), face = "bold"),
      plot.subtitle   = element_text(size = rel(0.85), colour = "grey35",
                                       margin = margin(b = 6)),
      plot.margin     = margin(8, 14, 6, 8)
    )
  p
}
