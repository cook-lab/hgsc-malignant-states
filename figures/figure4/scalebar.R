# Scale bar for the Figure 4 spatial panels (4b, 4c, 4h, 4j).
# Drawn in the plot's own data frame, so it works for both spatialCoords point
# plots (y positive) and cellSeg polygon plots (y negative): pass the limits the
# plot actually displays. Bottom-right corner, black bar with the length above it.
# `panel_h_in` is the height the panel is saved at; bar thickness and label size
# scale with it so they read the same once panels are sized for the page
# (~5.5 pt label for a 60 mm panel).
scalebar_layers <- function(xmin, xmax, ymin, ymax, um, panel_h_in) {
  W <- xmax - xmin; H <- ymax - ymin
  x1 <- xmax - 0.05 * W; x0 <- x1 - um; y0 <- ymin + 0.045 * H
  lab <- if (um >= 1000) paste0(um / 1000, " mm") else paste0(um, " µm")
  text_pt <- 2.33 * panel_h_in
  list(
    # white backing so the bar reads over dense tissue
    ggplot2::annotate("rect", xmin = x0 - 0.025 * W, xmax = x1 + 0.025 * W,
                      ymin = y0 - 0.02 * H, ymax = y0 + 0.068 * H,
                      fill = "white", alpha = 0.85, colour = NA),
    ggplot2::annotate("segment", x = x0, xend = x1, y = y0, yend = y0,
                      linewidth = 0.25 * panel_h_in, lineend = "butt", colour = "black"),
    ggplot2::annotate("text", x = (x0 + x1) / 2, y = y0 + 0.022 * H, label = lab,
                      size = text_pt / ggplot2::.pt, vjust = 0, colour = "black")
  )
}
