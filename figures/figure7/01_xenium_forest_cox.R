#!/usr/bin/env Rscript
# ============================================================================
# Figure 7A,7B — Forest plots, univariate Cox PH (Xenium TMA)
# ----------------------------------------------------------------------------
# PURPOSE
#   Side-by-side forest plots of univariate Cox PH hazard ratios for cell
#   densities + key secretory ratios / polarization, for 5-year OS (7A, n=97)
#   and PFS (7B, n=95). Features ordered by OS HR; ratio/polarization rows
#   highlighted.
#
# INPUTS
#   data_root/2026_final_xenium_analysis/output/10_clinical_v2/
#     cox_univariate_results.csv   (pre-computed Cox results)
#     per_patient_features_v2.csv  (feature reference)
#   Shared helpers: config/config.R, spatial/00_setup/00_setup.R (ref_palette, theme_lab).
#
# OUTPUTS
#   figures_dir/figure7/xenium_forest_cox_combined.{pdf,png,svg}
#
# MANUSCRIPT PANEL(S): Fig 7A (5-yr OS forest), Fig 7B (5-yr PFS forest)
# RUNTIME TIER: fast (reads pre-computed Cox CSV)
# ============================================================================

.here     <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(.here) || !nzchar(.here)) .here <- getwd()
source(file.path(.here, "..", "..", "config", "config.R"))
source(file.path(.here, "..", "..", "spatial", "00_setup", "00_setup.R"))
library(survival)

FIG_DIR <- cfg_path("figures_dir", "figure7")

# --- Font size constants (per style guide) ----------------------------------
FA <- 6
FK <- 5.5
FN <- 5

# --- Load Cox results --------------------------------------------------------
cox_raw <- read.csv(cfg_path("data_root", "2026_final_xenium_analysis", "output",
                             "10_clinical_v2", "cox_univariate_results.csv"),
                    stringsAsFactors = FALSE)

# --- Feature selection: cell densities + key secretory ratios ---------------
density_features <- grep("^dens_", unique(cox_raw$feature), value = TRUE)
ratio_features   <- c("prop_SecA_of_sec", "prop_SecB_of_sec",
                       "polar_mean")
keep_features    <- c(density_features, ratio_features)

cox_subset <- cox_raw[cox_raw$feature %in% keep_features, ]

# --- Feature label mapping (dens_X -> readable name) ------------------------
feature_label_map <- c(
  dens_Ciliated_epithelium          = "Ciliated epithelium",
  dens_SecA_epithelium              = "SecA epithelium",
  dens_SecB_epithelium              = "SecB epithelium",
  dens_Fibroblast                   = "Fibroblast",
  dens_Smooth_muscle                = "Smooth muscle",
  dens_Mesothelial                  = "Mesothelial",
  dens_Pericyte                     = "Pericyte",
  dens_Endothelial                  = "Endothelial",
  dens_T_cell                       = "T cell",
  dens_NK_cell                      = "NK cell",
  dens_B_cell                       = "B cell",
  dens_Plasma_cell                  = "Plasma cell",
  dens_Macrophage                   = "Macrophage",
  dens_Conventional_dendritic_cell  = "Conventional dendritic cell",
  dens_Plasmacytoid_dendritic_cell  = "Plasmacytoid dendritic cell",
  dens_Neutrophil                   = "Neutrophil",
  dens_Mast_cell                    = "Mast cell",
  prop_SecA_of_sec                  = "SecA proportion",
  prop_SecB_of_sec                  = "SecB proportion",
  polar_mean                        = "Polarization (mean UCell)"
)

cox_subset$feature_label <- feature_label_map[cox_subset$feature]
cox_subset <- cox_subset[!is.na(cox_subset$feature_label), ]

# --- Shared y-axis order (based on OS HR) -----------------------------------
os_data <- cox_subset[cox_subset$endpoint == "OS_5yr" & !is.na(cox_subset$HR), ]
os_data$is_ratio <- grepl("prop_Sec|polar", os_data$feature)
os_data <- os_data[order(!os_data$is_ratio, os_data$HR), ]
shared_levels <- feature_label_map[os_data$feature]

# --- Forest plot helper ------------------------------------------------------
plot_forest <- function(forest_df, endpoint_filter, endpoint_label,
                        y_levels, show_y_axis = TRUE) {
  fd <- forest_df[forest_df$endpoint == endpoint_filter & !is.na(forest_df$HR), ]
  fd$is_ratio <- grepl("prop_Sec|polar", fd$feature)
  fd$feature_label <- factor(fd$feature_label, levels = y_levels)

  fd$hr_text <- sprintf("%.2f (%.2f-%.2f)", fd$HR, fd$HR_lower, fd$HR_upper)
  fd$p_text  <- ifelse(fd$p_value < 0.001, "<.001",
                       sub("^0\\.", ".", formatC(fd$p_value, format = "f", digits = 3)))

  fd$pt_color <- ref_palette[as.character(fd$feature_label)]
  fd$pt_color[grepl("polar", fd$feature)]        <- "#C08E48"
  fd$pt_color[fd$feature == "prop_SecA_of_sec"]   <- "#E6A141"
  fd$pt_color[fd$feature == "prop_Intermediate_of_sec"] <- "#C08E48"
  fd$pt_color[fd$feature == "prop_SecB_of_sec"]   <- "#9A7D55"
  fd$pt_color[is.na(fd$pt_color)] <- "grey50"

  x_hr_pos <- 3.0
  x_p_pos  <- 10.0
  n_feat   <- length(y_levels)

  ratio_idx <- which(y_levels %in%
    fd$feature_label[fd$is_ratio])

  p <- ggplot(fd, aes(x = HR, y = feature_label))

  for (ri in ratio_idx) {
    p <- p +
      annotate("rect",
               xmin = 0.15, xmax = 18,
               ymin = ri - 0.45, ymax = ri + 0.45,
               fill = "#E6A141", alpha = 0.08)
  }

  p <- p +
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "grey50", linewidth = 0.3) +
    geom_errorbar(aes(xmin = HR_lower, xmax = HR_upper, color = pt_color),
                  width = 0.2, linewidth = 0.35, show.legend = FALSE,
                  orientation = "y") +
    geom_point(aes(fill = pt_color, color = pt_color), size = 2.2, shape = 21,
               stroke = 0.5, alpha = 0.7, show.legend = FALSE) +
    geom_point(aes(color = pt_color), size = 2.2, shape = 1,
               stroke = 0.5, show.legend = FALSE) +
    scale_color_identity() +
    scale_fill_identity() +
    geom_text(aes(x = x_hr_pos, label = hr_text),
              hjust = 0, size = FN / .pt, color = "black") +
    geom_text(aes(x = x_p_pos, label = p_text),
              hjust = 0, size = FN / .pt, color = "black") +
    annotate("text", x = x_hr_pos, y = n_feat + 0.8,
             label = "HR (95% CI)", hjust = 0,
             size = FN / .pt, fontface = "bold") +
    annotate("text", x = x_p_pos, y = n_feat + 0.8,
             label = "p", hjust = 0,
             size = FN / .pt, fontface = "bold") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2),
                  labels = c("0.25", "0.5", "1.0", "2.0")) +
    scale_y_discrete(drop = FALSE) +
    coord_cartesian(xlim = c(0.2, 18), clip = "off") +
    labs(x = paste0("Hazard Ratio - ", endpoint_label), y = NULL) +
    theme_lab(base_size = 6) +
    theme(
      plot.title         = element_blank(),
      axis.text.y        = element_text(size = FK),
      axis.text.x        = element_text(size = FK),
      axis.title.x       = element_text(size = FA),
      panel.grid.major.y = element_line(color = "grey93", linewidth = 0.15),
      plot.margin        = margin(12, 45, 4, 4)
    )

  if (!show_y_axis) {
    p <- p + theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank())
  }

  p
}

# --- Build side-by-side combined figure -------------------------------------
p_os  <- plot_forest(cox_subset, "OS_5yr",  "5-yr OS",
                     shared_levels, show_y_axis = TRUE)
p_pfs <- plot_forest(cox_subset, "PFS_5yr", "5-yr PFS",
                     shared_levels, show_y_axis = FALSE)

fig_forest <- p_os + p_pfs + plot_layout(widths = c(1, 1))

n_feats <- length(shared_levels)
w_in <- 210 / 25.4
h_in <- max(60, 12 + n_feats * 4) / 25.4

fname <- "xenium_forest_cox_combined"

pdf(file.path(FIG_DIR, paste0(fname, ".pdf")),
    width = w_in, height = h_in)
print(fig_forest)
invisible(dev.off())
ggsave(file.path(FIG_DIR, paste0(fname, ".png")),
       fig_forest, width = w_in, height = h_in, dpi = 450, bg = "white")
ggsave(file.path(FIG_DIR, paste0(fname, ".svg")),
       fig_forest, width = w_in, height = h_in, bg = "white")

message("  Saved: ", fname)
message("\nForest plots complete. Combined PDF + PNG in: ", FIG_DIR)
