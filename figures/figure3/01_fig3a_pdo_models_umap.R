#!/usr/bin/env Rscript
# ============================================================================
# Figure 3a — Baseline organoid (PDO) UMAP by model
# ----------------------------------------------------------------------------
# PURPOSE: UMAP of the eight PDO models (snRNA-seq), coloured by origin
#   (ascites / primary) and labelled by model, at the final letter-page
#   footprint (3.0 x 3.0 in) with final typography (8 pt text / 6 pt axis
#   ticks). Drawn at its true footprint (no tight-bbox crop) so 8 pt renders
#   true when the SVG is placed 1:1 in Affinity.
#
# INPUTS  : <organoids_root>/output/01_Data_Processing_and_QC/
#             seurat_untreated_baseline.rds
# OUTPUTS : figures_dir/figure3/organoids_models_umap_letter.{png,svg}
# MANUSCRIPT PANEL(S): Fig 3a. RUNTIME TIER: moderate.
# ============================================================================

.fig_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(.fig_dir) || !nzchar(.fig_dir)) .fig_dir <- "figures/figure3"
.config_path <- normalizePath(file.path(.fig_dir, "..", "..", "config", "config.R"),
                              mustWork = FALSE)
if (!file.exists(.config_path)) .config_path <- "config/config.R"
source(.config_path)

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(data.table); library(ggrepel)
})

set.seed(CFG$seed)

# --- Paths -------------------------------------------------------------------
INPUT   <- cfg_path("organoids_root",
                    "output/01_Data_Processing_and_QC/seurat_untreated_baseline.rds")
FIGDIR  <- cfg_path("figures_dir", "figure3")
if (!dir.exists(FIGDIR)) dir.create(FIGDIR, recursive = TRUE)
OUT_PNG <- file.path(FIGDIR, "organoids_models_umap_letter.png")
OUT_SVG <- sub("\\.png$", ".svg", OUT_PNG)

# --- Final letter-page typography: 8 pt text / 6 pt ticks / 6 pt labels ------
FA <- 8;  FK <- 6;  FN <- 6

MODEL_LEVELS   <- c("OCAD106", "OCAD93", "OCAD96", "OCAD97",
                    "OPTO112", "OPTO129", "OPTO98", "PDO66")
ASCITES_MODELS <- c("OCAD106", "OCAD93", "OCAD96", "OCAD97")
# Atlas-palette recolour: ascites = B cell (dark blue/indigo), primary = T/NK (light blue)
ASCITES_FILL   <- "#5665B6"   # atlas B cell
PRIMARY_FILL   <- "#87CEFA"   # atlas T/NK cell

# --- Load and extract --------------------------------------------------------
cat("Loading: ", INPUT, "\n", sep = "")
stopifnot(file.exists(INPUT))
s <- readRDS(INPUT)
cat("  cells: ", ncol(s), "  genes: ", nrow(s), "\n", sep = "")

stopifnot("umap"      %in% names(s@reductions))
stopifnot("pdo_model" %in% colnames(s@meta.data))

emb <- as.data.table(Embeddings(s, "umap"), keep.rownames = "cell")
setnames(emb, c("cell", "UMAP1", "UMAP2"))
emb[, pdo_model := factor(s@meta.data$pdo_model[match(cell, colnames(s))],
                          levels = MODEL_LEVELS)]
stopifnot(!any(is.na(emb$pdo_model)))

emb[, origin := factor(
  ifelse(as.character(pdo_model) %in% ASCITES_MODELS, "ascites", "primary"),
  levels = c("ascites", "primary"))]

emb <- emb[sample(.N)]

# Per-model centroids + an OUTWARD nudge (away from the global centroid) so each
# label is pushed off its cluster into nearby empty space; ggrepel then resolves
# any residual label/label overlap and adds a thin leader only when it must.
labels_dt <- emb[, .(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)), by = pdo_model]
g1 <- mean(labels_dt$UMAP1); g2 <- mean(labels_dt$UMAP2)
labels_dt[, rr := pmax(sqrt((UMAP1 - g1)^2 + (UMAP2 - g2)^2), 1e-6)]
NUDGE <- 3.4   # UMAP units to push each label radially outward
labels_dt[, nx := UMAP1 + NUDGE * (UMAP1 - g1) / rr]
labels_dt[, ny := UMAP2 + NUDGE * (UMAP2 - g2) / rr]

# --- Plot (no in-plot legend; user adds the colour key externally in Affinity) -
p <- ggplot(emb, aes(UMAP1, UMAP2, colour = origin)) +
  geom_point(size = 0.5, alpha = 0.85, shape = 16, stroke = 0) +
  scale_colour_manual(values = c(ascites = ASCITES_FILL, primary = PRIMARY_FILL),
                      guide = "none") +
  ggrepel::geom_text_repel(
    data = labels_dt, aes(x = UMAP1, y = UMAP2, label = pdo_model),
    inherit.aes = FALSE, size = FA / .pt, colour = "black",
    family = "Helvetica", fontface = "plain",
    nudge_x = labels_dt$nx - labels_dt$UMAP1,
    nudge_y = labels_dt$ny - labels_dt$UMAP2,
    box.padding = 0.5, point.padding = 0.6, min.segment.length = 0.4,
    segment.size = 0.2, segment.color = "grey45", max.overlaps = Inf,
    seed = CFG$seed) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  coord_fixed() +
  theme_classic(base_size = FA) +
  theme(
    text             = element_text(family = "Helvetica", colour = "black"),
    axis.line        = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.3, colour = "black"),
    axis.title       = element_text(size = FA, colour = "black"),
    axis.title.x     = element_text(margin = margin(t = 2)),
    axis.title.y     = element_text(margin = margin(r = 2)),
    axis.text        = element_text(size = FK, colour = "black"),
    legend.position  = "none",
    plot.margin      = margin(2, 2, 2, 2),
    plot.title       = element_blank())

# 3.0 x 3.0 in = 76.2 x 76.2 mm. ggsave writes the device at this exact size
# (no tight crop), so 8 pt text renders true when placed 1:1 in Affinity.
ggsave(OUT_PNG, p, width = 76.2, height = 76.2, units = "mm", dpi = 600)
ggsave(OUT_SVG, p, width = 76.2, height = 76.2, units = "mm")
cat("\nSaved:\n  ", basename(OUT_PNG), "\n  ", basename(OUT_SVG), "\n", sep = "")
