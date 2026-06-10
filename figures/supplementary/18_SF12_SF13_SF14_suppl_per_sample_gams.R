#!/usr/bin/env Rscript
# =============================================================================
# SF12 / SF13 / SF14 — Per-sample GAM curves
# -----------------------------------------------------------------------------
# Purpose:  Per-WT-sample GAM variability, three figure panels:
#   SF12 (env, 4 panels)  — epi density, vascular distance, matrix remodeling, hypoxia
#   SF13 (epi, ~14)       — UCell pathways (proliferation, myc, pi3k, wnt, DNA damage,
#                           epigenetic, NFkB, glycolysis) + genes (LDHA, CDH2, CTNNB1,
#                           ITGB5, ICAM1, MMP7, TGM2)
#   SF14 (mac, 8)         — glycolysis, NFkB + 6 macrophage genes
#   Each panel: thin per-sample GAMs + thick black dashed pooled GAM.
#   Dashed = chemo-treated, solid = treatment-naive.
#
# INPUTS  (all under output_root):
#   28_glandular_architecture/per_cell_architecture_wt.rds          (SF12)
#   22_vascular_proximity/vascular_distance_all_cells.csv           (SF12)
#   16b_niche_succession_gams_v2/neighborhood_features.rds          (SF12)
#   19d_gene_polarization_gams/epithelial_expression_polarization.rds (SF13)
#   19e_gene_gams_all_celltypes/tme_expression_polarization.rds     (SF14)
#   9b_scoring/pathway_gene_sets_v2.csv                             (SF13/SF14 UCell sets)
# OUTPUTS  (output_root/figures/supplementary/):
#   SF12_suppl_gam_env.{svg,pdf}
#   SF13_suppl_gam_epi.{svg,pdf}
#   SF14_suppl_gam_mac.{svg,pdf}
#   SF12_14_suppl_gam_legend.{svg,pdf}
# MANUSCRIPT PANEL(S):  SF12A-D, SF13A-S, SF14A-H.
# RUNTIME TIER:  heavy (many GAM fits over stratified subsamples).
# =============================================================================

# Config + shared spatial setup + helpers (2 levels deep -> ../../).
source(file.path("..", "..", "config", "config.R"))
source(file.path("..", "..", "spatial", "00_setup", "00_setup.R"))    # out_dir, theme_lab
source(file.path("..", "..", "spatial", "00_setup", "36_helpers.R"))  # compute_signature_score

suppressPackageStartupMessages({
  library(mgcv)
  library(data.table)
  library(patchwork)
  library(ggplot2)
})

OUT_DIR <- cfg_path("output_root", "figures", "supplementary")

FA <- 6; FK <- 5.5; FN <- 5

# =============================================================================
# SAMPLE METADATA  (published whole-tissue cohort = config.cohort.whole_tissue)
# =============================================================================

wt <- CFG$cohort$whole_tissue
sample_meta <- data.table(
  sample_id  = paste0("sfe_", wt),
  # treatment status per whole-tissue sample (chemo-treated marked with *)
  treatment  = c("naive", "chemo", "chemo", "naive",
                 "chemo", "naive", "naive", "naive")[match(
                   wt, c("OTB_2384", "OTB_2417", "OTB_2432", "OTB_2454",
                         "OTB_2457", "OTB_2461", "SP24_24824", "SP24_25573"))],
  short_name = c("2384", "2417*", "2432*", "2454",
                 "2457*", "2461", "24824", "25573")[match(
                   wt, c("OTB_2384", "OTB_2417", "OTB_2432", "OTB_2454",
                         "OTB_2457", "OTB_2461", "SP24_24824", "SP24_25573"))]
)

sample_colors <- c(
  "sfe_OTB_2384"   = "#E64B35",
  "sfe_OTB_2454"   = "#F39B7F",
  "sfe_OTB_2461"   = "#4DBBD5",
  "sfe_SP24_24824" = "#00A087",
  "sfe_SP24_25573" = "#3C5488",
  "sfe_OTB_2417"   = "#8491B4",
  "sfe_OTB_2432"   = "#91D1C2",
  "sfe_OTB_2457"   = "#B09C85"
)

sample_linetypes <- setNames(
  ifelse(sample_meta$treatment == "chemo", "dashed", "solid"),
  sample_meta$sample_id
)
sample_labels <- setNames(
  ifelse(sample_meta$treatment == "chemo",
         paste0(sample_meta$short_name, " (chemo)"),
         sample_meta$short_name),
  sample_meta$sample_id
)
sample_ids <- sample_meta$sample_id

# =============================================================================
# HELPERS
# =============================================================================

subsample_stratified <- function(df, x_col = "polarization_UCell",
                                  n_per = 5000, nbins = 20) {
  pol <- df[[x_col]]
  pol_noNA <- pol[!is.na(pol)]
  if (length(pol_noNA) < 100) return(df)
  breaks <- quantile(pol_noNA, probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE)
  breaks <- unique(breaks)
  df$pol_bin <- cut(pol, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  per_bin <- max(1, round(n_per / length(unique(df$pol_bin[!is.na(df$pol_bin)]))))
  idx_list <- lapply(split(seq_len(nrow(df)), df$pol_bin), function(ii) {
    if (length(ii) <= per_bin) ii else sample(ii, per_bin)
  })
  out <- df[unlist(idx_list), ]
  out$pol_bin <- NULL
  out
}

fit_and_predict <- function(y, x, k = 10, n_pred = 100) {
  ok <- !is.na(y) & !is.na(x)
  if (sum(ok) < 50) return(NULL)
  y <- y[ok]; x <- x[ok]
  tryCatch({
    fit <- gam(y ~ s(x, bs = "tp", k = k), method = "REML", family = gaussian())
    newx <- seq(min(x), max(x), length.out = n_pred)
    pred <- predict(fit, newdata = data.frame(x = newx), type = "response")
    data.frame(x_val = newx, predicted = as.numeric(pred))
  }, error = function(e) NULL)
}

gam_direction <- function(pred_df) {
  if (is.null(pred_df) || nrow(pred_df) < 2) return(NA_real_)
  sign(pred_df$predicted[nrow(pred_df)] - pred_df$predicted[1])
}

make_gam_panel <- function(data_sub, data_pooled, feature_col, x_col,
                           feature_label, sample_ids_use, k = 10,
                           colors = sample_colors, linetypes = sample_linetypes,
                           labels = sample_labels) {
  pred_list <- list()
  sample_directions <- c()
  for (sid in sample_ids_use) {
    chunk <- data_sub[data_sub$sample_id == sid, ]
    if (nrow(chunk) < 50) next
    y_vals <- chunk[[feature_col]]
    x_vals <- chunk[[x_col]]
    if (all(is.na(y_vals)) || all(is.na(x_vals))) next
    pred <- fit_and_predict(y_vals, x_vals, k = k)
    if (!is.null(pred)) {
      pred$sample_id <- sid
      pred_list[[sid]] <- pred
      sample_directions[sid] <- gam_direction(pred)
    }
  }
  if (length(pred_list) == 0) return(NULL)
  pred_df <- do.call(rbind, pred_list)

  pooled_pred <- fit_and_predict(data_pooled[[feature_col]], data_pooled[[x_col]], k = k)
  pooled_dir <- gam_direction(pooled_pred)

  n_tested <- length(sample_directions)
  n_agree  <- sum(sample_directions == pooled_dir, na.rm = TRUE)
  consist_label <- paste0(n_agree, "/", n_tested, " consistent")

  p <- ggplot(pred_df, aes(x = x_val, y = predicted,
                           color = sample_id, linetype = sample_id)) +
    geom_line(linewidth = 0.3, alpha = 0.85) +
    scale_color_manual(values = colors, labels = labels, name = "Sample") +
    scale_linetype_manual(values = linetypes, labels = labels, name = "Sample") +
    labs(title = feature_label, x = "Polarization (UCell)", y = feature_label) +
    theme_lab(base_size = 6) +
    theme(legend.position = "none",
          plot.title = element_text(size = FA, face = "bold"),
          axis.title = element_text(size = FK),
          axis.text = element_text(size = FN))

  if (!is.null(pooled_pred)) {
    p <- p + geom_line(data = pooled_pred, aes(x = x_val, y = predicted),
                       color = "black", linewidth = 0.5, linetype = "dashed",
                       inherit.aes = FALSE)
  }
  p <- p + annotate("text", x = Inf, y = Inf, label = consist_label,
                    hjust = 1.05, vjust = 1.3, size = FN / .pt,
                    color = "grey30", fontface = "italic")
  p
}

# =============================================================================
# LOAD PATHWAY GENE SETS
# =============================================================================

message("Loading pathway gene sets...")
pw_df <- read.csv(cfg_path("data_root", "2026_final_xenium_analysis", "output", "9b_scoring", "pathway_gene_sets_v2.csv"),
                  stringsAsFactors = FALSE)

get_gene_set <- function(pathway_name) pw_df$gene[pw_df$pathway == pathway_name]

# =============================================================================
# SF12: ENVIRONMENTAL GAMs (4 panels)
# =============================================================================

message("\n=== SF12: Environmental GAMs ===")

message("  Loading architecture data...")
arch <- readRDS(cfg_path("data_root", "2026_final_xenium_analysis", "output", "28_glandular_architecture", "per_cell_architecture_wt.rds"))
arch$polarization_UCell <- arch$polarization

message("  Loading vascular distance...")
# Vascular cache is read from output_root: it is REGENERATED by
# spatial/05_gradients_gams/04_vascular_proximity.R, which now carries a per-cell
# `cell_id` (SFE barcode) so it can be keyed to `arch`. The deposited copy under
# data_root lacks `cell_id` and a shared cell ordering, so the old positional
# join (stopifnot(nrow==) + column assignment) was invalid and is replaced with a
# keyed left-merge of dist_to_vascular onto arch by (sample_id, cell_id).
vasc <- read.csv(cfg_path("output_root", "22_vascular_proximity", "vascular_distance_all_cells.csv"),
                 stringsAsFactors = FALSE)
vasc_wt <- vasc[vasc$tissue == "whole_tissue", c("sample_id", "cell_id", "dist_to_vascular")]
stopifnot("cell_id" %in% colnames(vasc_wt), "cell_id" %in% colnames(arch))
.vd <- vasc_wt$dist_to_vascular[
  match(paste(arch$sample_id, arch$cell_id), paste(vasc_wt$sample_id, vasc_wt$cell_id))
]
arch$dist_to_vascular <- .vd
message(sprintf("  Vascular distance merged: %d of %d arch cells matched (%.1f%%), %d NA",
                sum(!is.na(.vd)), nrow(arch), 100 * mean(!is.na(.vd)), sum(is.na(.vd))))
rm(vasc, vasc_wt, .vd)

# Architecture sample_ids lack "sfe_" prefix — harmonize
arch_sid_map <- setNames(gsub("^sfe_", "", sample_ids), sample_ids)
arch_sample_ids_use <- intersect(arch_sid_map, unique(arch$sample_id))
arch_colors <- setNames(sample_colors, arch_sid_map[names(sample_colors)])
arch_linetypes <- setNames(sample_linetypes, arch_sid_map[names(sample_linetypes)])
arch_labels <- setNames(sample_labels, arch_sid_map[names(sample_labels)])

set.seed(CFG$seed)
message("  Subsampling arch data...")
arch_sub <- as.data.frame(do.call(rbind, lapply(arch_sample_ids_use, function(sid) {
  subsample_stratified(as.data.frame(arch[arch$sample_id == sid, ]), n_per = 5000)
})))
arch_pooled <- as.data.frame(subsample_stratified(as.data.frame(arch), n_per = 8000))

message("  Loading neighborhood features...")
nf <- readRDS(cfg_path("data_root", "2026_final_xenium_analysis", "output", "16b_niche_succession_gams_v2", "neighborhood_features.rds"))
nf_wt <- nf[nf$sample_id != "TMA", ]

set.seed(CFG$seed)
message("  Subsampling nf data...")
nf_sub <- as.data.frame(do.call(rbind, lapply(sample_ids, function(sid) {
  subsample_stratified(as.data.frame(nf_wt[nf_wt$sample_id == sid, ]), n_per = 5000)
})))
nf_pooled <- as.data.frame(subsample_stratified(as.data.frame(nf_wt), n_per = 8000))

message("  Fitting GAMs...")
plots_env <- list()
plots_env[[1]] <- make_gam_panel(
  arch_sub, arch_pooled, "epi_neighbors_50um", "polarization_UCell",
  "Epithelial density (50 µm)", arch_sample_ids_use, k = 10,
  colors = arch_colors, linetypes = arch_linetypes, labels = arch_labels)
plots_env[[2]] <- make_gam_panel(
  arch_sub, arch_pooled, "dist_to_vascular", "polarization_UCell",
  "Dist. to vascular cell (µm)", arch_sample_ids_use, k = 10,
  colors = arch_colors, linetypes = arch_linetypes, labels = arch_labels)
plots_env[[3]] <- make_gam_panel(
  nf_sub, nf_pooled, "nb_mean_pathway_matrix_remodeling", "polarization_UCell",
  "Matrix remodeling (UCell)", sample_ids, k = 10)
plots_env[[4]] <- make_gam_panel(
  nf_sub, nf_pooled, "nb_mean_pathway_hypoxia", "polarization_UCell",
  "Hypoxia (UCell)", sample_ids, k = 10)
plots_env <- plots_env[!sapply(plots_env, is.null)]
message("  Environmental panels: ", length(plots_env))

# =============================================================================
# SF13: EPITHELIAL-SPECIFIC GAMs (~14 panels)
# =============================================================================

message("\n=== SF13: Epithelial-specific GAMs ===")
message("  Loading epithelial expression data...")
epi <- readRDS(cfg_path("data_root", "2026_final_xenium_analysis", "output", "19d_gene_polarization_gams",
                        "epithelial_expression_polarization.rds"))
setDT(epi)

set.seed(CFG$seed)
message("  Subsampling epithelial cells...")
epi_sub <- do.call(rbind, lapply(sample_ids, function(sid) {
  subsample_stratified(epi[epi$sample_id == sid, ], n_per = 5000)
}))
epi_pooled <- subsample_stratified(epi, n_per = 8000)

message("  Computing UCell pathway scores...")
gene_cols <- setdiff(colnames(epi_sub),
                     c("cell_id", "polarization_UCell", "sample_id", "cell_label"))
logc_sub <- as.matrix(epi_sub[, ..gene_cols])
rownames(logc_sub) <- NULL
gene_cols_pooled <- setdiff(colnames(epi_pooled),
                            c("cell_id", "polarization_UCell", "sample_id", "cell_label"))
logc_pooled <- as.matrix(epi_pooled[, ..gene_cols_pooled])

epi_pathways <- list(
  "Proliferation" = get_gene_set("proliferation"),
  "Myc"           = get_gene_set("myc_activating"),
  "PI3K"          = get_gene_set("pi3k_activating"),
  "Wnt"           = get_gene_set("wnt_activating"),
  "DNA damage"    = get_gene_set("p53_activating"),
  "NFkB"          = get_gene_set("nfkb"),
  "Glycolysis"    = get_gene_set("glycolysis")
)

epigenetic_gs <- get_gene_set("epigenetic_remodeling")
if (length(epigenetic_gs) == 0) epigenetic_gs <- get_gene_set("epigenetic")
if (length(epigenetic_gs) > 0) {
  epi_pathways[["Epigenetic remodeling"]] <- epigenetic_gs
  message("    Epigenetic remodeling gene set found: ", length(epigenetic_gs), " genes")
} else {
  message("    NOTE: Epigenetic remodeling gene set not in pathway_gene_sets_v2.csv — skipping")
}

for (pw_name in names(epi_pathways)) {
  gs <- epi_pathways[[pw_name]]
  col_name <- paste0("pw_", gsub(" |\\(|\\)", "_", pw_name))
  epi_sub[[col_name]] <- compute_signature_score(t(logc_sub), gs)
  epi_pooled[[col_name]] <- compute_signature_score(t(logc_pooled), gs)
  message("    ", pw_name, " scored")
}

epi_genes <- c("TGFB1", "INHBA", "IL32", "CD55", "C7",
               "LDHA", "CDH2", "CTNNB1", "ITGB5", "ICAM1", "MMP7", "TGM2")

message("  Fitting epithelial GAMs...")
plots_epi <- list()
for (pw_name in names(epi_pathways)) {
  col_name <- paste0("pw_", gsub(" |\\(|\\)", "_", pw_name))
  p <- make_gam_panel(epi_sub, epi_pooled, col_name, "polarization_UCell",
                      pw_name, sample_ids, k = 10)
  if (!is.null(p)) plots_epi[[length(plots_epi) + 1]] <- p
}
for (gene in epi_genes) {
  p <- make_gam_panel(epi_sub, epi_pooled, gene, "polarization_UCell",
                      gene, sample_ids, k = 10)
  if (!is.null(p)) plots_epi[[length(plots_epi) + 1]] <- p
}
plots_epi <- plots_epi[!sapply(plots_epi, is.null)]
message("  Epithelial panels: ", length(plots_epi))

# =============================================================================
# SF14: MACROPHAGE-SPECIFIC GAMs (8 panels)
# =============================================================================

message("\n=== SF14: Macrophage-specific GAMs ===")
message("  Loading TME expression data...")
tme <- readRDS(cfg_path("data_root", "2026_final_xenium_analysis", "output", "19e_gene_gams_all_celltypes",
                        "tme_expression_polarization.rds"))
setDT(tme)

mac <- tme[tme$cell_label == "Macrophage", ]
message("  Macrophages: ", nrow(mac))

sid_map <- setNames(gsub("^sfe_", "", sample_ids), sample_ids)
mac_sample_ids_use <- intersect(sid_map, unique(mac$sample_id))
mac_colors <- setNames(sample_colors, sid_map[names(sample_colors)])
mac_linetypes <- setNames(sample_linetypes, sid_map[names(sample_linetypes)])
mac_labels <- setNames(sample_labels, sid_map[names(sample_labels)])

set.seed(CFG$seed)
message("  Subsampling macrophages...")
mac_sub <- do.call(rbind, lapply(mac_sample_ids_use, function(sid) {
  subsample_stratified(mac[mac$sample_id == sid, ],
                       x_col = "nearest_epi_polarization", n_per = 3000)
}))
mac_pooled <- subsample_stratified(mac, x_col = "nearest_epi_polarization", n_per = 6000)

message("  Computing macrophage pathway scores...")
mac_gene_cols <- setdiff(colnames(mac_sub),
                         c("cell_id", "nearest_epi_polarization", "sample_id", "cell_label"))
logc_mac_sub <- as.matrix(mac_sub[, ..mac_gene_cols])
logc_mac_pooled_cols <- setdiff(colnames(mac_pooled),
                                c("cell_id", "nearest_epi_polarization", "sample_id", "cell_label"))
logc_mac_pooled <- as.matrix(mac_pooled[, ..logc_mac_pooled_cols])

mac_pathways <- list(
  "Glycolysis" = get_gene_set("glycolysis"),
  "NFkB"       = get_gene_set("nfkb")
)
for (pw_name in names(mac_pathways)) {
  gs <- mac_pathways[[pw_name]]
  col_name <- paste0("pw_", gsub(" ", "_", pw_name))
  mac_sub[[col_name]] <- compute_signature_score(t(logc_mac_sub), gs)
  mac_pooled[[col_name]] <- compute_signature_score(t(logc_mac_pooled), gs)
  message("    ", pw_name, " scored")
}

mac_genes <- c("A2M", "C1QC", "CXCL10", "CXCL11", "CTSL", "ICAM1")

message("  Fitting macrophage GAMs...")
plots_mac <- list()
for (pw_name in names(mac_pathways)) {
  col_name <- paste0("pw_", gsub(" ", "_", pw_name))
  p <- make_gam_panel(mac_sub, mac_pooled, col_name, "nearest_epi_polarization",
                      paste0("Mac: ", pw_name), mac_sample_ids_use, k = 8,
                      colors = mac_colors, linetypes = mac_linetypes, labels = mac_labels)
  if (!is.null(p)) {
    p <- p + labs(x = "Nearest epi. polarization")
    plots_mac[[length(plots_mac) + 1]] <- p
  }
}
for (gene in mac_genes) {
  p <- make_gam_panel(mac_sub, mac_pooled, gene, "nearest_epi_polarization",
                      paste0("Mac: ", gene), mac_sample_ids_use, k = 8,
                      colors = mac_colors, linetypes = mac_linetypes, labels = mac_labels)
  if (!is.null(p)) {
    p <- p + labs(x = "Nearest epi. polarization")
    plots_mac[[length(plots_mac) + 1]] <- p
  }
}
message("  Macrophage panels: ", length(plots_mac))

# =============================================================================
# ASSEMBLE AND SAVE
# =============================================================================

message("\n=== Saving figures ===")
PG_W <- 8.5
NCOL_MAX <- 4

if (length(plots_env) > 0) {
  nrow_env <- ceiling(length(plots_env) / NCOL_MAX)
  fig_env <- wrap_plots(plots_env, ncol = NCOL_MAX)
  for (ext in c("svg", "pdf")) {
    ggsave(file.path(OUT_DIR, paste0("SF12_suppl_gam_env.", ext)),
           fig_env, width = PG_W, height = PG_W / NCOL_MAX * nrow_env, bg = "white")
  }
  message("  Saved: SF12_suppl_gam_env")
}

if (length(plots_epi) > 0) {
  nrow_epi <- ceiling(length(plots_epi) / NCOL_MAX)
  fig_epi <- wrap_plots(plots_epi, ncol = NCOL_MAX)
  for (ext in c("svg", "pdf")) {
    ggsave(file.path(OUT_DIR, paste0("SF13_suppl_gam_epi.", ext)),
           fig_epi, width = PG_W, height = PG_W / NCOL_MAX * nrow_epi, bg = "white")
  }
  message("  Saved: SF13_suppl_gam_epi")
}

if (length(plots_mac) > 0) {
  nrow_mac <- ceiling(length(plots_mac) / NCOL_MAX)
  fig_mac <- wrap_plots(plots_mac, ncol = NCOL_MAX)
  for (ext in c("svg", "pdf")) {
    ggsave(file.path(OUT_DIR, paste0("SF14_suppl_gam_mac.", ext)),
           fig_mac, width = PG_W, height = PG_W / NCOL_MAX * nrow_mac, bg = "white")
  }
  message("  Saved: SF14_suppl_gam_mac")
}

# Shared sample legend
legend_df <- data.frame(
  x = rep(seq(0, 1, length.out = 10), length(sample_ids)),
  y = rep(seq_along(sample_ids), each = 10),
  sample_id = rep(sample_ids, each = 10)
)
p_legend <- ggplot(legend_df, aes(x = x, y = y, color = sample_id, linetype = sample_id)) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = sample_colors, labels = sample_labels, name = "Sample") +
  scale_linetype_manual(values = sample_linetypes, labels = sample_labels, name = "Sample") +
  guides(color = guide_legend(override.aes = list(linewidth = 2)),
         linetype = guide_legend(override.aes = list(linewidth = 2))) +
  theme_void(base_size = 10) +
  theme(legend.position = "right", legend.key.width = unit(1.5, "cm"))

for (ext in c("svg", "pdf")) {
  ggsave(file.path(OUT_DIR, paste0("SF12_14_suppl_gam_legend.", ext)),
         p_legend, width = 3, height = 2.5, bg = "white")
}
message("  Saved: SF12_14_suppl_gam_legend")

message("\nDone.")
