# ============================================================================
# 01_secretory_spatial_autocorrelation.R
# ----------------------------------------------------------------------------
# PURPOSE: Spatial autocorrelation of SecA / SecB UCell scores: univariate global+local Moran's I (LISA), bivariate Lee's L, and BiLISA HH/HL/LH/LL regime classification (segregation vs co-occurrence).
#
# INPUTS:
#   - SFEs (load_sfe) with cell_label, SecA_UCell/SecB_UCell/polarization_UCell, coords; 06f override
#
# OUTPUTS:
#   - output/44_spatial_autocorrelation/interpretation_summary.csv
#   - tma_patient_level_lee.csv, bilisa_vs_label_crosstab.csv
#
# MANUSCRIPT PANEL(S): Fig 4E, Fig 4F.
# RUNTIME TIER: heavy
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

library(spdep)
library(dbscan)
library(data.table)

set.seed(CFG$seed)

message("\n", strrep("=", 70))
message("Script 44: Bivariate Spatial Autocorrelation of Secretory UCell Scores")
message(strrep("=", 70))

# ============================================================================
# CONFIGURATION
# ============================================================================

out_path <- file.path(out_dir, "44_spatial_autocorrelation")
fig_path <- file.path(fig_dir, "44_spatial_autocorrelation")
lisa_path <- file.path(out_path, "lisa_results")
bilisa_path <- file.path(out_path, "bilisa_results")
for (d in c(out_path, fig_path, lisa_path, bilisa_path)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

sfe_names <- c("sfe_tma_filtered", sfe_names_wt)

RADIUS       <- 50      # µm — fixed-radius neighborhood (consistent with 09, 16a, 23)
NSIM         <- 999     # permutations for local inference
FDR_ALPHA    <- 0.05
MIN_SEC      <- 200     # minimum secretory cells per sample
N_SUBSAMPLE  <- 15000   # max secretory cells per sample for local statistics

# Secretory cell types
sec_types <- c("SecA epithelium", "Intermediate epithelium", "SecB epithelium")

# UCell score columns (stored in colData by 06d/06f)
score_cols <- c("SecA_UCell", "SecB_UCell", "polarization_UCell")

# BiLISA regime palette
bilisa_pal <- c(
  "HH" = "#D7191C",   # High SecA — High SecB (co-localized / transitioning zones)
  "HL" = "#E6A141",   # High SecA — Low SecB  (SecA-dominant clusters)
  "LH" = "#9A7D55",   # Low SecA  — High SecB (SecB-dominant clusters)
  "LL" = "#2C7BB6",   # Low both  (depleted / non-secretory border)
  "NS" = "#D9D9D9"    # Not significant
)

# LISA cluster palette
lisa_pal <- c(
  "High-High" = "#D7191C",
  "Low-Low"   = "#2C7BB6",
  "High-Low"  = "#FDAE61",
  "Low-High"  = "#ABD9E9",
  "NS"        = "#D9D9D9"
)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Build row-standardized spatial weights from fixed-radius neighbors
#' @param xy Nx2 coordinate matrix
#' @param radius Neighborhood radius in coordinate units (µm)
#' @return listw object (row-standardized)
build_weights <- function(xy, radius) {
  nb <- spdep::dnearneigh(xy, 0, radius)
  spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
}

#' Global Moran's I with permutation test
#' @param x Numeric vector
#' @param listw listw object
#' @param nsim Number of permutations
#' @return list with I, E_I, variance, z, p_value
global_moran <- function(x, listw, nsim = 999) {
  mt <- spdep::moran.mc(x, listw, nsim = nsim, zero.policy = TRUE)
  list(
    I       = mt$statistic,
    E_I     = mean(mt$res[-length(mt$res)]),
    p_value = mt$p.value,
    nsim    = nsim
  )
}

#' Lee's bivariate spatial association statistic
#' Measures spatial correlation between x and spatially-lagged y.
#' L = (n / S0) * (x' W y) / (sqrt(x'x) * sqrt(y'y))
#' where x, y are centered.
#' @param x Numeric vector (SecA_UCell)
#' @param y Numeric vector (SecB_UCell)
#' @param listw listw object (row-standardized)
#' @param nsim Number of permutations for inference
#' @return list with L, p_value
lee_bivariate <- function(x, y, listw, nsim = 999) {
  n <- length(x)
  xc <- x - mean(x)
  yc <- y - mean(y)

  W <- spdep::listw2mat(listw)
  Wy <- as.numeric(W %*% yc)

  S0 <- sum(W)
  L_obs <- (n / S0) * sum(xc * Wy) / (sqrt(sum(xc^2)) * sqrt(sum(yc^2)))

  L_perm <- numeric(nsim)
  for (i in seq_len(nsim)) {
    perm_idx <- sample.int(n)
    yc_perm <- yc[perm_idx]
    Wy_perm <- as.numeric(W %*% yc_perm)
    L_perm[i] <- (n / S0) * sum(xc * Wy_perm) /
      (sqrt(sum(xc^2)) * sqrt(sum(yc_perm^2)))
  }

  p_val <- (sum(abs(L_perm) >= abs(L_obs)) + 1) / (nsim + 1)

  list(L = L_obs, p_value = p_val, nsim = nsim,
       L_perm_mean = mean(L_perm), L_perm_sd = sd(L_perm))
}

#' Local Moran's I (LISA) with FDR-adjusted significance
#' @param x Numeric vector
#' @param listw listw object
#' @param nsim Number of conditional permutations
#' @param alpha FDR threshold
#' @return data.table with Ii, p_adj, quadrant, cluster
local_moran_lisa <- function(x, listw, nsim = 999, alpha = 0.05) {
  lm <- spdep::localmoran_perm(x, listw, nsim = nsim,
                                zero.policy = TRUE, iseed = 42)

  xz <- scale(x)[, 1]
  W <- spdep::listw2mat(listw)
  Wxz <- as.numeric(W %*% xz)

  p_col <- if ("Pr(z != E(Ii)) Sim" %in% colnames(lm)) {
    "Pr(z != E(Ii)) Sim"
  } else {
    grep("^Pr\\(", colnames(lm), value = TRUE)[1]
  }

  p_raw <- lm[, p_col]
  p_adj <- p.adjust(p_raw, method = "BH")

  quadrant <- ifelse(xz >= 0 & Wxz >= 0, "High-High",
              ifelse(xz >= 0 & Wxz < 0,  "High-Low",
              ifelse(xz < 0  & Wxz >= 0, "Low-High",
                                          "Low-Low")))

  cluster <- ifelse(p_adj < alpha, quadrant, "NS")

  data.table(
    Ii       = lm[, "Ii"],
    z_score  = if ("Z.Ii" %in% colnames(lm)) lm[, "Z.Ii"] else NA_real_,
    p_raw    = p_raw,
    p_adj    = p_adj,
    quadrant = quadrant,
    cluster  = cluster
  )
}

#' Local bivariate indicators (BiLISA): x vs spatially-lagged y
#' @param x Numeric vector (SecA_UCell)
#' @param y Numeric vector (SecB_UCell)
#' @param listw listw object
#' @param nsim Permutations
#' @param alpha FDR threshold
#' @return data.table with Li, p_adj, regime
local_bivariate <- function(x, y, listw, nsim = 999, alpha = 0.05) {
  n <- length(x)
  xz <- scale(x)[, 1]
  yz <- scale(y)[, 1]

  W <- spdep::listw2mat(listw)
  Wyz <- as.numeric(W %*% yz)

  Li_obs <- xz * Wyz

  Li_perm <- matrix(NA_real_, nrow = n, ncol = nsim)
  for (s in seq_len(nsim)) {
    perm_idx <- sample.int(n)
    yz_perm <- yz[perm_idx]
    Wyz_perm <- as.numeric(W %*% yz_perm)
    Li_perm[, s] <- xz * Wyz_perm
  }

  p_raw <- numeric(n)
  for (i in seq_len(n)) {
    p_raw[i] <- (sum(abs(Li_perm[i, ]) >= abs(Li_obs[i])) + 1) / (nsim + 1)
  }
  p_adj <- p.adjust(p_raw, method = "BH")

  regime <- ifelse(xz >= 0 & Wyz >= 0, "HH",
            ifelse(xz >= 0 & Wyz < 0,  "HL",
            ifelse(xz < 0  & Wyz >= 0, "LH",
                                        "LL")))
  regime_sig <- ifelse(p_adj < alpha, regime, "NS")

  data.table(
    Li      = Li_obs,
    p_raw   = p_raw,
    p_adj   = p_adj,
    regime  = regime,
    regime_sig = regime_sig
  )
}

# ============================================================================
# PART 1: Per-sample global statistics
# ============================================================================

message("\n", strrep("=", 70))
message("PART 1: Global Moran's I + Lee's bivariate L per sample")
message(strrep("=", 70))

global_results <- list()

for (sname in sfe_names) {
  t0 <- Sys.time()
  message("\n--- Processing ", sname, " ---")

  sfe <- load_sfe(sname)
  cd <- as.data.frame(colData(sfe))

  # Reconcile legacy SFE label with refactor naming (idempotent; no-op if absent)
  cd$cell_label[cd$cell_label == "Transitioning epithelium"] <- "Intermediate epithelium"

  is_sec <- cd$cell_label %in% sec_types
  n_sec <- sum(is_sec)
  message("  Secretory cells: ", format(n_sec, big.mark = ","))

  if (n_sec < MIN_SEC) {
    message("  SKIPPING: fewer than ", MIN_SEC, " secretory cells.")
    rm(sfe, cd); gc(verbose = FALSE)
    next
  }

  # Extract coordinates and scores for secretory cells
  xy_all <- spatialCoords(sfe)
  xy_sec <- xy_all[is_sec, , drop = FALSE]
  secA <- cd$SecA_UCell[is_sec]
  secB <- cd$SecB_UCell[is_sec]
  pol  <- cd$polarization_UCell[is_sec]

  rm(sfe); gc(verbose = FALSE)

  # For global statistics, subsample if very large
  if (n_sec > N_SUBSAMPLE) {
    message("  Subsampling to ", N_SUBSAMPLE, " cells for global stats")
    idx <- sort(sample.int(n_sec, N_SUBSAMPLE))
    xy_g <- xy_sec[idx, , drop = FALSE]
    secA_g <- secA[idx]
    secB_g <- secB[idx]
    pol_g  <- pol[idx]
  } else {
    xy_g <- xy_sec
    secA_g <- secA
    secB_g <- secB
    pol_g  <- pol
  }

  message("  Building spatial weights (r = ", RADIUS, " µm)...")
  listw <- build_weights(xy_g, RADIUS)

  # Global Moran's I for each score
  message("  Computing global Moran's I (SecA_UCell)...")
  gm_secA <- global_moran(secA_g, listw, nsim = NSIM)

  message("  Computing global Moran's I (SecB_UCell)...")
  gm_secB <- global_moran(secB_g, listw, nsim = NSIM)

  message("  Computing global Moran's I (polarization_UCell)...")
  gm_pol <- global_moran(pol_g, listw, nsim = NSIM)

  # Lee's bivariate L: SecA vs SecB
  message("  Computing Lee's bivariate L (SecA vs SecB)...")
  lee_res <- lee_bivariate(secA_g, secB_g, listw, nsim = NSIM)

  global_results[[sname]] <- data.table(
    sample        = sname,
    n_sec         = n_sec,
    n_used        = length(secA_g),
    moranI_SecA   = gm_secA$I,
    moranI_SecA_p = gm_secA$p_value,
    moranI_SecB   = gm_secB$I,
    moranI_SecB_p = gm_secB$p_value,
    moranI_pol    = gm_pol$I,
    moranI_pol_p  = gm_pol$p_value,
    lee_L         = lee_res$L,
    lee_L_p       = lee_res$p_value,
    lee_L_perm_sd = lee_res$L_perm_sd
  )

  message(sprintf("  Moran's I: SecA=%.3f (p=%s), SecB=%.3f (p=%s), pol=%.3f (p=%s)",
                  gm_secA$I, format(gm_secA$p_value, digits = 3),
                  gm_secB$I, format(gm_secB$p_value, digits = 3),
                  gm_pol$I, format(gm_pol$p_value, digits = 3)))
  message(sprintf("  Lee's L (SecA vs SecB) = %.4f (p=%s)",
                  lee_res$L, format(lee_res$p_value, digits = 3)))

  elapsed <- difftime(Sys.time(), t0, units = "mins")
  message(sprintf("  Done in %.1f min", as.numeric(elapsed)))

  rm(xy_all, xy_sec, secA, secB, pol, xy_g, secA_g, secB_g, pol_g, listw, cd)
  gc(verbose = FALSE)
}

# Combine and save
global_dt <- rbindlist(global_results)
fwrite(global_dt, file.path(out_path, "global_moran_summary.csv"))
message("\nSaved: global_moran_summary.csv (", nrow(global_dt), " samples)")

# ============================================================================
# PART 2: Local statistics (LISA + BiLISA) per sample
# ============================================================================

message("\n", strrep("=", 70))
message("PART 2: Local Moran's I (LISA) + Local Bivariate (BiLISA)")
message(strrep("=", 70))

lisa_summaries <- list()
bilisa_summaries <- list()

for (sname in sfe_names) {
  t0 <- Sys.time()
  message("\n--- LISA/BiLISA for ", sname, " ---")

  sfe <- load_sfe(sname)
  cd <- as.data.frame(colData(sfe))

  # Reconcile legacy SFE label with refactor naming (idempotent; no-op if absent)
  cd$cell_label[cd$cell_label == "Transitioning epithelium"] <- "Intermediate epithelium"

  is_sec <- cd$cell_label %in% sec_types
  n_sec <- sum(is_sec)

  if (n_sec < MIN_SEC) {
    message("  SKIPPING: fewer than ", MIN_SEC, " secretory cells.")
    rm(sfe, cd); gc(verbose = FALSE)
    next
  }

  xy_all <- spatialCoords(sfe)
  xy_sec <- xy_all[is_sec, , drop = FALSE]
  labels_sec <- cd$cell_label[is_sec]
  secA <- cd$SecA_UCell[is_sec]
  secB <- cd$SecB_UCell[is_sec]

  rm(sfe); gc(verbose = FALSE)

  # Subsample for local analysis (permutation-intensive)
  if (n_sec > N_SUBSAMPLE) {
    message("  Subsampling to ", N_SUBSAMPLE, " for local stats")
    idx <- sort(sample.int(n_sec, N_SUBSAMPLE))
  } else {
    idx <- seq_len(n_sec)
  }

  xy_l <- xy_sec[idx, , drop = FALSE]
  secA_l <- secA[idx]
  secB_l <- secB[idx]
  labels_l <- labels_sec[idx]

  message("  Building spatial weights...")
  listw <- build_weights(xy_l, RADIUS)

  # LISA for polarization (SecB - SecA)
  pol_l <- secB_l - secA_l
  message("  Computing LISA (polarization)...")
  lisa_pol <- local_moran_lisa(pol_l, listw, nsim = NSIM, alpha = FDR_ALPHA)
  lisa_pol[, `:=`(x = xy_l[, 1], y = xy_l[, 2], cell_label = labels_l,
                  sample = sname)]

  fwrite(lisa_pol, file.path(lisa_path, paste0("lisa_pol_", sname, ".csv")))

  lisa_tab <- lisa_pol[, .N, by = cluster]
  lisa_tab[, pct := round(100 * N / sum(N), 1)]
  lisa_summaries[[sname]] <- copy(lisa_tab)[, sample := sname]

  message("  LISA polarization clusters:")
  print(lisa_tab)

  # BiLISA: SecA vs SecB
  message("  Computing BiLISA (SecA vs SecB)...")
  bilisa <- local_bivariate(secA_l, secB_l, listw, nsim = NSIM,
                             alpha = FDR_ALPHA)
  bilisa[, `:=`(x = xy_l[, 1], y = xy_l[, 2], cell_label = labels_l,
                sample = sname)]

  fwrite(bilisa, file.path(bilisa_path, paste0("bilisa_", sname, ".csv")))

  bilisa_tab <- bilisa[, .N, by = regime_sig]
  bilisa_tab[, pct := round(100 * N / sum(N), 1)]
  bilisa_summaries[[sname]] <- copy(bilisa_tab)[, sample := sname]

  message("  BiLISA regimes:")
  print(bilisa_tab)

  elapsed <- difftime(Sys.time(), t0, units = "mins")
  message(sprintf("  Done in %.1f min", as.numeric(elapsed)))

  rm(xy_all, xy_sec, secA, secB, labels_sec, xy_l, secA_l, secB_l, pol_l,
     labels_l, listw, lisa_pol, bilisa, cd)
  gc(verbose = FALSE)
}

# Combine summaries
lisa_summary_dt <- rbindlist(lisa_summaries)
fwrite(lisa_summary_dt, file.path(out_path, "lisa_cluster_summary.csv"))

bilisa_summary_dt <- rbindlist(bilisa_summaries)
fwrite(bilisa_summary_dt, file.path(out_path, "bilisa_regime_summary.csv"))

# ============================================================================
# PART 3: Figures
# ============================================================================

message("\n", strrep("=", 70))
message("PART 3: Figures")
message(strrep("=", 70))

# --- 3a. Global Moran's I barplot -------------------------------------------

message("  [3a] Global Moran's I barplot")

gm_long <- melt(global_dt,
                 id.vars = c("sample", "n_sec", "n_used"),
                 measure.vars = c("moranI_SecA", "moranI_SecB", "moranI_pol"),
                 variable.name = "score", value.name = "I")
gm_long[, score := gsub("moranI_", "", score)]
gm_long[, score := factor(score, levels = c("SecA", "SecB", "pol"),
                           labels = c("SecA UCell", "SecB UCell", "Polarization"))]

# Add significance stars
p_long <- melt(global_dt,
               id.vars = "sample",
               measure.vars = c("moranI_SecA_p", "moranI_SecB_p", "moranI_pol_p"),
               variable.name = "score", value.name = "p")
p_long[, score := gsub("moranI_|_p", "", score)]
p_long[, score := factor(score, levels = c("SecA", "SecB", "pol"),
                          labels = c("SecA UCell", "SecB UCell", "Polarization"))]
gm_long <- merge(gm_long, p_long, by = c("sample", "score"))
gm_long[, sig := ifelse(p < 0.001, "***",
                 ifelse(p < 0.01, "**",
                 ifelse(p < 0.05, "*", "")))]

p_moran <- ggplot(gm_long, aes(x = reorder(sample, -I), y = I, fill = score)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = sig), position = position_dodge(0.8),
            vjust = -0.3, size = 2.5) +
  scale_fill_manual(values = c("SecA UCell" = "#E6A141",
                                "SecB UCell" = "#9A7D55",
                                "Polarization" = "#4A4A4A")) +
  labs(title = "Global Moran's I — Secretory UCell Scores",
       subtitle = paste0("Permutation test (n=", NSIM, "), ",
                          RADIUS, " µm radius"),
       x = NULL, y = "Moran's I", fill = NULL) +
  theme_lab() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

ggsave(file.path(fig_path, "global_moran_barplot.pdf"), p_moran,
       width = 8, height = 4.5)

# --- 3b. Lee's bivariate L barplot ------------------------------------------

message("  [3b] Lee's bivariate L barplot")

global_dt[, lee_sig := ifelse(lee_L_p < 0.001, "***",
                       ifelse(lee_L_p < 0.01, "**",
                       ifelse(lee_L_p < 0.05, "*", "")))]

p_lee <- ggplot(global_dt, aes(x = reorder(sample, -lee_L), y = lee_L)) +
  geom_col(fill = "#7570B3", width = 0.65) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text(aes(label = lee_sig), vjust = ifelse(global_dt$lee_L >= 0, -0.3, 1.3),
            size = 3) +
  labs(title = "Lee's Bivariate L — SecA vs SecB UCell Scores",
       subtitle = paste0("Spatial correlation between SecA_UCell and neighbor SecB_UCell\n",
                          NSIM, " permutations, ", RADIUS, " µm radius"),
       x = NULL,
       y = "Lee's L (negative = spatial segregation)") +
  theme_lab() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

ggsave(file.path(fig_path, "lee_bivariate_barplot.pdf"), p_lee,
       width = 7, height = 4.5)

# --- 3c. BiLISA spatial maps (whole tissue samples) -------------------------

message("  [3c] BiLISA spatial maps")

wt_samples <- grep("^sfe_OTB|^sfe_SP24", sfe_names, value = TRUE)
bilisa_plots <- list()

for (sname in wt_samples) {
  fpath <- file.path(bilisa_path, paste0("bilisa_", sname, ".csv"))
  if (!file.exists(fpath)) next

  bl <- fread(fpath)

  p <- ggplot(bl, aes(x = x, y = y, color = regime_sig)) +
    geom_point(shape = 16, size = 0.1, alpha = 0.8) +
    scale_color_manual(values = bilisa_pal, name = "BiLISA regime") +
    coord_fixed() +
    labs(title = gsub("sfe_", "", sname),
         subtitle = sprintf("SecA vs SecB BiLISA  (n=%s secretory)",
                            format(nrow(bl), big.mark = ","))) +
    theme_lab() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), axis.line = element_blank())

  bilisa_plots[[sname]] <- p

  ggsave(file.path(fig_path, paste0("bilisa_spatial_", sname, ".pdf")), p,
         width = 7, height = 6)
}

# --- 3d. LISA polarization spatial maps (whole tissue) ----------------------

message("  [3d] LISA polarization spatial maps")

for (sname in wt_samples) {
  fpath <- file.path(lisa_path, paste0("lisa_pol_", sname, ".csv"))
  if (!file.exists(fpath)) next

  ls <- fread(fpath)

  p <- ggplot(ls, aes(x = x, y = y, color = cluster)) +
    geom_point(shape = 16, size = 0.1, alpha = 0.8) +
    scale_color_manual(values = lisa_pal, name = "LISA cluster") +
    coord_fixed() +
    labs(title = gsub("sfe_", "", sname),
         subtitle = "Local Moran's I — Polarization UCell") +
    theme_lab() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), axis.line = element_blank())

  ggsave(file.path(fig_path, paste0("lisa_pol_spatial_", sname, ".pdf")), p,
         width = 7, height = 6)
}

# --- 3e. BiLISA regime composition vs discrete cell_label -------------------

message("  [3e] BiLISA vs discrete label agreement")

all_bilisa <- rbindlist(lapply(
  list.files(bilisa_path, pattern = "^bilisa_.*\\.csv$", full.names = TRUE),
  fread
))

if (nrow(all_bilisa) > 0) {
  cross_tab <- all_bilisa[regime_sig != "NS",
                           .N, by = .(cell_label, regime_sig)]
  cross_tab[, total := sum(N), by = cell_label]
  cross_tab[, pct := round(100 * N / total, 1)]

  p_cross <- ggplot(cross_tab,
                     aes(x = cell_label, y = pct, fill = regime_sig)) +
    geom_col(position = "stack", width = 0.7) +
    scale_fill_manual(values = bilisa_pal, name = "BiLISA regime") +
    labs(title = "BiLISA Regime Composition by Discrete Label",
         subtitle = "Significant cells only (FDR < 0.05)",
         x = NULL, y = "% of significant cells") +
    theme_lab() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(fig_path, "bilisa_composition_by_label.pdf"), p_cross,
         width = 6, height = 4.5)

  fwrite(cross_tab, file.path(out_path, "bilisa_vs_label_crosstab.csv"))
}

# --- 3f. Summary interpretation table ---------------------------------------

message("  [3f] Summary interpretation")

if (nrow(global_dt) > 0) {
  interpretation <- global_dt[, .(
    sample,
    n_sec,
    SecA_spatial_clustering = sprintf("I=%.3f (p=%s)", moranI_SecA,
                                      format(moranI_SecA_p, digits = 3)),
    SecB_spatial_clustering = sprintf("I=%.3f (p=%s)", moranI_SecB,
                                      format(moranI_SecB_p, digits = 3)),
    SecA_SecB_segregation   = sprintf("L=%.4f (p=%s)", lee_L,
                                      format(lee_L_p, digits = 3)),
    interpretation = fifelse(
      lee_L < 0 & lee_L_p < 0.05,
      "SecA and SecB are spatially SEGREGATED",
      fifelse(
        lee_L > 0 & lee_L_p < 0.05,
        "SecA and SecB spatially CO-LOCALIZE",
        "No significant bivariate spatial association"
      )
    )
  )]

  fwrite(interpretation, file.path(out_path, "interpretation_summary.csv"))
  message("\n  Interpretation summary:")
  print(interpretation[, .(sample, SecA_SecB_segregation, interpretation)])
}

# ============================================================================
# PART 4: Bivariate Lee's L specifically within TMA cores
# ============================================================================

message("\n", strrep("=", 70))
message("PART 4: Per-core Lee's L (TMA)")
message(strrep("=", 70))

tma_checkpoint <- file.path(out_path, "tma_per_core_lee.csv")

if (file.exists(tma_checkpoint)) {
  message("  Loading existing per-core results")
  core_lee_dt <- fread(tma_checkpoint)
} else {
  sfe <- load_sfe("sfe_tma_filtered")
  cd <- as.data.frame(colData(sfe))
  xy <- spatialCoords(sfe)

  # Reconcile legacy SFE label with refactor naming (idempotent; no-op if absent)
  cd$cell_label[cd$cell_label == "Transitioning epithelium"] <- "Intermediate epithelium"

  is_sec <- cd$cell_label %in% sec_types
  core_ids <- unique(cd$core_id[is_sec & !is.na(cd$core_id)])
  message("  Processing ", length(core_ids), " TMA cores")

  core_lee_list <- list()

  for (cid in core_ids) {
    idx <- which(cd$core_id == cid & is_sec)

    if (length(idx) < 50) next

    xy_c <- xy[idx, , drop = FALSE]
    secA_c <- cd$SecA_UCell[idx]
    secB_c <- cd$SecB_UCell[idx]

    listw_c <- tryCatch(build_weights(xy_c, RADIUS), error = function(e) NULL)
    if (is.null(listw_c)) next

    lee_c <- tryCatch(
      lee_bivariate(secA_c, secB_c, listw_c, nsim = 499),
      error = function(e) NULL
    )
    if (is.null(lee_c)) next

    core_lee_list[[cid]] <- data.table(
      core_id    = cid,
      patient_id = cd$patient_id[idx[1]],
      n_sec      = length(idx),
      lee_L      = lee_c$L,
      lee_L_p    = lee_c$p_value
    )
  }

  core_lee_dt <- rbindlist(core_lee_list)
  fwrite(core_lee_dt, tma_checkpoint)

  rm(sfe, cd, xy); gc(verbose = FALSE)
  message("  Computed Lee's L for ", nrow(core_lee_dt), " cores")
}

# Per-core distribution figure
if (nrow(core_lee_dt) > 0) {
  core_lee_dt[, sig := lee_L_p < 0.05]

  p_core <- ggplot(core_lee_dt, aes(x = lee_L, fill = sig)) +
    geom_histogram(bins = 40, color = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    scale_fill_manual(values = c("TRUE" = "#7570B3", "FALSE" = "#D9D9D9"),
                      labels = c("TRUE" = "p < 0.05", "FALSE" = "NS"),
                      name = NULL) +
    labs(title = "Per-Core Lee's L (SecA vs SecB UCell)",
         subtitle = sprintf("TMA cores (n=%d), %d µm radius, %d permutations",
                            nrow(core_lee_dt), RADIUS, 499),
         x = "Lee's L (negative = spatial segregation)",
         y = "Number of cores") +
    theme_lab()

  ggsave(file.path(fig_path, "tma_per_core_lee_histogram.pdf"), p_core,
         width = 6, height = 4)

  # Patient-level summary (mean across replicate cores)
  patient_lee <- core_lee_dt[, .(
    mean_lee_L = mean(lee_L),
    n_cores    = .N,
    n_sig      = sum(sig)
  ), by = patient_id]

  fwrite(patient_lee, file.path(out_path, "tma_patient_level_lee.csv"))
  message("  Patient-level Lee's L: ",
          sprintf("mean=%.4f, %d/%d patients with ≥1 sig core",
                  mean(patient_lee$mean_lee_L),
                  sum(patient_lee$n_sig > 0), nrow(patient_lee)))
}

# ============================================================================
# FINISH
# ============================================================================

message("\n", strrep("=", 70))
message("Script 44 complete.")
message(strrep("=", 70))

message("\nOutput directory: ", out_path)
message("Figure directory: ", fig_path)
message("\nKey outputs:")
message("  global_moran_summary.csv      — Global I for SecA/SecB/polarization per sample")
message("  bivariate_lee_summary.csv     — Lee's L (SecA vs SecB) global per sample")
message("  lisa_results/*.csv            — Per-cell LISA clusters (polarization)")
message("  bilisa_results/*.csv          — Per-cell BiLISA regimes (SecA vs SecB)")
message("  interpretation_summary.csv    — Plain-language per-sample interpretation")
message("  tma_per_core_lee.csv          — Per-core Lee's L for TMA")
message("  tma_patient_level_lee.csv     — Patient-averaged Lee's L")
message("\nInterpretation:")
message("  Moran's I > 0: spatial clustering of that score (high near high)")
message("  Lee's L < 0: SecA and SecB are spatially SEGREGATED (high-SecA neighborhoods")
message("    have low-SecB neighbors and vice versa)")
message("  Lee's L > 0: SecA and SecB spatially co-localize (transitioning zones)")
message("  BiLISA HH: Transitioning zones (high both), HL/LH: subtype-pure regions")

log_session()
