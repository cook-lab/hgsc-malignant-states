#!/usr/bin/env Rscript
# ============================================================================
# Fig 4d -- global bivariate Moran's I between SecA and SecB scores (SingleR)
# ----------------------------------------------------------------------------
# PURPOSE
#   One number per whole tissue (and per TMA core -> patient) for the spatial association
#   between a secretory cell's SecA score and its neighbours' SecB score: Wartenberg's global
#   bivariate Moran's I. Fig 4e (BiLISA) is the local form of the same statistic, so this run
#   uses the SAME weights as 4e: binary 50 um fixed-radius neighbours (style "B").
#
#   I = (n / S0) * sum_i sum_j w_ij xc_i yc_j / (|xc| |yc|),  x = SecA_UCell, y = SecB_UCell,
#   w_ij = 1 for secretory cells j within 50 um of i (j != i); S0 = sum of weights.
#   Per tissue: conditional permutation (SecB scores shuffled across cells, SecA fixed; as the
#   BiLISA test), 999 permutations; two-sided. With 10^4-10^5 cells this p is at its floor, so the
#   inference is across units: exact sign test over the 8 tissues; one-sample Wilcoxon (and a sign
#   test) over TMA patients. Within-cell SecA-SecB Pearson r is reported alongside.
#
# INPUTS
#   SFE objects (cfg sfe dir; SecA_UCell / SecB_UCell in colData; spatialCoords)
#   SingleR labels via spatial/00_setup/singler_labels.R (cell_label_singler)
#   output_root/singler_recompute/<tma>_singler_cells.csv.gz  (patient_id, core_id, sample_type)
#   the final TMA core set (spatial/00_setup/final_tma_cores.R), cores with >= 500 secretory cells
#
# OUTPUTS (output_root/singler_recompute/)
#   bivmoran_wt_SINGLER.csv          tissue, n_cells, mean_neighbours, I, p_perm, within_cell_r
#   bivmoran_tma_percore_SINGLER.csv tma, core, patient_unit, n_cells, I, p_perm (per-core record of the TMA arm)
#   bivmoran_tma_patient_SINGLER.csv patient_unit, n_cores, mean_I, n_sig_cores
#   bivmoran_tests_SINGLER.csv       cohort, n_units, n_negative, median_I, sign-test p, Wilcoxon p
#
# MANUSCRIPT PANEL(S): Fig 4d (plot: 05_fig4d_bivariate_moran_plot.py)
# RUNTIME TIER: medium (~20-40 min; 999 permutations on up to 3e5 cells per tissue)
# ============================================================================
source("spatial/00_setup/00_setup.R")
source("spatial/00_setup/singler_labels.R")   # singler_labels(): SingleR labels (cell_label_singler)
source("spatial/00_setup/final_tma_cores.R")
suppressWarnings(suppressMessages({
  library(SpatialFeatureExperiment); library(SummarizedExperiment)
  library(spdep); library(data.table)}))
set.seed(CFG$seed)
R <- cfg_path("output_root"); OUT <- file.path(R, "singler_recompute")
TMA_LABEL <- c(sfe_tma = "A", sfe_COEUR_TMA_G = "B", sfe_COEUR_TMA_H = "B")
RADIUS <- 50; NSIM_WT <- 999; NSIM_TMA <- 199; MIN_SEC <- 200; MIN_CORE <- 500
EPI <- c("SecA epithelium", "Intermediate epithelium", "SecB epithelium")
WT  <- paste0("sfe_", CFG$cohort$whole_tissue)
TMA <- names(TMA_LABEL)

# Global bivariate Moran's I, computed sparsely (lag.listw gives W %*% y without forming W).
bivariate_moran <- function(x, y, listw, nsim) {
  n <- length(x); xc <- x - mean(x); yc <- y - mean(y)
  S0 <- sum(unlist(listw$weights)); sx <- sqrt(sum(xc^2)); sy <- sqrt(sum(yc^2))
  I_obs <- (n / S0) * sum(xc * spdep::lag.listw(listw, yc, zero.policy = TRUE)) / (sx * sy)
  Ip <- vapply(seq_len(nsim), function(i) {
    ycp <- yc[sample.int(n)]
    (n / S0) * sum(xc * spdep::lag.listw(listw, ycp, zero.policy = TRUE)) / (sx * sy)
  }, numeric(1))
  list(I = I_obs, p = (sum(abs(Ip) >= abs(I_obs)) + 1) / (nsim + 1), S0 = S0)
}
get_cells <- function(nm) {
  sfe <- load_sfe(nm)
  L <- singler_labels(sfe, nm)
  cd <- data.table(barcode = colnames(sfe), lab = L$lab,
                   secA = as.numeric(colData(sfe)$SecA_UCell), secB = as.numeric(colData(sfe)$SecB_UCell),
                   x = spatialCoords(sfe)[, 1], y = spatialCoords(sfe)[, 2])
  rm(sfe); gc(verbose = FALSE)
  cd[lab %in% EPI & is.finite(secA) & is.finite(secB)]
}
binary_w <- function(xy) nb2listw(dnearneigh(xy, 0, RADIUS), style = "B", zero.policy = TRUE)

# ---- whole tissues ----
wt <- list()
for (nm in WT) {
  t0 <- Sys.time(); d <- get_cells(nm)
  if (nrow(d) < MIN_SEC) { message(nm, ": <", MIN_SEC, " secretory cells, skipped"); next }
  lw <- binary_w(as.matrix(d[, .(x, y)]))
  r <- bivariate_moran(d$secA, d$secB, lw, NSIM_WT)
  wt[[nm]] <- data.table(tissue = sub("^sfe_", "", nm), n_cells = nrow(d), mean_neighbours = r$S0 / nrow(d),
                         I = r$I, p_perm = r$p, within_cell_r = cor(d$secA, d$secB))
  message(sprintf("[WT] %s n=%d I=%.4f p=%.4f r=%.3f (%.1f min)", nm, nrow(d), r$I, r$p,
                  cor(d$secA, d$secB), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(d, lw); gc(verbose = FALSE)
}
wt <- rbindlist(wt)
fwrite(wt, file.path(OUT, "bivmoran_wt_SINGLER.csv"))

# ---- TMA cores -> patients (patients keyed by TMA: CHUM and COEUR reuse numbers) ----
tc <- list()
for (nm in TMA) {
  d <- get_cells(nm)
  s <- fread(file.path(OUT, paste0(sub("^sfe_", "", nm), "_singler_cells.csv.gz")),
             select = c("barcode", "patient_id", "core_id", "sample_type"),
             colClasses = list(character = c("patient_id", "core_id")))
  d[s, on = "barcode", `:=`(pat = i.patient_id, core = i.core_id, st = i.sample_type)]
  keep <- final_core_ids(nm)                       # final TMA cores (tumour; QC; HGSC; treatment-naive)
  d <- d[!is.na(pat) & pat != "" & !is.na(core) & core %in% keep]
  for (cid in unique(d$core)) {
    cc <- d[core == cid]
    if (nrow(cc) < MIN_CORE) next
    lw <- tryCatch(binary_w(as.matrix(cc[, .(x, y)])), error = function(e) NULL)
    if (is.null(lw)) next
    r <- tryCatch(bivariate_moran(cc$secA, cc$secB, lw, NSIM_TMA), error = function(e) NULL)
    if (is.null(r)) next
    tc[[length(tc) + 1]] <- data.table(tma = sub("^sfe_", "", nm), core = cid,
                                       patient_unit = paste0(TMA_LABEL[[nm]], "|", cc$pat[1]),
                                       n_cells = nrow(cc), I = r$I, p_perm = r$p)
  }
  message("[TMA] ", nm, ": cores scored ", sum(vapply(tc, function(z) z$tma == sub("^sfe_", "", nm), logical(1))))
  rm(d); gc(verbose = FALSE)
}
tc <- rbindlist(tc)
fwrite(tc, file.path(OUT, "bivmoran_tma_percore_SINGLER.csv"))
pt <- tc[, .(n_cores = .N, mean_I = mean(I), n_sig_cores = sum(p_perm < 0.05)), by = patient_unit]
fwrite(pt, file.path(OUT, "bivmoran_tma_patient_SINGLER.csv"))

# ---- across-unit tests ----
sign_p <- function(v) binom.test(sum(v < 0), length(v), 0.5)$p.value
tests <- rbind(
  data.table(cohort = "whole tissue", unit = "tissue", n_units = nrow(wt), n_negative = sum(wt$I < 0),
             median_I = median(wt$I), min_I = min(wt$I), max_I = max(wt$I),
             sign_test_p = sign_p(wt$I), wilcoxon_p = wilcox.test(wt$I, mu = 0, exact = TRUE)$p.value),
  data.table(cohort = "TMA", unit = "patient (mean over cores)", n_units = nrow(pt), n_negative = sum(pt$mean_I < 0),
             median_I = median(pt$mean_I), min_I = min(pt$mean_I), max_I = max(pt$mean_I),
             sign_test_p = sign_p(pt$mean_I), wilcoxon_p = wilcox.test(pt$mean_I, mu = 0)$p.value))
fwrite(tests, file.path(OUT, "bivmoran_tests_SINGLER.csv"))
print(wt); print(tests)
message("DONE")
