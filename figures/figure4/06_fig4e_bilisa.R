#!/usr/bin/env Rscript
# ============================================================================
# Fig 4e -- BiLISA neighbourhoods of secretory cells in the 8 whole tissues (SingleR)
# ----------------------------------------------------------------------------
# PURPOSE
#   Bivariate local indicator of spatial association (BiLISA) between each secretory
#   cell's SecA score and its neighbours' SecB score; regimes HH / HL / LH / LL, kept
#   where the permutation p is significant after BH within tissue; then the % of each
#   secretory state in each regime (the 4e bars).
#     cells   SingleR cell_label_singler (SecA / Intermediate / SecB)
#     scores  SecA_UCell / SecB_UCell from colData (UCell on raw counts)
#     weights binary 50 um dnearneigh (style "B"), 499 permutations, BH alpha 0.05
#   The whole-tissue loop also computes a global row-standardised bivariate statistic
#   (lee_bivariate(), "Lee's L" in earlier code). No panel uses it; see the note at the
#   function for why it is still computed.
#
# INPUTS
#   SFE objects of the 8 whole tissues (SecA_UCell / SecB_UCell in colData; spatialCoords)
#   SingleR labels via spatial/00_setup/singler_labels.R
#
# OUTPUTS (output_root/singler_recompute/)
#   bilisa_percell_SINGLER.csv.gz   sample, cell_label, regime, regime_gated (per-cell record)
#   bilisa_crosstab_SINGLER.csv     cell_label, regime_sig, N, total, pct (read by 07_fig4e_bilisa_plot.py)
#
# MANUSCRIPT PANEL(S): Fig 4e (plot: 07_fig4e_bilisa_plot.py); "89% of SecB cells" in the Results
# RUNTIME TIER: medium (permutations on all secretory cells of each tissue)
# ============================================================================
source("spatial/00_setup/00_setup.R")
source("spatial/00_setup/singler_labels.R")   # singler_labels(): SingleR labels (cell_label_singler)
suppressWarnings(suppressMessages({
  library(SpatialFeatureExperiment); library(SummarizedExperiment)
  library(spdep); library(data.table)}))
set.seed(CFG$seed)
R <- cfg_path("output_root"); OUT <- file.path(R,"singler_recompute")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
RADIUS <- 50; NSIM_LEE <- 999; NSIM_BIL <- 499; ALPHA <- 0.05
MIN_SEC <- 200; MIN_CORE <- 500
EPI <- c("SecA epithelium","Intermediate epithelium","SecB epithelium")
WT  <- paste0("sfe_", CFG$cohort$whole_tissue)

# Global bivariate statistic with row-standardised weights, computed sparsely:
# lag.listw() gives W %*% yc without forming W, so it runs on all secretory cells.
# S0 = total weight = sum of the per-row weight sums (row-standardised, so this
# equals the number of non-isolated cells under zero.policy).
# NOT USED BY ANY PANEL. It is still computed in the whole-tissue loop because its
# 999 permutations draw from the random-number stream (set.seed above) before each
# tissue's BiLISA permutations; removing the call would change the BiLISA p-values
# and so the Fig 4e regimes. Its result is not written.
lee_bivariate <- function(x, y, listw, nsim = NSIM_LEE) {
  n <- length(x); xc <- x-mean(x); yc <- y-mean(y)
  S0 <- sum(unlist(listw$weights))
  sxx <- sqrt(sum(xc^2))
  Wy <- spdep::lag.listw(listw, yc, zero.policy=TRUE)
  L_obs <- (n/S0)*sum(xc*Wy)/(sxx*sqrt(sum(yc^2)))
  Lp <- numeric(nsim)
  for (i in seq_len(nsim)) {
    ycp <- yc[sample.int(n)]
    Wyp <- spdep::lag.listw(listw, ycp, zero.policy=TRUE)
    Lp[i] <- (n/S0)*sum(xc*Wyp)/(sxx*sqrt(sum(ycp^2)))
  }
  list(L=L_obs, p=(sum(abs(Lp)>=abs(L_obs))+1)/(nsim+1))
}
local_bivariate <- function(x, y, listw, nsim=NSIM_BIL, alpha=ALPHA) {
  n <- length(x); xz <- (x-mean(x))/sd(x); yz <- (y-mean(y))/sd(y)
  Wyz <- spdep::lag.listw(listw, yz, zero.policy=TRUE); Li <- xz*Wyz; aLi <- abs(Li)
  ge <- integer(n)
  for (s in seq_len(nsim)) {
    p <- sample.int(n)
    ge <- ge + (abs(xz*spdep::lag.listw(listw, yz[p], zero.policy=TRUE)) >= aLi)
  }
  p_adj <- p.adjust((ge+1)/(nsim+1), "BH")
  regime <- ifelse(xz>=0 & Wyz>=0,"HH", ifelse(xz>=0 & Wyz<0,"HL",
                   ifelse(xz<0 & Wyz>=0,"LH","LL")))
  data.table(regime=regime, regime_gated=ifelse(p_adj<alpha, regime, "NS"))
}
get_cells <- function(nm) {
  sfe <- load_sfe(nm)
  L <- singler_labels(sfe, nm)
  cd <- data.table(barcode=colnames(sfe), lab=L$lab,
                   secA=as.numeric(colData(sfe)$SecA_UCell),
                   secB=as.numeric(colData(sfe)$SecB_UCell),
                   x=spatialCoords(sfe)[,1], y=spatialCoords(sfe)[,2])
  ex <- intersect(c("core_id","sample_type"), colnames(colData(sfe)))
  if (length(ex)) cd <- cbind(cd, as.data.table(as.data.frame(colData(sfe)[, ex, drop=FALSE])))
  rm(sfe); gc(verbose=FALSE)
  cd[lab %in% EPI & is.finite(secA) & is.finite(secB)]
}
wt_lee <- list(); bil <- list()
for (nm in WT) {
  message("[WT] ", nm)
  d <- get_cells(nm)
  if (nrow(d) < MIN_SEC) { message("   <",MIN_SEC," secretory -- skip"); next }
  xy <- as.matrix(d[, .(x,y)])
  # row-standardised global statistic: not used, kept for the random-number stream (see lee_bivariate)
  lw_w <- nb2listw(dnearneigh(xy,0,RADIUS), style="W", zero.policy=TRUE)
  r <- lee_bivariate(d$secA, d$secB, lw_w)
  wt_lee[[nm]] <- data.table(sample=nm, lee_L=r$L, p=r$p, n=nrow(d))
  message(sprintf("   n=%d  Lee L=%.4f  p=%.4f", nrow(d), r$L, r$p))
  lw_b <- nb2listw(dnearneigh(xy,0,RADIUS), style="B", zero.policy=TRUE)
  b <- local_bivariate(d$secA, d$secB, lw_b)
  bil[[nm]] <- data.table(sample=nm, cell_label=d$lab, regime=b$regime, regime_gated=b$regime_gated)
  rm(d, lw_w, lw_b); gc(verbose=FALSE)
}
allb <- rbindlist(bil)
fwrite(allb, file.path(OUT,"bilisa_percell_SINGLER.csv.gz"))
sig <- allb[regime_gated!="NS"]
g <- sig[, .N, by=.(cell_label, regime_gated)]
tot <- sig[, .(total=.N), by=cell_label]
g <- merge(g, tot, by="cell_label"); g[, pct := round(100*N/total,1)]
setnames(g,"regime_gated","regime_sig"); setcolorder(g, c("cell_label","regime_sig","N","total","pct"))
fwrite(g[order(cell_label, regime_sig)], file.path(OUT,"bilisa_crosstab_SINGLER.csv"))
cat("\n=== BiLISA (SingleR), % of BH-significant cells ===\n")
print(dcast(g, cell_label ~ regime_sig, value.var="pct", fill=0))
message("DONE")
