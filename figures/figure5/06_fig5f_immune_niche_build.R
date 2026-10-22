#!/usr/bin/env Rscript
# ============================================================================
# Fig 5f -- immune composition of the local niche vs local epithelial hypoxia (x_hyp)
# ----------------------------------------------------------------------------
#   * cell identity = SingleR cell_label_singler (via spatial/00_setup/singler_labels.R)
#   * local epithelial hypoxia uses the SFE `pathway_hypoxia` column written by
#     spatial/04_pathway_scoring/04c_pathway_scoring_allcells.R (UCell, full
#     477-gene panel, maxRank=477, all cells)
#   * coordinates are spatialCoords (cellSeg is incomplete for sfe_tma)
#   * 50um dnearneigh, BINARY weights, x = mean hypoxia over EPITHELIAL
#     neighbours (sum over epithelial neighbours / count of them), z-scored within
#     sample; rows are the non-epithelial (TME) cells.
#   * x_pol (mean polarization of the epithelial neighbours) is also written; Fig 5f
#     uses x_hyp only (07_fig5f_immune_niche_fit.R with FIG5C_AXIS=hyp).
# Rows: the 8 whole tissues and the 3 TMA arrays. The panel shows whole tissues only,
# but the TMA rows enter the 1st-99th percentile grid of x in 07, so they are kept.
# OUTPUT: output_root/split_recompute/fig5c_tme_percell_presplit.csv (folder and file
#   names are historical); read by 07 and 11_fig5_compose_panels.py.
# MANUSCRIPT PANEL(S): Fig 5f. RUNTIME TIER: heavy (loads 11 SFEs).
# ============================================================================
source("spatial/00_setup/00_setup.R")
source("spatial/00_setup/singler_labels.R")   # presplit_labels() = singler_labels(): SingleR labels
suppressWarnings(suppressMessages({
  library(SpatialFeatureExperiment); library(SummarizedExperiment)
  library(sf); library(spdep); library(data.table)
}))
set.seed(CFG$seed)
R <- cfg_path("output_root"); OUT <- cfg_path("output_root","split_recompute")
RADIUS <- 50
WT  <- paste0("sfe_", CFG$cohort$whole_tissue)
TMA <- c("sfe_tma","sfe_COEUR_TMA_G","sfe_COEUR_TMA_H")
zt <- function(v){ s<-sd(v,na.rm=TRUE); if(!is.finite(s)||s==0) return(v*NA_real_); (v-mean(v,na.rm=TRUE))/s }

res <- list()
for (nm in c(WT, TMA)) {
  t0 <- Sys.time(); message("[", format(t0), "] ", nm)
  sfe <- load_sfe(nm)
  L   <- presplit_labels(sfe, nm)
  bc  <- colnames(sfe)
  lab <- L$lab
  hyp <- as.numeric(colData(sfe)[["pathway_hypoxia"]])
  polv <- L$polarization
  # spatialCoords, NOT cellSeg centroids: sfe_tma carries only 309,797 cellSeg
  # polygons for 671,453 cells, so a geometry-based frame is incomplete there.
  # Neighbour finding is invariant to the spatialCoords/cellSeg y-mirroring.
  xy  <- spatialCoords(sfe)
  rm(sfe); gc(verbose=FALSE)

  keep <- !is.na(lab) & nzchar(lab)
  lab <- lab[keep]; hyp <- hyp[keep]; polv <- polv[keep]; xy <- xy[keep, , drop=FALSE]
  isEpi <- grepl("epithelium", lab, ignore.case=TRUE)
  message("   SPLIT-kept ", length(lab), " ; epithelial ", sum(isEpi))

  nb <- dnearneigh(xy, 0, RADIUS); lw <- nb2listw(nb, style="B", zero.policy=TRUE)
  num <- lag.listw(lw, ifelse(isEpi & is.finite(hyp), hyp, 0), zero.policy=TRUE)
  den <- lag.listw(lw, as.numeric(isEpi), zero.policy=TRUE)
  nb_hyp <- ifelse(den > 0, num/den, NA_real_)
  # local epithelial polarization (x_pol), same definition as 08_fig5i: mean
  # polarization_UCell over EPITHELIAL neighbours (not used by Fig 5f).
  npol <- lag.listw(lw, ifelse(isEpi & is.finite(polv), polv, 0), zero.policy=TRUE)
  nb_pol <- ifelse(den > 0, npol/den, NA_real_)

  dt <- data.table(sample=nm, cohort=ifelse(nm %in% WT, "WT", "TMA"),
                   cell_label=lab, x_hyp=zt(nb_hyp), x_pol=zt(nb_pol))
  dt <- dt[!isEpi & !is.na(cell_label) & (is.finite(x_hyp) | is.finite(x_pol))]
  res[[nm]] <- dt
  message(sprintf("   TME cells %d  done in %.1f min", nrow(dt),
                  as.numeric(difftime(Sys.time(), t0, units="mins"))))
  rm(nb, lw); gc(verbose=FALSE)
}
all <- rbindlist(res)
fwrite(all, file.path(OUT, "fig5c_tme_percell_presplit.csv"))
message("\nrows: ", nrow(all), "  samples: ", uniqueN(all$sample))
print(all[, .N, by=.(cohort)][order(cohort)])
message("DONE")
