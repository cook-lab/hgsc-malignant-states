#!/usr/bin/env Rscript
# PURPOSE: Fig 4f per-cell GAM inputs for all SingleR secretory cells (SecA / Intermediate / SecB) of the 8 whole
#   tissues: distance to the nearest pericyte or endothelial cell; number of cells within 50 um per mm^2 of the
#   circle (cell_density; cell_density_epi counts secretory neighbours only); mean pathway_hypoxia and
#   pathway_matrix_remodeling of the cells within 50 um.
# INPUTS: SFE pathway_hypoxia / pathway_matrix_remodeling and spatialCoords (8 WT);
#   singler_recompute/<sample>_singler_cells.csv.gz (cell_label_singler, polarization_UCell).
# OUTPUTS: singler_recompute/singler_gam_percell.csv (1,235,350 secretory cells; the column cell_label_split holds
#   the SingleR label). Read by 09_fig4f_gam_fit.R and by figures/supplementary/21_figS16-18_per_sample_gams.R.
# MANUSCRIPT PANEL(S): Fig 4f; Figure S16. RUNTIME TIER: moderate.
source("spatial/00_setup/00_setup.R")
suppressPackageStartupMessages({
  library(SpatialFeatureExperiment); library(SummarizedExperiment)
  library(data.table); library(dbscan); library(RANN)
})
R <- cfg_path("output_root")
SGL <- file.path(R,"singler_recompute")
OUT <- SGL
RADIUS <- 50
CORE8 <- paste0("sfe_", c("WT01","WT02","WT03","WT04",
                          "WT05","WT06","WT07","WT08"))
VASC <- c("Pericyte","Endothelial")
EPI  <- c("SecA epithelium","Intermediate epithelium","SecB epithelium")

res <- list()
for (nm in CORE8) {
  t0 <- Sys.time(); message("[", format(t0), "] ", nm)
  s <- sub("^sfe_","",nm)
  sgl <- fread(file.path(SGL, paste0(s,"_singler_cells.csv.gz")),
               select=c("barcode","cell_label_singler","polarization_UCell"))
  sfe <- load_sfe(nm)
  cd  <- as.data.table(as.data.frame(colData(sfe)[, c("pathway_hypoxia","pathway_matrix_remodeling")]))
  xy  <- spatialCoords(sfe); bc <- colnames(sfe); rm(sfe); gc()
  lab <- sgl$cell_label_singler[match(bc, sgl$barcode)]
  polv<- sgl$polarization_UCell[match(bc, sgl$barcode)]
  keep <- !is.na(lab) & lab!=""
  dt <- data.table(barcode=bc, lab=lab, x=xy[,1], y=xy[,2],
                   hyp=cd$pathway_hypoxia, mat=cd$pathway_matrix_remodeling, pol=polv)[keep]
  message("   SingleR-labeled ", nrow(dt))
  iv <- which(dt$lab %in% VASC)
  dt[, dist_to_vascular := NA_real_]
  sec <- which(dt$lab %in% EPI)
  if (length(iv) > 0)
    dt[sec, dist_to_vascular := RANN::nn2(as.matrix(dt[iv, .(x,y)]),
                                          as.matrix(dt[sec, .(x,y)]), k=1)$nn.dists[,1]]
  fr <- dbscan::frNN(as.matrix(dt[, .(x,y)]), eps = RADIUS)
  is_epi <- dt$lab %in% EPI
  hyp <- dt$hyp; mat <- dt$mat
  AREA <- pi * (RADIUS/1000)^2
  nb_hyp <- nb_mat <- rep(NA_real_, nrow(dt))
  dens_all <- dens_epi <- rep(0L, nrow(dt))
  for (i in sec) {
    id <- fr$id[[i]]
    dens_all[i] <- length(id); dens_epi[i] <- sum(is_epi[id])
    if (length(id)) { nb_hyp[i] <- mean(hyp[id], na.rm=TRUE); nb_mat[i] <- mean(mat[id], na.rm=TRUE) }
  }
  out <- dt[sec, .(sample=nm, barcode, cell_label_split=lab, polarization=pol, dist_to_vascular)]
  out[, `:=`(cell_density = dens_all[sec]/AREA, cell_density_epi = dens_epi[sec]/AREA,
             nb_hypoxia = nb_hyp[sec], nb_matrix = nb_mat[sec])]
  res[[nm]] <- out
  message(sprintf("   secretory %d  done in %.1f min", nrow(out),
                  as.numeric(difftime(Sys.time(), t0, units="mins"))))
  rm(dt, fr); gc()
}
all <- rbindlist(res)
setcolorder(all, c("sample","barcode","polarization","dist_to_vascular","nb_hypoxia","cell_density","nb_matrix"))
fwrite(all, file.path(OUT, "singler_gam_percell.csv"))
message("\nrows: ", nrow(all), "  -> singler_gam_percell.csv"); message("DONE")
