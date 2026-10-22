#!/usr/bin/env Rscript
# ============================================================================
# Fig 5i -- MACROPHAGE-INTRINSIC programs vs the local epithelial niche
# ----------------------------------------------------------------------------
# Every score is computed on the macrophage ITSELF (SingleR cell_label_singler ==
# "Macrophage", via spatial/00_setup/singler_labels.R), using UCell on the FULL
# 477-gene panel, maxRank = 477 -- the same settings as the SFE pathway_* columns.
# Gene sets are exactly those named in the results text:
#   Glycolytic         LDHA, ENO1, PGK1, ALDOA
#   NF-kB activated    ICAM1, NFKBIA, BIRC3
#   Tissue-resident    A2M, C1QC
#   Chemokine          CXCL10, CXCL11
# Two x-axes are computed (both 50um, binary weights, over EPITHELIAL neighbours,
# z-scored within sample):
#   nb_pol  local epithelial polarization (the Fig 5i and Figure S18 axis)
#   nb_hyp  local epithelial hypoxia      (fitted by 09 as well; not shown)
# OUTPUT: output_root/split_recompute/fig5e_macrophage_percell_presplit.csv (folder and
#   file names are historical); read by 09, 11_fig5_compose_panels.py and
#   figures/supplementary/21_figS16-18_per_sample_gams.R.
# MANUSCRIPT PANEL(S): Fig 5i; Figure S18. RUNTIME TIER: heavy (loads the 8 whole-tissue SFEs).
# ============================================================================
source("spatial/00_setup/00_setup.R")
source("spatial/00_setup/singler_labels.R")   # presplit_labels() = singler_labels(): SingleR labels
suppressWarnings(suppressMessages({
  library(SpatialFeatureExperiment); library(SummarizedExperiment)
  library(spdep); library(data.table); library(Seurat); library(UCell)
}))
set.seed(CFG$seed)
R <- cfg_path("output_root"); OUT <- cfg_path("output_root","split_recompute")
RADIUS <- 50
SETS <- list(glycolytic      = c("LDHA","ENO1","PGK1","ALDOA"),
             nfkb_activated  = c("ICAM1","NFKBIA","BIRC3"),
             tissue_resident = c("A2M","C1QC"),
             chemokine       = c("CXCL10","CXCL11"))
WT <- paste0("sfe_", CFG$cohort$whole_tissue)
zt <- function(v){ s<-sd(v,na.rm=TRUE); if(!is.finite(s)||s==0) return(v*NA_real_); (v-mean(v,na.rm=TRUE))/s }

res <- list()
for (nm in WT) {
  t0 <- Sys.time(); message("[", format(t0), "] ", nm)
  sfe <- load_sfe(nm)
  L <- presplit_labels(sfe, nm); bc <- colnames(sfe); lab <- L$lab
  hyp <- as.numeric(colData(sfe)[["pathway_hypoxia"]])
  pgly <- as.numeric(colData(sfe)[["pathway_glycolysis"]])
  pnfkb <- as.numeric(colData(sfe)[["pathway_nfkb"]])
  polv <- L$polarization
  xy <- spatialCoords(sfe)
  keep <- !is.na(lab) & nzchar(lab)
  M <- counts(sfe)[, keep, drop=FALSE]; M <- as(M,"CsparseMatrix")
  rm(sfe); gc(verbose=FALSE)
  lab<-lab[keep]; hyp<-hyp[keep]; polv<-polv[keep]; xy<-xy[keep,,drop=FALSE]
  pgly<-pgly[keep]; pnfkb<-pnfkb[keep]; bck<-bc[keep]

  isEpi <- grepl("epithelium", lab, ignore.case=TRUE)
  nb <- dnearneigh(xy, 0, RADIUS); lw <- nb2listw(nb, style="B", zero.policy=TRUE)
  den  <- lag.listw(lw, as.numeric(isEpi), zero.policy=TRUE)
  nh   <- lag.listw(lw, ifelse(isEpi & is.finite(hyp),  hyp,  0), zero.policy=TRUE)
  np   <- lag.listw(lw, ifelse(isEpi & is.finite(polv), polv, 0), zero.policy=TRUE)
  nb_hyp <- zt(ifelse(den>0, nh/den, NA_real_))
  nb_pol <- zt(ifelse(den>0, np/den, NA_real_))
  rm(nb, lw); gc(verbose=FALSE)

  so <- CreateSeuratObject(counts=M)
  so <- AddModuleScore_UCell(so, features=SETS, maxRank=nrow(M), name=NULL)
  md <- so@meta.data; rm(so, M); gc(verbose=FALSE)

  ix <- which(lab == "Macrophage" & (is.finite(nb_hyp) | is.finite(nb_pol)))
  dt <- data.table(sample=nm, barcode=bck[ix], nb_hyp=nb_hyp[ix], nb_pol=nb_pol[ix],
                   glycolytic=md$glycolytic[ix], nfkb_activated=md$nfkb_activated[ix],
                   tissue_resident=md$tissue_resident[ix], chemokine=md$chemokine[ix],
                   pathway_glycolysis=pgly[ix], pathway_nfkb=pnfkb[ix])
  res[[nm]] <- dt
  message(sprintf("   macrophages %d  done in %.1f min", nrow(dt),
                  as.numeric(difftime(Sys.time(), t0, units="mins"))))
}
all <- rbindlist(res)
fwrite(all, file.path(OUT,"fig5e_macrophage_percell_presplit.csv"))
message("\nmacrophages total: ", nrow(all), "  samples: ", uniqueN(all$sample))
message("DONE")
