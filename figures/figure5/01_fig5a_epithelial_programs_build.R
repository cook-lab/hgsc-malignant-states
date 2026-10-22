#!/usr/bin/env Rscript
# ============================================================================
# Fig 5a -- SECRETORY-EPITHELIUM-INTRINSIC program scores vs polarization
# ----------------------------------------------------------------------------
# Strictly epithelium-intrinsic: every score is computed on the secretory
# epithelial cell ITSELF (SecA / Intermediate / SecB), not on its neighbourhood.
# Cells: all SingleR secretory cells with a polarization score in the 8 whole
# tissues (labels via spatial/00_setup/singler_labels.R).
#
# Curves (8), matching the results text:
#   UP    NF-kB, Glycolysis, LDHA, TNFA
#   DOWN  Oncogenic growth, DNA damage repair, Epigenetic silencing, Proliferation
#
# Sources:
#   * SFE colData pathway_*  -- written by spatial/04_pathway_scoring/04c_pathway_scoring_allcells.R
#       (UCell, full 477-gene panel, maxRank=477, all cells)
#       nfkb, glycolysis, proliferation, and the four growth arms
#   * computed here with identical UCell settings (not among the 37 modules):
#       DNA damage repair  = BRCA1, BRCA2, ATM, CHEK2
#       Epigenetic silencing = DNMT3A, HDAC1, KDM5A
#   * logcounts, z-scored: LDHA and TNF (single genes, as in the results text)
#
#   Oncogenic growth = mean of the z-scored Myc / PI3K / Wnt / Notch activating
#   pathway scores (the four growth arms named in the text).
#
# OUTPUT: output_root/split_recompute/fig5a_epi_program_scores_presplit.csv (1,235,350
#   secretory cells; also read by figures/supplementary/21_figS16-18_per_sample_gams.R).
#   The folder and file names are historical; the column cell_label_split holds the SingleR label.
# MANUSCRIPT PANEL(S): Fig 5a (fit: 02_fig5a_epithelial_programs_fit.R); Figure S17.
# RUNTIME TIER: heavy (loads the 8 whole-tissue SFEs; UCell on raw counts).
# ============================================================================
source("spatial/00_setup/00_setup.R")
source("spatial/00_setup/singler_labels.R")   # presplit_labels() = singler_labels(): SingleR labels
suppressPackageStartupMessages({
  library(SpatialFeatureExperiment); library(SummarizedExperiment)
  library(data.table); library(Seurat); library(UCell)
})
R   <- cfg_path("output_root")
OUT <- cfg_path("output_root", "split_recompute")
EPI <- c("SecA epithelium","Intermediate epithelium","SecB epithelium")
CORE8 <- paste0("sfe_", CFG$cohort$whole_tissue)
stopifnot(length(CORE8) == 8)

DDR  <- c("BRCA1","BRCA2","ATM","CHEK2")
EPIS <- c("DNMT3A","HDAC1","KDM5A")
GROWTH <- paste0("pathway_", c("myc_activating","pi3k_activating",
                               "wnt_activating","notch_activating"))
FROM_SFE <- c("pathway_nfkb","pathway_glycolysis","pathway_proliferation", GROWTH)



res <- list()
for (nm in CORE8) {
  t0 <- Sys.time(); message("[", format(t0), "] ", nm)
  sfe <- load_sfe(nm)
  L   <- presplit_labels(sfe, nm)
  ps  <- L[lab %in% EPI & is.finite(polarization)]
  idx <- match(ps$barcode, colnames(sfe))
  keep <- !is.na(idx); ps <- ps[keep]; idx <- idx[keep]
  message("   secretory epithelial cells: ", length(idx))

  cd <- as.data.table(as.data.frame(colData(sfe)[idx, FROM_SFE, drop=FALSE]))
  stopifnot(!anyNA(cd$pathway_nfkb))

  # DDR + epigenetic silencing: same UCell settings as the SFE pathway_* columns
  M <- counts(sfe)[, idx, drop=FALSE]; M <- as(M, "CsparseMatrix")
  sets <- list(ddr = intersect(DDR, rownames(M)), epi_silencing = intersect(EPIS, rownames(M)))
  message("   ddr ", length(sets$ddr), "/", length(DDR),
          " ; epi_silencing ", length(sets$epi_silencing), "/", length(EPIS),
          " ; maxRank ", nrow(M))
  so <- CreateSeuratObject(counts = M)
  so <- AddModuleScore_UCell(so, features = sets, maxRank = nrow(M), name = NULL)

  lg <- logcounts(sfe)[c("LDHA","TNF"), idx, drop=FALSE]

  out <- data.table(sample = nm, barcode = ps$barcode,
                    cell_label_split = ps$lab,
                    polarization = ps$polarization)
  out[, nfkb          := cd$pathway_nfkb]
  out[, glycolysis    := cd$pathway_glycolysis]
  out[, proliferation := cd$pathway_proliferation]
  out[, ddr           := so@meta.data$ddr]
  out[, epi_silencing := so@meta.data$epi_silencing]
  out[, LDHA          := as.numeric(lg["LDHA", ])]
  out[, TNFA          := as.numeric(lg["TNF", ])]
  # Oncogenic growth = mean of the four z-scored growth arms (within sample)
  gz <- as.data.table(lapply(GROWTH, function(g){ v<-cd[[g]]; (v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE) }))
  out[, oncogenic_growth := rowMeans(as.matrix(gz), na.rm = TRUE)]
  res[[nm]] <- out
  message(sprintf("   done in %.1f min", as.numeric(difftime(Sys.time(), t0, units="mins"))))
  rm(sfe, M, so, lg, cd, gz); gc(verbose = FALSE)
}
all <- rbindlist(res)
fwrite(all, file.path(OUT, "fig5a_epi_program_scores_presplit.csv"))
message("\nrows: ", nrow(all), "  samples: ", uniqueN(all$sample))
print(all[, .(n=.N), by=cell_label_split])
message("DONE")
