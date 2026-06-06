#!/usr/bin/env Rscript
# Bootstrap the R environment for the HGSC epitype DOWNSTREAM analyses.
# R 4.5.2. Prefers CRAN/Bioconductor BINARY packages (macOS) to avoid source compilation.
# Run:  Rscript renv_bootstrap.R   (then renv::snapshot() to freeze renv.lock)
#
# Spatial stack (SpatialFeatureExperiment/Voyager/sf/spdep/spatstat) underpins all Xenium work;
# survival/mgcv/BayesPrism/ConsensusOV/SingleR/UCell cover atlas + clinical.

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 600)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# Bioconductor release compatible with R 4.5.x
BiocManager::install(version = "3.21", ask = FALSE, update = FALSE)

cran_pkgs <- c(
  "remotes", "renv",
  "tidyverse", "data.table", "arrow",          # data wrangling + parquet
  "ggplot2", "patchwork", "cowplot", "ggrepel", "scales", "RColorBrewer", "viridis",
  "Matrix",
  "mgcv",                                       # GAMs (Fig 4G/5A/5G/6B/6G/6H, SF12-14)
  "survival", "survminer",                      # Cox + Kaplan-Meier (Fig 7)
  "spdep", "spatstat", "sf", "dbscan", "RANN",  # spatial: Lee's L/BiLISA, frNN, nn2, morphometry
  "cluster",                                    # silhouette (k selection)
  "pheatmap", "circlize"                        # heatmaps; ComplexHeatmap below
)

bioc_pkgs <- c(
  "SingleCellExperiment", "SummarizedExperiment", "S4Vectors",
  "HDF5Array", "DelayedArray", "rhdf5",         # HDF5-backed SFE objects
  "SpatialFeatureExperiment", "Voyager",        # Xenium spatial framework
  "SingleR", "celldex",                         # spatial cell-type annotation
  "UCell",                                       # signature scoring (atlas 18b + xenium 06d/9b)
  "scran", "scater", "scuttle",                  # normalization/QC
  "ComplexHeatmap",
  "consensusOV"                                  # TCGA molecular subtyping (Fig 3H) — note lowercase 'c'
)

github_pkgs <- c(
  "navinlabcode/copykat",                        # CNV (Fig 1J, SF4C, SF7)
  "Danko-Lab/BayesPrism/BayesPrism",             # bulk deconvolution (Fig 7 epithelial fraction)
  "dmcable/spacexr"                              # RCTD/SPLIT (06g_clean_split)
)

install.packages(cran_pkgs)
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)
for (p in github_pkgs) try(remotes::install_github(p, upgrade = "never"), silent = FALSE)

cat("\n--- verification ---\n")
need <- c("SpatialFeatureExperiment","Voyager","SingleR","UCell","mgcv","spdep","spatstat",
          "survival","sf","dbscan","RANN","arrow","ConsensusOV","copykat","BayesPrism","spacexr")
for (p in need) cat(sprintf("%-28s %s\n", p, requireNamespace(p, quietly = TRUE)))
