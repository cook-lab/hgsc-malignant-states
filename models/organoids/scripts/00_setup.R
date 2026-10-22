# Trimmed copy (2026-10-21) of 202603_organoids/scripts/00_setup.R, the script that ran. Kept: Seurat, the project
# paths and the SecA / SecB gene lists that 04, 08 and 10a use (Table S6 parses the two lists). Removed: the plot
# theme, palettes, plotting packages and unused paths. The untrimmed script is in the lab project folder.

# 00_setup.R
# Environment, packages, and shared paths for the organoids project

# --- Packages ---
library(Seurat)

# --- Paths ---
project_dir <- "/Volumes/CookLab/Sarah/organoids"
output_dir  <- file.path(project_dir, "output")

# --- SecA / SecB gene signatures ---
secB_genes <- c("KRT17", "KRT19", "KRT7", "LCN2", "PRSS22", "SLPI", "TACSTD2")
secA_genes <- c("BCAM", "MECOM", "FBXO21", "LGR5", "LPAR3", "PBX1", "SOX17", "RCN2")

cat("Setup loaded.\n")
