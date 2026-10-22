# 10_markers_config.R
# Annotation marker panel for GSE299200 (primary HGSC tumor cells + autologous PBMCs).
# Cell-type scheme matches Plesselova et al. (bioRxiv 2025.06.23.661098):
#   epithelial (HGSC), CAF, CD4 T, CD8 T, NK, B, macrophage/monocyte (+ cycling).
# Sourced by 10a_geo_cluster_annotate.R.  Sarah Nersesian, 2026-07-14.
#
# NOTE ON KERATINS: KRT19/KRT7/KRT17 are epithelial identity markers AND SecB
# signature genes. To avoid circularity, the epithelial compartment is anchored on
# EPCAM + PAX8 + WT1 (HGSC/Mullerian), NOT keratins. Keratins are reserved for the
# downstream SecA/SecB UCell scoring (step 10b).

annotation_markers <- list(

  # --- Epithelial / HGSC tumor (anchor on EPCAM/PAX8/WT1, not keratins) ---
  Epithelial = c("EPCAM", "PAX8", "WT1", "MUC16", "MSLN", "CD24", "KRT8", "KRT18", "CDH1"),

  # --- Cancer-associated fibroblast ---
  CAF = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "PDGFRB", "ACTA2", "FAP", "THY1", "POSTN"),

  # --- Pan-immune gate ---
  Immune_pan = c("PTPRC"),

  # --- T cells (pan + CD4 / CD8 / Treg) ---
  T_pan  = c("CD3D", "CD3E", "CD3G", "TRAC"),
  T_CD4  = c("CD4", "IL7R"),
  T_reg  = c("FOXP3", "IL2RA"),
  T_CD8  = c("CD8A", "CD8B", "GZMK", "GZMB"),

  # --- NK (CD3-negative distinguishes from T) ---
  NK = c("NCAM1", "NKG7", "GNLY", "KLRD1", "KLRF1"),

  # --- B / plasma ---
  B = c("CD19", "MS4A1", "CD79A", "CD79B", "MZB1", "SDC1"),

  # --- Myeloid / macrophage / monocyte ---
  Myeloid = c("CD14", "LYZ", "CD68", "CD163", "FCGR3A", "C1QA", "ITGAM", "AIF1"),

  # --- Cross-cutting: proliferation ---
  Cycling = c("MKI67", "TOP2A", "PCNA")
)

# Final cell-type labels to resolve (paper scheme, matched exactly):
target_celltypes <- c("Epithelial", "CAF", "CD4_T", "CD8_T", "NK", "B", "Myeloid", "Cycling")
