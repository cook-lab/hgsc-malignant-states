#!/usr/bin/env Rscript
# ============================================================================
# Fig 5c, 5d (WT01 ROI 06) and 5e, 5g, 5h (TMA A core 2) -- export cell polygons + values
# ----------------------------------------------------------------------------
#   * labels  = SingleR cell_label_singler (via spatial/00_setup/singler_labels.R)
#   * scores  = SFE pathway_* (04c: UCell, full 477-gene panel, ALL cells);
#               LDHA = logcounts
#   * ROI selection in the spatialCoords frame (y POSITIVE); the cellSeg polygons
#     are stored with y NEGATIVE (the two frames are mirrored in y).
# Output (output_root/figures/figure5/roi_data/): polygons as vertex tables (exterior
# ring) + per-cell attributes, WT01_roi06_{poly.csv.gz,attr.csv} and
# tma_core2_{poly.csv.gz,attr.csv}, drawn by 11_fig5_compose_panels.py.
# MANUSCRIPT PANEL(S): Fig 5c, 5d, 5e, 5g, 5h. RUNTIME TIER: moderate (loads 2 SFEs).
# ============================================================================
source("spatial/00_setup/00_setup.R")
source("spatial/00_setup/singler_labels.R")   # singler_labels(): SingleR labels
suppressWarnings(suppressMessages({
  library(SpatialFeatureExperiment); library(SummarizedExperiment); library(sf); library(data.table)
}))
OUT <- cfg_path("output_root","figures","figure5","roi_data"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

export <- function(sfe, keep, stem, extra) {
  bc <- colnames(sfe)[keep]
  cs <- colGeometry(sfe, "cellSeg"); cs <- cs[rownames(cs) %in% bc, ]
  message(sprintf("  %s: %d cells selected, %d with polygons", stem, length(bc), nrow(cs)))
  g <- sf::st_cast(sf::st_geometry(cs), "MULTIPOLYGON")
  xy <- sf::st_coordinates(g)                     # X Y L1(ring) L2(part) L3(feature)
  v <- data.table(cell = rownames(cs)[xy[, "L3"]], part = xy[, "L2"], ring = xy[, "L1"],
                  x = xy[, "X"], y = xy[, "Y"])[ring == 1]
  v[, ring := NULL]
  fwrite(v, file.path(OUT, paste0(stem, "_poly.csv.gz")))
  a <- extra[barcode %in% rownames(cs)]
  fwrite(a, file.path(OUT, paste0(stem, "_attr.csv")))
}

# ---- 5c, 5d: WT01 ROI 06 (800 x 800 um) ---------------------------------
sfe <- load_sfe("sfe_WT01"); L <- singler_labels(sfe, "sfe_WT01")
co <- spatialCoords(sfe)
keep <- co[,1] >= 1499 & co[,1] <= 2299 & co[,2] >= 7133 & co[,2] <= 7933
cd <- colData(sfe)
attr <- data.table(barcode = colnames(sfe), lab = L$lab, polarization = L$polarization,
                   proliferation = cd$pathway_proliferation, nfkb = cd$pathway_nfkb,
                   glycolysis = cd$pathway_glycolysis,
                   LDHA = as.numeric(logcounts(sfe)["LDHA", ]))[keep]
export(sfe, keep, "WT01_roi06", attr)
print(attr[, .N, by = lab][order(-N)][1:8])
rm(sfe); gc(verbose = FALSE)

# ---- 5e, 5g, 5h: TMA A core 2 (TMA-P003) ----------------------------------
sfe <- load_sfe("sfe_tma"); L <- singler_labels(sfe, "sfe_tma")
s <- fread(file.path(cfg_path("output_root"),"singler_recompute","tma_singler_cells.csv.gz"),
           select = c("barcode","core_id","patient_id"))
m <- match(colnames(sfe), s$barcode)
keep <- as.character(s$core_id[m]) == "2" & !is.na(m)
attr <- data.table(barcode = colnames(sfe), lab = L$lab, polarization = L$polarization,
                   hypoxia = colData(sfe)$pathway_hypoxia,
                   patient_id = s$patient_id[m])[keep]
export(sfe, keep, "tma_core2", attr)
print(attr[, .N, by = lab][order(-N)][1:8])
message("DONE")
