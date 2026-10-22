# Trimmed copy (2026-10-22) of ~/shift_and_merge.R, the script that ran (2025-10-14, interactively, in the home folder).
# Merges the two TMA A slides into one object: slide 2 shifted 10,000 um in x and its cell ids prefixed "tma2_" (the
# prefix 05_rebuild_existing_cohort_rerun.R uses to match the labels). Removed: the cell-type palette and one core plot
# at the end. The merged object was saved as xenium_tma_final_sfe_hdf5 in a later interactive step (the copy in
# 2026_final_xenium_analysis/data/previously processed/ is dated 2025-10-17).


#shift and merge
library(SingleR)
library(HDF5Array)
library(SingleCellExperiment)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(viridis)
library(sf)

tma_1 <- loadHDF5SummarizedExperiment(dir = "/Users/sarahnersesian/Desktop/tma_1_sfe_hdf5")
tma_2 <- loadHDF5SummarizedExperiment(dir = "/Users/sarahnersesian/Desktop/tma_2_sfe_hdf5")

#### SHIFT (10,000 only in the x direction) ####
geom_list <- colGeometries(tma_2)
geom_bbox <- lapply(geom_list, function(g) st_bbox(g))
lapply(geom_list, st_bbox)

shift_amount <- c(10000, 0)
geom_shifted <- lapply(geom_list, function(g) {
  st_geometry(g) <- st_geometry(g) + shift_amount
  return(g)
})

colGeometries(tma_2) <- SimpleList(geom_shifted)
lapply(colGeometries(tma_2), st_bbox)
centroids_shifted <- colGeometries(tma_2)$centroids
coords <- st_coordinates(centroids_shifted)
head(coords)
colData(tma_2)$x <- coords[,1]
colData(tma_2)$y <- coords[,2]

new_ids <- paste0("tma2_", colnames(tma_2))
colnames(tma_2) <- new_ids
rownames(colData(tma_2)) <- new_ids

cells_1 <- colnames(tma_1)
cells_2 <- colnames(tma_2)

cat("Number of cells in tma_1:", length(cells_1), "\n")
cat("Number of cells in tma_2:", length(cells_2), "\n")

overlap <- intersect(cells_1, cells_2)
cat("Overlapping cell IDs:\n")
print(overlap)  

for(a in assayNames(tma_1)) {
  cat("tma_1 assay:", a, " - dims:", dim(assay(tma_1, a)), "\n")
}
for(a in assayNames(tma_2)) {
  cat("tma_2 assay:", a, " - dims:", dim(assay(tma_2, a)), "\n")
}

cat("tma_1 colData rows:", nrow(colData(tma_1)), "\n")
cat("tma_2 colData rows:", nrow(colData(tma_2)), "\n")

cat("tma_1 colGeometries names and types:\n")
lapply(colGeometries(tma_1), class)

cat("tma_2 colGeometries names and types:\n")
lapply(colGeometries(tma_2), class)

cat("tma_1 bounding boxes:\n")
lapply(colGeometries(tma_1), st_bbox)

cat("tma_2 bounding boxes:\n")
lapply(colGeometries(tma_2), st_bbox)

#### MERGE by building a new sfe ####

library(SingleCellExperiment)
library(SpatialFeatureExperiment)
library(S4Vectors)

assay_names <- assayNames(tma_1)
combined_assays <- lapply(assay_names, function(a) {
  cbind(assay(tma_1, a), assay(tma_2, a))
})
names(combined_assays) <- assay_names

combined_colData <- rbind(colData(tma_1), colData(tma_2))
combined_colData <- DataFrame(combined_colData)

combined_centroids <- rbind(
  colGeometries(tma_1)[["centroids"]],
  colGeometries(tma_2)[["centroids"]]
)

cellSeg_1 <- st_as_sf(colGeometries(tma_1)[["cellSeg"]])
cellSeg_2 <- st_as_sf(colGeometries(tma_2)[["cellSeg"]])
combined_cellSeg <- rbind(cellSeg_1, cellSeg_2)

combined_geom_list <- SimpleList(
  centroids = combined_centroids,
  cellSeg = combined_cellSeg
)

combined_sfe <- SpatialFeatureExperiment(
  assays = combined_assays,
  colData = combined_colData,
  rowData = rowData(tma_1),  # assuming identical rowData
  colGeometries = combined_geom_list
)
