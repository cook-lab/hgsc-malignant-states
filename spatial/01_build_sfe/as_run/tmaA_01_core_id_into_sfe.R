# Copy (2026-10-22) of ~/Core_ID_into_SFE.R, the script that ran (2025-09-19, interactively, in the home folder). TMA A
# core and patient labels: each cell of slide 1 and slide 2 is assigned to the core whose bounding box (from the Xenium
# core map, extended by 5 um) contains it; cells in no box are "Off core". The objects tma_1 / tma_2 and readxl were
# already loaded in the session. The code is unchanged.
# Input not on disk: Desktop/TMA_xenium_map.xlsx (per core: TMA_id, core_id, patient_id and the x / y min / max of its
# box). The labels survive in previously processed/xenium_tma_final_sfe_hdf5/se.rds, which
# spatial/01_build_sfe/05_rebuild_existing_cohort_rerun.R reads (one patient per core there, checked 2026-10-22).

#to incorporate the core_id regions into the metadata of the SFE
library(dplyr)
library(sf)

TMA_xenium_map <- read_excel("~/Desktop/TMA_xenium_map.xlsx")
TMA_xenium_map_bbox <- TMA_xenium_map %>%
  dplyr::mutate(
    xmin_ext = `x min` - 5,
    xmax_ext = `x max` + 5,
    ymin_ext = `y min` - 5, 
    ymax_ext = `y max` + 5)

core_sf <- TMA_xenium_map_bbox %>%
  filter(TMA_id == 1,   
         !is.na(xmin_ext) & !is.na(xmax_ext) &
           !is.na(ymin_ext) & !is.na(ymax_ext)) %>%
  rowwise() %>%
  mutate(geometry = st_sfc(
    st_polygon(list(matrix(
      c(xmin_ext, ymin_ext,
        xmax_ext, ymin_ext,
        xmax_ext, ymax_ext,
        xmin_ext, ymax_ext,
        xmin_ext, ymin_ext),
      ncol = 2, byrow = TRUE)))
  )) %>%
  st_as_sf(crs = NA_crs_)

coords <- spatialCoords(tma_1)
cells_sf <- st_as_sf(
  data.frame(cell_id = colnames(tma_1),
             x = coords[,1],
             y = coords[,2]),
  coords = c("x","y"),
  crs = st_crs(core_sf))
cells_with_core <- st_join(cells_sf,
  core_sf %>% select(core_id, TMA_id, patient_id))

cells_with_core$core_id <- ifelse(
  is.na(cells_with_core$core_id), "Off core", cells_with_core$core_id)
colData(tma_1)$core_id <- cells_with_core$core_id[match(colnames(tma_1), cells_with_core$cell_id)]
colData(tma_1)$patient_id <- cells_with_core$patient_id[match(colnames(tma_1), cells_with_core$cell_id)]

saveHDF5SummarizedExperiment(tma_1, dir = "tma_1_sfe_hdf5", replace = TRUE)

core_sf <- TMA_xenium_map_bbox %>%
  filter(TMA_id == 2,   
         !is.na(xmin_ext) & !is.na(xmax_ext) &
           !is.na(ymin_ext) & !is.na(ymax_ext)) %>%
  rowwise() %>%
  mutate(geometry = st_sfc(
    st_polygon(list(matrix(
      c(xmin_ext, ymin_ext,
        xmax_ext, ymin_ext,
        xmax_ext, ymax_ext,
        xmin_ext, ymax_ext,
        xmin_ext, ymin_ext),
      ncol = 2, byrow = TRUE))))) %>%
  st_as_sf(crs = NA_crs_)
coords <- spatialCoords(tma_2)
cells_sf <- st_as_sf(
  data.frame(cell_id = colnames(tma_2),
             x = coords[,1],
             y = coords[,2]),
  coords = c("x","y"),
  crs = st_crs(core_sf))
cells_with_core <- st_join(cells_sf,
  core_sf %>% select(core_id, TMA_id, patient_id))
cells_with_core$core_id <- ifelse(is.na(cells_with_core$core_id), "Off core", cells_with_core$core_id)
colData(tma_2)$core_id <- cells_with_core$core_id[match(colnames(tma_2), cells_with_core$cell_id)]
colData(tma_2)$patient_id <- cells_with_core$patient_id[match(colnames(tma_2), cells_with_core$cell_id)]


saveHDF5SummarizedExperiment(tma_2, dir = "tma_2_sfe_hdf5", replace = TRUE)

