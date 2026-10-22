# ============================================================================
# final_tma_cores.R: the final TMA core set for R figure and table scripts
# ----------------------------------------------------------------------------
# PURPOSE: every TMA analysis uses the cores of spatial/02_qc/06b_final_tma_cores.py (author decision 2026-10-19):
#   tumour cores that pass per-array core QC on SingleR labels, from patients with primary HGSC who were
#   treatment-naive at sampling and have >= 500 cells on those cores. Source this file instead of reading a
#   core-QC table.
# INPUTS: $FINAL_TMA_CORES, default $OUTPUT_ROOT/07_core_qc/final_tma_cores.csv (raw ids; lab-internal)
# USAGE:
#   source(".../spatial/00_setup/final_tma_cores.R")
#   keep <- final_core_ids("sfe_COEUR_TMA_G")      # character vector of raw core ids
#   d <- d[as.character(core_id) %in% keep]
#   Arrays: "tma"/"sfe_tma" (TMA A), "COEUR_TMA_G", "COEUR_TMA_H" (TMA B blocks G and H).
# ============================================================================
FINAL_TMA_CORES_CSV <- path.expand(Sys.getenv("FINAL_TMA_CORES",
                                              file.path(Sys.getenv("OUTPUT_ROOT", "output"), "07_core_qc/final_tma_cores.csv")))

final_tma_cores <- function() {
  data.table::fread(FINAL_TMA_CORES_CSV, colClasses = list(character = c("core_id_raw", "patient_id_raw")))
}

final_tma_patients <- function() {
  data.table::fread(file.path(dirname(FINAL_TMA_CORES_CSV), "final_tma_patients.csv"),
                    colClasses = list(character = "patient_id_raw"))
}

final_core_ids <- function(array) {
  a <- sub("^sfe_", "", array)
  tb <- switch(a, tma = c("A", "A"), COEUR_TMA_G = c("B", "G"), COEUR_TMA_H = c("B", "H"),
               stop("unknown TMA array: ", array))
  x <- final_tma_cores()
  x[final == TRUE & tma == tb[1] & block == tb[2], core_id_raw]
}
