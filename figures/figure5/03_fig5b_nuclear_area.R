#!/usr/bin/env Rscript
# ============================================================================
# Fig 5b -- secretory epithelial nuclear morphology (SingleR labels)
# ----------------------------------------------------------------------------
# Cell identity is the SingleR label (cell_label_singler, via
# spatial/00_setup/singler_labels.R), secretory cells split SecA / Intermediate /
# SecB at the frozen polarization thresholds.
# Measures, per cell: Xenium nucleus_area (um^2) and N:C ratio = nucleus_area/cell_area
# (the panel and the Results use nucleus_area; nc_ratio is kept in the table).
# Units: WT = section (core 8); TMA = patient, on the final TMA core set
# (spatial/00_setup/final_tma_cores.R: per-array core QC on SingleR labels, primary HGSC,
# treatment-naive, >= 500 cells per patient).
# OUTPUT: output_root/split_recompute/fig5b_morphometrics_presplit.csv (per-unit medians per state;
#   folder and file names are historical); read by 11_fig5_compose_panels.py (5b) and
#   04_fig5b_nuclear_area_medians.py (medians in the Results).
# MANUSCRIPT PANEL(S): Fig 5b. RUNTIME TIER: heavy (loads 11 SFEs).
# ============================================================================
source("spatial/00_setup/00_setup.R")
source("spatial/00_setup/singler_labels.R")   # presplit_labels() = singler_labels(): SingleR labels
source("spatial/00_setup/final_tma_cores.R")
suppressPackageStartupMessages({
  library(SpatialFeatureExperiment); library(SummarizedExperiment); library(data.table)
})
R   <- cfg_path("output_root")
OUT <- cfg_path("output_root", "split_recompute")
EPI <- c("SecA epithelium","Intermediate epithelium","SecB epithelium")
MIN_PATIENT_CELLS <- 500   # cohort rule (as Fig 4)
CORE8 <- paste0("sfe_", CFG$cohort$whole_tissue)


percell <- function(nm){
  sfe <- load_sfe(nm)
  cd <- as.data.table(as.data.frame(colData(sfe)[, c("cell_area","nucleus_area")]))
  cd[, barcode := colnames(sfe)]
  extra <- intersect(c("core_id","sample_type"), colnames(colData(sfe)))
  if (length(extra)) cd <- cbind(cd, as.data.table(as.data.frame(colData(sfe)[, extra, drop=FALSE])))
  L <- presplit_labels(sfe, nm)
  rm(sfe); gc(verbose=FALSE)
  cd[L, on="barcode", lab := i.lab]
  # keep ALL cells: the >=500 cohort rule is a patient TOTAL across every cell
  # type, so the epithelial restriction must come AFTER that filter.
  cd <- cd[!is.na(nucleus_area) & !is.na(cell_area) & cell_area > 0]
  cd[, nc_ratio := nucleus_area/cell_area]
  cd[, is_epi := lab %in% EPI]
  cd
}
rows <- list()
# ---- WT: unit = section ----
for (nm in CORE8) {
  message("[WT] ", nm)
  cd <- percell(nm)
  g <- cd[is_epi == TRUE, .(n=.N, nucleus_area=median(nucleus_area), nc_ratio=median(nc_ratio)), by=lab]
  # WT: all 8 sections far exceed the 500-cell floor
  rows[[nm]] <- data.table(unit=sub("^sfe_","",nm), source="WT",
                           cell_label_split=g$lab, n=g$n,
                           nucleus_area=g$nucleus_area, nc_ratio=g$nc_ratio)
  message("   ", paste(sprintf("%s n=%d", g$lab, g$n), collapse=" | "))
}
# ---- TMA: unit = patient ----
# patient_id / core_id / sample_type come from the SingleR tables themselves.
# Cells of the final TMA core set only (core QC, primary HGSC, treatment-naive, patient >= 500 cells on
# those cores), so the patient-total >= 500 check below always passes; kept as a guard.
specs <- list(list(nm="sfe_tma"), list(nm="sfe_COEUR_TMA_G"), list(nm="sfe_COEUR_TMA_H"))
tma <- list()
for (sp in specs) {
  nm <- sp$nm; message("[TMA] ", nm)
  cd <- percell(nm)                       # cell_area, nucleus_area, lab (SingleR)
  stem <- sub("^sfe_","",nm)
  sgl <- fread(file.path(R,"singler_recompute",paste0(stem,"_singler_cells.csv.gz")),
               select=c("barcode","patient_id","core_id","sample_type"))
  cd[sgl, on="barcode", `:=`(patient=i.patient_id, core=as.character(i.core_id),
                             stype=as.character(i.sample_type))]
  cd <- cd[!is.na(core) & core %in% final_core_ids(nm)]   # final TMA cores only
  cd <- cd[!is.na(patient)]
  cd[, patient := if (grepl("COEUR", nm)) paste0("COEUR-", patient) else as.character(patient)]
  message("   cells ", nrow(cd), "  patients ", uniqueN(cd$patient))
  tma[[nm]] <- cd[, .(barcode, lab, nucleus_area, nc_ratio, patient, is_epi)]
}
td <- rbindlist(tma)
ptot <- td[, .(n_total=.N), by=patient][n_total >= MIN_PATIENT_CELLS]   # ALL cell types
message("\nTMA patients passing >=", MIN_PATIENT_CELLS, " cells: ", nrow(ptot),
        " of ", uniqueN(td$patient))
td <- td[patient %in% ptot$patient]
tg <- td[is_epi == TRUE, .(n=.N, nucleus_area=median(nucleus_area), nc_ratio=median(nc_ratio)),
         by=.(patient, lab)]
out <- rbind(rbindlist(rows),
             data.table(unit=tg$patient, source="TMA", cell_label_split=tg$lab,
                        n=tg$n, nucleus_area=tg$nucleus_area, nc_ratio=tg$nc_ratio))
fwrite(out, file.path(OUT, "fig5b_morphometrics_presplit.csv"))
message("\nWT units: ", uniqueN(out[source=="WT", unit]),
        " | TMA units: ", uniqueN(out[source=="TMA", unit]))
message("DONE")
