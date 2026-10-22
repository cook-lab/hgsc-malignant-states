# ============================================================================
# 06_extract_singler.R — SingleR label table per sample
# ============================================================================
# PURPOSE: Write one per-cell table per SFE with the manuscript's 18 Xenium cell
#   types. cell_label_singler = singler_label, with SingleR "Secretory
#   epithelium" split on polarization_UCell at the frozen thresholds of
#   04_reclassification_polarization.R: SecB if >= 0.210501, SecA if <= 0.109429,
#   Intermediate otherwise (NA without a score). "Transitioning epithelium" is
#   read as "Intermediate epithelium". Every Xenium figure and table script,
#   and spatial/02_qc/06b_final_tma_cores.py, take cell identity from these
#   tables (R scripts through spatial/00_setup/singler_labels.R).
#
# USAGE: Rscript spatial/03_annotation_polarization/06_extract_singler.R <sfe_name>
#   once per SFE: sfe_WT01 … sfe_WT08, sfe_tma, sfe_COEUR_TMA_G, sfe_COEUR_TMA_H
#
# INPUTS:
#   - SFE <sfe_name> (load_sfe): singler_label, singler_pruned, polarization_UCell,
#     spatialCoords, patient_id, core_id, sample_type, cell_area
#
# OUTPUTS:
#   - OUT/<sample>_singler_cells.csv.gz (<sample> = sfe_name without "sfe_"):
#     barcode, x, y, singler_label, singler_pruned, polarization_UCell,
#     cell_label_singler, patient_id, core_id, sample_type, cell_area.
#     OUT = $OUTPUT_ROOT/singler_recompute.
#
# MANUSCRIPT PANEL(S): the cell types of every Xenium panel (Figs 4–7, Figures
#   S11–S18, Tables S7–S11).
#
# RUNTIME TIER: moderate (loads one SFE per call)
# ============================================================================
source("spatial/00_setup/00_setup.R")
suppressPackageStartupMessages({library(SpatialFeatureExperiment); library(data.table)})
args <- commandArgs(trailingOnly=TRUE); nm <- args[1]
OUT <- file.path(Sys.getenv("OUTPUT_ROOT", "output"), "singler_recompute"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
T_LOW <- 0.109429; T_HIGH <- 0.210501
s <- sub("^sfe_","",nm)
message("[",Sys.time(),"] ",nm)
sfe <- load_sfe(nm)
cd <- as.data.frame(colData(sfe))
xy <- spatialCoords(sfe)
sl <- as.character(cd$singler_label)
slp <- if("singler_pruned" %in% colnames(cd)) as.character(cd$singler_pruned) else sl
pol <- if("polarization_UCell" %in% colnames(cd)) cd$polarization_UCell else rep(NA_real_, nrow(cd))
# derive cell_label from SingleR + polarization split
lab <- sl
lab[lab=="Transitioning epithelium"] <- "Intermediate epithelium"
is_sec <- !is.na(lab) & lab=="Secretory epithelium"
sub <- ifelse(is.na(pol), NA_character_,
              ifelse(pol >= T_HIGH, "SecB epithelium",
              ifelse(pol <= T_LOW, "SecA epithelium", "Intermediate epithelium")))
lab[is_sec] <- sub[is_sec]
getcol <- function(cn) if(cn %in% colnames(cd)) as.character(cd[[cn]]) else rep(NA_character_, nrow(cd))
getnum <- function(cn) if(cn %in% colnames(cd)) as.numeric(cd[[cn]]) else rep(NA_real_, nrow(cd))
dt <- data.table(barcode=colnames(sfe), x=xy[,1], y=xy[,2],
                 singler_label=sl, singler_pruned=slp, polarization_UCell=pol,
                 cell_label_singler=lab,
                 patient_id=getcol("patient_id"), core_id=getcol("core_id"),
                 sample_type=getcol("sample_type"), cell_area=getnum("cell_area"))
fwrite(dt, file.path(OUT, paste0(s,"_singler_cells.csv.gz")))
# report
cat("  ",s,": ",nrow(dt)," cells\n",sep="")
cat("   singler secretory (pre-split):", sum(is_sec), "\n")
cat("   epithelial after split:", sum(lab %in% c("SecA epithelium","Intermediate epithelium","SecB epithelium"), na.rm=TRUE), "\n")
cat("   ciliated:", sum(lab=="Ciliated epithelium", na.rm=TRUE), "\n")
cat("   polarization non-NA:", sum(!is.na(pol)), "\n")
cat("EXTRACT_DONE ",s,"\n")
