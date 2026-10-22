# ============================================================================
# singler_labels.R: SingleR cell labels for the Xenium figure scripts
# ----------------------------------------------------------------------------
# PURPOSE: singler_labels(sfe, name) returns one row per SFE cell (barcode, lab,
#   polarization) from the SingleR label table of that sample. Sourced by the
#   Fig 4d / 4e builders and the Fig 5 builders.
#
# INPUTS: <output_root>/singler_recompute/<sample>_singler_cells.csv.gz
#   (spatial/03_annotation_polarization/06_extract_singler.R):
#   cell_label_singler = singler_label with "Secretory epithelium" split into
#     SecA / Intermediate / SecB on polarization_UCell at the frozen thresholds
#     (>= 0.210501 SecB; <= 0.109429 SecA; else Intermediate)
#   polarization_UCell = SecB - SecA UCell score
#
# OUTPUTS: none (defines singler_labels() and its alias presplit_labels()).
#
# NOTE: colData(sfe)$cell_label is not used: that SFE column was last written by
# the RCTD/SPLIT step (spatial/03_annotation_polarization/05_clean_split_rctd.R),
# whose labels the manuscript does not use. Every SFE cell is returned; cells
# missing from the table get NA.
# ============================================================================
singler_labels <- function(sfe, nm) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  stem <- sub("^sfe_", "", nm)
  f <- file.path(cfg_path("output_root"), "singler_recompute",
                 paste0(stem, "_singler_cells.csv.gz"))
  stopifnot("SingleR table missing" = file.exists(f))
  s <- data.table::fread(f, select = c("barcode","cell_label_singler","polarization_UCell"))
  s[cell_label_singler == "Transitioning epithelium",
    cell_label_singler := "Intermediate epithelium"]
  bc <- colnames(sfe)
  m  <- match(bc, s$barcode)
  lab <- s$cell_label_singler[m]
  pol <- s$polarization_UCell[m]
  message(sprintf("   SingleR labels: %d SFE cells; matched %d; labelled %d",
                  length(bc), sum(!is.na(m)), sum(!is.na(lab) & nzchar(lab))))
  data.table::data.table(barcode = bc, lab = lab, polarization = pol)
}
# alias called by the Fig 5 builders (figures/figure5/01, 03, 06, 08)
presplit_labels <- function(sfe, nm) singler_labels(sfe, nm)
