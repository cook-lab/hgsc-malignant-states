# ============================================================================
# 04_reclassification_polarization.R — Atlas-aligned polarization reclassification
# ============================================================================
# PURPOSE: Derive the frozen polarization thresholds that split SingleR
#   secretory cells into SecA / Intermediate / SecB, and apply them to the SFE
#   cell_label. On the atlas secretory cells scored by 03_ucell_scoring_noBCAM.R
#   (polarization_UCell = SecB_UCell - SecA_UCell, on the same 457 genes as
#   Xenium), with the atlas states taken from the NMF labels:
#     t_low  = 75th percentile of polarization in atlas SecA cells   (0.109429)
#     t_high = 25th percentile of polarization in atlas SecB cells   (0.210501)
#   SFE secretory cells are then SecA if polarization < t_low, Intermediate if
#   t_low <= polarization < t_high, and SecB if polarization >= t_high.
#
#   The thresholds are written to threshold_summary.csv. Their values are typed
#   into 06_extract_singler.R and spatial/02_qc/05b_core_qc.R (which use
#   SecA <= t_low); 05_clean_split_rctd.R reads the file.
#
# SAMPLES: the SFE writeback covers every SFE in sfe_names_all.
#
# NAMING: atlas NMF labels "Transitioning epithelium" are read as
#   "Intermediate epithelium".
#
# INPUTS:
#   - <output_root>/06d_annotation_noBCAM/atlas_ucell_scores.csv
#   - <sfe_dir>/sfe_*  (cell_label from 02_adaptive_secretory_noBCAM.R;
#       polarization_UCell from 03_ucell_scoring_noBCAM.R)
#
# OUTPUTS:
#   - <output_root>/06f_reclassification_polarization/threshold_summary.csv
#       (t_low, t_high and summary rows: atlas AUC, accuracy, macro-F1, class counts)
#   - SFEs: cell_label overwritten for the secretory cells (SecA / Intermediate /
#       SecB by polarization; other cells unchanged); the previous cell_label is
#       kept as cell_label_06e
#   - <output_root>/06f_reclassification_polarization/.sfe_update_complete
#       (checkpoint: when present, the SFE writeback is skipped)
#
# MANUSCRIPT PANEL(S): the SecA / Intermediate / SecB thresholds of every Xenium
#   panel (Figs 4–7, Figures S11–S18, Tables S7–S11), through 06_extract_singler.R.
#
# RUNTIME TIER: heavy (SFE writeback over all samples)
# ============================================================================

source("spatial/00_setup/00_setup.R")

if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(pROC)

message("\n=== Atlas-aligned Polarization Reclassification ===")
message("[", Sys.time(), "] Starting...")

# --- Paths + constants ------------------------------------------------------

step_dir  <- file.path(out_dir, "06f_reclassification_polarization")
fig_out   <- file.path(step_dir, "figures")
flag_file <- file.path(step_dir, ".sfe_update_complete")

for (d in c(step_dir, fig_out)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

sec_classes <- c("SecA epithelium", "Intermediate epithelium", "SecB epithelium")

# Writeback scope: every SFE in sfe_names_all.
sfe_names <- sfe_names_all

# Standardize atlas NMF label naming (Transitioning -> Intermediate).
remap_nmf_label <- function(x) {
  x <- as.character(x)
  x[x == "Transitioning epithelium"] <- "Intermediate epithelium"
  x
}

# Known valid cell_label values (used for SFE write validation)
valid_cell_labels_extra <- c(
  "Ciliated epithelium", "Mesothelial", "Fibroblast", "Smooth muscle",
  "Pericyte", "Endothelial", "T cell", "NK cell", "B cell", "Plasma cell",
  "Macrophage", "Conventional dendritic cell", "Plasmacytoid dendritic cell",
  "Neutrophil", "Mast cell"
)

# --- Helpers ----------------------------------------------------------------

classify_polarization <- function(pol, t_low, t_high) {
  fifelse(
    pol < t_low, "SecA epithelium",
    fifelse(pol < t_high, "Intermediate epithelium", "SecB epithelium")
  )
}

# ============================================================================
# PART 1: Threshold Derivation on Atlas
# ============================================================================

message("\n[", Sys.time(), "] === PART 1: Atlas Threshold Derivation ===")

atlas_dt <- fread(file.path(out_dir, "06d_annotation_noBCAM", "atlas_ucell_scores.csv"))
if ("polarization" %in% names(atlas_dt)) {
  setnames(atlas_dt, "polarization", "polarization_UCell")
}
atlas_dt[, celltype_nmf := remap_nmf_label(celltype_nmf)]
message("  atlas cells: ", format(nrow(atlas_dt), big.mark = ","))
message("  NMF classes: ", paste(unique(atlas_dt$celltype_nmf), collapse = ", "))

atlas_sec <- atlas_dt[celltype_nmf %in% sec_classes]
atlas_sec[, celltype_nmf := factor(celltype_nmf, levels = sec_classes)]
message("  atlas secretory cells: ", format(nrow(atlas_sec), big.mark = ","))
print(atlas_sec[, .N, by = celltype_nmf])

# Thresholds: SecA p75 + SecB p25 on polarization
t_low  <- unname(quantile(atlas_sec[celltype_nmf == "SecA epithelium", polarization_UCell], 0.75,
                           na.rm = TRUE))
t_high <- unname(quantile(atlas_sec[celltype_nmf == "SecB epithelium", polarization_UCell], 0.25,
                           na.rm = TRUE))

message(sprintf("  t_low  (atlas SecA polarization p75): %.4f", t_low))
message(sprintf("  t_high (atlas SecB polarization p25): %.4f", t_high))

if (!(t_low < t_high)) {
  stop(sprintf(
    "Inverted band: t_low (%.4f) >= t_high (%.4f). Atlas SecA/SecB polarization distributions overlap too heavily — Strategy c is not viable. Investigate before proceeding.",
    t_low, t_high
  ))
}

# Binary AUC on polarization (SecA vs SecB) for reference
atlas_binary <- atlas_sec[celltype_nmf %in% c("SecA epithelium", "SecB epithelium")]
atlas_binary[, label := fifelse(celltype_nmf == "SecB epithelium", 1L, 0L)]
roc_obj <- pROC::roc(response = atlas_binary$label,
                      predictor = atlas_binary$polarization_UCell,
                      direction = "<", levels = c(0, 1), quiet = TRUE)
auc_val <- as.numeric(pROC::auc(roc_obj))
message(sprintf("  AUC (polarization, SecA vs SecB): %.4f  (06e SecB-only was 0.9365)", auc_val))

atlas_sec[, threshold_class := factor(classify_polarization(polarization_UCell, t_low, t_high),
                                        levels = sec_classes)]

conf_atlas <- atlas_sec[, .N, by = .(celltype_nmf, threshold_class)]
conf_atlas_wide <- dcast(conf_atlas, celltype_nmf ~ threshold_class,
                          value.var = "N", fill = 0)
message("  Atlas confusion (rows=NMF, cols=threshold):")
print(conf_atlas_wide)

# Proper per-class F1
cm <- table(atlas_sec$celltype_nmf, atlas_sec$threshold_class)
cm <- cm[sec_classes, sec_classes]
per_class <- lapply(seq_along(sec_classes), function(i) {
  tp <- cm[i, i]
  fn <- sum(cm[i, ]) - tp
  fp <- sum(cm[, i]) - tp
  prec <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  rec  <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  f1   <- if (is.na(prec) || is.na(rec) || (prec + rec) == 0) NA_real_ else 2 * prec * rec / (prec + rec)
  data.table(class = sec_classes[i], n = sum(cm[i, ]),
              precision = prec, recall = rec, f1 = f1)
})
per_class_dt <- rbindlist(per_class)
macro_f1   <- mean(per_class_dt$f1, na.rm = TRUE)
overall_accuracy <- sum(diag(cm)) / sum(cm)

message(sprintf("  atlas overall accuracy: %.3f  macro-F1: %.3f", overall_accuracy, macro_f1))
print(per_class_dt)

# --- Save threshold summary -------------------------------------------------

thresh_summary <- data.table(
  metric = c("strategy", "t_low_SecA_p75", "t_high_SecB_p25",
             "AUC_polarization_SecA_vs_SecB", "AUC_06e_SecB_only_ref",
             "atlas_overall_accuracy", "atlas_macro_F1",
             "n_atlas_SecA", "n_atlas_Intermediate", "n_atlas_SecB"),
  value = c(
    "SecA_p75 + SecB_p25 on polarization_UCell",
    sprintf("%.6f", t_low), sprintf("%.6f", t_high),
    sprintf("%.6f", auc_val), "0.9365",
    sprintf("%.6f", overall_accuracy), sprintf("%.6f", macro_f1),
    sum(atlas_sec$celltype_nmf == "SecA epithelium"),
    sum(atlas_sec$celltype_nmf == "Intermediate epithelium"),
    sum(atlas_sec$celltype_nmf == "SecB epithelium")
  )
)
fwrite(thresh_summary, file.path(step_dir, "threshold_summary.csv"))


# ============================================================================
# PART 3: SFE Writeback
# ============================================================================

message("\n[", Sys.time(), "] === PART 3: SFE Writeback ===")

if (file.exists(flag_file)) {
  message("  Checkpoint flag found — SFEs already updated. Skipping.")
} else {
  for (nm in sfe_names) {
    t0 <- Sys.time()
    message("[", Sys.time(), "] Updating ", nm, "...")
    sfe <- load_sfe(nm)
    cd  <- colData(sfe)

    if (!"polarization_UCell" %in% colnames(cd)) {
      rm(sfe); gc(verbose = FALSE)
      stop(sprintf("SFE %s lacks polarization_UCell — rerun 03_ucell_scoring_noBCAM.R first.", nm))
    }

    # Keep the current cell_label (from 02_adaptive_secretory_noBCAM.R) as cell_label_06e
    cd$cell_label_06e <- as.character(cd$cell_label)

    cur_label <- as.character(cd$cell_label)
    pol       <- as.numeric(cd$polarization_UCell)
    is_sec    <- cur_label %in% sec_classes
    has_pol   <- !is.na(pol)
    do_update <- is_sec & has_pol

    new_label <- cur_label
    new_label[do_update] <- classify_polarization(pol[do_update], t_low, t_high)

    cd$cell_label <- new_label
    colData(sfe)  <- cd

    allowed <- c(sec_classes, valid_cell_labels_extra)
    bad <- setdiff(unique(new_label), allowed)
    if (length(bad) > 0) {
      warning(sprintf("SFE %s has unexpected cell_label values after update: %s",
                       nm, paste(bad, collapse = ", ")))
    }

    n_upd     <- sum(do_update)
    n_changed_sfe <- sum(do_update & (cd$cell_label_06e[do_update] != new_label[do_update]))
    message(sprintf("  %s: %s secretory cells reclassified  (%s changed, %s unchanged, %s non-secretory preserved)",
                    nm,
                    format(n_upd, big.mark = ","),
                    format(n_changed_sfe, big.mark = ","),
                    format(n_upd - n_changed_sfe, big.mark = ","),
                    format(sum(!is_sec), big.mark = ",")))

    save_sfe(sfe, nm)
    rm(sfe); gc(verbose = FALSE)
    message(sprintf("  [%s] %s update complete in %.1f s",
                     Sys.time(), nm, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }

  writeLines(as.character(Sys.time()), flag_file)
  message("[", Sys.time(), "] All ", length(sfe_names), " SFEs updated. Checkpoint flag written.")
}

# --- Session info -----------------------------------------------------------

log_session()
message("\n[", Sys.time(), "] 06f reclassification complete.")
