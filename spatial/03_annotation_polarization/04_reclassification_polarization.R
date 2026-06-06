# ============================================================================
# 04_reclassification_polarization.R — Atlas-aligned polarization reclassification
# ============================================================================
# PURPOSE: Reclassify Xenium secretory cells into SecA / Intermediate / SecB
#   using atlas-calibrated thresholds on the bivariate score
#   polarization_UCell = SecB_UCell - SecA_UCell (this captures both axes and
#   matches the atlas NMF Intermediate definition, unlike a SecB-only cut).
#
# APPROACH (Strategy c):
#   t_low  = 75th percentile of atlas SecA-class polarization
#   t_high = 25th percentile of atlas SecB-class polarization
#   SecA          <- polarization < t_low
#   Intermediate  <- t_low <= polarization < t_high
#   SecB          <- polarization >= t_high
#   These thresholds are written to threshold_summary.csv and FROZEN for use by
#   05_clean_split_rctd.R.
#
# COHORT PIN: writes back to sfe_tma, sfe_tma_filtered, and the 8 published
#   whole tissues; FTE whole-tissue samples excluded.
#
# DETERMINISM: the 50k atlas scatter sub-sample is seeded from config
#   (set.seed(CFG$seed)) so figures/derived stats are reproducible.
#
# NAMING: atlas NMF labels remapped on read "Transitioning epithelium" ->
#   "Intermediate epithelium".
#
# INPUTS:
#   - <output_root>/06d_annotation_noBCAM/atlas_ucell_scores.csv
#   - <output_root>/06d_annotation_noBCAM/xenium_ucell_scores.csv
#   - <sfe_dir>/sfe_* (polarization_UCell already in colData)
#
# OUTPUTS:
#   - <output_root>/06f_reclassification_polarization/threshold_summary.csv
#       (FROZEN thresholds), comparison_vs_atlas_nmf.csv, atlas_per_class_metrics.csv,
#       comparison_vs_06e.csv, reclassified_xenium_scores.csv,
#       biological_plausibility_checks_06f.csv, report HTML + figures
#   - SFEs: cell_label overwritten for secretory cells; cell_label_06e preserved
#
# MANUSCRIPT PANEL(S): the reclassified_xenium_scores.csv override is consumed
#   by nearly every Fig 4–6 ROI/composition panel.
#
# RUNTIME TIER: heavy (SFE writeback over all samples)
# ============================================================================

source("spatial/00_setup/00_setup.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(base64enc)
  library(ragg)
  library(knitr)
  library(patchwork)
})

if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(pROC)

message("\n=== Atlas-aligned Polarization Reclassification ===")
message("[", Sys.time(), "] Starting...")

# --- Paths + constants ------------------------------------------------------

step_dir  <- file.path(out_dir, "06f_reclassification_polarization")
fig_out   <- file.path(step_dir, "figures")
html_out  <- file.path(step_dir, "06f_reclassification_report.html")
flag_file <- file.path(step_dir, ".sfe_update_complete")

for (d in c(step_dir, fig_out)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

sec_classes <- c("SecA epithelium", "Intermediate epithelium", "SecB epithelium")
class_pal   <- ref_palette[sec_classes]

# Writeback scope: merged TMA, filtered TMA, and the 8 published whole tissues.
sfe_names <- c("sfe_tma", "sfe_tma_filtered", sfe_names_wt)

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

save_and_embed <- function(p, name, width = 1200, height = 600, res = 150) {
  pdf_path <- file.path(fig_out, paste0(name, ".pdf"))
  ggsave(pdf_path, plot = p, width = width / res, height = height / res,
         device = cairo_pdf)
  tmp <- tempfile(fileext = ".png")
  ragg::agg_png(tmp, width = width, height = height, res = res)
  print(p)
  dev.off()
  b64 <- base64enc::base64encode(tmp); unlink(tmp)
  sprintf('<img src="data:image/png;base64,%s" style="max-width:100%%;">', b64)
}

make_html_table <- function(df, digits = 3) {
  for (j in seq_along(df)) {
    if (is.numeric(df[[j]])) df[[j]] <- round(df[[j]], digits)
  }
  knitr::kable(df, format = "html", row.names = FALSE,
               table.attr = 'class="styled-table"')
}

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

# Strategy c: SecA p75 + SecB p25 on polarization
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

# --- Save threshold summary + atlas confusion -------------------------------

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
fwrite(conf_atlas_wide, file.path(step_dir, "comparison_vs_atlas_nmf.csv"))
fwrite(per_class_dt, file.path(step_dir, "atlas_per_class_metrics.csv"))
message("  Saved: threshold_summary.csv, comparison_vs_atlas_nmf.csv, atlas_per_class_metrics.csv")


# ============================================================================
# PART 2: Xenium reclassification (CSV-driven, for reporting stats)
# ============================================================================

message("\n[", Sys.time(), "] === PART 2: Xenium Reclassification (reporting) ===")

xen_dt <- fread(file.path(out_dir, "06d_annotation_noBCAM", "xenium_ucell_scores.csv"))
message("  xenium secretory cells (06d CSV scope): ",
        format(nrow(xen_dt), big.mark = ","))

# Load current cell_label (pre-06f) for the SFEs covered by the 06d CSV.
# The 06d CSV was produced before sfe_tma_filtered existed.
csv_sfe_names <- setdiff(sfe_names, "sfe_tma_filtered")

label_list <- lapply(csv_sfe_names, function(nm) {
  sfe <- load_sfe(nm)
  cd  <- colData(sfe)
  dt  <- data.table(barcode_orig = colnames(sfe),
                     cell_label   = as.character(cd$cell_label),
                     sample       = nm)
  dt[, barcode_unique := paste0(nm, "_", barcode_orig)]
  rm(sfe); gc(verbose = FALSE)
  dt
})
labels_dt <- rbindlist(label_list)
message("  total cells with current cell_label in CSV scope: ",
        format(nrow(labels_dt), big.mark = ","))

xen_dt <- merge(xen_dt, labels_dt[, .(barcode_unique, cell_label)],
                 by = "barcode_unique", all.x = TRUE)
setnames(xen_dt, "cell_label", "cell_label_06e")

n_na_pol <- sum(is.na(xen_dt$polarization_UCell))
if (n_na_pol > 0) {
  warning(sprintf("  %d cells with NA polarization_UCell (they will retain cell_label_06e).", n_na_pol))
  fwrite(xen_dt[is.na(polarization_UCell)],
         file.path(step_dir, "NA_polarization_cells.csv"))
}

xen_dt[, cell_label_06f := classify_polarization(polarization_UCell, t_low, t_high)]
xen_dt[is.na(polarization_UCell), cell_label_06f := cell_label_06e]
xen_dt[, cell_label_06f := factor(cell_label_06f, levels = sec_classes)]
xen_dt[, cell_label_06e := factor(cell_label_06e, levels = sec_classes)]

conf_06 <- xen_dt[, .N, by = .(cell_label_06e, cell_label_06f)]
conf_06_wide <- dcast(conf_06, cell_label_06e ~ cell_label_06f,
                       value.var = "N", fill = 0)
fwrite(conf_06_wide, file.path(step_dir, "comparison_vs_06e.csv"))

n_changed <- xen_dt[as.character(cell_label_06e) != as.character(cell_label_06f), .N]
n_total   <- nrow(xen_dt)
pct_changed <- round(n_changed / n_total * 100, 1)
message(sprintf("  Cells changed prior -> 06f: %s / %s (%.1f%%)",
                format(n_changed, big.mark = ","),
                format(n_total,   big.mark = ","), pct_changed))

xen_dt[, changed_vs_06e := as.character(cell_label_06e) != as.character(cell_label_06f)]
fwrite(xen_dt, file.path(step_dir, "reclassified_xenium_scores.csv"))
message("  Saved: comparison_vs_06e.csv, reclassified_xenium_scores.csv")

# Biological plausibility per sample
biop <- xen_dt[, .(
  n_secretory      = .N,
  n_SecA           = sum(cell_label_06f == "SecA epithelium"),
  n_Intermediate   = sum(cell_label_06f == "Intermediate epithelium"),
  n_SecB           = sum(cell_label_06f == "SecB epithelium"),
  pct_SecA_of_sec  = round(100 * sum(cell_label_06f == "SecA epithelium") / .N, 2),
  pct_Int_of_sec   = round(100 * sum(cell_label_06f == "Intermediate epithelium") / .N, 2),
  pct_SecB_of_sec  = round(100 * sum(cell_label_06f == "SecB epithelium") / .N, 2),
  median_pol_SecA  = round(median(polarization_UCell[cell_label_06f == "SecA epithelium"], na.rm = TRUE), 4),
  median_pol_Int   = round(median(polarization_UCell[cell_label_06f == "Intermediate epithelium"], na.rm = TRUE), 4),
  median_pol_SecB  = round(median(polarization_UCell[cell_label_06f == "SecB epithelium"], na.rm = TRUE), 4)
), by = sample]
biop[, flag_low_SecB   := pct_SecB_of_sec < 5]
biop[, flag_high_SecB  := pct_SecB_of_sec > 30]
biop[, flag_pol_SecA_sign := median_pol_SecA >= 0]   # should be negative
biop[, flag_pol_SecB_sign := median_pol_SecB <= 0]   # should be positive
fwrite(biop, file.path(step_dir, "biological_plausibility_checks_06f.csv"))
message("  Saved: biological_plausibility_checks_06f.csv")
print(biop)


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

    # Preserve current labels for audit
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


# ============================================================================
# PART 4: Figures + HTML Report
# ============================================================================

message("\n[", Sys.time(), "] === PART 4: Figures + HTML ===")

# --- Fig A: ROC curve (binary) ----------------------------------------------

roc_df <- data.table(fpr = 1 - roc_obj$specificities, tpr = roc_obj$sensitivities)
p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_line(color = "#0066cc", linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  annotate("text", x = 0.6, y = 0.25,
           label = sprintf("AUC = %.4f (polarization)\n       %.4f (06e SecB-only)",
                            auc_val, 0.9365),
           size = 3.2, color = "#0066cc", hjust = 0) +
  labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)",
       title = "ROC: polarization_UCell discriminating atlas SecA vs SecB",
       subtitle = sprintf("Strategy c thresholds: t_low=%.4f (SecA p75), t_high=%.4f (SecB p25)",
                           t_low, t_high)) +
  coord_fixed() +
  theme_lab(base_size = 9)

img_roc <- save_and_embed(p_roc, "roc_curve_polarization",
                            width = 800, height = 800, res = 150)

# --- Fig B: atlas polarization density by NMF class with threshold lines ----

p_dens <- ggplot(atlas_sec, aes(x = polarization_UCell, fill = celltype_nmf)) +
  geom_density(alpha = 0.35, color = NA) +
  geom_vline(xintercept = t_low,  linetype = "dashed",  color = "#cc0000", linewidth = 0.8) +
  geom_vline(xintercept = t_high, linetype = "dotted", color = "#0066cc", linewidth = 0.8) +
  annotate("text", x = t_low,  y = Inf, label = sprintf("t_low = %.3f",  t_low),
           vjust = 2, hjust = 1.05, size = 2.8, color = "#cc0000") +
  annotate("text", x = t_high, y = Inf, label = sprintf("t_high = %.3f", t_high),
           vjust = 2, hjust = -0.05, size = 2.8, color = "#0066cc") +
  scale_fill_manual(values = class_pal, name = "Atlas NMF class") +
  labs(x = "Polarization UCell (SecB - SecA)", y = "Density",
       title = "Atlas polarization distributions with 06f thresholds") +
  theme_lab(base_size = 9) +
  theme(legend.position = "bottom")

img_dens <- save_and_embed(p_dens, "polarization_density_by_atlas_class",
                            width = 1200, height = 600, res = 150)

# --- Fig C: atlas scatter SecA x SecB with polarization diagonals ------------

set.seed(CFG$seed)
atlas_samp <- atlas_sec[sample(.N, min(.N, 50000))]
p_scatter <- ggplot(atlas_samp, aes(x = SecA_UCell, y = SecB_UCell, color = celltype_nmf)) +
  geom_point(size = 0.25, alpha = 0.4) +
  geom_abline(slope = 1, intercept = t_low,  color = "#cc0000", linetype = "dashed",  linewidth = 0.6) +
  geom_abline(slope = 1, intercept = t_high, color = "#0066cc", linetype = "dotted", linewidth = 0.6) +
  scale_color_manual(values = class_pal, name = "Atlas NMF class") +
  labs(title = "Atlas SecA vs SecB UCell (sampled 50k) with polarization thresholds",
       subtitle = sprintf("Dashed red: SecB = SecA + %.3f (SecA/Intermediate border)\nDotted blue: SecB = SecA + %.3f (Intermediate/SecB border)",
                           t_low, t_high),
       x = "SecA UCell", y = "SecB UCell") +
  coord_fixed() +
  theme_lab(base_size = 9) +
  theme(legend.position = "bottom")

img_scatter <- save_and_embed(p_scatter, "threshold_visualization",
                                width = 900, height = 900, res = 150)

# --- Fig D: per-sample composition prior vs 06f ------------------------------

old_props <- xen_dt[, .N, by = .(sample, cell_label_06e)]
old_props[, pct := N / sum(N) * 100, by = sample]; old_props[, scheme := "prior"]
setnames(old_props, "cell_label_06e", "class")

new_props <- xen_dt[, .N, by = .(sample, cell_label_06f)]
new_props[, pct := N / sum(N) * 100, by = sample]; new_props[, scheme := "06f"]
setnames(new_props, "cell_label_06f", "class")

prop_dt <- rbind(old_props, new_props)
prop_dt[, class := factor(class, levels = sec_classes)]
prop_dt[, scheme := factor(scheme, levels = c("prior", "06f"))]

p_persample <- ggplot(prop_dt, aes(x = scheme, y = pct, fill = class)) +
  geom_col(width = 0.75, position = "stack") +
  scale_fill_manual(values = class_pal, name = "Class") +
  facet_wrap(~ sample, nrow = 1) +
  labs(x = NULL, y = "% of secretory cells",
       title = "Secretory composition per sample: prior vs 06f") +
  theme_lab(base_size = 7) +
  theme(axis.text.x = element_text(angle = 0, size = rel(0.9)),
        strip.text = element_text(size = rel(0.75)),
        legend.position = "bottom")

img_persample <- save_and_embed(p_persample, "xenium_reclassification_comparison",
                                  width = 1800, height = 550, res = 150)

# --- Fig E: per-sample % changed bar ---------------------------------------

stay_change <- xen_dt[, .(
  stayed  = sum(as.character(cell_label_06e) == as.character(cell_label_06f)),
  changed = sum(as.character(cell_label_06e) != as.character(cell_label_06f))
), by = sample]
stay_long <- melt(stay_change, id.vars = "sample",
                   variable.name = "status", value.name = "n")
stay_long[, pct := n / sum(n) * 100, by = sample]

p_change <- ggplot(stay_long, aes(x = sample, y = pct, fill = status)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(stayed = "grey70", changed = "#cc4444"),
                     name = "Status") +
  labs(x = NULL, y = "% of secretory cells",
       title = "Proportion of cells that changed classification (prior -> 06f)") +
  theme_lab(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

img_change <- save_and_embed(p_change, "per_sample_change_bar",
                              width = 1200, height = 500, res = 150)

# --- Fig F: prior vs 06f confusion flow heatmap ------------------------------

conf_flow <- xen_dt[, .N, by = .(cell_label_06e, cell_label_06f)]
conf_flow[, pct := 100 * N / sum(N)]
conf_flow[, lab := paste0(format(N, big.mark = ","), "\n(", round(pct, 1), "%)")]

p_conf <- ggplot(conf_flow,
                  aes(x = cell_label_06f, y = cell_label_06e, fill = N)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = lab), size = 3, lineheight = 0.85) +
  scale_fill_gradient(low = "#f7f7f7", high = "#0066cc", name = "n cells") +
  labs(x = "06f (polarization-based)",
       y = "prior (SecB-only)",
       title = "prior -> 06f cell flow",
       subtitle = sprintf("%s / %s cells changed (%.1f%%)",
                           format(n_changed, big.mark = ","),
                           format(n_total, big.mark = ","), pct_changed)) +
  coord_fixed() +
  theme_lab(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank())

img_conf <- save_and_embed(p_conf, "confusion_flow_06e_vs_06f",
                            width = 900, height = 750, res = 150)

# --- Fig G: atlas NMF agreement heatmap --------------------------------------

nmf_agree <- atlas_sec[, .N, by = .(celltype_nmf, threshold_class)]
nmf_agree[, pct := 100 * N / sum(N), by = celltype_nmf]
nmf_agree[, lab := paste0(format(N, big.mark = ","), "\n(", round(pct, 1), "%)")]

p_nmf <- ggplot(nmf_agree,
                 aes(x = threshold_class, y = celltype_nmf, fill = pct)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = lab), size = 3, lineheight = 0.85) +
  scale_fill_gradient(low = "#f7f7f7", high = "#228B22", name = "% of row") +
  labs(x = "06f threshold class (polarization)", y = "Atlas NMF class",
       title = "Atlas NMF vs 06f polarization-thresholded class",
       subtitle = sprintf("Overall accuracy: %.3f  |  macro-F1: %.3f",
                           overall_accuracy, macro_f1)) +
  coord_fixed() +
  theme_lab(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank())

img_nmf <- save_and_embed(p_nmf, "atlas_nmf_agreement_heatmap",
                            width = 900, height = 750, res = 150)

# --- Fig H: per-class score distributions (atlas dashed vs Xenium solid) ----

score_cols   <- c("SecA_UCell", "SecB_UCell", "polarization_UCell")
score_labels <- c("SecA_UCell" = "SecA UCell", "SecB_UCell" = "SecB UCell",
                   "polarization_UCell" = "Polarization UCell")

xen_long <- melt(xen_dt, id.vars = "cell_label_06f",
                  measure.vars = score_cols,
                  variable.name = "score", value.name = "value")
xen_long[, source := "Xenium (06f)"]; setnames(xen_long, "cell_label_06f", "class")

atl_long <- melt(atlas_sec[, c("celltype_nmf", score_cols), with = FALSE],
                  id.vars = "celltype_nmf",
                  measure.vars = score_cols,
                  variable.name = "score", value.name = "value")
atl_long[, source := "Atlas NMF"]; setnames(atl_long, "celltype_nmf", "class")

both <- rbind(xen_long, atl_long)
both[, class := factor(class, levels = sec_classes)]
both[, score := factor(score, levels = score_cols, labels = unname(score_labels))]

p_dist_list <- lapply(unname(score_labels), function(sc) {
  ggplot(both[score == sc], aes(x = value, fill = class, linetype = source)) +
    geom_density(alpha = 0.22, color = "grey30", linewidth = 0.4) +
    scale_fill_manual(values = class_pal, name = "Class") +
    scale_linetype_manual(values = c("Atlas NMF" = "dashed",
                                       "Xenium (06f)" = "solid"),
                            name = "Source") +
    facet_wrap(~ class, nrow = 1) +
    labs(x = sc, y = "Density", title = sc) +
    theme_lab(base_size = 7) +
    theme(legend.position = "bottom")
})
p_dist <- wrap_plots(p_dist_list, ncol = 1) + plot_layout(guides = "collect") &
           theme(legend.position = "bottom")

img_dist <- save_and_embed(p_dist, "score_distributions_06f",
                            width = 1500, height = 1200, res = 150)

# --- HTML tables -----------------------------------------------------------

tbl_thresh  <- make_html_table(as.data.frame(thresh_summary),   digits = 6)
tbl_atlas   <- make_html_table(as.data.frame(conf_atlas_wide))
tbl_per_cls <- make_html_table(as.data.frame(per_class_dt),     digits = 3)
tbl_conf06  <- make_html_table(as.data.frame(conf_06_wide))
tbl_biop    <- make_html_table(as.data.frame(biop),             digits = 2)

# --- HTML report ------------------------------------------------------------

html_css <- '
<style>
  body {font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;
    max-width:1200px;margin:0 auto;padding:30px 20px;background:white;color:#333;line-height:1.6;}
  h1 {border-bottom:2px solid #333;padding-bottom:10px;font-size:1.8em;}
  h2 {color:#444;border-bottom:1px solid #ddd;padding-bottom:5px;margin-top:40px;}
  h3 {color:#555;margin-top:25px;}
  .date {color:#888;font-size:0.9em;margin-bottom:20px;}
  .meta {background:#f8f8f8;padding:12px 18px;border-radius:5px;margin:15px 0;font-size:0.9em;}
  .toc  {background:#f8f8f8;padding:15px 25px;border-radius:5px;margin:20px 0;}
  .toc a{color:#0066cc;text-decoration:none;}
  .toc ul{list-style:none;padding-left:0;} .toc li{margin:6px 0;}
  img {border:1px solid #eee;border-radius:3px;margin:10px 0;}
  .note {background:#f0f7ff;border-left:3px solid #0066cc;padding:10px 15px;margin:15px 0;font-size:0.9em;}
  .result {background:#f0fff0;border-left:3px solid #228B22;padding:10px 15px;margin:15px 0;font-size:0.9em;}
  .styled-table{border-collapse:collapse;width:100%%;font-size:0.85em;margin:15px 0;}
  .styled-table th{background:#f4f4f4;text-align:left;padding:8px 12px;border-bottom:2px solid #ddd;}
  .styled-table td{padding:6px 12px;border-bottom:1px solid #eee;}
</style>
'

html_toc <- '
<div class="toc"><strong>Table of Contents</strong>
<ul>
  <li><a href="#s1">1. Thresholds + atlas ROC (polarization)</a></li>
  <li><a href="#s2">2. Atlas NMF agreement</a></li>
  <li><a href="#s3">3. Threshold visualization</a></li>
  <li><a href="#s4">4. Xenium reclassification (prior -> 06f)</a></li>
  <li><a href="#s5">5. Score distributions by new class</a></li>
  <li><a href="#s6">6. Per-sample change rates</a></li>
  <li><a href="#s7">7. Biological plausibility</a></li>
</ul></div>
'

meta_block <- paste0(
  '<div class="meta">',
  '<strong>Strategy (c):</strong> t_low = atlas SecA polarization p75 = ', sprintf("%.4f", t_low),
  '  |  t_high = atlas SecB polarization p25 = ', sprintf("%.4f", t_high), '<br>',
  '<strong>Classification rule:</strong> SecA (pol &lt; ', sprintf("%.4f", t_low),
  '), Intermediate (', sprintf("%.4f", t_low), ' &le; pol &lt; ', sprintf("%.4f", t_high),
  '), SecB (pol &ge; ', sprintf("%.4f", t_high), ')<br>',
  '<strong>AUC:</strong> polarization = ', sprintf("%.4f", auc_val),
  '  (06e SecB-only baseline = 0.9365)<br>',
  '<strong>Atlas agreement:</strong> accuracy = ', sprintf("%.3f", overall_accuracy),
  ', macro-F1 = ', sprintf("%.3f", macro_f1), '<br>',
  '<strong>Xenium cells (CSV scope) changed prior -> 06f:</strong> ',
  format(n_changed, big.mark = ","), ' / ', format(n_total, big.mark = ","),
  '  (', pct_changed, '%)<br>',
  '<strong>SFE writeback:</strong> ', length(sfe_names),
  ' SFEs updated (cell_label overwritten; cell_label_06e preserves prior labels).',
  '</div>'
)

html_body <- paste0(
  '<!DOCTYPE html>\n<html>\n<head>\n<meta charset="utf-8">\n',
  '<title>06f: Atlas-aligned Polarization Reclassification</title>\n',
  html_css,
  '\n</head>\n<body>\n',
  '<h1>06f: Atlas-aligned Polarization Reclassification</h1>\n',
  '<p class="date">Generated: ', format(Sys.time(), "%Y-%m-%d %H:%M"), '</p>\n',
  '<p>06f replaces the SecB-only decision axis with the bivariate ',
  '<code>polarization_UCell = SecB_UCell - SecA_UCell</code> score, and calibrates ',
  'thresholds against the atlas three-class NMF labels using the symmetric ',
  'SecA p75 + SecB p25 anchors (Strategy c).</p>\n',
  meta_block, '\n', html_toc,

  '\n<h2 id="s1">1. Thresholds + Atlas ROC</h2>\n',
  '<p>ROC of polarization_UCell on the atlas binary SecA vs SecB set. ',
  'Strategy c anchors are shown in the threshold summary.</p>\n',
  '<h3>Threshold summary</h3>\n', tbl_thresh, '\n',
  img_roc, '\n',

  '\n<h2 id="s2">2. Atlas NMF Agreement</h2>\n',
  '<p>Applying the 06f polarization rule to atlas secretory cells and ',
  'comparing to their NMF labels.</p>\n',
  '<h3>Confusion matrix (rows = atlas NMF, cols = 06f class)</h3>\n', tbl_atlas, '\n',
  '<h3>Per-class precision / recall / F1</h3>\n', tbl_per_cls, '\n',
  img_nmf, '\n',

  '\n<h2 id="s3">3. Threshold Visualization</h2>\n',
  '<p>Atlas polarization distributions coloured by NMF class, with the two ',
  'Strategy c thresholds overlaid. The scatter on the right shows the same ',
  'thresholds as parallel diagonals in the SecA-SecB UCell plane.</p>\n',
  img_dens, '\n', img_scatter, '\n',

  '\n<h2 id="s4">4. Xenium Reclassification (prior -> 06f)</h2>\n',
  '<p>Per-sample stacked bars show secretory composition under the prior and 06f ',
  'schemes; the confusion heatmap shows how cells flow between classes.</p>\n',
  img_persample, '\n',
  '<h3>Confusion matrix (overall)</h3>\n', tbl_conf06, '\n',
  img_conf, '\n',

  '\n<h2 id="s5">5. Score Distributions by New Class</h2>\n',
  '<p>Each facet column is a 06f class; dashed = atlas NMF reference, solid = Xenium (06f). ',
  'Good alignment means dashed + solid curves overlap.</p>\n',
  img_dist, '\n',

  '\n<h2 id="s6">6. Per-Sample Change Rates</h2>\n',
  img_change, '\n',

  '\n<h2 id="s7">7. Biological Plausibility</h2>\n',
  '<p>Flags: <code>flag_low_SecB</code> if SecB% &lt; 5 of secretory; ',
  '<code>flag_high_SecB</code> if &gt; 30; polarization sign flags indicate ',
  'whether SecA and SecB class medians lie on the expected sides of 0.</p>\n',
  tbl_biop, '\n',

  '\n<hr>\n<p style="color:#888;font-size:0.8em;">',
  'Report generated by spatial/03_annotation_polarization/04_reclassification_polarization.R on ',
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '</p>\n',
  '</body>\n</html>'
)

writeLines(html_body, html_out)
message("[", Sys.time(), "] HTML report written: ", html_out)

# --- Session info -----------------------------------------------------------

log_session()
message("\n[", Sys.time(), "] 06f reclassification complete.")
