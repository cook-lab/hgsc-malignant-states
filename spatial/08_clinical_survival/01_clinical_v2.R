# ============================================================================
# 01_clinical_v2.R
# ----------------------------------------------------------------------------
# PURPOSE: Per-patient clinical feature construction (epithelial densities, SecA/SecB proportions, ratios) and Cox univariate survival association.
#
# INPUTS:
#   - SFEs (load_sfe) with cell_label + clinical_data_clean.csv
#
# OUTPUTS:
#   - output/10_clinical_v2/per_patient_features_v2.csv
#   - output/10_clinical_v2/cox_univariate_results.csv
#
# MANUSCRIPT PANEL(S): Fig 4J, Fig 7A, Fig 7B.
# RUNTIME TIER: heavy
#
# Migrated from 2026_final_xenium_analysis/scripts/. Analytical logic preserved;
# paths routed through central config, seed from CFG$seed, epithelial label
# "Transitioning" -> "Intermediate", SecA/SecB from shared/signatures.yml.
# ============================================================================

# --- Config + shared setup (replaces hardcoded /Volumes/CookLab/Sarah paths) ---
here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
source(file.path(here, "..", "..", "config", "config.R"))   # CFG, cfg_obj, cfg_path
source(file.path(here, "..", "00_setup", "00_setup.R"))      # load_sfe, save_sfe, theme_lab, nb_names, palettes
set.seed(CFG$seed)

library(survival)
library(broom)
library(base64enc)
library(ragg)
library(circlize)

message("\n=== 10_clinical_v2: Per-Patient Clinical Associations (Reclassified) ===")
message("[", Sys.time(), "] Starting...")

# --- Paths ------------------------------------------------------------------

step_dir <- file.path(out_dir, "10_clinical_v2")
fig_out  <- file.path(step_dir, "figures")
html_out <- file.path(step_dir, "10_clinical_v2_report.html")

for (d in c(step_dir, fig_out)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# --- Constants ---------------------------------------------------------------

epi_types <- c("Ciliated epithelium", "SecA epithelium",
               "Intermediate epithelium", "SecB epithelium")
immune_types <- c("T cell", "NK cell", "B cell", "Plasma cell",
                  "Macrophage", "Conventional dendritic cell",
                  "Plasmacytoid dendritic cell", "Neutrophil", "Mast cell")
stromal_types <- c("Fibroblast", "Smooth muscle", "Mesothelial",
                   "Pericyte", "Endothelial")
celltype_order <- c(epi_types, stromal_types, immune_types)

sec_classes <- c("SecA epithelium", "Intermediate epithelium", "SecB epithelium")
sec_pal <- c("SecA" = "#E6A141", "Intermediate" = "#C08E48", "SecB" = "#9A7D55")
km_pal <- c("Low" = "#0072B2", "High" = "#D55E00")
chemo_pal <- c("Chemo-treated" = "#D55E00", "Treatment-naive" = "#0072B2")
sample_type_pal <- c("tumour" = "#D55E00", "FT" = "#56B4E9")

# Chemo-treated samples
chemo_samples <- c("sfe_OTB_2417", "sfe_OTB_2432", "sfe_OTB_2457")

# SFE names: sfe_tma_filtered + 8 WT
sfe_wt_names <- sfe_names_wt

# Pathway column prefix
pathway_prefix <- "pathway_"

# --- Helper: save plot as PNG + PDF, return base64 img tag ------------------

save_and_embed <- function(p, name, width = 1200, height = 600, res = 150) {
  # PDF
  pdf_path <- file.path(fig_out, paste0(name, ".pdf"))
  ggsave(pdf_path, plot = p, width = width / res, height = height / res,
         device = cairo_pdf)
  message("  Saved PDF: ", pdf_path)

  # PNG for HTML embedding
  tmp <- tempfile(fileext = ".png")
  ragg::agg_png(tmp, width = width, height = height, res = res)
  print(p)
  dev.off()
  b64 <- base64enc::base64encode(tmp)
  unlink(tmp)
  sprintf('<img src="data:image/png;base64,%s" style="max-width:100%%;">', b64)
}

# Helper: embed ComplexHeatmap
save_and_embed_heatmap <- function(hm, name, width = 1200, height = 800, res = 150) {
  pdf_path <- file.path(fig_out, paste0(name, ".pdf"))
  pdf(pdf_path, width = width / res, height = height / res)
  draw(hm)
  dev.off()
  message("  Saved PDF: ", pdf_path)

  tmp <- tempfile(fileext = ".png")
  ragg::agg_png(tmp, width = width, height = height, res = res)
  draw(hm)
  dev.off()
  b64 <- base64enc::base64encode(tmp)
  unlink(tmp)
  sprintf('<img src="data:image/png;base64,%s" style="max-width:100%%;">', b64)
}

# Helper: styled HTML table
make_html_table <- function(df, digits = 3) {
  for (j in seq_along(df)) {
    if (is.numeric(df[[j]])) df[[j]] <- round(df[[j]], digits)
  }
  knitr::kable(df, format = "html", row.names = FALSE,
               table.attr = 'class="styled-table"')
}


# ============================================================================
# PART 1: Per-Patient Feature Extraction
# ============================================================================

message("\n[", Sys.time(), "] === PART 1: Per-Patient Feature Extraction ===")

# --- 1a. Load neighborhood assignments -------------------------------------

message("[", Sys.time(), "] Loading neighborhood assignments...")
nb_assign_all <- fread(file.path(out_dir, "09_neighborhood",
                                  "neighborhood_assignments.csv"))
nb_levels <- paste0("nb_", 1:10)

# --- 1b. Extract features from TMA -----------------------------------------

message("[", Sys.time(), "] Processing sfe_tma_filtered...")
sfe_tma <- load_sfe("sfe_tma_filtered")
cd_tma <- as.data.table(as.data.frame(colData(sfe_tma)))
cd_tma[, cell_id := colnames(sfe_tma)]
# Idempotent rename: SFE colData still carries the legacy epithelial label
# "Transitioning epithelium"; standardize to "Intermediate epithelium" before
# any cell_label match (epi_types/sec_classes) so the Intermediate epitype is
# not silently dropped from densities/proportions.
cd_tma[cell_label == "Transitioning epithelium", cell_label := "Intermediate epithelium"]
rm(sfe_tma); gc(verbose = FALSE)

# Identify pathway columns and polarization_UCell
all_cd_cols <- names(cd_tma)
pathway_cols <- grep(paste0("^", pathway_prefix), all_cd_cols, value = TRUE)
has_polarization <- "polarization_UCell" %in% all_cd_cols

message("  Pathway columns found: ", length(pathway_cols))
message("  polarization_UCell present: ", has_polarization)

# Keep tumour cores with valid patient_id
cd_tma[, patient_id := as.character(patient_id)]
cd_tma[, core_id := as.character(core_id)]
tma_tumour <- cd_tma[sample_type == "tumour" & !is.na(core_id) & !is.na(patient_id)]
message("  Tumour cells: ", format(nrow(tma_tumour), big.mark = ","),
        " | Patients: ", length(unique(tma_tumour$patient_id)))

# Merge neighborhoods (drop existing if present to avoid .x/.y collision)
nb_tma <- nb_assign_all[sample_id == "sfe_tma_filtered"]
if ("neighborhood" %in% names(tma_tumour)) tma_tumour[, neighborhood := NULL]
tma_tumour <- merge(tma_tumour, nb_tma[, .(cell_id, neighborhood)],
                    by = "cell_id", all.x = TRUE)

# --- Per-core cell counts and densities ---
core_area <- tma_tumour[, .(total_area_mm2 = sum(cell_area, na.rm = TRUE) / 1e6),
                        by = core_id]

# Cell type counts per core
core_counts <- dcast(tma_tumour, core_id ~ cell_label,
                     fun.aggregate = length, value.var = "cell_id", fill = 0)

# Density per core per type
core_density <- merge(core_counts, core_area, by = "core_id")
dens_cols_present <- intersect(celltype_order, names(core_density))
for (ct in dens_cols_present) {
  set(core_density, j = ct,
      value = core_density[[ct]] / core_density$total_area_mm2)
}

# Map core -> patient
core_patient <- unique(tma_tumour[, .(core_id, patient_id)])
core_density <- merge(core_density, core_patient, by = "core_id")

# Average across cores per patient
patient_density_tma <- core_density[,
  lapply(.SD, mean, na.rm = TRUE),
  by = patient_id,
  .SDcols = dens_cols_present
]
setnames(patient_density_tma, dens_cols_present,
         paste0("dens_", gsub(" ", "_", dens_cols_present)))

# --- Per-core neighborhood proportions ---
core_n <- tma_tumour[, .(n_total = .N), by = core_id]
nb_dt <- tma_tumour[!is.na(neighborhood) & neighborhood != "nb_unassigned"]
core_nb <- dcast(nb_dt, core_id ~ neighborhood,
                 fun.aggregate = length, value.var = "cell_id", fill = 0)
core_nb <- merge(core_nb, core_n, by = "core_id")

nb_present <- intersect(nb_levels, names(core_nb))
for (nb in nb_present) {
  set(core_nb, j = nb, value = core_nb[[nb]] / core_nb$n_total)
}
core_nb <- merge(core_nb, core_patient, by = "core_id")
patient_nb_tma <- core_nb[,
  lapply(.SD, mean, na.rm = TRUE),
  by = patient_id,
  .SDcols = nb_present
]
setnames(patient_nb_tma, nb_present, paste0("prop_", nb_present))

# --- Per-patient polarization_UCell stats (TMA) ---
polar_tma <- NULL
if (has_polarization) {
  polar_tma <- tma_tumour[!is.na(polarization_UCell), .(
    polar_mean   = mean(polarization_UCell),
    polar_median = median(polarization_UCell),
    polar_sd     = sd(polarization_UCell),
    polar_iqr    = IQR(polarization_UCell),
    polar_n      = .N
  ), by = patient_id]
}

# --- Per-patient pathway score means (TMA) ---
pw_tma <- NULL
if (length(pathway_cols) > 0) {
  pw_tma <- tma_tumour[, lapply(.SD, mean, na.rm = TRUE),
                       by = patient_id, .SDcols = pathway_cols]
}

# --- Per-patient composition proportions (TMA) ---
tma_comp <- tma_tumour[, .N, by = .(patient_id, cell_label)]
tma_total <- tma_tumour[, .(total = .N), by = patient_id]
tma_comp <- merge(tma_comp, tma_total, by = "patient_id")
tma_comp[, prop := N / total]

# Secretory proportions
sec_tma <- tma_tumour[cell_label %in% sec_classes, .N, by = .(patient_id, cell_label)]
sec_total_tma <- tma_tumour[cell_label %in% sec_classes, .(sec_total = .N), by = patient_id]
sec_tma <- merge(sec_tma, sec_total_tma, by = "patient_id")
sec_tma[, sec_prop := N / sec_total]
sec_wide_tma <- dcast(sec_tma, patient_id ~ cell_label, value.var = "sec_prop", fill = 0)
setnames(sec_wide_tma, sec_classes,
         c("prop_SecA_of_sec", "prop_Int_of_sec", "prop_SecB_of_sec"))

# Total counts and lineage proportions
lineage_tma <- tma_tumour[, .(
  total_cells = .N,
  n_epi       = sum(cell_label %in% epi_types),
  n_immune    = sum(cell_label %in% immune_types),
  n_stromal   = sum(cell_label %in% stromal_types)
), by = patient_id]
lineage_tma[, `:=`(
  prop_epi     = n_epi / total_cells,
  prop_immune  = n_immune / total_cells,
  prop_stromal = n_stromal / total_cells
)]

# Merge TMA features
tma_features <- patient_density_tma
tma_features <- merge(tma_features, patient_nb_tma, by = "patient_id", all.x = TRUE)
tma_features <- merge(tma_features, sec_wide_tma, by = "patient_id", all.x = TRUE)
tma_features <- merge(tma_features, lineage_tma[, .(patient_id, total_cells,
                                                     prop_epi, prop_immune, prop_stromal)],
                      by = "patient_id", all.x = TRUE)
if (!is.null(polar_tma)) {
  tma_features <- merge(tma_features, polar_tma, by = "patient_id", all.x = TRUE)
}
if (!is.null(pw_tma)) {
  tma_features <- merge(tma_features, pw_tma, by = "patient_id", all.x = TRUE)
}

tma_features[, source := "TMA"]
message("  TMA: ", nrow(tma_features), " patients extracted")

rm(cd_tma, tma_tumour, nb_dt, core_density, core_nb); gc(verbose = FALSE)

# --- 1c. Extract features from whole-tissue samples ------------------------

message("[", Sys.time(), "] Processing whole-tissue samples...")

wt_features_list <- list()

for (sname in sfe_wt_names) {
  message("  Loading ", sname, "...")
  sfe <- load_sfe(sname)
  cd <- as.data.table(as.data.frame(colData(sfe)))
  cd[, cell_id := colnames(sfe)]
  # Idempotent rename: SFE colData still carries the legacy epithelial label
  # "Transitioning epithelium"; standardize to "Intermediate epithelium" before
  # any cell_label match (epi_types/sec_classes) so the Intermediate epitype is
  # not silently dropped from densities/proportions.
  cd[cell_label == "Transitioning epithelium", cell_label := "Intermediate epithelium"]
  rm(sfe); gc(verbose = FALSE)

  # Detect available columns for this sample
  wt_pathway_cols <- grep(paste0("^", pathway_prefix), names(cd), value = TRUE)
  wt_has_polar <- "polarization_UCell" %in% names(cd)

  # Patient ID = sample name (strip sfe_ prefix)
  pid <- gsub("^sfe_", "", sname)
  cd[, patient_id := pid]

  # Merge neighborhoods (drop existing if present)
  nb_wt <- nb_assign_all[sample_id == sname]
  if ("neighborhood" %in% names(cd)) cd[, neighborhood := NULL]
  cd <- merge(cd, nb_wt[, .(cell_id, neighborhood)],
              by = "cell_id", all.x = TRUE)

  n_cells <- nrow(cd)

  # --- Cell type densities (cells/mm2 using cell_area) ---
  total_area_mm2 <- sum(cd$cell_area, na.rm = TRUE) / 1e6
  ct_counts <- cd[, .N, by = cell_label]

  dens_dt <- data.table(patient_id = pid)
  for (ct in celltype_order) {
    ct_n <- ct_counts[cell_label == ct, N]
    if (length(ct_n) == 0) ct_n <- 0
    col_nm <- paste0("dens_", gsub(" ", "_", ct))
    set(dens_dt, j = col_nm, value = ct_n / total_area_mm2)
  }

  # --- Neighborhood proportions ---
  nb_wt_counts <- cd[!is.na(neighborhood) & neighborhood != "nb_unassigned",
                     .N, by = neighborhood]
  nb_total <- sum(nb_wt_counts$N)
  for (nb in nb_levels) {
    nb_n <- nb_wt_counts[neighborhood == nb, N]
    if (length(nb_n) == 0) nb_n <- 0
    set(dens_dt, j = paste0("prop_", nb),
        value = if (nb_total > 0) nb_n / nb_total else 0)
  }

  # --- Secretory proportions ---
  sec_cd <- cd[cell_label %in% sec_classes]
  n_sec <- nrow(sec_cd)
  for (sc in sec_classes) {
    nm_short <- c("SecA epithelium" = "prop_SecA_of_sec",
                  "Intermediate epithelium" = "prop_Int_of_sec",
                  "SecB epithelium" = "prop_SecB_of_sec")[sc]
    set(dens_dt, j = nm_short,
        value = if (n_sec > 0) sum(sec_cd$cell_label == sc) / n_sec else 0)
  }

  # --- Lineage proportions ---
  set(dens_dt, j = "total_cells", value = n_cells)
  set(dens_dt, j = "prop_epi",
      value = sum(cd$cell_label %in% epi_types) / n_cells)
  set(dens_dt, j = "prop_immune",
      value = sum(cd$cell_label %in% immune_types) / n_cells)
  set(dens_dt, j = "prop_stromal",
      value = sum(cd$cell_label %in% stromal_types) / n_cells)

  # --- Polarization UCell ---
  if (wt_has_polar) {
    polar_vals <- cd[!is.na(polarization_UCell), polarization_UCell]
    if (length(polar_vals) > 0) {
      set(dens_dt, j = "polar_mean",   value = mean(polar_vals))
      set(dens_dt, j = "polar_median", value = median(polar_vals))
      set(dens_dt, j = "polar_sd",     value = sd(polar_vals))
      set(dens_dt, j = "polar_iqr",    value = IQR(polar_vals))
      set(dens_dt, j = "polar_n",      value = length(polar_vals))
    }
  }

  # --- Pathway means ---
  if (length(wt_pathway_cols) > 0) {
    for (pc in wt_pathway_cols) {
      vals <- cd[[pc]]
      set(dens_dt, j = pc, value = mean(vals, na.rm = TRUE))
    }
  }

  # Treatment status
  dens_dt[, treatment_status := fifelse(sname %in% chemo_samples,
                                        "Chemo-treated", "Treatment-naive")]
  dens_dt[, source := "WT"]

  wt_features_list[[sname]] <- dens_dt

  rm(cd, sec_cd); gc(verbose = FALSE)
  message("    ", pid, ": ", format(n_cells, big.mark = ","), " cells")
}

wt_features <- rbindlist(wt_features_list, fill = TRUE)
message("  WT: ", nrow(wt_features), " samples extracted")

# --- 1d. Load clinical data and merge with TMA features ---------------------

message("[", Sys.time(), "] Merging with clinical data...")
clinical <- fread(file.path(data_dir, "clinical_data_clean.csv"))
clinical[, patient_id := as.character(patient_id)]

# TMA patients: merge with clinical
tma_merged <- merge(tma_features, clinical, by = "patient_id", all.x = TRUE)

# Derive treatment_status from clinical data for TMA patients
tma_merged[, treatment_status := fifelse(
  !is.na(treatment_status) & treatment_status != "",
  treatment_status,
  "Treatment-naive"  # default for TMA if not annotated
)]

# 5-year censoring
tma_merged[!is.na(survival_months), `:=`(
  os_time_5y  = pmin(survival_months, 60),
  os_event_5y = fifelse(survival_months > 60, 0L, as.integer(survival_outcome))
)]
tma_merged[!is.na(progression_months), `:=`(
  pfs_time_5y  = pmin(progression_months, 60),
  pfs_event_5y = fifelse(progression_months > 60, 0L, as.integer(progression_outcome))
)]

# --- 1e. Compute derived features ---

# SecB:SecA ratio
eps_ratio <- 0.1

compute_ratio <- function(dt) {
  dens_seca <- dt[["dens_SecA_epithelium"]]
  dens_secb <- dt[["dens_SecB_epithelium"]]
  if (is.null(dens_seca) || is.null(dens_secb)) return(dt)
  dt[, ratio_SecA_SecB := (dens_SecA_epithelium + eps_ratio) /
                           (dens_SecB_epithelium + eps_ratio)]
  dt[, log2_ratio := log2(ratio_SecA_SecB)]
  dt
}

tma_merged <- compute_ratio(tma_merged)
wt_features <- compute_ratio(wt_features)

# --- 1f. Combine all patients ------------------------------------------------

per_patient <- rbindlist(list(tma_merged, wt_features), fill = TRUE)

message("  Combined: ", nrow(per_patient), " patients (",
        sum(per_patient$source == "TMA"), " TMA + ",
        sum(per_patient$source == "WT"), " WT)")

# Save
fwrite(per_patient, file.path(step_dir, "per_patient_features_v2.csv"))
message("  Saved: per_patient_features_v2.csv")


# ============================================================================
# PART 2: Composition Analysis
# ============================================================================

message("\n[", Sys.time(), "] === PART 2: Composition Analysis ===")

html_sections <- list()

# --- 2a. Stacked bar: cell type proportions per patient ---------------------

message("[", Sys.time(), "] Generating composition bar plot...")

dens_cols_all <- grep("^dens_", names(per_patient), value = TRUE)

# Compute proportions from densities
prop_mat <- per_patient[, ..dens_cols_all]
prop_mat <- prop_mat / rowSums(prop_mat, na.rm = TRUE)
prop_mat[, patient_id := per_patient$patient_id]
prop_mat[, source := per_patient$source]

# Order patients by SecB proportion
prop_mat[, secb_prop := `dens_SecB_epithelium` / rowSums(.SD, na.rm = TRUE),
         .SDcols = dens_cols_all]
setorder(prop_mat, -secb_prop)
pat_order <- prop_mat$patient_id

bar_long <- melt(prop_mat, id.vars = c("patient_id", "source", "secb_prop"),
                 variable.name = "cell_type", value.name = "proportion")
bar_long[, cell_type := gsub("^dens_", "", cell_type)]
bar_long[, cell_type := gsub("_", " ", cell_type)]
bar_long[, patient_id := factor(patient_id, levels = pat_order)]
bar_long[, cell_type := factor(cell_type, levels = rev(celltype_order))]

p_bar <- ggplot(bar_long[!is.na(proportion)],
                aes(x = patient_id, y = proportion, fill = cell_type)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = ref_palette, name = "Cell type") +
  labs(x = "Patient", y = "Proportion", title = "Cell type composition per patient") +
  theme_lab(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
        legend.key.size = unit(0.5, "lines"))

img_bar <- save_and_embed(p_bar, "composition_stacked_bar",
                          width = 1600, height = 600)

# --- 2b. Heatmap: patient x cell type proportions (z-scored) ---------------

message("[", Sys.time(), "] Generating composition heatmap...")

hm_dt <- per_patient[, c("patient_id", dens_cols_all), with = FALSE]
hm_mat <- as.matrix(hm_dt[, -1])
rownames(hm_mat) <- hm_dt$patient_id

# Convert to proportions then z-score by column
hm_props <- sweep(hm_mat, 1, rowSums(hm_mat, na.rm = TRUE), "/")
hm_props[is.nan(hm_props)] <- 0
hm_z <- scale(hm_props)
hm_z[is.nan(hm_z)] <- 0

# Clean column names
colnames(hm_z) <- gsub("^dens_", "", colnames(hm_z))
colnames(hm_z) <- gsub("_", " ", colnames(hm_z))

# Annotation
ha_row <- ComplexHeatmap::rowAnnotation(
  Source = per_patient$source,
  Treatment = per_patient$treatment_status,
  col = list(
    Source    = c("TMA" = "#999999", "WT" = "#333333"),
    Treatment = c("Chemo-treated" = "#D55E00", "Treatment-naive" = "#0072B2",
                  "naive" = "#0072B2", "post-chemotherapy" = "#D55E00")
  ),
  show_legend = TRUE,
  simple_anno_size = unit(3, "mm")
)

col_fun <- colorRamp2(c(-2, 0, 2), c("#0072B2", "white", "#D55E00"))

hm <- ComplexHeatmap::Heatmap(
  hm_z,
  name = "z-score",
  col = col_fun,
  row_split = per_patient$source,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 7),
  left_annotation = ha_row,
  use_raster = FALSE,
  row_title_gp = gpar(fontsize = 8),
  column_title = "Cell type composition (z-scored proportions)"
)

img_hm <- save_and_embed_heatmap(hm, "composition_heatmap",
                                  width = 1400, height = 1000)

# Save composition CSV
fwrite(as.data.table(hm_props, keep.rownames = "patient_id"),
       file.path(step_dir, "composition_by_patient.csv"))

html_sections[["composition"]] <- paste0(
  '<h2>2. Cell Type Composition</h2>\n',
  '<h3>2a. Stacked bar (ordered by SecB proportion)</h3>\n',
  img_bar, '\n',
  '<h3>2b. Z-scored composition heatmap with patient clustering</h3>\n',
  img_hm, '\n'
)

message("  Composition analysis complete")


# ============================================================================
# PART 3: Clinical Associations (TMA patients with clinical data)
# ============================================================================

message("\n[", Sys.time(), "] === PART 3: Clinical Associations ===")

tma_clin <- per_patient[source == "TMA" & !is.na(age)]
message("  TMA patients with clinical data: ", nrow(tma_clin))

html_clin <- ""

# --- 3a. Tumour vs FT comparison -------------------------------------------
# Check if both sample types exist (FT cores vs tumour cores)
# In this analysis, all TMA patients are already filtered to tumour cores.
# If colData has sample_type, we could compare. For now, skip unless data present.

html_tumour_ft <- '<h3>3a. Tumour vs FT comparison</h3>\n'
if ("sample_type" %in% names(per_patient) &&
    length(unique(per_patient$sample_type[!is.na(per_patient$sample_type)])) > 1) {
  html_tumour_ft <- paste0(html_tumour_ft,
    '<p>Multiple sample types detected; comparison available in per_patient_features_v2.csv.</p>\n')
} else {
  html_tumour_ft <- paste0(html_tumour_ft,
    '<p>Only tumour cores included in TMA extraction. No FT comparison available.</p>\n')
}

# --- 3b. SecB proportion by treatment status --------------------------------

html_chemo <- '<h3>3b. Composition by treatment status</h3>\n'

# For WT samples, compare chemo vs naive
wt_clin <- per_patient[source == "WT"]
if (nrow(wt_clin) > 2) {
  p_chemo_secb <- ggplot(wt_clin, aes(x = treatment_status, y = prop_SecB_of_sec,
                                       fill = treatment_status)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.3, alpha = 0.7) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
    scale_fill_manual(values = chemo_pal, guide = "none") +
    labs(x = NULL, y = "SecB proportion (of secretory)",
         title = "SecB proportion: Chemo-treated vs Treatment-naive (WT samples)") +
    theme_lab()

  # Wilcoxon (small N, so report with caveat)
  wt_chemo <- wt_clin[treatment_status == "Chemo-treated", prop_SecB_of_sec]
  wt_naive <- wt_clin[treatment_status == "Treatment-naive", prop_SecB_of_sec]
  wt_pval <- tryCatch(
    wilcox.test(wt_chemo, wt_naive, exact = FALSE)$p.value,
    error = function(e) NA
  )
  p_chemo_secb <- p_chemo_secb +
    labs(subtitle = sprintf("Wilcoxon p = %s (N = %d vs %d)",
                            ifelse(is.na(wt_pval), "NA", signif(wt_pval, 3)),
                            length(wt_chemo), length(wt_naive)))

  img_chemo <- save_and_embed(p_chemo_secb, "chemo_secb_wt",
                               width = 600, height = 500)
  html_chemo <- paste0(html_chemo, img_chemo, '\n')
}

# --- 3c. Univariate tests: features vs sample_type for TMA -----------------

message("[", Sys.time(), "] Running univariate Wilcoxon tests (TMA)...")

# Define features for testing
dens_feats <- grep("^dens_", names(tma_clin), value = TRUE)
nb_feats <- grep("^prop_nb_", names(tma_clin), value = TRUE)
ratio_feats <- intersect(c("log2_ratio", "prop_SecA_of_sec",
                            "prop_SecB_of_sec", "prop_Int_of_sec",
                            "prop_epi", "prop_immune", "prop_stromal"),
                         names(tma_clin))
polar_feats <- intersect(c("polar_mean", "polar_median", "polar_sd", "polar_iqr"),
                         names(tma_clin))
test_features <- c(dens_feats, nb_feats, ratio_feats, polar_feats)

# Binary outcomes for testing
binary_outcomes <- list()
if ("residual_binary" %in% names(tma_clin) &&
    sum(!is.na(tma_clin$residual_binary)) >= 10) {
  binary_outcomes[["residual_binary"]] <- "Debulking (optimal vs suboptimal)"
}
if ("chemo_response_12mo" %in% names(tma_clin) &&
    sum(!is.na(tma_clin$chemo_response_12mo)) >= 10) {
  binary_outcomes[["chemo_response_12mo"]] <- "Chemo response 12mo (1=sens, 0=resist)"
}

univar_results <- list()

for (outcome_name in names(binary_outcomes)) {
  for (feat in test_features) {
    dat_test <- tma_clin[!is.na(get(feat)) & !is.na(get(outcome_name))]
    n_total <- nrow(dat_test)

    if (n_total < 10) next

    groups <- split(dat_test[[feat]], dat_test[[outcome_name]])
    if (length(groups) < 2 || any(sapply(groups, length) < 3)) next

    wt <- tryCatch(
      wilcox.test(dat_test[[feat]] ~ dat_test[[outcome_name]], exact = FALSE),
      error = function(e) NULL
    )

    if (!is.null(wt)) {
      univar_results[[length(univar_results) + 1]] <- data.table(
        outcome     = outcome_name,
        outcome_label = binary_outcomes[[outcome_name]],
        feature     = feat,
        n           = n_total,
        p_value     = wt$p.value,
        statistic   = wt$statistic
      )
    }
  }
}

univar_dt <- rbindlist(univar_results)
if (nrow(univar_dt) > 0) {
  univar_dt[, padj := p.adjust(p_value, method = "BH"), by = outcome]
  setorder(univar_dt, outcome, p_value)
  fwrite(univar_dt, file.path(step_dir, "univariate_tests.csv"))
  message("  Saved: univariate_tests.csv (",
          nrow(univar_dt[p_value < 0.05]), " nominal p < 0.05)")
}

# --- 3d. Univariate Cox PH (if survival data available) --------------------

message("[", Sys.time(), "] Running univariate Cox PH...")

cox_results <- list()
surv_endpoints <- list()

if (sum(!is.na(tma_clin$os_time_5y) & !is.na(tma_clin$os_event_5y)) >= 15) {
  surv_endpoints[["OS_5yr"]] <- c("os_time_5y", "os_event_5y")
}
if (sum(!is.na(tma_clin$pfs_time_5y) & !is.na(tma_clin$pfs_event_5y)) >= 15) {
  surv_endpoints[["PFS_5yr"]] <- c("pfs_time_5y", "pfs_event_5y")
}

for (ep_name in names(surv_endpoints)) {
  time_var  <- surv_endpoints[[ep_name]][1]
  event_var <- surv_endpoints[[ep_name]][2]

  for (feat in test_features) {
    dat_surv <- tma_clin[!is.na(get(feat)) & !is.na(get(time_var)) &
                          !is.na(get(event_var)) & get(time_var) > 0]
    if (nrow(dat_surv) < 15) next

    fit <- tryCatch({
      dat_surv[, x_scaled := as.numeric(scale(get(feat)))]
      coxph(Surv(get(time_var), get(event_var)) ~ x_scaled, data = dat_surv)
    }, error = function(e) NULL)

    if (!is.null(fit)) {
      s <- summary(fit)
      cox_results[[length(cox_results) + 1]] <- data.table(
        endpoint  = ep_name,
        feature   = feat,
        n         = s$n,
        n_events  = s$nevent,
        HR        = s$conf.int[1, 1],
        HR_lower  = s$conf.int[1, 3],
        HR_upper  = s$conf.int[1, 4],
        p_value   = s$coefficients[1, 5]
      )
    }
  }
}

cox_dt <- rbindlist(cox_results)
if (nrow(cox_dt) > 0) {
  cox_dt[, padj := p.adjust(p_value, method = "BH"), by = endpoint]
  setorder(cox_dt, endpoint, p_value)
  fwrite(cox_dt, file.path(step_dir, "cox_univariate_results.csv"))
  message("  Saved: cox_univariate_results.csv (",
          nrow(cox_dt[p_value < 0.05]), " nominal p < 0.05)")
}

# Build HTML for clinical section
html_univar <- '<h3>3c. Univariate Wilcoxon tests (TMA)</h3>\n'
if (nrow(univar_dt) > 0) {
  top_hits <- univar_dt[p_value < 0.05][order(p_value)][1:min(.N, 20)]
  if (nrow(top_hits) > 0) {
    html_univar <- paste0(html_univar,
      '<p>Features with p < 0.05 (top 20):</p>\n',
      make_html_table(as.data.frame(top_hits[, .(outcome, feature, n, p_value, padj)])),
      '\n')
  } else {
    html_univar <- paste0(html_univar,
      '<p>No features reached p < 0.05.</p>\n')
  }
}

html_cox <- '<h3>3d. Univariate Cox PH (TMA, 5-year censored)</h3>\n'
if (nrow(cox_dt) > 0) {
  top_cox <- cox_dt[p_value < 0.05][order(p_value)][1:min(.N, 20)]
  if (nrow(top_cox) > 0) {
    html_cox <- paste0(html_cox,
      '<p>Features with p < 0.05 (top 20):</p>\n',
      make_html_table(as.data.frame(top_cox[, .(endpoint, feature, n, n_events,
                                                 HR, HR_lower, HR_upper, p_value, padj)])),
      '\n')
  } else {
    html_cox <- paste0(html_cox,
      '<p>No features reached p < 0.05 in Cox models.</p>\n')
  }
}

html_sections[["clinical"]] <- paste0(
  '<h2>3. Clinical Associations (TMA)</h2>\n',
  html_tumour_ft,
  html_chemo,
  html_univar,
  html_cox
)


# ============================================================================
# PART 4: UCell-Based Associations
# ============================================================================

message("\n[", Sys.time(), "] === PART 4: UCell-Based Associations ===")

html_ucell <- '<h2>4. Polarization UCell Metrics</h2>\n'

ucell_patients <- per_patient[!is.na(polar_mean)]
message("  Patients with polarization data: ", nrow(ucell_patients))

if (nrow(ucell_patients) >= 3) {

  # --- 4a. Per-patient polarization summary table ---
  ucell_summary <- ucell_patients[, .(patient_id, source, treatment_status,
                                       polar_mean, polar_median, polar_sd,
                                       polar_iqr, polar_n)]
  fwrite(ucell_summary, file.path(step_dir, "ucell_summary_by_patient.csv"))

  # --- 4b. Polarization by treatment status ---
  p_polar_treat <- ggplot(ucell_patients, aes(x = treatment_status, y = polar_mean,
                                               fill = treatment_status)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.3, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
    scale_fill_manual(values = c("Chemo-treated" = "#D55E00",
                                 "Treatment-naive" = "#0072B2",
                                 "naive" = "#0072B2",
                                 "post-chemotherapy" = "#D55E00"),
                      guide = "none") +
    labs(x = NULL, y = "Mean polarization_UCell",
         title = "Per-patient mean polarization by treatment status") +
    theme_lab()

  img_polar_treat <- save_and_embed(p_polar_treat, "polar_by_treatment",
                                     width = 700, height = 500)

  # --- 4c. Polarization variability (SD) by treatment ---
  p_polar_sd <- ggplot(ucell_patients, aes(x = treatment_status, y = polar_sd,
                                            fill = treatment_status)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.3, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
    scale_fill_manual(values = c("Chemo-treated" = "#D55E00",
                                 "Treatment-naive" = "#0072B2",
                                 "naive" = "#0072B2",
                                 "post-chemotherapy" = "#D55E00"),
                      guide = "none") +
    labs(x = NULL, y = "SD polarization_UCell",
         title = "Per-patient polarization variability") +
    theme_lab()

  img_polar_sd <- save_and_embed(p_polar_sd, "polar_sd_by_treatment",
                                  width = 700, height = 500)

  # --- 4d. Correlation: polarization vs immune proportion ---
  p_cor_immune <- ggplot(ucell_patients, aes(x = polar_mean, y = prop_immune)) +
    geom_point(aes(color = source), size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.5) +
    scale_color_manual(values = c("TMA" = "#999999", "WT" = "#333333")) +
    labs(x = "Mean polarization_UCell", y = "Immune proportion",
         title = "Polarization vs immune composition") +
    theme_lab()

  # Spearman correlation
  cor_test <- tryCatch(
    cor.test(ucell_patients$polar_mean, ucell_patients$prop_immune,
             method = "spearman", exact = FALSE),
    error = function(e) NULL
  )
  if (!is.null(cor_test)) {
    p_cor_immune <- p_cor_immune +
      labs(subtitle = sprintf("Spearman rho = %.3f, p = %s",
                              cor_test$estimate, signif(cor_test$p.value, 3)))
  }

  img_cor_immune <- save_and_embed(p_cor_immune, "polar_vs_immune",
                                    width = 700, height = 500)

  # --- 4e. Correlation: polarization vs SecB proportion ---
  p_cor_secb <- ggplot(ucell_patients, aes(x = polar_mean, y = prop_SecB_of_sec)) +
    geom_point(aes(color = source), size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.5) +
    scale_color_manual(values = c("TMA" = "#999999", "WT" = "#333333")) +
    labs(x = "Mean polarization_UCell", y = "SecB proportion (of secretory)",
         title = "Polarization vs SecB enrichment") +
    theme_lab()

  cor_secb <- tryCatch(
    cor.test(ucell_patients$polar_mean, ucell_patients$prop_SecB_of_sec,
             method = "spearman", exact = FALSE),
    error = function(e) NULL
  )
  if (!is.null(cor_secb)) {
    p_cor_secb <- p_cor_secb +
      labs(subtitle = sprintf("Spearman rho = %.3f, p = %s",
                              cor_secb$estimate, signif(cor_secb$p.value, 3)))
  }

  img_cor_secb <- save_and_embed(p_cor_secb, "polar_vs_secb",
                                  width = 700, height = 500)

  # --- 4f. UCell clinical tests ---
  ucell_test_results <- list()

  # Wilcoxon tests for polarization features against clinical outcomes
  ucell_feats <- intersect(c("polar_mean", "polar_median", "polar_sd", "polar_iqr"),
                           names(tma_clin))
  for (outcome_name in names(binary_outcomes)) {
    for (feat in ucell_feats) {
      dat_u <- tma_clin[!is.na(get(feat)) & !is.na(get(outcome_name))]
      if (nrow(dat_u) < 10) next

      groups <- split(dat_u[[feat]], dat_u[[outcome_name]])
      if (length(groups) < 2 || any(sapply(groups, length) < 3)) next

      wt_u <- tryCatch(
        wilcox.test(dat_u[[feat]] ~ dat_u[[outcome_name]], exact = FALSE),
        error = function(e) NULL
      )
      if (!is.null(wt_u)) {
        ucell_test_results[[length(ucell_test_results) + 1]] <- data.table(
          outcome = outcome_name,
          feature = feat,
          n       = nrow(dat_u),
          p_value = wt_u$p.value
        )
      }
    }
  }

  ucell_tests_dt <- rbindlist(ucell_test_results)
  if (nrow(ucell_tests_dt) > 0) {
    ucell_tests_dt[, padj := p.adjust(p_value, method = "BH")]
    fwrite(ucell_tests_dt, file.path(step_dir, "ucell_clinical_tests.csv"))
    message("  Saved: ucell_clinical_tests.csv")
  }

  html_ucell <- paste0(html_ucell,
    '<h3>4a. Polarization by treatment status</h3>\n',
    img_polar_treat, '\n',
    '<h3>4b. Polarization variability by treatment status</h3>\n',
    img_polar_sd, '\n',
    '<h3>4c. Polarization vs immune composition</h3>\n',
    img_cor_immune, '\n',
    '<h3>4d. Polarization vs SecB enrichment</h3>\n',
    img_cor_secb, '\n'
  )

  if (nrow(ucell_tests_dt) > 0) {
    html_ucell <- paste0(html_ucell,
      '<h3>4e. UCell feature clinical tests (TMA)</h3>\n',
      make_html_table(as.data.frame(ucell_tests_dt)), '\n')
  }

} else {
  html_ucell <- paste0(html_ucell,
    '<p>Insufficient patients with polarization_UCell data for analysis.</p>\n')
}

html_sections[["ucell"]] <- html_ucell

# Save pathway summary
pw_cols_patient <- grep(paste0("^", pathway_prefix), names(per_patient), value = TRUE)
if (length(pw_cols_patient) > 0) {
  pw_summary <- per_patient[, c("patient_id", "source", "treatment_status",
                                 pw_cols_patient), with = FALSE]
  fwrite(pw_summary, file.path(step_dir, "pathway_summary_by_patient.csv"))
  message("  Saved: pathway_summary_by_patient.csv")
}


# ============================================================================
# PART 5: HTML Report
# ============================================================================

message("\n[", Sys.time(), "] === PART 5: Generating HTML Report ===")

# --- Section 1: Feature table summary ---
feat_summary <- per_patient[, .(patient_id, source, treatment_status,
                                 total_cells, prop_epi, prop_immune, prop_stromal,
                                 prop_SecA_of_sec, prop_Int_of_sec, prop_SecB_of_sec)]
if ("log2_ratio" %in% names(per_patient)) {
  feat_summary[, log2_ratio := per_patient$log2_ratio]
}
if ("polar_mean" %in% names(per_patient)) {
  feat_summary[, polar_mean := per_patient$polar_mean]
  feat_summary[, polar_sd := per_patient$polar_sd]
}

html_feat_table <- paste0(
  '<h2>1. Per-Patient Feature Summary</h2>\n',
  '<p>Total patients: ', nrow(per_patient),
  ' (TMA: ', sum(per_patient$source == "TMA"),
  ', WT: ', sum(per_patient$source == "WT"), ')</p>\n',
  make_html_table(as.data.frame(feat_summary)),
  '\n'
)

# --- Section 6: Key findings ---
html_findings <- '<h2>5. Key Findings Summary</h2>\n<ul>\n'

html_findings <- paste0(html_findings,
  '<li>Extracted per-patient features from ', nrow(per_patient),
  ' patients using atlas-calibrated (06e) reclassified cell labels.</li>\n')

if (nrow(univar_dt) > 0) {
  n_sig <- nrow(univar_dt[p_value < 0.05])
  n_adj <- nrow(univar_dt[padj < 0.05])
  html_findings <- paste0(html_findings,
    '<li>Univariate Wilcoxon tests: ', n_sig, ' features with nominal p < 0.05 (',
    n_adj, ' after BH correction).</li>\n')
}

if (nrow(cox_dt) > 0) {
  n_cox_sig <- nrow(cox_dt[p_value < 0.05])
  html_findings <- paste0(html_findings,
    '<li>Univariate Cox PH: ', n_cox_sig,
    ' features with nominal p < 0.05 for survival.</li>\n')
}

if (nrow(ucell_patients) > 0) {
  html_findings <- paste0(html_findings,
    '<li>Polarization_UCell metrics computed for ', nrow(ucell_patients),
    ' patients (mean, median, SD, IQR).</li>\n')
  if (!is.null(cor_secb)) {
    html_findings <- paste0(html_findings,
      '<li>Polarization-SecB correlation: rho = ',
      round(cor_secb$estimate, 3), ', p = ',
      signif(cor_secb$p.value, 3), '.</li>\n')
  }
}

html_findings <- paste0(html_findings,
  '<li>All intermediate CSVs saved to output/10_clinical_v2/ for reproducibility.</li>\n',
  '</ul>\n')

# --- Assemble HTML ----------------------------------------------------------

html_css <- '
<style>
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
       max-width: 1200px; margin: 40px auto; padding: 0 20px; color: #333; line-height: 1.6; }
h1 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
h2 { color: #34495e; margin-top: 40px; border-bottom: 1px solid #bdc3c7; padding-bottom: 5px; }
h3 { color: #7f8c8d; }
img { max-width: 100%; height: auto; margin: 15px 0; }
.styled-table { border-collapse: collapse; margin: 15px 0; font-size: 0.85em; width: 100%;
                box-shadow: 0 0 10px rgba(0,0,0,0.1); }
.styled-table th { background-color: #34495e; color: white; padding: 8px 12px; text-align: left; }
.styled-table td { padding: 6px 12px; border-bottom: 1px solid #ecf0f1; }
.styled-table tr:nth-child(even) { background-color: #f8f9fa; }
.styled-table tr:hover { background-color: #e8f4fd; }
.callout { background: #fff3e0; border-left: 4px solid #FF9800; padding: 12px; margin: 15px 0; }
.callout-blue { background: #e3f2fd; border-left: 4px solid #2196F3; padding: 12px; margin: 15px 0; }
</style>
'

html_body <- paste0(
  '<!DOCTYPE html>\n<html lang="en">\n<head>\n',
  '<meta charset="UTF-8">\n',
  '<title>10_clinical_v2: Clinical Associations (Atlas-Calibrated Labels)</title>\n',
  html_css, '\n</head>\n<body>\n',
  '<h1>Phase 10 v2: Per-Patient Clinical Associations</h1>\n',
  '<p>Atlas-calibrated UCell reclassified cell labels (06e). Generated: ',
  format(Sys.time(), "%Y-%m-%d %H:%M"), '</p>\n',
  '<div class="callout-blue"><strong>Analysis update</strong>: This version uses ',
  'reclassified cell_label from 06e (atlas-derived UCell thresholds) and adds ',
  'polarization_UCell metrics and pathway score summaries not present in the ',
  'original Phase 10 analysis.</div>\n',
  '\n',
  html_feat_table, '\n',
  html_sections[["composition"]], '\n',
  html_sections[["clinical"]], '\n',
  html_sections[["ucell"]], '\n',
  html_findings, '\n',
  '\n<hr>\n<p style="color: #999; font-size: 0.8em;">',
  'Generated by 10_clinical_v2.R | ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  '</p>\n</body>\n</html>'
)

writeLines(html_body, html_out)
message("  Saved: ", html_out)

# --- Session info -----------------------------------------------------------

message("\n=== 10_clinical_v2 Complete ===")
message("  per_patient_features_v2.csv: ", nrow(per_patient), " patients")
message("  Report: ", html_out)
log_session()
