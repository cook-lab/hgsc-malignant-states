# Figure 4: SecB-polarized secretory cells co-localize within avascular and hypoxic niches

Xenium spatial transcriptomics of 8 whole HGSC tissues (WT01–WT08), the TMA cohort (378 patients) and 15 fallopian
tube epithelium (FTE) cores. The figure shows tissue composition, WT01 cell-type maps, the spatial segregation of SecA
and SecB scores (bivariate Moran's I, BiLISA), and how distance to vasculature, local cell density, matrix remodeling
and hypoxia change along the polarization axis.

## Run order

**Requires** (backend, see [`spatial/README.md`](../../spatial/README.md)):
- SFEs in `$OUTPUT_ROOT/sfe` with SingleR labels and UCell scores: `spatial/01_build_sfe/` →
  `spatial/02_qc/04_qc_filter_normalize.R` → `spatial/03_annotation_polarization/01_annotation_singler.R` …
  `05_clean_split_rctd.R`. SingleR reads the probe-QC exclusion list `$OUTPUT_ROOT/05_probe_qc/genes_exclude_singler.txt`
  (from `spatial/02_qc/as_run/explore_probe_qc.Rmd`, copied by hand; see [`spatial/02_qc/README.md`](../../spatial/02_qc/README.md)).
- SingleR label tables: `spatial/03_annotation_polarization/06_extract_singler.R` →
  `$OUTPUT_ROOT/singler_recompute/<sample>_singler_cells.csv.gz` (all panels).
- Final TMA core set: `spatial/02_qc/05b_core_qc.R` (×3) → `spatial/02_qc/06b_final_tma_cores.py` →
  `$OUTPUT_ROOT/07_core_qc/final_tma_cores.csv` (4a, 4d, 4g, 4i).
- Pathway scores: `spatial/04_pathway_scoring/04_pathway_scoring.R` (or `04b_pathway_gene_sets.R`) →
  `$OUTPUT_ROOT/9b_scoring/pathway_gene_sets_v2.csv` → `spatial/04_pathway_scoring/04c_pathway_scoring_allcells.R` →
  SFE `pathway_*` (4f, 4j).

Outputs go to `$OUTPUT_ROOT/singler_recompute/` (tables) and `$OUTPUT_ROOT/figures/figure4/` (panels).

| Step | Script | Reads | Writes | Panel(s) |
|---|---|---|---|---|
| 1 | `01_fig4a_composition_table.py` | label tables of WT01–WT08, `tma`, `COEUR_TMA_G`, `COEUR_TMA_H`; final core set | `singler_recompute/fig4a_composition_by_patient_SINGLER.csv` | a |
| 2 | `02_fig4a_composition_plot.py` | step 1 table (`COMP_CSV` overrides) | `figures/figure4/xenium_composition_by_tissue_dendro_panelA.{png,svg}` | a |
| 3 | `03_fig4bc_wt01_celltype_maps.R` (+ `scalebar.R`) | `sfe_WT01`; WT01 label table | `figures/figure4/fig4C_WT01_split_celltype.png`, `fig4D_ROI_C_split_celltype.png` | b, c |
| 4 | `04_fig4d_bivariate_moran.R` | SFE `SecA_UCell`, `SecB_UCell`, coordinates (8 WT, `sfe_tma`, `sfe_COEUR_TMA_G`, `sfe_COEUR_TMA_H`); labels via `spatial/00_setup/singler_labels.R`; TMA label tables; final core set via `spatial/00_setup/final_tma_cores.R` | `singler_recompute/bivmoran_{wt,tma_percore,tma_patient,tests}_SINGLER.csv` | d |
| 5 | `05_fig4d_bivariate_moran_plot.py` | step 4 tables (except `tma_percore`); `config/study_codes.py` | `figures/figure4/Fig4d_bivariate_moran.{png,svg}`, `Fig4d_bivariate_moran_{wholetissue,tests}.csv` | d |
| 6 | `06_fig4e_bilisa.R` | SFE `SecA_UCell`, `SecB_UCell`, coordinates of the 8 WT; labels via `singler_labels.R` | `singler_recompute/bilisa_crosstab_SINGLER.csv`, `bilisa_percell_SINGLER.csv.gz` (per-cell record) | e |
| 7 | `07_fig4e_bilisa_plot.py` | `bilisa_crosstab_SINGLER.csv` (`BILISA_CROSSTAB` overrides) | `figures/figure4/bilisa_regime_by_label.{png,svg}` | e |
| 8 | `08_fig4f_gam_features.R` | SFE `pathway_hypoxia`, `pathway_matrix_remodeling`, coordinates (8 WT); label tables | `singler_recompute/singler_gam_percell.csv` (1,235,350 secretory cells) | f (also Figure S16) |
| 9 | `09_fig4f_gam_fit.R` | `singler_gam_percell.csv` | `figures/figure4/gam_curves.csv`, `gam_labels.csv` | f |
| 10 | `10_fig4f_gam_plot.py` | `gam_curves.csv`, `gam_labels.csv` | `figures/figure4/gam_microenvironment_polarization.{png,svg}` | f |
| 11 | `11_fig4g_vascular_distance.py` | label tables; final core set | `singler_recompute/vascular_distance_summary_SINGLER.csv` | g |
| 12 | `12_fig4g_vascular_distance_plot.py` | `vascular_distance_summary_SINGLER.csv` (`VASC_SUMMARY_CSV` overrides) | `figures/figure4/xenium_vascular_distance_paired.{png,svg}` | g |
| 13 | `13_fig4h_wt01_roi_vascular_distance.R` (+ `scalebar.R`) | `sfe_WT01`; WT01 label table | `figures/figure4/xenium_roi_WT01_vascular_distance.{png,svg}` | h |
| 14 | `14_fig4i_local_density.py` | label tables; final core set | `singler_recompute/per_patient_local_density_SINGLER.csv` | i |
| 15 | `15_fig4i_local_density_plot.py` | `per_patient_local_density_SINGLER.csv` (`PER_PATIENT_CSV` overrides) | `figures/figure4/xenium_epi_density_paired.{png,svg}`, `_stats.csv`, `_units.csv` | i |
| 16 | `16_fig4j_wt01_roi_hypoxia.R` (+ `scalebar.R`) | `sfe_WT01` `pathway_hypoxia`; `9b_scoring/pathway_gene_sets_v2.csv` (gene list printed for the log) | `figures/figure4/xenium_roi_WT01_hypoxia.{png,svg}` | j |

`03_fig4bc_wt01_celltype_maps.R` writes `fig4C_WT01_split_celltype.png`, which is panel b, and
`fig4D_ROI_C_split_celltype.png`, which is panel c. "4C", "4D" and "split" in these names are not the published
letters or labels; the labels drawn are SingleR. Script headers use the published letters.

## Environment

| Variable | Value in the final run | Step(s) |
|---|---|---|
| `DATA_ROOT`, `OUTPUT_ROOT` | deposit root; output root (default `./output`) | all |
| `LANG`, `LC_ALL` | `en_US.UTF-8` | R steps (3, 4, 6, 8, 9, 13, 16) |
| `RERUN_SFE_DIR` | `$OUTPUT_ROOT/sfe` | R steps that load SFEs (3, 4, 6, 8, 13, 16) |
| `FINAL_TMA_CORES` | default, `$OUTPUT_ROOT/07_core_qc/final_tma_cores.csv` | 1, 4, 11, 14 |
| `COMP_CSV` | default, `$OUTPUT_ROOT/singler_recompute/fig4a_composition_by_patient_SINGLER.csv` | 2 |
| `PANEL_A_ONLY` | default, `1` (the Fig 4a layout: bars, dendrogram and keys) | 2 |
| `BILISA_CROSSTAB` | default, `$OUTPUT_ROOT/singler_recompute/bilisa_crosstab_SINGLER.csv` | 7 |
| `VASC_SUMMARY_CSV` | default, `$OUTPUT_ROOT/singler_recompute/vascular_distance_summary_SINGLER.csv` | 12 |
| `PER_PATIENT_CSV` | default, `$OUTPUT_ROOT/singler_recompute/per_patient_local_density_SINGLER.csv` | 15 |
| `STUDY_CODE_KEY` | private key (lab data only; coded names pass through) | 5 |

**Defaults.** Steps 2, 7, 12 and 15 read the table that the step before them writes; the variables only point them at
another copy. If the table is missing, steps 2, 7 and 12 stop (7 and 12 with a `NameError`); step 15 ignores a
`PER_PATIENT_CSV` that does not exist and reads its default. `PANEL_A_ONLY` must not be set to an empty string: the
Fig 4a layout is the only one the script draws.

```bash
export RERUN_SFE_DIR=$OUTPUT_ROOT/sfe LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
python figures/figure4/01_fig4a_composition_table.py && python figures/figure4/02_fig4a_composition_plot.py
Rscript figures/figure4/03_fig4bc_wt01_celltype_maps.R
Rscript figures/figure4/04_fig4d_bivariate_moran.R && python figures/figure4/05_fig4d_bivariate_moran_plot.py
Rscript figures/figure4/06_fig4e_bilisa.R && python figures/figure4/07_fig4e_bilisa_plot.py
Rscript figures/figure4/08_fig4f_gam_features.R && Rscript figures/figure4/09_fig4f_gam_fit.R && python figures/figure4/10_fig4f_gam_plot.py
python figures/figure4/11_fig4g_vascular_distance.py && python figures/figure4/12_fig4g_vascular_distance_plot.py
Rscript figures/figure4/13_fig4h_wt01_roi_vascular_distance.R
python figures/figure4/14_fig4i_local_density.py && python figures/figure4/15_fig4i_local_density_plot.py
Rscript figures/figure4/16_fig4j_wt01_roi_hypoxia.R
```

## Panel → script

- **4a** Per-sample composition, hierarchical clustering (ward.D2). Drawn by `02_fig4a_composition_plot.py`; table by
  `01_fig4a_composition_table.py`: FTE cores (TMA A cores with sample type "fallopian", 15 cores, 34,281 cells), 8 whole
  tissues (1,898,772 cells) and TMA patients of the final core set (378); samples need ≥ 500 cells.
- **4b** WT01 whole-tissue cell-type map, 1 mm scale bar: `03_fig4bc_wt01_celltype_maps.R`.
- **4c** WT01 ROI cell-type map, 200 µm scale bar: `03_fig4bc_wt01_celltype_maps.R`. The ROI is a fixed box
  (x 2800–4000, y 6600–7800 µm in the SFE `spatialCoords` frame), shared with 4h and 4j.
- **4d** Global bivariate Moran's I (SecA score vs neighbours' SecB score; binary 50 µm weights; 999 permutations per
  tissue, 199 per TMA core; TMA cores with ≥ 500 secretory cells, averaged per patient). Computed by
  `04_fig4d_bivariate_moran.R`, which also runs the tests: exact sign test over the 8 tissues, one-sample Wilcoxon over
  TMA patients (n = 334). Drawn by `05_fig4d_bivariate_moran_plot.py`.
- **4e** BiLISA neighbourhoods of secretory cells in the 8 whole tissues (binary 50 µm weights, 499 permutations, BH
  within tissue): `06_fig4e_bilisa.R`; drawn by `07_fig4e_bilisa_plot.py`. In its whole-tissue loop step 6 also
  computes, and does not write, a global row-standardised bivariate statistic: its permutations come before each
  tissue's BiLISA permutations in the random-number stream, so it is kept to reproduce the BiLISA p-values exactly.
- **4f** GAMs along focal polarization. Features by `08_fig4f_gam_features.R` (per secretory cell: distance to the
  nearest pericyte or endothelial cell; number of cells within 50 µm per mm²; mean `pathway_matrix_remodeling` and
  `pathway_hypoxia` of the cells within 50 µm). Fit by `09_fig4f_gam_fit.R` (all cells, equal weight per tissue,
  weighted z-score, `mgcv::bam(y ~ s(x, k = 10))`, curves over the 1st–99th percentile of polarization). Drawn by
  `10_fig4f_gam_plot.py`. The per-tissue tests of these four features (top vs bottom 20%, paired Wilcoxon, n = 8) are
  computed by `figures/supplementary/21_figS16-18_per_sample_gams.R` (Figure S16), not here.
- **4g** Distance to the nearest pericyte or endothelial cell for SecA, Intermediate and SecB cells. Per-unit medians
  by `11_fig4g_vascular_distance.py` (WT per section; TMA per core, then median per patient on final cores). Drawn and
  tested by `12_fig4g_vascular_distance_plot.py` (paired Wilcoxon SecA vs SecB; WT n = 8, TMA patients n = 361).
- **4h** WT01 ROI coloured by distance to vasculature: `13_fig4h_wt01_roi_vascular_distance.R`.
- **4i** Local cell density around SecA, Intermediate and SecB cells (all cells within 50 µm, per mm² of the circle,
  median per unit): `14_fig4i_local_density.py` (TMA neighbours searched within the core; a patient's final cores
  pooled). Drawn and tested by `15_fig4i_local_density_plot.py` (paired Wilcoxon SecA vs SecB; WT n = 8, TMA n = 361;
  tests use the median neighbour counts). The `*_core50` columns are an edge check that is not plotted.
- **4j** WT01 ROI coloured by the hypoxia UCell score (2nd–98th percentile colour limits):
  `16_fig4j_wt01_roi_hypoxia.R`.

## Not in this folder

- Backend steps (see Requires) and the Figure S11–S16 scripts in [`figures/supplementary/`](../supplementary/README.md),
  which share them.
- **Manual:** the ROI box on the whole-tissue map in 4b is drawn by hand during assembly; panel letters and the
  cell-type key are also placed by hand. Panels are copied and renamed into the figure and assembled in Affinity
  Designer. The probe-QC exclusion list is copied by hand into `$OUTPUT_ROOT/05_probe_qc/` (see Requires).
- Source-data tables: step 5 writes study codes through `config/study_codes.py`. The per-unit table of step 15 keeps
  the sample names of its input.
- Upstream of the TMA arms (4a, 4d, 4g, 4i): the TMA A and TMA B de-array code is in `spatial/01_build_sfe/as_run/`;
  the TMA A core map it read is **[not on disk]** (see [`spatial/README.md`](../../spatial/README.md)).
