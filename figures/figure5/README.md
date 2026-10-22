# Figure 5: SecB cells and macrophages adapt to hypoxic niches through a shared glycolytic survival program

Xenium whole tissues (WT01–WT08) and TMA patients. The figure shows epithelium-intrinsic programs along the
polarization axis, nuclear size, WT01 and TMA ROI maps of cell types and pathway scores, immune composition along local
epithelial hypoxia, macrophage-intrinsic programs along local epithelial polarization, and macrophage size in SecA
versus SecB niches.

## Run order

**Requires** (backend, see [`spatial/README.md`](../../spatial/README.md)):
- SFEs in `$OUTPUT_ROOT/sfe`: `spatial/01_build_sfe/` → `spatial/02_qc/04_qc_filter_normalize.R` →
  `spatial/03_annotation_polarization/01_annotation_singler.R` … `05_clean_split_rctd.R` (SingleR reads the probe-QC
  exclusion list `$OUTPUT_ROOT/05_probe_qc/genes_exclude_singler.txt`, copied by hand; see
  [`spatial/02_qc/README.md`](../../spatial/02_qc/README.md)).
- SingleR label tables: `spatial/03_annotation_polarization/06_extract_singler.R` →
  `$OUTPUT_ROOT/singler_recompute/<sample>_singler_cells.csv.gz`, read through `spatial/00_setup/singler_labels.R`
  (all R steps).
- Pathway scores: `spatial/04_pathway_scoring/04c_pathway_scoring_allcells.R` → SFE `pathway_*` (5a, c, d, f, g, i).
- Final TMA core set: `spatial/02_qc/05b_core_qc.R` (×3) → `spatial/02_qc/06b_final_tma_cores.py`, read through
  `spatial/00_setup/final_tma_cores.R` (5b, 5j).

Tables go to `$OUTPUT_ROOT/split_recompute/` and `$OUTPUT_ROOT/figures/figure5/`; panels to
`$OUTPUT_ROOT/figures/figure5/panels/`. The folder name `split_recompute`, the suffix `_presplit` and the column name
`cell_label_split` do not mean a different labelling: every table here is built on the SingleR labels and raw counts.

| Step | Script | Reads | Writes | Panel(s) |
|---|---|---|---|---|
| 1 | `01_fig5a_epithelial_programs_build.R` | SFEs of the 8 WT (`pathway_*`, `counts`, `logcounts`); labels | `split_recompute/fig5a_epi_program_scores_presplit.csv` (1,235,350 secretory cells) | a (also Figure S17) |
| 2 | `02_fig5a_epithelial_programs_fit.R` | step 1 | `figures/figure5/gam_curves_5a.csv` | a |
| 3 | `03_fig5b_nuclear_area.R` | SFE `nucleus_area`, `cell_area` (8 WT, `sfe_tma`, `sfe_COEUR_TMA_G`, `sfe_COEUR_TMA_H`); labels; TMA label tables; final core set | `split_recompute/fig5b_morphometrics_presplit.csv` (per-unit medians) | b |
| 4 | `04_fig5b_nuclear_area_medians.py` | step 3 | per-state medians, n and paired p printed (no plot) | b (medians in Results) |
| 5 | `05_fig5cdegh_roi_export.R` | `sfe_WT01`, `sfe_tma` (`cellSeg` polygons, `pathway_*`, `logcounts` of LDHA); labels; `tma_singler_cells.csv.gz` | `figures/figure5/roi_data/{WT01_roi06,tma_core2}_{poly.csv.gz,attr.csv}` | c, d, e, g, h |
| 6 | `06_fig5f_immune_niche_build.R` | SFE `pathway_hypoxia`, coordinates (8 WT and the 3 TMA arrays); labels | `split_recompute/fig5c_tme_percell_presplit.csv` | f |
| 7 | `07_fig5f_immune_niche_fit.R` (`FIG5C_AXIS=hyp`, the default) | step 6 | `figures/figure5/fig5c_gam_persample.csv` | f |
| 8 | `08_fig5i_macrophage_programs_build.R` | SFE `counts`, `pathway_*`, coordinates (8 WT); labels | `split_recompute/fig5e_macrophage_percell_presplit.csv` | i (also Figure S18) |
| 9 | `09_fig5i_macrophage_programs_fit.R` | step 8 | `figures/figure5/fig5e_mac_gam_curves.csv` | i |
| 10 | `10_fig5j_macrophage_area.R` | SFE `cellSeg` / `nucSeg` polygons and coordinates (8 WT and the 3 TMA arrays); labels; TMA label tables; final core set | `split_recompute/fig5e_mac_morphology_{percell_presplit.csv.gz, stats_presplit.csv, unitmedians_presplit.csv}` | j |
| 11 | `11_fig5_compose_panels.py` | outputs of steps 2, 3, 5, 6, 7, 8, 9, 10; `config/study_codes.py` | `figures/figure5/panels/Fig5*.{svg,png}` (one file per subpanel; ROI maps PNG only); `Figure5_composite_stats_eh.csv` | a–j |

`11_fig5_compose_panels.py` draws every panel on one page grid, which sizes the subpanels, and cuts one file per
subpanel; the page itself is not saved. Its output names carry other letters than the published panels:

| File written | Panel |
|---|---|
| `Fig5a_epithelial_programs` | a |
| `Fig5b_epithelial_nuclear_area` | b |
| `Fig5c_ROI_celltype`, `Fig5c_celltype_legend` | c (the legend also serves e) |
| `Fig5c_ROI_glycolysis`, `Fig5c_ROI_proliferation`, `Fig5c_ROI_LDHA`, `Fig5c_ROI_nfkb` | d |
| `Fig5d_TMA_core_celltype` | e |
| `Fig5e_immune_composition_hypoxia` | f |
| `Fig5f_TMA_core_hypoxia` | g |
| `Fig5g_TMA_core_highlights` | h |
| `Fig5h_macrophage_programs_polarization` | i |
| `Fig5i_macrophage_cell_area` | j |

In `Figure5_composite_stats_eh.csv`, rows "e" are panel f and rows "h" are panel i. The build tables also carry other
letters: `fig5c_*` is panel f, `fig5e_macrophage*` / `fig5e_mac_gam_curves` panel i, `fig5e_mac_morphology_*` panel j.
Script headers use the published letters.

## Environment

| Variable | Value in the final run | Step(s) |
|---|---|---|
| `DATA_ROOT`, `OUTPUT_ROOT` | deposit root; output root (default `./output`) | all |
| `LANG`, `LC_ALL` | `en_US.UTF-8` | R steps |
| `RERUN_SFE_DIR` | `$OUTPUT_ROOT/sfe` | R steps that load SFEs (1, 3, 5, 6, 8, 10) |
| `FINAL_TMA_CORES` | default, `$OUTPUT_ROOT/07_core_qc/final_tma_cores.csv` | 3, 10 |
| `FIG5C_AXIS` | `hyp` (default; `pol` writes `_pol` files that no panel uses) | 7 |
| `STUDY_CODE_KEY` | private key (lab data only; coded names pass through) | 11 (title of panel e) |

```bash
export RERUN_SFE_DIR=$OUTPUT_ROOT/sfe LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
Rscript figures/figure5/01_fig5a_epithelial_programs_build.R && Rscript figures/figure5/02_fig5a_epithelial_programs_fit.R
Rscript figures/figure5/03_fig5b_nuclear_area.R && python figures/figure5/04_fig5b_nuclear_area_medians.py
Rscript figures/figure5/05_fig5cdegh_roi_export.R
Rscript figures/figure5/06_fig5f_immune_niche_build.R && Rscript figures/figure5/07_fig5f_immune_niche_fit.R
Rscript figures/figure5/08_fig5i_macrophage_programs_build.R && Rscript figures/figure5/09_fig5i_macrophage_programs_fit.R
Rscript figures/figure5/10_fig5j_macrophage_area.R
python figures/figure5/11_fig5_compose_panels.py
```

## Panel → script

- **5a** Epithelium-intrinsic program scores along focal polarization (8 WT, all secretory cells). Scores by
  `01_fig5a_epithelial_programs_build.R`: NF-κB, glycolysis, proliferation and the Myc, PI3K, Wnt and Notch activating
  arms from SFE `pathway_*`; DNA damage repair (BRCA1, BRCA2, ATM, CHEK2) and epigenetic silencing (DNMT3A, HDAC1,
  KDM5A) scored inline with UCell (maxRank 477); LDHA and TNFA as single-gene `logcounts` (the gene is TNF); oncogenic
  growth as the mean of the four z-scored growth arms. Fit by `02_fig5a_epithelial_programs_fit.R` (equal weight per
  tissue, `mgcv::bam(y ~ s(x, k = 10))`, 1st–99th percentile). Drawn by `11`.
- **5b** Nuclear area (Xenium `nucleus_area`) of SecA, Intermediate and SecB cells. Unit medians by
  `03_fig5b_nuclear_area.R` (WT per section; TMA per patient on final cores, ≥ 500 cells). Drawn and tested by `11`
  (paired Wilcoxon SecA vs SecB: WT n = 8; TMA patients with all three states, n = 361). The medians quoted in the
  Results (SecB 39.5 vs SecA 47.7 µm² in the whole tissues) are printed by `04_fig5b_nuclear_area_medians.py`.
- **5c** WT01 ROI (800 × 800 µm; x 1499–2299, y 7133–7933 µm in `spatialCoords`) coloured by cell type: exported by
  `05_fig5cdegh_roi_export.R`, drawn by `11` (100 µm scale bar).
- **5d** Same ROI: UCell `pathway_glycolysis`, `pathway_proliferation`, `pathway_nfkb` and LDHA `logcounts`. `11`
  colours SecA, Intermediate, SecB and ciliated epithelial cells (2nd–98th percentile colour limits) and greys the
  rest. Exported by `05`, drawn by `11`.
- **5e** TMA A core "2" (patient TMA-P003) coloured by cell type: exported by `05`; drawn by `11`, whose panel title
  takes the study code from `config/study_codes.py`.
- **5f** Immune cell share of the local niche along local epithelial hypoxia. `06_fig5f_immune_niche_build.R`: one row
  per non-epithelial cell, x = mean `pathway_hypoxia` of epithelial cells within 50 µm, z-scored per sample.
  `07_fig5f_immune_niche_fit.R`: binomial GAM per tissue and cell type (k = 6; fitted when the tissue has ≥ 50 cells of
  that type, hence B cells in 5 and plasma cells in 7 tissues). `11` draws the mean of the per-tissue curves and runs
  the per-tissue test (lowest vs highest 20% of x, paired Wilcoxon, n = 8), written to `Figure5_composite_stats_eh.csv`.
  The curve grid spans the 1st–99th percentile of x pooled over the whole-tissue **and** TMA rows of step 6, so the
  quoted curve endpoints (macrophages 13.4% → 32.6%) depend on the TMA rows although the panel shows whole tissues
  only.
- **5g** Hypoxia UCell score over the TMA core: exported by `05`; drawn by `11`, which also computes the box (5th–95th
  percentile of SecB cell centroids ± 25 µm).
- **5h** The TMA core with SecB cells, macrophages, T, NK, B and plasma cells highlighted (same box): `05` → `11`.
- **5i** Macrophage-intrinsic UCell programs along local epithelial polarization. `08_fig5i_macrophage_programs_build.R`
  scores the four sets inline (glycolytic LDHA, ENO1, PGK1, ALDOA; NF-κB activated ICAM1, NFKBIA, BIRC3;
  tissue-resident A2M, C1QC; chemokine CXCL10, CXCL11; maxRank 477) and the x axis (`nb_pol`: summed polarization of
  epithelial neighbours within 50 µm divided by their number, ciliated neighbours counting as 0; z-scored per tissue).
  `09_fig5i_macrophage_programs_fit.R` fits all 100,372 macrophages (same GAM spec as 5a). `11` draws the panel and
  runs the per-tissue test (stats rows "h").
- **5j** Median macrophage cell area in the lowest vs highest tertile of local secretory polarization.
  `10_fig5j_macrophage_area.R`: niche = mean polarization of secretory neighbours within 50 µm (≥ 5 neighbours);
  tertiles within unit; units with ≥ 10 macrophages in both tertiles; areas from `cellSeg` / `nucSeg` polygons; TMA
  units are patients on final cores. Drawn and tested by `11` (paired Wilcoxon on unit medians; WT n = 8, TMA n = 297).

## Not in this folder

- Backend steps (see Requires). Figures S16–S18 (per-tissue GAMs, including the per-tissue tests of the Figure 5a and
  5i features): `figures/supplementary/21_figS16-18_per_sample_gams.R` and
  `figures/supplementary/22_figS16-18_per_sample_gams_plot.py`.
- **Manual:** panel letters and the order of the four 5d maps are set during assembly in Affinity Designer; the panel
  files are copied and renamed into the figure. The probe-QC exclusion list is copied by hand into
  `$OUTPUT_ROOT/05_probe_qc/` (see Requires).
- Upstream of the TMA panels (5b, e, g, h, j and the TMA rows of 5f): the TMA A and TMA B de-array code is in
  `spatial/01_build_sfe/as_run/`; the TMA A core map it read is **[not on disk]** (see
  [`spatial/README.md`](../../spatial/README.md)).
