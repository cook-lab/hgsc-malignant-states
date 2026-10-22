# `spatial/04_pathway_scoring/`: the 37 UCell pathway modules

Defines the 37 pathway modules of Table S9 (325 gene entries) and scores them with UCell on every cell of every SFE
used in the manuscript, writing one `pathway_<module>` column per module into the SFEs. These scores feed the
hypoxia, matrix-remodeling, glycolysis, NF-κB and proliferation curves and maps of Figures 4 and 5 and Figures
S16–S18.

## Run order

| Step | Script | Reads | Writes | Used by |
|---|---|---|---|---|
| 1 | `04_pathway_scoring.R` (the final run's writer) | none (gene sets defined inline) | `$OUTPUT_ROOT/9b_scoring/pathway_gene_sets_v2.csv` (pathway, gene) | step 2; `figures/figure4/16_fig4j_wt01_roi_hypoxia.R` (prints the hypoxia genes); `tables/07_tableS7_xenium_gene_panel.py`, `tables/09_tableS9_ucell_gene_sets.py`; a hand copy at `$OUTPUT_ROOT/fig6_singler/inputs/pathway_gene_sets_v2.csv` is read by `figures/figure6/05_fig6fh_xenium_reciprocal_effects.R` |
| 1 (alt) | `04b_pathway_gene_sets.R` | none (the same 37 sets, inline) | the same `pathway_gene_sets_v2.csv` | `tables/09_tableS9_ucell_gene_sets.py` parses its inline sets and asserts they match the file |
| 2 | `04c_pathway_scoring_allcells.R` | SFEs `sfe_WT01`–`sfe_WT08`, `sfe_tma`, `sfe_COEUR_TMA_G`, `sfe_COEUR_TMA_H` (raw `counts`, full 477-gene panel, all cells); `pathway_gene_sets_v2.csv`; `$DATA_ROOT/2026_final_xenium_analysis/output/05_probe_qc/flag_qc.csv` | SFE `pathway_<module>` (UCell, maxRank = 477) saved to `$OUTPUT_ROOT/sfe`; `$OUTPUT_ROOT/9b_scoring/pathway_coverage_allcells.csv`, `scoring_summary_allcells.csv` | `figures/figure4/08_fig4f_gam_features.R`, `16_fig4j_wt01_roi_hypoxia.R`; `figures/figure5/01_fig5a_epithelial_programs_build.R`, `05_fig5cdegh_roi_export.R`, `06_fig5f_immune_niche_build.R`, `08_fig5i_macrophage_programs_build.R`; `figures/supplementary/21_figS16-18_per_sample_gams.R`; `tables/09_tableS9_ucell_gene_sets.py` (the two CSVs) |

Steps 1 and 1 (alt) write byte-identical files in seconds; run either. Step 2 is heavy (about 1.5 h over the 11
SFEs).

```bash
Rscript -e 'source("spatial/04_pathway_scoring/04_pathway_scoring.R")'   # or: Rscript spatial/04_pathway_scoring/04b_pathway_gene_sets.R
Rscript spatial/04_pathway_scoring/04c_pathway_scoring_allcells.R
mkdir -p $OUTPUT_ROOT/fig6_singler/inputs
cp $OUTPUT_ROOT/9b_scoring/pathway_gene_sets_v2.csv $OUTPUT_ROOT/fig6_singler/inputs/   # read by figures/figure6/05
```

`04_pathway_scoring.R` finds `config/` relative to its own folder, which R only knows when the file is `source()`d: run
it with `Rscript -e 'source(...)'` from the repository root, as above, or from inside its folder.
`04b_pathway_gene_sets.R` runs with `Rscript` directly.

## Notes

- **Substrate.** Scores use the raw counts of the full 477-gene panel on all cells, with maxRank fixed at 477 so that
  scores compare across samples. The 28 probe-QC exclusions apply to annotation only, not to scoring. Any `pathway_*`
  columns already in an SFE are dropped before scoring.
- **Probe-QC flags.** `pathway_coverage_allcells.csv` lists, per module and sample, the genes found and any probe-QC
  FAIL or WARN genes. Modules with a FAIL gene: `matrix_remodeling` (MMP7), `adc_targets` (ERBB2, VTCN1),
  `apoptosis_pro` (BAX), `ecm_adhesion` (CLDN3). A missing `flag_qc.csv` silently leaves these columns empty; scores
  are unaffected.
- **Gene sets defined in figure scripts.** A few sets are not among the 37 and are scored inside figure scripts: DNA
  damage repair and epigenetic silencing (`figures/figure5/01_fig5a_epithelial_programs_build.R`), the four macrophage
  programs (`figures/figure5/08_fig5i_macrophage_programs_build.R`), and the glycolytic and resident macrophage programs
  of Figure 6g (`figures/figure6/06_fig6g_macrophage_secB_neighbours.py`). Table S9 parses them from those scripts.

## Environment

| Variable | Final run | Step |
|---|---|---|
| `RERUN_SFE_DIR` | `$OUTPUT_ROOT/sfe` | 2 (`load_sfe()` falls back to `sfe_dir`) |
| `LANG`, `LC_ALL` | `en_US.UTF-8` | all |
