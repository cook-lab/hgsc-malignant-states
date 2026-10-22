# figures/_prep: shared loaders and per-figure prep

The code here does two jobs, and none of it draws a panel:
- **Loaders**, which figure scripts import.
- **Prep scripts**, which run once before a figure script and write a cache.

The outputs go to `$OUTPUT_ROOT/_prep_caches/` or to `$OUTPUT_ROOT/figure2/…`. `__init__.py` makes the folder a
package, so that `from figures._prep.<module> import …` works. Run every script from the repository root.

## Run order

Paths without a root are under `$DATA_ROOT/2026_final_atlas/`.

| Step | Script | Reads | Writes | Feeds |
|---|---|---|---|---|
| 1 | `nath2021_fix.py` (loader; builds its cache the first time the fix is applied, or run it directly to rebuild and print the old → new sample map) | `pre-processing and integration/processed/nath_2021.h5ad` (obs only); `pre-processing and integration/atlas_metadata.xlsx`; `pre-processing and integration/raw/nath_2021/GSE158722.cell_annotations.txt.gz`; `hgsc_atlas_final.h5ad` obs | `$OUTPUT_ROOT/_prep_caches/nath2021_cell_correction.parquet`, `nath2021_sample_map.csv` | step 2; Table S1 (`tables/01_tableS1_atlas_studies.py`); Fig 6 atlas pseudobulk (`figures/figure6/08_fig6fh_atlas_macrophage_pseudobulk.py`) |
| 2 | `atlas_final_tables.py` (loader: `load_obs()`, `load_epithelial()`; builds its caches on first use, or run it directly) | `hgsc_atlas_final.h5ad` (obs and `obsm/X_umap`, through h5py; X is not loaded); `output/11d_epithelial_nmf/11d_nmf_usage.csv`; `output/fig_secretory_polarization/data/meta.parquet` (UMAP1/UMAP2 only); step 1 | `$OUTPUT_ROOT/_prep_caches/atlas_final_obs.parquet` (1,980,703 cells), `atlas_final_epithelial.parquet` (575,366 cells, with Factor_1–10 and `state`) | the Figure 1, Figure 2 and Figure S2–S9 scripts, directly or through steps 3–9 |
| 3 | `fig1g_polarization_pc_alignment.py` | `celltype_h5ad/hgsc_atlas_final_epithelial.h5ad` (`layers/counts`); `load_epithelial()`; `output/19_cnv/per_sample/<sample>/copykat_{prediction,subclones,subclone_qc}.csv` | `$OUTPUT_ROOT/_prep_caches/fig1g_polarization_pc_alignment.csv` (one row per sample, 250 samples; `pol_top_pc`, plus the PC-basis variance shares `pc_*` that the paper does not use) | Fig 1g (`figures/figure1/04_fig1g_polarization_pc_alignment.py`) |
| 4 | `fig2a_scfea_all_cells.py` | epithelial h5ad (log-normalised X); `tools/scFEA/data/{module_gene_m168,cmMat_c70_m168,cName_c70_m168}.csv`; `load_epithelial()` | `$OUTPUT_ROOT/figure2/fig2a_scfea_all/{21_flux_per_cell.parquet, 21_balance_per_cell.parquet, run_log.json}` (575,366 cells × 168 modules; only the flux is read downstream) | step 5; Figure S8 |
| 5 | `fig2a_functional_state_summaries.py` | `output/21_epitype_functional_characterization/21_{progeny,hallmark,dorothea}_per_cell.parquet`; step 4; `load_epithelial()` | `$OUTPUT_ROOT/figure2/fig2a_functional/21_{progeny,hallmark,dorothea,flux}_{zscored,means}.csv`, `coverage.csv` | Fig 2a (`figures/figure2/01_fig2a_radar_characterization.py`); Figure S8 (`figures/supplementary/14_figS8_functional_heatmap.py`) |
| 6 | `fig2b_cell_cycle_phase.py` | epithelial h5ad (`layers/counts`); `load_epithelial()` | `$OUTPUT_ROOT/figure2/fig2b_cell_cycle/{cell_cycle_per_cell.parquet, cell_cycle_phase_proportions.csv, cell_cycle_phase_counts.csv}` | Fig 2b–e (`figures/figure2/02_fig2bcde_cellcycle_composition_paired.py`, which needs this output even for c–e) |
| 7 | `fig2f_consensusov.py` (calls `fig2f_consensusov.R` through `Rscript`) | `hgsc_atlas_final.h5ad` (`layers/counts` rows of the selected samples); `load_obs()`, `load_epithelial()` | `$OUTPUT_ROOT/figure2/fig2f_consensusov/{pseudobulk_counts.tsv.gz, pseudobulk_metadata.csv, consensusov_calls.csv, gene_id_mapping.csv, fig2f_per_sample_joined.csv}` (the R script writes the calls and the gene map) | step 8; Fig 2f |
| 8 | `fig2f_consensusov_multiseed.R` (run by hand; command below) | `pseudobulk_counts.tsv.gz` from step 7 | `$OUTPUT_ROOT/figure2/fig2f_consensusov/consensusov_calls_multiseed.csv` and `consensusov_calls_multiseed_per_seed.csv` | Fig 2f (`figures/figure2/03_fig2f_subtype_state_boxes.py`) |
| 9 | `text_state_coexistence.py` (run after `figures/supplementary/13_figS7_cnv_clone_validation.py`) | `load_epithelial()`; `$OUTPUT_ROOT/_prep_caches/sf7_cnv/within_clone_coexistence.csv` (written by the Figure S7 script) | `$OUTPUT_ROOT/_prep_caches/text_state_coexistence.csv` | Discussion: "93.5% of patients with ≥ 100 epithelial cells contain SecA, Intermediate, and SecB cells" (144 of 154 patients) |
| – | `final_tma_cores.py` (loader: `final_core_ids()`, `final_cores()`, `final_patients()`) | `$FINAL_TMA_CORES` (default `$OUTPUT_ROOT/07_core_qc/final_tma_cores.csv`) and `final_tma_patients.csv` next to it, both written by `spatial/02_qc/06b_final_tma_cores.py` | – | the TMA arm of Fig 4a, 4g and 4i (`figures/figure4/01_fig4a_composition_table.py`, `figures/figure4/11_fig4g_vascular_distance.py`, `figures/figure4/14_fig4i_local_density.py`), Fig 6g (`figures/figure6/06_fig6g_macrophage_secB_neighbours.py`), Figures S14 and S15 (`figures/supplementary/19_figS14_annotation_validation_plot.py`, `figures/supplementary/20_figS15_polarization_distribution.py`), Tables S8 and S11. The R twin is `spatial/00_setup/final_tma_cores.R`. |

## What the loaders do

**`nath2021_fix.py`.**
- **The problem it fixes.** In the deposited atlas, each Nath 2021 (GSE158722) "sample" is a whole GEO patient that
  carries another sample's name, patient and treatment.
- **The fix.** The loader reassigns every Nath cell from its GEO cell id to the curated sample (patient, time point,
  treatment, TP53, BRCA). It does this when the data are loaded. The deposited object is not changed, and no other
  study is touched.
- **The result.** With the fix, the atlas has 158 patients and 362 samples.
- **`sample_id_deposited`.** `apply_nath2021_fix(obs)` adds this column, the sample id before the fix. The deposited
  CopyKAT folders are named by it, so Fig 1g, Fig 1l,m, Figure S4c and Figure S7 join CopyKAT on it.

**`atlas_final_tables.py`.** This is the single source of cell metadata for every atlas panel.
- **Inputs.** Metadata come from the final atlas, after the Nath 2021 fix. The epithelial UMAP coordinates come from
  `meta.parquet`, and NMF usage from the deposited `11d` run. (Table S6 and the Fig 6h atlas row read the `schema_nmf`
  state column of `meta.parquet` directly.)
- **Epithelial state.** Among non-ciliated epithelial cells:
  - Factor 2 usage below p50 is SecA.
  - p50 to p75 is Intermediate.
  - p75 and above is SecB.
- **Asserted counts.** 284,344 SecA, 142,172 Intermediate, 142,172 SecB and 6,678 ciliated cells; 575,366 epithelial
  cells in total.

## Run commands

Steps 1–9 are atlas steps. The final run used one output root for them and for the Figure 1, Figure 2 and
Figure S1–S9 scripts, and another for the Xenium steps, which is where `final_tma_cores.py` reads its table. The
scripts here also run with a single `OUTPUT_ROOT` for both.

```bash
export DATA_ROOT=/path/to/deposit OUTPUT_ROOT=/path/to/outputs
python figures/_prep/nath2021_fix.py
python figures/_prep/atlas_final_tables.py
python figures/_prep/fig1g_polarization_pc_alignment.py
python figures/_prep/fig2a_scfea_all_cells.py              # PyTorch on CPU (--threads, default 12); torch seed 16 is set in the script
python figures/_prep/fig2a_functional_state_summaries.py
python figures/_prep/fig2b_cell_cycle_phase.py
python figures/_prep/fig2f_consensusov.py                  # needs Rscript (E1: consensusOV, edgeR, org.Hs.eg.db)
Rscript figures/_prep/fig2f_consensusov_multiseed.R \
    "$OUTPUT_ROOT/figure2/fig2f_consensusov/pseudobulk_counts.tsv.gz" \
    "$OUTPUT_ROOT/figure2/fig2f_consensusov/consensusov_calls_multiseed.csv" 25
python figures/supplementary/13_figS7_cnv_clone_validation.py
python figures/_prep/text_state_coexistence.py
```

## Environment

| Variable | Used by | Value in the final run |
|---|---|---|
| `DATA_ROOT` | steps 1–9 (step 9 only when the loader rebuilds its caches) | deposit root (holds `2026_final_atlas/`) |
| `OUTPUT_ROOT` | all | steps 1–9: the atlas output root; `final_tma_cores.py`: the Xenium output root |
| `FINAL_TMA_CORES` | `final_tma_cores.py` | not set; the default is `$OUTPUT_ROOT/07_core_qc/final_tma_cores.csv` |

**Caches are reused.**
- `atlas_final_tables.py` reuses `atlas_final_{obs,epithelial}.parquet` in `$OUTPUT_ROOT/_prep_caches/` whenever they
  exist and carry the Nath 2021 flag.
- `nath2021_fix.py` reuses its own cache when it is imported; run directly, it deletes and rebuilds it.
- Neither checks whether the atlas or NMF inputs have changed. Delete the caches to force a rebuild.

**Options that the paper does not use.**
- `fig1g_polarization_pc_alignment.py --samples <ids>`: a debug subset, written to
  `fig1g_polarization_pc_alignment_subset.csv`.
- `fig2a_scfea_all_cells.py --validate`: re-runs the 20,000 cells of the deposited
  `output/21_epitype_functional_characterization/21_flux_per_cell.parquet` and writes `validate_vs_2026-05_subset.csv`
  instead of the flux.
- `fig2f_consensusov.py --cohort all_tumours`: every non-healthy sample, written to `fig2f_consensusov_all/`.

**Software.** Python steps run in E2 (Python 3.12.9). Step 7's R call and step 8 run in E1 (R 4.4.3; consensusOV
1.25.1, edgeR 4.4.2, org.Hs.eg.db 3.20.0).

## Figure 2f calls

The Fig 2f subtype of each sample is the 25-seed consensus of step 8: seeds 1–25, mean random-forest probability,
highest mean wins. Step 7 also makes one call with seed 42 (`fig2f_consensusov.R`); it fills the subtype column of
`fig2f_per_sample_joined.csv`, which `figures/figure2/03_fig2f_subtype_state_boxes.py` replaces with the consensus
call on read, so the paper uses only the consensus.
