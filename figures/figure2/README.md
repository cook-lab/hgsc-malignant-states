# Figure 2: programs, cell cycle and distribution of SecA and SecB epithelium

Figure 2 compares SecA and SecB epithelium in the atlas:
- **a:** pathway, gene-set, transcription-factor and metabolic activity
- **b:** cell-cycle phase
- **c:** composition by anatomic site and by treatment
- **d, e:** patient-matched shifts across sites and after chemotherapy
- **f:** composition by consensusOV subtype
- **g:** pseudobulk differential expression, with the 7 + 7 SecA / SecB signature genes highlighted

## Run order

**Requires.**
- Deposited atlas outputs under `$DATA_ROOT/2026_final_atlas/` (see [`../../atlas/README.md`](../../atlas/README.md)):
  - `hgsc_atlas_final.h5ad` (`layers/counts` for f)
  - `celltype_h5ad/hgsc_atlas_final_epithelial.h5ad` (X for a; counts for b and g)
  - `output/11d_epithelial_nmf/11d_nmf_usage.csv` and `output/fig_secretory_polarization/data/meta.parquet` (through
    the loader)
  - `output/21_epitype_functional_characterization/21_{progeny,hallmark,dorothea}_per_cell.parquet`, written by
    `atlas/04_functional/as_run/21_epitype_functional_characterization.py`
  - `tools/scFEA/data/` (the scFEA human M168 model files)
- The loaders: `figures/_prep/nath2021_fix.py` → `figures/_prep/atlas_final_tables.py`.
- `shared/signatures.yml`, for the 7 + 7 genes in g.

Renders go to `$OUTPUT_ROOT/figures/figure2/` as `.svg` and `.png`. The statistics CSVs of b–f are written next to
them.

| Step | Script | Reads | Writes | Panel(s) |
|---|---|---|---|---|
| 1 | `figures/_prep/fig2a_scfea_all_cells.py` | epithelial h5ad X; scFEA M168 model files; `load_epithelial()` | `$OUTPUT_ROOT/figure2/fig2a_scfea_all/21_flux_per_cell.parquet` (+ `21_balance_per_cell.parquet`, not read; `run_log.json`) | a (numbers); Figure S8 |
| 2 | `figures/_prep/fig2a_functional_state_summaries.py` | the deposited per-cell PROGENy / Hallmark / DoRothEA parquets; step 1; `load_epithelial()` | `$OUTPUT_ROOT/figure2/fig2a_functional/21_{progeny,hallmark,dorothea,flux}_{zscored,means}.csv`, `coverage.csv` | a (numbers); Figure S8 |
| 3 | `01_fig2a_radar_characterization.py` | step 2 `21_*_zscored.csv` | `atlas_radar_characterization_1x4` | a |
| 4 | `figures/_prep/fig2b_cell_cycle_phase.py` | epithelial h5ad `layers/counts`; `load_epithelial()` | `$OUTPUT_ROOT/figure2/fig2b_cell_cycle/{cell_cycle_per_cell.parquet, cell_cycle_phase_proportions.csv, cell_cycle_phase_counts.csv}` | b (numbers). Must exist before step 5, even for c–e. |
| 5 | `02_fig2bcde_cellcycle_composition_paired.py` | `load_epithelial()`; `cell_cycle_phase_proportions.csv` | `row2_b_cell_cycle` (b), `row2_c_composition` (c), `row2_d_paired_sites` (d), `row2_e_paired_treatment_secb_pct` (e); `row2_stats.csv` (n, percentages and tests of b–e); `row2_e_variant_secb_pct_stats.csv` (the panel e test) | b, c, d, e |
| 6 | `figures/_prep/fig2f_consensusov.py` (calls `figures/_prep/fig2f_consensusov.R`) | `hgsc_atlas_final.h5ad` `layers/counts` (rows of the 96 samples); `load_obs()`, `load_epithelial()` | `$OUTPUT_ROOT/figure2/fig2f_consensusov/{pseudobulk_counts.tsv.gz, pseudobulk_metadata.csv, consensusov_calls.csv, gene_id_mapping.csv, fig2f_per_sample_joined.csv}` (single-seed calls; see f below) | f (numbers) |
| 7 | `Rscript figures/_prep/fig2f_consensusov_multiseed.R <…/pseudobulk_counts.tsv.gz> <…/consensusov_calls_multiseed.csv> 25` (run by hand) | `pseudobulk_counts.tsv.gz` | `$OUTPUT_ROOT/figure2/fig2f_consensusov/consensusov_calls_multiseed.csv` (+ `_per_seed.csv`) | f (numbers) |
| 8 | `03_fig2f_subtype_state_boxes.py --no-intermediate` | `fig2f_per_sample_joined.csv`; `consensusov_calls_multiseed.csv` | `atlas_subtype_epitype_grouped_boxes_SecA_SecB`; `atlas_subtype_epitype_grouped_boxes_stats_SecA_SecB.csv` (per state: Kruskal–Wallis over samples and over one sample per patient; medians and n per subtype) | f |
| 9 | `04_fig2g_pseudobulk_volcano.py` | epithelial h5ad `layers/counts`; `load_epithelial()`; `shared/signatures.yml` | `atlas_volcano_secA_secB_pseudobulk`; `$OUTPUT_ROOT/figure2/atlas_volcano_secA_secB_pseudobulk_deseq2.csv` (full DE table) | g |

Steps 1–3, 4–5, 6–8 and 9 are independent chains. Where a file name carries a panel letter, it is the published one.

**Options that the paper does not use.**
- Without `--no-intermediate`, step 8 also draws Intermediate boxes and drops the `_SecA_SecB` suffix. Panel f shows
  SecA and SecB only; the statistics table tests all three states either way.
- `--cohort all_tumours` (steps 6 and 8) runs on every non-healthy sample and writes to `fig2f_consensusov_all/`.
- `figures/_prep/fig2a_scfea_all_cells.py --validate` re-runs the 20,000 cells of the deposited
  `21_epitype_functional_characterization/21_flux_per_cell.parquet` and writes a comparison table instead of the flux.

**Used downstream.** The DE table of step 9 is also read from `$OUTPUT_ROOT/figure2/` by Table S6
(`tables/06_tableS6_secA_secB_signatures.py`) and by the Fig 6h atlas row (`figures/figure6/11_fig6h_concordance.py`,
`figures/figure6/13_fig6_draw_panels.py`). Table S6 also reads the copy in `output/figure2/` inside the repository and
asserts that the two are identical. The final run used one output root for the atlas figures and another for the
Xenium figures, and the table was copied by hand from the first to the second; with a single `OUTPUT_ROOT` no copy is
needed.

## Environment

| Variable | Value in the final run |
|---|---|
| `DATA_ROOT` | the deposit root, which holds `2026_final_atlas/` |
| `OUTPUT_ROOT` | a separate atlas output root, shared with Figure 1 and Figures S1–S9 |

- **R.** `Rscript` must be on `PATH` with E1 (R 4.4.3; consensusOV 1.25.1, edgeR 4.4.2, org.Hs.eg.db 3.20.0) for
  steps 6–7.
- **Python.** Python steps run in E2. Step 1 uses PyTorch on the CPU (`--threads`, default 12), and step 9 uses
  pydeseq2.
- **Seeds.**
  - scFEA uses torch seed 16, set in the script.
  - consensusOV uses seeds 1–25 for the consensus (step 7, which also records the seed-42 call) and seed 42 for the
    single call of step 6.
  - DESeq2 and the point jitter use the `config.yml` seed (42).

## Panel → script

| Panel | What it shows | Drawn by | Numbers computed by |
|---|---|---|---|
| **a** | radar plots, SecA vs SecB, of PROGENy, Hallmark, DoRothEA and scFEA flux (six features per resource, chosen in the plot script) | `01_fig2a_radar_characterization.py` | per cell: `atlas/04_functional/as_run/21_epitype_functional_characterization.py` (PROGENy and DoRothEA by decoupler MLM; Hallmark by scanpy `score_genes`) and `figures/_prep/fig2a_scfea_all_cells.py` (scFEA M168, all 575,366 cells). Per state: `figures/_prep/fig2a_functional_state_summaries.py` (means, z-scored across SecA, Intermediate, SecB and ciliated). |
| **b** | cell-cycle phase of SecA and SecB cells. Scored with scanpy `score_genes_cell_cycle` and the Tirosh S and G2/M lists. | `02_fig2bcde_cellcycle_composition_paired.py` | `figures/_prep/fig2b_cell_cycle_phase.py`; `row2_stats.csv` (panel b) |
| **c** | composition of the epithelium by site (pre-treatment cells pooled per site): primary 99 samples / 81 patients, ascites 40 / 40, metastasis 131 / 89. Composition of primary tumours, pre-treatment vs post-chemotherapy (99 / 81 and 12 / 12). | `02_fig2bcde_cellcycle_composition_paired.py` | the same script; `row2_stats.csv` (panel c) |
| **d** | patient-matched SecB % across sites (n = 14): pre-treatment patients with ascites and at least one other site, ≥ 100 epithelial cells per site; two-sided paired Wilcoxon | `02_fig2bcde_cellcycle_composition_paired.py` | the same script; `row2_stats.csv` (panel d) |
| **e** | patient-matched SecB % before and after chemotherapy (n = 18; ≥ 100 epithelial cells per time point; samples pooled per time point, not site-matched; paired Wilcoxon) | `02_fig2bcde_cellcycle_composition_paired.py` | the same script; `row2_e_variant_secb_pct_stats.csv` and `row2_stats.csv` (panel e) |
| **f** | SecA and SecB % per consensusOV subtype in 96 pre-treatment adnexal samples from 77 patients (≥ 500 cells and ≥ 100 epithelial cells per sample); Kruskal–Wallis per state | `03_fig2f_subtype_state_boxes.py --no-intermediate` | `figures/_prep/fig2f_consensusov.py` (whole-sample raw-count pseudobulk → edgeR TMM → log2(CPM + 1) → Entrez → consensusOV) and `figures/_prep/fig2f_consensusov_multiseed.R` (the final call is the 25-seed consensus) |
| **g** | volcano plot of SecB vs SecA, from pseudobulks paired within tumours (522 profiles from 261 tumours). PyDESeq2, design `~sample + epitype`, Wald test, BH. The 7 + 7 signature genes are highlighted. | `04_fig2g_pseudobulk_volcano.py` | the same script. A sample × state pseudobulk needs ≥ 10 cells and a tumour needs both states; genes need a total count ≥ 10. The 261 tumours are asserted. |

**Numbers in the text from this figure.**

| Text | Value | Where it is computed |
|---|---|---|
| Results | SecA 37.8% S/G2M; SecB 80.6% G1 | `row2_stats.csv`, panel b |
| Results | ascites 68% SecB vs primary 14% and metastasis 17% | `row2_stats.csv`, panel c |
| Results | n = 14 matched patients; median SecB 60.7% in ascites vs 10.2% in primary | `row2_stats.csv`, panel d |
| Results | after chemotherapy, SecA 60% → 53% and SecB 14% → 24% | `row2_stats.csv`, panel c |
| Results | SecB rose in 12 of 18 patients, 58% → 69%, p = 0.11 | `row2_e_variant_secb_pct_stats.csv` |
| Results | SecA medians 76%, 75%, 38% and 37% in the proliferative, immunoreactive, differentiated and mesenchymal subtypes | `atlas_subtype_epitype_grouped_boxes_stats_SecA_SecB.csv` and the per-sample table |

## Not in this folder

- **Prep and loaders.** These are in `figures/_prep/`: `fig2a_scfea_all_cells.py`,
  `fig2a_functional_state_summaries.py`, `fig2b_cell_cycle_phase.py`, `fig2f_consensusov.py`, `fig2f_consensusov.R`,
  `fig2f_consensusov_multiseed.R`, `atlas_final_tables.py` and `nath2021_fix.py`.
  See [`../_prep/README.md`](../_prep/README.md).
- **Per-cell functional scores.** These come from `atlas/04_functional/as_run/21_epitype_functional_characterization.py`
  and are deposited.
- **Manual steps.**
  - Step 7 is run by hand; no script calls it.
  - Each render was copied by hand into the figure layout.
  - The page was assembled in Affinity Designer; no code does this.
- **The signature genes in g.** The 7 + 7 genes highlighted in g (and used for polarization scoring throughout the
  paper) are fixed in `shared/signatures.yml`. They were chosen from the NMF loadings and the DE results, as the Methods
  describe; no script makes the selection.
- **Missing code.** None for this figure. Missing upstream code is listed in [`../../atlas/README.md`](../../atlas/README.md).
