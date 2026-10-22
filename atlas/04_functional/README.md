# atlas/04_functional: per-cell pathway, gene-set and transcription-factor scores

This stage scores every epithelial cell of the atlas with:
- PROGENy pathway activity (decoupler, multivariate linear model)
- Hallmark gene-set scores (scanpy `score_genes`)
- DoRothEA A/B/C transcription-factor activity (decoupler, multivariate linear model)

Fig 2a and Figure S8 average these per epithelial state.

## Run order

Paths without a root are under `$DATA_ROOT/2026_final_atlas/`.

| Step | Script | Reads | Writes (deposited) | Used by |
|---|---|---|---|---|
| 1 | `as_run/21_epitype_functional_characterization.py` | `celltype_h5ad/hgsc_atlas_final_epithelial.h5ad`; `output/11d_epithelial_nmf/11d_nmf_usage.csv` (Factor 2 → state); `tools/scFEA/data/` | in `output/21_epitype_functional_characterization/`: `21_{progeny,hallmark,dorothea}_per_cell.parquet` and the per-state `21_{progeny,hallmark,dorothea}_zscored.csv`; scFEA flux of a 20,000-cell subset, `21_flux_per_cell.parquet` and `21_flux_zscored.csv` | the per-cell PROGENy, Hallmark and DoRothEA scores: `figures/_prep/fig2a_functional_state_summaries.py` (Fig 2a) and `figures/supplementary/14_figS8_functional_heatmap.py` (Figure S8) |

**What the figures use.**
- Fig 2a and Figure S8 read only the per-cell PROGENy, Hallmark and DoRothEA parquets. They recompute the per-state
  summaries from them, using the final-atlas states (`figures/_prep/atlas_final_tables.py`). The per-state
  `21_*_zscored.csv` tables of this step are not read.
- The flux in Fig 2a and Figure S8 comes from `figures/_prep/fig2a_scfea_all_cells.py`, which runs scFEA (human M168
  model) on all 575,366 epithelial cells with this step's network, loss and hyperparameters. The 20,000-cell flux of
  this step is read only by that script's optional `--validate` check.

**Software.** Python (E2): decoupler for PROGENy, DoRothEA and the Hallmark gene sets; torch for scFEA.

## `as_run/`

The copy is trimmed from the script that ran in the atlas project folder `202601_final_atlas`, and keeps that
project's paths. The code that wrote the per-cell parquets, the per-state `21_*_zscored.csv` tables and the
20,000-cell scFEA flux is kept as it ran. The heatmaps, violins and radars, the means, radar, balance and summary
outputs and the closing report were removed; the copy opens with a note on what was removed. The untrimmed script
stays in `202601_final_atlas`.

## Missing code

None for this stage.
