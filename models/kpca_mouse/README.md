# models/kpca_mouse: KPCA.B orthotopic syngeneic mouse tumour

Trimmed copies of the scripts that ran, from the mouse project folder `202608_kpca_secb_mouse` (see
[`../README.md`](../README.md)). The code that wrote the outputs the paper uses is kept as it ran; plots and analyses
not in the paper were removed. The untrimmed scripts stay in the project folder. `18_infercnv.R` and
`33_consolidate_handoff.py` had nothing to remove and are unchanged copies. The originals are in `scripts/`;
`41b_export_1p8in.py` sat in the project's `output/41_fig3_umaps/` folder.

Seven steps (`18a`, `21a`, `21`, `22`, `22b`, `27b`, `27c`) ran as inline code in the analysis session and were
never saved as files. They are recovered from the session log; each opens with a note on when it ran and what it
printed. `upstream/` records the Cell Ranger multi configuration and the CellBender command. Re-run on 2026-10-22
from the step-01 object (15b with `--z 0.75`, then 18a, 19, 21a, 21, 22, 22b; and 27b, 27c from `kpca_annotated.h5ad`),
they reproduce `output/21_working/kpca_working.h5ad` and the epithelial columns of the handoff object exactly (same
cells, identical clusters, UMAP and scores).

The scripts take one KPCA.B tumour profiled with 10x Flex from Cell Ranger / CellBender output through QC, doublet
removal, inferCNV and cell typing (20,598 cells, 12,429 tumour cells), de novo NMF of the tumour cells and its
alignment to the human NMF factors. The result is one handoff object.

The repository reads these outputs from the deposited project folder, `$DATA_ROOT/202608_kpca_secb_mouse/output/`:
- the handoff object `33_handoff/kpca_final_annotated.h5ad`, read by `figures/figure3/04_…` (Fig 3j),
  `figures/figure3/05_…` (Fig 3k; both take `$KPCA_H5AD` as an override) and
  `figures/supplementary/16_figS10_kpca_mouse_vs_human_programs.py` (Figure S10);
- the orthologue map `02_orthologs/ortholog_map.csv` (Fig 3j, override `$KPCA_ORTHO`; Figure S10) and
  `02_orthologs/02_coverage_summary.csv` (Figure S10);
- the step-31 and step-32 tables in `31_fulltumour_nmf/` and `32_secb_align/` (Figure S10).

Fig 3i is drawn here, by `41b_export_1p8in.py`.

Paths below are relative to the mouse project folder. `<atlas>` is the atlas project folder `202601_final_atlas`.

## Run order

| Step | Script | Reads | Writes | Used by |
|---|---|---|---|---|
| – | `00_setup.py` | – | (imported by every step) paths, human gene lists, `SEED = 42`, `N_HVG = 3000`, NMF k = 4–12 × 20 replicates, `N_PERM = 10,000`, ortholog helpers | all; step 02 uses its gene lists for the coverage report |
| 1 | `upstream/cellranger_multi_config.csv` (Cell Ranger multi 9.0.0, 10x Cloud), then `upstream/cellbender_command.sh` (CellBender 0.3.2) | Flex pool `FlexPool2` FASTQs | KPCA.B raw-data folder `2025_kpca/KPCA-B`: `sample_raw_feature_bc_matrix.h5`, `sample_filtered_feature_bc_matrix.h5`, `cell_types.csv`, `cellbender/KPCA-B_background_removed_filtered.h5` | 01 |
| 2 | `01_load_qc_kpca.py` | the raw-data folder `2025_kpca/KPCA-B` (step 1) | `output/01_kpca_qc/{kpca_all.h5ad, kpca_epithelial.h5ad, 01_*.csv}` | 15b; 02 |
| 3 | `02_orthologs.py` | MGI `HOM_MouseHumanSequence.rpt` (downloaded once to `data/external/`); `output/01_kpca_qc/kpca_epithelial.h5ad` (gene space); `<atlas>/output/11e_nmf_characterization/11e_gene_classification.csv`; `<atlas>/output/11d_epithelial_nmf/11d_nmf_loadings.csv` | `output/02_orthologs/{ortholog_map.csv, 02_coverage_summary.csv}` | 22, 31, 32; Fig 3j (`$KPCA_ORTHO`); Figure S10a, c |
| 4 | `15b_apply_doublet_cut.py --z 0.75 --genes 2500 --min-lineages 2` | `output/01_kpca_qc/kpca_all.h5ad` | `output/15_doublet_cleanup/{kpca_all_clean.h5ad, kpca_epithelial_clean.h5ad}` (22,437 cells; no longer on disk) | 18a, 18, 19, 21a, 21 |
| 5 | `18a_infercnv_gene_positions.R` | the step-15b object (Ensembl IDs); `org.Mm.eg.db` | `output/18_infercnv/gene_positions.tsv` (19,123 genes) | 18, 19 |
| 6 | `18_infercnv.R` | the step-15b object; `gene_positions.tsv` | `output/18_infercnv/{run/, 18_cnv_signal_per_cell.csv}` (no longer on disk) | the tumour-cluster call (no later step reads it) |
| 7 | `19_percell_identity_confidence.py` | the step-15b object; `gene_positions.tsv` | `output/19_identity_confidence/19_percell.csv` (per-cell `cnv_score`; no longer on disk) | 21a |
| 8 | `21a_add_cnv_score.py` | the step-15b object; `19_percell.csv` | the step-15b object, in place, with `cnv_score` | 21 |
| 9 | `21_working_object.py` | the step-15b object | `output/21_working/kpca_working.h5ad`: 20,634 cells, 16 Leiden clusters, `X_umap`, `cnv_score`, `compartment` (on disk) | 22, 24, 27 |
| 10 | `22_secb_axis_export.py`, then `22b_secb_axis_ucell.R` | `kpca_working.h5ad`; `ortholog_map.csv` | `output/22_secb_axis/{seca_genes.csv, secb_genes.csv, kpca_tumour.h5ad, tumour_X.mtx, genes.csv, cells.csv, 22_ucell_scores.csv}` (not on disk) | 24 |
| 11 | `24_secb_whole_umap.py` | `output/21_working/kpca_working.h5ad`; `output/22_secb_axis/*` | `output/24_whole_umap/24_allcell_ucell.csv` | 27 |
| 12 | `27_final_annotation_split.py` | `output/21_working/kpca_working.h5ad`; `output/24_whole_umap/24_allcell_ucell.csv` | `output/26_final_annotation/kpca_annotated.h5ad` (20,598 cells) | 27b, 31, 32, 33 |
| 13 | `27b_epithelial_umap.py`, then `27c_epithelial_split.py` | the tumour cells of `kpca_annotated.h5ad` (12,429) | `epi_umap.h5ad`, then `epi_split.h5ad`: epithelial UMAP (`X_umap_epi`), epithelial Leiden (`epi_leiden`), SecA / SecB UCell from 5 + 4 genes. Both were written to the session's temporary folder and are gone | 33 |
| 14 | `31_fulltumour_nmf.py` (imports `04_mouse_nmf.py`) | `output/26_final_annotation/kpca_annotated.h5ad` (`is_tumour`); `ortholog_map.csv` | `output/31_fulltumour_nmf/{31_program_loadings.csv, 31_program_usage.csv, 31_program_usage_raw.csv, 31_program_annotation.csv, 31_program_summary.csv, kpca_fulltumour_programs.h5ad, 04_k_selection.csv}` | 32, 33; Figure S10c, d |
| 15 | `32_identify_secb_align.py` (imports `05_crossspecies_alignment.py`) | `output/31_fulltumour_nmf/*`; `kpca_annotated.h5ad` (`polarization`); `<atlas>/output/11d_epithelial_nmf/11d_nmf_loadings.csv`; `ortholog_map.csv` | `output/32_secb_align/{32_alignment_stats.csv, 32_verdict.csv}` | Figure S10c |
| 16 | `33_consolidate_handoff.py` | `kpca_annotated.h5ad`; `epi_split.h5ad`; `output/31_fulltumour_nmf/31_program_usage.csv` | `output/33_handoff/{kpca_final_annotated.h5ad, kpca_final_metadata.csv}` (deposited) | 41b; Fig 3j, k; Figure S10 |
| 17 | `41b_export_1p8in.py` | `output/33_handoff/kpca_final_annotated.h5ad` | `output/41_fig3_umaps/41_mouse_umap_compartment_1p8in.{svg,png}`, which is Fig 3i | Fig 3i |

`04_mouse_nmf.py` and `05_crossspecies_alignment.py` are helpers. Steps 31 and 32 import their functions; the trimmed
copies keep only those functions, not their own step-04 and step-05 runs.

Order constraints:
- Step 02 needs only step 01, but must run before 22, 31 and 32, which read `ortholog_map.csv`.
- Step 15b ran with `--z 0.75`. The script's default (z 1.0) keeps 23,024 cells, not the 22,437 that step 21 started
  from.
- 21a writes into the step-15b object, so it runs after 19 and before 21.
- `27` writes into the folder `output/26_final_annotation/`; it replaced an earlier step 26, which is not included.

## What each step does

- **01.** Takes the CellBender counts, restricted to the 28,973 Cell Ranger cell calls. Keeps cells within 5 MADs of
  the median number of genes (log scale, at least 200 genes). Then clusters (Leiden) and labels compartments with
  mouse markers.
  The Cell Ranger automated annotation (`cell_types.csv`) is only attached to the object (`cr_coarse_cell_type`,
  `cr_fine_cell_type`); it is not used for the labels.
- **02.** Builds the human → mouse orthologue map from MGI, keeping the orthologue class. It falls back to case
  conversion where MGI has no entry, and reports signature coverage. Figure 3j uses the 1:1 orthologues of the 7 + 7
  signature genes.
- **15b.** Removes doublets by a rule: at least 2 lineage programs (epithelial, endothelial, fibroblast, myeloid,
  lymphoid) each with z > 0.75 (run with `--z 0.75`; the script's default is 1.0), together with ≥ 2,500 detected
  genes.
- **18a.** Maps the Ensembl IDs of the step-15b object to mouse chromosome positions with `org.Mm.eg.db`, keeping
  chromosomes 1–19, X and Y.
- **18.**
  - Runs inferCNV with immune clusters as the diploid reference and up to 400 cells per group.
  - Compares candidate tumour clusters with epithelial and fibroblast anchor clusters. Cluster numbers in the script
    are those of the step-15b object.
  - Writes a per-cell CNV signal. This is the basis of "inferCNV-confirmed" tumour cells. No later script reads the
    file; the per-cell `cnv_score` of the objects comes from step 19.
- **19.** A per-cell CNV score for all cells: expression residuals against the immune clusters, ordered by gene
  position and smoothed per chromosome (window 100), correlated with the profile of the epithelial anchor clusters.
  The script also scores competing identities per cell; only `cnv_score` is carried forward.
- **21a, 21.** 21a adds `cnv_score` to the step-15b object. 21 removes the old Leiden clusters 6 and 15 (the doublet
  band; 22,437 → 20,634 cells), then recomputes variable genes, PCA, neighbours, the UMAP of Fig 3i and Leiden
  (resolution 1.0, 16 clusters).
- **22, 22b.** Re-embed the tumour clusters of step 21 (8, 9 and 11), map the human 8 + 7 SecA / SecB genes to mouse
  orthologues and score them with R UCell. Step 24 takes its gene lists from here.
- **24.**
  - Scores SecA and SecB with a Python UCell (maxRank 1500, log-normalised matrix) on every cell, and checks it
    against the saved step-22 tumour scores. `figures/figure3/04_…` uses the same UCell implementation for Fig 3j.
  - `27` copies the scores into the object as `SecA_UCell`, `SecB_UCell` and `polarization`; step 32 uses
    `polarization`.
- **27 (final cell types).**
  - Drops Leiden cluster 1 (36 cycling granulocytes), which takes 20,634 cells to 20,598.
  - Sets `celltype`. Tumour is Leiden clusters 5, 6, 8, 9, 11, 12 and 15: 12,429 cells, `is_tumour`. The other types
    are Fibroblast, Endothelial, Macrophage, Monocyte, Neutrophil, Dendritic, T/NK and B.
- **31 (de novo NMF).**
  - Fits the consensus NMF of `04_mouse_nmf.py` on all 12,429 tumour cells: 3,000 variable genes, k = 4–12, and 20
    replicates per k on 80% cell subsamples.
  - Chooses k by stability; this gives k = 4, programs mP1–mP4. Per-cell usage is computed by NNLS.
  - Flags each program for cycling, ambient and stress genes. Heavy.
- **32 (alignment to human).**
  - Extends the program loadings to all genes by NNLS with the usage fixed.
  - Correlates them (Spearman, cosine) with each of the 10 human 11d factors and with the human polarization loading,
    over shared orthologues.
  - p values come from a 10,000-permutation gene-label null, with Benjamini–Hochberg correction across the
    4 × 11 = 44 tests (`perm_null_pvalues` in `05_crossspecies_alignment.py`). These are the FDR values of Figure S10c.
- **27b, 27c (epithelial split).** 27b computes the epithelial UMAP of the 12,429 tumour cells (3,000 variable
  genes, 30 PCs, 15 neighbours, seed 0). 27c clusters the same cells (Leiden, resolution 0.5) and scores SecA (Bcam,
  Mecom, Pbx1, Sox17, Rcn2) and SecB (Krt19, Krt7, Lcn2, Prss22) with the Python UCell of step 24.
- **33 (handoff object).** Adds to the step-27 object:
  - `epi_umap_1/2` and `epi_subcluster`, the epithelial UMAP and clusters of Fig 3j, 3k, S10b and S10d;
  - `SecB_pole` (epithelial cluster 5) and `SecB_high_top25`;
  - `SecA_UCell_epi`, `SecB_UCell_epi` and `polarization_epi` (the 5 + 4-gene scores);
  - `NMF_mP1`–`NMF_mP4`, and `final_identity`.

  `X` is log-normalised. Non-tumour cells have NA in the epithelial score, UMAP and NMF columns, and
  `epi_subcluster = "non-epithelial"`.
- **41b (Fig 3i).**
  - Draws the all-cell UMAP (`X_umap`) as a 1.8 in square, coloured by `celltype`, with Tumour shown as
    "Epithelial". Despite the file name, it does not use `compartment`.
  - Has no in-plot legend and uses seed 0 for the point order.
  - Fig 3j and 3k are drawn by `figures/figure3/04_…` and `05_…`.

## Paths in the scripts

Edit these to run outside the lab file server. All other steps take their paths from `00_setup.py` (`S.OUTPUT`).

| Script | Path constants |
|---|---|
| `00_setup.py` | `BASE` (mouse project folder), `KPCA_DIR` (raw-data folder `2025_kpca/KPCA-B`), `ATLAS` (atlas project folder) |
| `18_infercnv.R` | `PROJ` |
| `18a_…`, `21a_…`, `21_…`, `22_…`, `22b_…` | paths relative to the mouse project folder: run them from it |
| `27b_…`, `27c_…` | the temporary session folder for `epi_umap.h5ad` and `epi_split.h5ad` (the same folder as `SCRATCH` in 33); the other paths are relative to the mouse project folder |
| `33_consolidate_handoff.py` | `SCRATCH`, the temporary session folder that held `epi_split.h5ad` |
| `41b_export_1p8in.py` | `H5AD`, `OUTDIR` |

## Panels and Methods numbers fed by this folder

| Item | Steps |
|---|---|
| Fig 3i | upstream → 01 → 15b → 18a → 18, 19 → 21a → 21 → 02 → 22, 22b → 24 → 27 → 27b, 27c → 31 → 33 → 41b |
| Fig 3j | as Fig 3i up to 33, plus 02 → `figures/figure3/04_fig3j_kpca_polarization_7plus7.py` |
| Fig 3k | as Fig 3i up to 33 → `figures/figure3/05_fig3k_kpca_krt7_umap.py` |
| Figure S10a (orthologue expression) | 02, 33 (+ Table S6 % expressing from `tables/06_…`) → `figures/supplementary/16_figS10_…` |
| Figure S10b (epithelial UMAPs of SecA / SecB genes) | 33 (`epi_umap_1/2`, `X`) → `16_figS10_…` |
| Figure S10c (loading correlations, FDR) | 02 → 31 → 32 → `16_figS10_…` (recomputes the correlations and asserts that they match `32_alignment_stats.csv`) |
| Figure S10d (program usage on the epithelial UMAP) | 31 (`31_program_usage.csv`) + 33 → `16_figS10_…` |
| Methods: 20,598 cells after ambient-RNA removal and doublet filtering | upstream (CellBender) → 01 → 15b → 21a → 21 → 27 |
| Methods: 12,429 inferCNV-confirmed tumour cells | 18, 19 → 21a → 21 → 27 (Tumour clusters) → 33; asserted in `figures/figure3/04_…`, `05_…` and the S10 script |
| Methods: tumour cells scored with mouse orthologues of the polarization signature | 02 → `figures/figure3/04_…` |
| Methods: NMF programs compared with human factors by loading correlation, permutation null; BH within the S10 loading correlations | 31 (+ 04) → 32 (+ 05) |

## Missing code

None. The steps that ran as inline code (18a, 21a, 21, 22, 22b, 27b, 27c) are recovered from the analysis session
log, and the Cell Ranger multi configuration and CellBender command are recorded in `upstream/` (from the run's web
summary and the CellBender log in `2025_kpca/KPCA-B`).
