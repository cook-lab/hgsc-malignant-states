# Figure 1: an integrated HGSC atlas and the SecA–SecB axis of secretory polarization

Figure 1 introduces the integrated scRNA-seq atlas (1,980,703 cells from 13 studies) and its epithelial compartment
(575,366 cells). NMF finds a SecA program and a SecB program. The axis between them is the leading axis of variation
within most tumours, and SecA, Intermediate and SecB cells occur within individual CopyKAT clones.

## Run order

**Requires.**
- Deposited atlas outputs under `$DATA_ROOT/2026_final_atlas/` (see [`../../atlas/README.md`](../../atlas/README.md)).
  They were written by the `as_run/` copies in `atlas/`:
  - `hgsc_atlas_final.h5ad`: the final filter of `atlas/02_annotation/as_run/09d_publication_figures.py`. Its
    `X_umap` (panel a) is the UMAP of `atlas/01_preprocess_qc/as_run/03b_refilter_umap.py`.
  - `celltype_h5ad/hgsc_atlas_final_epithelial.h5ad` (`atlas/02_annotation/as_run/9a_atlas_celltypelevel1_v2.py`;
    its counts are read for g)
  - `output/11d_epithelial_nmf/11d_nmf_usage.csv` (`atlas/03_epithelial_nmf/as_run/11d_epithelial_nmf.py`)
  - `output/fig_secretory_polarization/data/meta.parquet`, for the epithelial UMAP1/UMAP2
    (`atlas/03_epithelial_nmf/as_run/fig_secretory_polarization_00_prepare_data.py`)
  - `output/19_cnv/per_sample/` (CopyKAT, `atlas/05_cnv/as_run/`)
  - the inputs of the Nath 2021 fix in `pre-processing and integration/` (see [`../_prep/README.md`](../_prep/README.md))
- `figures/_prep/nath2021_fix.py` → `figures/_prep/atlas_final_tables.py`, which write
  `$OUTPUT_ROOT/_prep_caches/atlas_final_{obs,epithelial}.parquet`. Both build their caches the first time a script
  below loads them.

Renders go to `$OUTPUT_ROOT/figures/figure1/` as `.svg` and `.png`.

| Step | Script | Reads | Writes | Panel(s) |
|---|---|---|---|---|
| 1 | `01_fig1a_atlas_level1_umap.py` | `load_obs()`: UMAP1/UMAP2, `celltype_level1` | `atlas_level1_umap` | a |
| 2 | `02_fig1bcdek_atlas_study_contribution.py` | `load_obs()`: study, patient, metastatic site, treatment, level-1 type; `load_epithelial()`: state | `atlas_study_contribution_with_epi_composition` (five stacked panels); `$OUTPUT_ROOT/figure1/atlas_study_contribution_values.csv` (the plotted values per study) | b, c, d, e, k |
| 3 | `03_fig1f_epithelial_secretory_ciliated_umap.py` | `load_epithelial()`: epithelial UMAP, state | `atlas_epithelial_secretory_ciliated_umap` | f |
| 4 | `figures/_prep/fig1g_polarization_pc_alignment.py` | epithelial h5ad `layers/counts`; `load_epithelial()`; CopyKAT `copykat_{prediction,subclones,subclone_qc}.csv` per sample | `$OUTPUT_ROOT/_prep_caches/fig1g_polarization_pc_alignment.csv` | g (numbers) |
| 5 | `04_fig1g_polarization_pc_alignment.py` | the step 4 cache | `atlas_polarization_pc_alignment`; `$OUTPUT_ROOT/figure1/atlas_polarization_pc_alignment_data.csv` (rank → number of tumours) | g |
| 6 | `05_fig1h_secA_nmf_factor_umap.py` | `load_epithelial()`: Factor_3 | `atlas_SecA_nmf_factor_umap` | h |
| 7 | `06_fig1i_secB_nmf_factor_umap.py` | `load_epithelial()`: Factor_2 | `atlas_SecB_nmf_factor_umap` | i |
| 8 | `07_fig1j_program_factor_umaps.py` | `load_epithelial()`: Factor_5, 6, 4, 9 | `atlas_program_{g2m,ap1,ifn,emt}_umap` | j |
| 9 | `08_fig1lm_clone_composition_and_polarization.py` | `load_epithelial()`: state, Factor_2; CopyKAT `copykat_subclones.csv`, `barcodes.csv`, `DONE.txt` | `atlas_clone_polarization_and_composition_lm` (l and m stacked on one clone axis, as published); `$OUTPUT_ROOT/figure1/atlas_copykat_coverage.csv` (the barcode-coverage check behind an assertion) | l, m |

The steps can run in any order, except that step 4 must run before step 5.

The output file names carry no panel letters. Step 9 writes l and m as one stacked render, which the figure uses for
both panels; the page was assembled by hand from these renders.

## Environment

| Variable | Value in the final run |
|---|---|
| `DATA_ROOT` | the deposit root, which holds `2026_final_atlas/` |
| `OUTPUT_ROOT` | a separate atlas output root, shared with Figure 2 and Figures S1–S9 |

No other variables are read.
- **Caches.** The loaders reuse the caches in `$OUTPUT_ROOT/_prep_caches/` when they exist. Delete them after
  changing any atlas input.
- **Software.** Python E2 (3.12.9).
- **Seed.** `config.yml` (42).

## Panel → script

| Panel | What it shows | Drawn by | Numbers computed by |
|---|---|---|---|
| **a** | UMAP of 1,980,703 cells coloured by the 13 level-1 lineages | `01_fig1a_atlas_level1_umap.py` | `figures/_prep/atlas_final_tables.py` (cells); the UMAP is `X_umap` of the deposited atlas (`atlas/01_preprocess_qc/as_run/03b_refilter_umap.py`) |
| **b** | cells per study (log scale), with the number of patients per study | `02_fig1bcdek_atlas_study_contribution.py` | the same script. It prints the cells, studies and patients and writes `atlas_study_contribution_values.csv`. |
| **c** | metastatic-site composition per study (primary, metastasis, ascites, healthy) | `02_fig1bcdek_atlas_study_contribution.py` | the same script; the Vázquez-García ascites samples carry the site written into the atlas by `atlas/02_annotation/as_run/fix_metastatic_site_vazquez.py` |
| **d** | treatment-status composition per study | `02_fig1bcdek_atlas_study_contribution.py` | the same script, on the labels after the Zhang 2022 fix (`atlas/02_annotation/as_run/12a_fix_zhang2022_treatment_status.py`) and the Nath 2021 fix. The three post-chemotherapy maintenance groups are drawn as one "Other" segment. The script's treatment categories do not include `post-chemotherapy_niraparib_pembro`, the status of one Nath 2021 sample after the fix, so those cells are not counted in d. |
| **e** | level-1 cell-type composition per study | `02_fig1bcdek_atlas_study_contribution.py` | the same script |
| **f** | epithelial UMAP (575,366 cells), secretory vs ciliated (6,678) | `03_fig1f_epithelial_secretory_ciliated_umap.py` | the counts are asserted in `figures/_prep/atlas_final_tables.py` |
| **g** | per tumour, the rank of the principal component most correlated with the polarization score (250 samples) | `04_fig1g_polarization_pc_alignment.py` | `figures/_prep/fig1g_polarization_pc_alignment.py` (see the note below) |
| **h** | epithelial UMAP coloured by Factor 3 (SecA) usage | `05_fig1h_secA_nmf_factor_umap.py` | deposited NMF (`atlas/03_epithelial_nmf/as_run/11d_epithelial_nmf.py`) |
| **i** | epithelial UMAP coloured by Factor 2 (SecB) usage | `06_fig1i_secB_nmf_factor_umap.py` | deposited NMF |
| **j** | four UMAPs: Factor 5 (G2M / cell cycle), Factor 6 (AP-1 / early stress), Factor 4 (IFN), Factor 9 (EMT) | `07_fig1j_program_factor_umaps.py` | deposited NMF |
| **k** | SecA / Intermediate / SecB / ciliated proportion per study | `02_fig1bcdek_atlas_study_contribution.py` | state rule in `figures/_prep/atlas_final_tables.py` |
| **l** | state composition of each CopyKAT subclone (≥ 20 cells) in three representative tumours, one monoclonal and two polyclonal | `08_fig1lm_clone_composition_and_polarization.py` | the same script (CopyKAT calls from `atlas/05_cnv/as_run/`) |
| **m** | Factor 2 usage per subclone (violins) for the same tumours; dashed lines at the atlas Factor 2 p50 (Intermediate) and p75 (SecB) thresholds | `08_fig1lm_clone_composition_and_polarization.py` | the same script |

**How panel g is computed.** For each sample with at least 150 CopyKAT-aneuploid epithelial cells:
- normalise and log-transform the counts; score polarization as scanpy `score_genes` SecB minus SecA, with the
  organoid SecA and SecB gene sets listed in the script;
- take the top 2,000 HVGs (seurat flavour), excluding those gene sets, the NMF marker genes and the cell-cycle genes;
- compute 30 principal components;
- record the PC that correlates best with the polarization score (`pol_top_pc`).

Panel g plots the distribution of that rank. The cache also holds the shares of the PC variance explained by
polarization, library size, cell cycle and clone (`pc_*` columns); the paper does not use them.

**Numbers in the text from this figure.**

| Number | Where it is computed |
|---|---|
| 1,980,703 cells from 158 patients in 13 studies | printed by `02_fig1bcdek_atlas_study_contribution.py` |
| 362 samples and 12 anatomic sites | distinct values in `load_obs()`; also the Total row of Table S1 (`tables/01_tableS1_atlas_studies.py`) |
| 13 lineages and 51 cell types | `load_obs()`; the 51 is asserted in `tables/02_tableS2-4_annotation_markers.py` |
| 575,366 epithelial cells and 6,678 ciliated cells (1.2%) | asserted in `figures/_prep/atlas_final_tables.py` |
| PC1 in 54% of samples, PC1–3 in 82%, 250 samples | printed by `04_fig1g_polarization_pc_alignment.py` and `figures/_prep/fig1g_polarization_pc_alignment.py` |

## Not in this folder

- **Prep and loaders.** `figures/_prep/fig1g_polarization_pc_alignment.py`, `figures/_prep/atlas_final_tables.py` and
  `figures/_prep/nath2021_fix.py`.
- **Upstream atlas steps.** Annotation, NMF and CopyKAT ran as the `as_run/` scripts of `atlas/`; their outputs are
  deposited. `atlas/03_epithelial_nmf/01_epithelial_nmf.py` re-implements the NMF step and is checked against the
  deposited NMF by `atlas/03_epithelial_nmf/01b_verify_against_deposited.py`; the figures read the deposited NMF.
- **The epithelial UMAP.** The coordinates used in f and h–j were computed by
  `atlas/02_annotation/as_run/9e_celltypelevel2_umaps.py` (stored as `X_umap_local`) and saved as UMAP1/UMAP2 in
  `meta.parquet`. The epithelial object on disk no longer holds `X_umap_local`, so the figures read the coordinates
  from `meta.parquet`.
- **Manual steps.** Each render was copied by hand from `$OUTPUT_ROOT/figures/figure1/` into the figure layout. The page
  was assembled in Affinity Designer; no code does this.
- **A file another script depends on.** `tables/05_tableS5_epithelial_nmf_programs.py` reads the factor labels
  (`PANELS`) from `07_fig1j_program_factor_umaps.py` by path. Keep that file's name.
- **Missing code.** None for this figure. Missing upstream code is listed in [`../../atlas/README.md`](../../atlas/README.md).
