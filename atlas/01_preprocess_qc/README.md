# atlas/01_preprocess_qc: data preparation, QC and integration

This stage turns the raw count matrices of the 13 public studies into the integrated atlas object
`hgsc_atlas_scanvi.h5ad`. The steps are:
1. prepare each study
2. concatenate the studies on 14,470 shared genes, then cell QC and Scrublet doublet removal
3. HVG selection
4. CellAssign
5. scVI → scANVI integration
6. a second doublet pass and the atlas UMAP

The stage is heavy: steps 4–6 of the run order need a GPU, and steps 3 and 7 a high-memory machine. The outputs are
deposited, so no figure needs this stage to be re-run.

## Run order

Paths without a root are under `$DATA_ROOT/2026_final_atlas/`. `…/integration/` stands for
`$OUTPUT_ROOT/01_preprocess_qc/integration/`.

| Step | Script | Reads | Writes | Used by |
|---|---|---|---|---|
| 1 | `as_run/01_dataprep.ipynb` | `raw/<study>/` (per-study count files); `atlas_metadata.xlsx` (curated sample sheet) | `processed/<study>.h5ad`, one per study (13) | step 2; `figures/_prep/nath2021_fix.py` (`nath_2021.h5ad`) |
| 2 | `as_run/02_concat_qc_doublets.ipynb` | `processed/<study>.h5ad` | `processed/atlas_concat_counts_only_X.h5ad` (2,731,632 cells, before QC); `processed/atlas_concatenated_filtered.h5ad` (2,326,532 cells, after QC and Scrublet) | steps 3, 4, 7; Figure S1 (`figures/supplementary/00_figS1_extract_atlas_obs.py`) |
| 3 | `03_preprocess_hvg.py` | `processed/atlas_concatenated_filtered.h5ad` | `…/integration/anndata/preprocessed.h5ad` (4,000 seurat_v3 HVGs, batch = sample) | steps 5, 6 |
| 4 | `04_cellassign.py` (+ `shared/cellassign_markers_v3.csv`: 82 genes × 16 types) | `processed/atlas_concatenated_filtered.h5ad` | `…/integration/cellassign/{predictions,probabilities}.csv` | steps 6, 7 |
| 5 | `05a_integrate_scvi.py` | `preprocessed.h5ad` | `…/integration/models/scvi/`; `…/integration/embeddings/scvi/embedding.npz` (not read downstream) | step 6 (scANVI starts from this scVI model) |
| 6 | `05b_integrate_scanvi.py` | `preprocessed.h5ad`; CellAssign predictions (the scANVI seed labels); the scVI model | `…/integration/models/scanvi/`; `…/integration/embeddings/scanvi/embedding.npz` | step 7 |
| 7 | `07_finalize.py --method scanvi` | `processed/atlas_concatenated_filtered.h5ad`; CellAssign predictions; the scANVI embedding | `…/integration/anndata/integrated_scanvi.h5ad`, with `obsm["X_scanvi"]` and obs `celltype_pred`, `celltype_probability` and `celltype_confident` from CellAssign | step 8 |
| 8 | `as_run/03b_refilter_umap.py` | `pre-processing and integration/20260213_integration/output/anndata/integrated_scanvi.h5ad` (the deposited output of steps 3–7) | `hgsc_atlas_scanvi.h5ad` (config `atlas_scanvi`) | `atlas/02_annotation` |

**QC rules.**
- Step 2 keeps cells with a study label from samples whose `histological_subtype` is Serous. It removes cells with
  fewer than 500 UMIs or fewer than 300 genes, then samples with fewer than 500 cells. Scrublet then runs per sample
  (samples with at least 1,000 cells; expected doublet rate 0.8% per 1,000 cells; seed 42). Cells with a score of 0.3
  or more are removed, and samples left with fewer than 500 cells are removed again.
- Step 8 removes the cells Scrublet called doublets and cells with a Scrublet score of 0.25 or more. It then
  recomputes neighbours (n = 10) and the UMAP (min_dist 0.2) on `X_scanvi`. The 0.25 cut is the line drawn in
  Figure S1c.
- Artifact clusters are removed later, at level-1 annotation (`../02_annotation`).

**Cell-type labels.** `celltype_pred` and `celltype_probability` in the final atlas are the CellAssign label and
probability that step 7 copies from step 4. Figure S4a plots them.

**Paths.**
- Steps 1, 2 and 8 are `as_run/` copies and use the paths of the atlas project folder (`…/2026_final_atlas/raw/`,
  `processed/` and `atlas_metadata.xlsx`; the Hornburg 2021 cell reads `…/2026_final_atlas/hornburg/EGAF*/`). In the
  deposit these sit under `pre-processing and integration/` (Hornburg 2021 in `raw/hornburg_2021/`).
- Steps 3, 4 and 7 read `$DATA_ROOT/2026_final_atlas/processed/atlas_concatenated_filtered.h5ad`. The deposited copy
  is `pre-processing and integration/processed/atlas_concatenated_filtered.h5ad`; link it to `processed/` first.
- Steps 3–7 write under `$OUTPUT_ROOT/01_preprocess_qc/integration/`. The run that produced the deposit is
  `pre-processing and integration/20260213_integration/output/`; step 8 reads its `anndata/integrated_scanvi.h5ad`.

**Integration in the paper.** The paper uses scVI → scANVI with CellAssign seed labels. `07_finalize.py` ran with
`--method scanvi`; its other `--method` choices are not used.

**Environments.**
- Steps 1–2 ran in E3 (Python 3.12.3). They use Scrublet as built into scanpy, and mygene, which needs network access.
- Steps 3–7 ran in E4, a GPU cluster with scvi-tools 1.4.1. Steps 4–7 set their seed from `config.yml` (42).
- Step 8 ran in E2 (Python 3.12.9).

GPU training is not bit-identical across hardware, so the deposited objects are the reference for everything
downstream.

## `as_run/` and the cluster scripts

`as_run/` holds trimmed copies of what ran in the atlas project folder `202601_final_atlas`. The cells and code that
wrote the deposited files are kept as they ran. QC plots, check and inspection cells, empty cells and two cells whose
files a later cell overwrote (the first Nath 2021 loader and the first Vázquez-García write) were removed, and the
notebook outputs were cleared. Each copy opens with a note on what was removed. The untrimmed notebooks and script
stay in `202601_final_atlas`. The copies keep that project's paths.

In `01_dataprep.ipynb` the "new nath" cell stopped with an error after the concatenation. The next cell repeats its
metadata merge and save, and that save wrote the deposited `nath_2021.h5ad`. The Hornburg 2021 cell has no execution
count in the notebook, but it is the only code that writes `hornburg_2021.h5ad`.

Steps 3–7 have no `as_run/` copy: the files in this folder are the scripts that ran on the GPU cluster, with their
paths resolved through `config/` and their QC plots, metadata JSON and UMAP previews removed. They have not been
re-run from this repository.

## Missing code

None for this stage. The 13 public datasets were downloaded from their accessions (Table S1) into
`pre-processing and integration/raw/<study>/`; no script does this.
