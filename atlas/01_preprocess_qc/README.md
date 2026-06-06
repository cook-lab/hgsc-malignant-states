# atlas/01_preprocess_qc — preprocessing, QC, integration (raw → integrated atlas)

The **full preprocessing chain from raw per-study data through the integrated atlas object**
that is the entry point for all downstream analysis. The integration sub-sequence
(steps 02–07) consists of the **original cluster scripts** (`atlas_00..05`), path-centralised
into this repo — NOT reconstructions. The scVI/scANVI steps are GPU-heavy and **optional to
re-run**: the integrated atlas object is deposited as a config entry object, so downstream work
can start from it.

## Ordered scripts

| # | Script | Source (authoritative) | Inputs | Outputs |
|---|--------|------------------------|--------|---------|
| 01 | `01_dataprep.py` | (repo; pre-existing) | per-study raw matrices (`DATA_ROOT/2026_final_atlas/raw/<study>/`), `atlas_metadata.xlsx` | harmonised per-study h5ad → `processed/<study>.h5ad` |
| 02 | `02_aggregate.py` | `atlas_00_aggregate.py` | `raw_datasets/{vazquez_garcia_2022,luo_2024,zheng_2023,ovarian_cancer_aggregate}.h5ad` | `…/integration/ovca_atlas_raw.h5ad` (common genes, concat, <500-cell-sample filter) |
| 03 | `03_preprocess_hvg.py` | `atlas_01_preprocess.py` | `ovca_atlas_raw.h5ad` | `…/integration/ovca_atlas_preprocess.h5ad` (HVG 4000, seurat_v3, batch=sample_id, subset) |
| 04 | `04_scvi.py` | `atlas_03_scvi.py` | `ovca_atlas_preprocess.h5ad` | `…/integration/ovca_atlas_integrated.h5ad` (`obsm["X_scVI"]`), `scvi_model_hgsc/`. **GPU; optional** |
| 04b | `04b_harmony_comparison.py` | `atlas_02_harmony.py` | `ovca_atlas_preprocess.h5ad` (needs `obsm["X_pca"]`) | `…/integration/ovca_atlas_harmony.h5ad` (`obsm["Harmony"]`). **COMPARISON / optional** |
| 05 | `05_cellassign.py` | `atlas_04a_cellassign.py` | `ovca_atlas_raw.h5ad`, `shared/cellassign_markers.csv` | `ovca_atlas_raw.h5ad` (in place: + `obs["celltype_pred"]`). **GPU; optional** |
| 06 | `06_scanvi.py` | `atlas_04b_scanvi.py` | `ovca_atlas_integrated.h5ad`, `ovca_atlas_raw.h5ad` (labels), `scvi_model_hgsc/` | `ovca_atlas_integrated.h5ad` (+ `obsm["X_scANVI"]`), `X_scANVI_hvg.npz`. **GPU; optional** |
| 07 | `07_process.py` | `atlas_05_process.py` | `ovca_atlas_raw.h5ad`, `X_scANVI_hvg.npz` | `…/integration/ovca_atlas_final.h5ad` (normalize/log1p, neighbours/UMAP/Leiden) |
| 02* | `02_concat_qc_doublets.py` | (repo; pre-existing) | 13 per-study h5ad | `processed/atlas_concatenated_filtered.h5ad` (per-cell QC + Scrublet<0.3). **See QC FLAG below** |
| 08 | `08_refilter_umap.py` | (repo; was `03b`) | `integrated_scanvi.h5ad` (20260213 run) | `obj("atlas_scanvi")` = `hgsc_atlas_scanvi.h5ad` (post-integration Scrublet<0.25, UMAP) |
| 09 | `09_atlas_qc_review.py` | (repo; was `04`) | `obj("atlas_scanvi")` | per-study QC summary CSVs + QC histograms/UMAPs |
| 10 | `10_umap_suite.py` | (repo; was `05`) | `obj("atlas_scanvi")` | 9 publication metadata UMAP SVGs |

The canonical integration chain is **02 → 03 → 04 → 05 → 06 → 07**. `04b` (Harmony) is a
comparison method, off the canonical path. Integration-chain intermediates
(`ovca_atlas_raw/preprocess/integrated/final.h5ad`, `scvi_model_hgsc/`, `X_scANVI_hvg.npz`) are
**regenerated** and written under `output_root/01_preprocess_qc/integration/` — they are not
deposited primary inputs. The QC/post-integration steps (`02_concat_qc_doublets`, `08`–`10`)
write under `output_root/01_preprocess_qc/`.

## These are the original scripts (optional GPU re-run)

Steps 02–07 are the **authoritative original cluster scripts** (`atlas_00_aggregate.py` …
`atlas_05_process.py`), migrated verbatim with only the hardcoded cluster paths
(`/project/6090753/dcook/...`, `/global/scratch/hpc3837/...`) replaced by central config and
seeds sourced from `config.SEED`. **Analytical parameters were not altered.** Re-running the
scVI/scANVI/CellAssign steps needs a GPU and is computationally expensive; it is **not
required** to reproduce downstream results because the integrated atlas object is **deposited as
a config entry object** (`obj("atlas_scanvi")`). The full code is included so every step from
raw data is open to scrutiny and can be independently re-executed.

## Superseded / excluded originals

- `atlas_04_annotation.py` — **SUPERSEDED.** It performed CellAssign + scANVI in one file using
  `adata.uns["all_counts"]` and `adata.raw`. The cleaner two-file split
  (`atlas_04a_cellassign` → `atlas_04b_scanvi`, i.e. steps 05 → 06) replaces it: same
  algorithm, but memory-safe (size factors on `.X`, labels round-tripped via the raw object and
  an NPZ for a CPU merge). Not migrated.
- `atlas_integration.py` / `atlas_integration2.py` — **SUPERSEDED.** Older monolithic all-in-one
  versions (aggregate → preprocess → Harmony → scVI → CellAssign → scANVI in one script). The
  00–05 split supersedes them. They use `n_top_genes=3000`, `target_sum=1e4`+`log1p` before HVG,
  inline `sc.tl.pca`, and `batch_size=64` — superseded by the split's parameters. The only thing
  the split's `04b_harmony_comparison` lacks is the inline PCA the monolith ran before Harmony
  (see Harmony FLAG below). Not migrated.

## Marker-matrix reconciliation

The authoritative `05_cellassign.py` (atlas_04a) references `cellassign_markers.csv`. That file
(now `shared/cellassign_markers.csv`, copied from the original
`/Users/dpcook/Projects/ovca_states/data/cellassign_markers.csv`) is the **53-gene × 11-type**
matrix: `Epithelial, Mesothelial, Fibroblasts, Endothelial, T/NK, B cells, Macrophage, DC,
Plasma cells, Mast, Other`.

This is **NOT** the same as `shared/cellassign_markers_v3.csv` (**81-gene × 16-type**:
adds Smooth_Muscle, Pericyte, splits T_cell/NK_cell, Neutrophil, Erythrocyte, etc.), which is
the matrix used by the later **20260213** production integration run whose output
(`integrated_scanvi.h5ad`) feeds step 08 and yields the deposited `hgsc_atlas_scanvi.h5ad`. The
two marker matrices are different vintages and are **not interchangeable**. They are both kept
in `shared/`; each script points at the one its own pipeline version used. **FLAG (see below).**

> The 81-gene/16-type expectation belongs to the `cellassign_markers_v3.csv` lineage
> (16 columns incl. Other = 16 types; 81 marker rows). The migrated `atlas_04a` originals
> instead use the 53-gene/11-type `cellassign_markers.csv`.

## RECONCILIATION FLAGS (need user input — not silently resolved)

1. **QC / Scrublet placement.** The authoritative `atlas_00..05` chain (steps 02–07) has **no
   per-cell UMI/gene QC and no Scrublet** — only the `<500-cell-sample` filter in `02_aggregate`.
   But the Methods describe pre-integration Scrublet < 0.3 (2,398,571 → 2,326,532 cells), and the
   repo's `02_concat_qc_doublets.py` implements exactly that (total_counts≥500, n_genes≥300,
   Scrublet<0.3). **Where does QC/Scrublet actually belong in the authoritative chain?** Options:
   (a) it lives in the `02_concat_qc_doublets` lineage and the `atlas_00` aggregate is an
   alternate/earlier entry that skipped it; (b) it is a missing step that should sit between
   aggregate and preprocess; (c) per-study QC happened in `01_dataprep`. These two lineages
   (`atlas_00` aggregate of pre-merged study h5ads vs. `02_concat_qc_doublets` concat of 13
   per-study h5ads) do not obviously share the same QC. **Not guessed — please resolve.**

2. **Object-naming / pipeline-version bridge.** The authoritative chain ends at
   `ovca_atlas_final.h5ad` (step 07, via `obsm["X_scANVI"]`, normalize/log1p, UMAP, Leiden
   res=0.15). The repo's `08_refilter_umap.py` instead reads
   `integrated_scanvi.h5ad` (from the **20260213** integration run) and writes
   `hgsc_atlas_scanvi.h5ad` (post-integration Scrublet<0.25 + UMAP min_dist=0.2). **How does
   `ovca_atlas_final.h5ad` map to `integrated_scanvi.h5ad` / `hgsc_atlas_scanvi.h5ad`?**
   These appear to be **two integration vintages**: the migrated `atlas_00..05` scripts (4-study
   aggregate, 4000 HVG, `X_scANVI`, 53-gene markers) vs. the 20260213 run (13-study concat,
   `X_scanvi` lowercase, 81-gene v3 markers) that produced the deposited object. Is step 08 a
   later add-on stage on top of the 20260213 run (not on top of `atlas_05`)? Is
   `ovca_atlas_final` an earlier object superseded by the deposited `hgsc_atlas_scanvi.h5ad`?
   **Not guessed — please resolve.** (`obsm` key casing differs too: `X_scANVI` here vs.
   `X_scanvi` in step 08's input.)

3. **CellAssign marker matrix** (see section above): authoritative 05_cellassign uses the
   53-gene/11-type `cellassign_markers.csv`; the deposited atlas used the 81-gene/16-type
   `cellassign_markers_v3.csv`. Both kept in `shared/`; confirm which lineage is canonical for
   the manuscript.

4. **Harmony PCA dependency.** `04b_harmony_comparison.py` reads `obsm["X_pca"]`, but
   `03_preprocess_hvg` (atlas_01) does not compute PCA — the monolith ran `sc.tl.pca` inline
   before Harmony. To run the Harmony comparison, a PCA pass must be added to the preprocess
   object. Logic preserved as-is; flagged here.

## Conventions applied

Central config for all paths (`raw_datasets` config key for per-study raw inputs;
`output_root/01_preprocess_qc/integration/` for regenerated intermediates); header docstrings
(purpose / INPUTS / OUTPUTS / runtime tier / manuscript role); seeds from `config.SEED` where
the originals set them (Harmony random_state, Scrublet, UMAP). Analytical parameters of the
authoritative scripts preserved EXACTLY. The original `atlas_03_scvi.py` output-filename typo
(`ovca_atlas_integatrated.h5ad`) was corrected to `ovca_atlas_integrated.h5ad` so the chain is
coherent.

## Figures supported

QC/overview substrate only. SF1 (per-study QC violins), SF2 (integration UMAPs), Fig 1A (atlas
UMAP) and SF3 (metadata UMAPs) are rendered by per-panel scripts in `figures/` from
`obj("atlas_scanvi")` / `obj("atlas_final")`.
