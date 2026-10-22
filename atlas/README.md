# atlas/: the scRNA-seq atlas backend

This backend builds the integrated HGSC scRNA-seq atlas: 13 public studies, 1,980,703 cells, 362 samples and 158
patients. It also produces the atlas results that the figures read:
- cell-type annotation
- epithelial NMF programs, which define the SecA / Intermediate / SecB states
- functional scores
- CopyKAT clones
- TCGA-OV signature scores
- the SingleR reference used to annotate Xenium

You do not need to run this backend to reproduce a figure. Its outputs are deposited (see
[`../data/README.md`](../data/README.md)), and every atlas figure starts from them. Paths below without a root are
relative to `$DATA_ROOT/2026_final_atlas/`, the deposited atlas project folder. In the lab this is the project folder
`202601_final_atlas`; the scripts' own paths use both names.

## Stages, in run order

There is no stage 06.

| Stage | What it does | Main deposited outputs | Read by |
|---|---|---|---|
| [`01_preprocess_qc/`](01_preprocess_qc/README.md) | Prepares each study, concatenates them on 14,470 shared genes, runs cell QC and Scrublet, selects HVGs, runs CellAssign and scVI → scANVI, then a second doublet pass and the atlas UMAP | `pre-processing and integration/processed/*.h5ad`; `pre-processing and integration/20260213_integration/output/anndata/integrated_scanvi.h5ad`; `hgsc_atlas_scanvi.h5ad` | stage 02; Figure S1; Figure S4a (the CellAssign labels carried into the final atlas) |
| [`02_annotation/`](02_annotation/README.md) | Annotates level 1 (13 lineages) and level 2 (51 types), applies the final filter, fixes metadata, splits the lineage subsets and computes the per-lineage UMAPs | `hgsc_atlas_final.h5ad`; `celltype_h5ad/hgsc_atlas_final_<lineage>.h5ad`; `output/08_celltype_level2/markers/` | every atlas figure; Tables S1–S3 |
| [`03_epithelial_nmf/`](03_epithelial_nmf/README.md) | Fits the k = 10 NMF of the epithelial compartment; classifies genes; runs per-patient consensus NMF; exports the atlas secretory cells for UCell; writes the epithelial UMAP table | `output/11d_epithelial_nmf/`; `output/11e_nmf_characterization/`; `output/11f_consensus_nmf/`; `output/11g_secAB_consensus_nmf/`; `output/18_ucell_atlas/`; `output/fig_secretory_polarization/data/meta.parquet` | Figs 1–2; Fig 3 (atlas reference); Fig 6h; Figures S5–S10; Tables S5–S6; Xenium polarization thresholds; Fig 7g–j |
| [`04_functional/`](04_functional/README.md) | Scores PROGENy, Hallmark and DoRothEA per epithelial cell | `output/21_epitype_functional_characterization/` | Fig 2a; Figure S8 |
| [`05_cnv/`](05_cnv/README.md) | Runs CopyKAT per sample and calls subclones | `output/19_cnv/per_sample/` | Fig 1g, 1l, 1m; Figures S4c and S7; Discussion (93.5%) |
| [`07_deconvolution_survival/`](07_deconvolution_survival/README.md) | Estimates the BayesPrism epithelial fraction of TCGA-OV tumours and scores them with the 177 + 177 NMF-specific genes | `output/22_tcga_deconvolution/22d_signature_scores.csv` | Fig 7g–j |
| [`08_xenium_reference/`](08_xenium_reference/README.md) | Builds the 16-type SingleR reference for Xenium | `output/xenium_celltype/xenium_celltype_downsampled.h5ad` | `spatial/03_annotation_polarization/`; Figure S14; Table S7 |

## Entry points

The figure and table scripts, and the model pipelines in `models/`, read these deposited files directly:

| File | Content used | Read by |
|---|---|---|
| `hgsc_atlas_final.h5ad` (config `atlas_final`) | obs, `obsm/X_umap`, `X`, `layers/counts` | every atlas figure, through `figures/_prep/atlas_final_tables.py`; Fig 2f (counts; `figures/_prep/fig2f_consensusov.py`); Figure S4b (`X`; `figures/supplementary/08_figS4b_marker_stats.py`); Tables S1–S3 (obs) |
| `celltype_h5ad/hgsc_atlas_final_epithelial.h5ad` | epithelial counts and log-normalised `X` | Fig 1g prep; Fig 2a (scFEA); Fig 2b; Fig 2g; Fig 6h (`figures/figure6/10_fig6h_atlas_epithelial_pct.py`); Figure S6 (when config `atlas_epithelial` is absent); Table S6 |
| `celltype_h5ad/hgsc_atlas_final_macrophage.h5ad` | macrophage counts and level-2 labels | Fig 6f,h (`figures/figure6/08_fig6fh_atlas_macrophage_pseudobulk.py`, which reads a copy at `$OUTPUT_ROOT/atlas_local/hgsc_atlas_final_macrophage.h5ad` unless a path is given) |
| `output/11d_epithelial_nmf/11d_nmf_usage.csv`, `11d_nmf_loadings.csv`, `11d_factor_scores.csv` | per-cell usage of the 10 factors (Factor 2 defines the states); gene loadings; factor scores | every epithelial panel (usage); Tables S5–S6; Figure S10 and the KPCA.B alignment (loadings; `models/kpca_mouse/`) |
| `output/11e_nmf_characterization/11e_gene_classification.csv` | SecA-specific / SecB-specific / shared genes | Tables S5–S6; `models/kpca_mouse/scripts/02_orthologs.py` |
| `output/11f_consensus_nmf/`, `output/11g_secAB_consensus_nmf/` | consensus metaprograms | Table S5 |
| `output/fig_secretory_polarization/data/meta.parquet` | `UMAP1` / `UMAP2` (the epithelial UMAP) and `schema_nmf` (the epithelial states) | every epithelial UMAP, through `figures/_prep/atlas_final_tables.py` (`UMAP1`/`UMAP2` only); Table S6 and Fig 6h (`schema_nmf`; `tables/06_tableS6_secA_secB_signatures.py`, `figures/figure6/10_fig6h_atlas_epithelial_pct.py`) |
| `output/19_cnv/per_sample/` | CopyKAT calls and subclones | Fig 1g, 1l, 1m; Figures S4c and S7 |
| `output/21_epitype_functional_characterization/` | per-cell PROGENy, Hallmark and DoRothEA scores | Fig 2a; Figure S8 |
| `output/18_ucell_atlas/` | the atlas secretory-cell export (counts, barcodes, genes, states); atlas UCell scores; gene space shared with the organoids | Xenium polarization thresholds (`spatial/03_annotation_polarization/03_ucell_scoring_noBCAM.R`); Fig 3 and the organoid scripts; Table S6 |
| `output/22_tcga_deconvolution/22d_signature_scores.csv` | TCGA-OV scores and covariates | Fig 7g–j |
| `data/cibersort_data_prev/tcga_ecotyper.txt` | TCGA-OV expression (also the input of stage 07) | Table S6 (gene presence) |
| `output/xenium_celltype/xenium_celltype_downsampled.h5ad` (config `xenium_ref`) | SingleR reference, 1,000 cells per type | Xenium SingleR; Figure S14; Table S7 |
| `pre-processing and integration/`: `processed/`, `atlas_metadata.xlsx`, `raw/` | pre- and post-QC concatenated objects and `nath_2021.h5ad`; the curated sample sheet; the raw study files (Nath 2021 GEO cell annotations; file names that carry the GEO sample accessions) | Figure S1; `figures/_prep/nath2021_fix.py`; Table S1 |
| `output/06_celltype_level1/markers/`, `output/08_celltype_level2/markers/`, `output/supplemental/` | level-1 and level-2 Wilcoxon markers; the per-cell metadata export and the ciliated-vs-secretory table | Tables S1–S4 |
| `tools/scFEA/data/` | scFEA human M168 model files | Fig 2a; Figure S8 |
| `output/spectrum_epithelial_crosstab/` | published per-cell labels of Vázquez-García (GSE180661) and Geistlinger | Figure S9f, S9g |

## How the figures read the atlas

Every atlas panel reads cell metadata through `figures/_prep/atlas_final_tables.py`. This covers Figs 1–2 and Figures
S2–S9. The loader applies `figures/_prep/nath2021_fix.py`, which reassigns each Nath 2021 cell to its GEO patient and
time point while the caches are built. The deposited object is not changed. With the fix applied, the atlas has 158
patients and 362 samples.

The loader also defines the epithelial states. Among non-ciliated epithelial cells:
- SecA: Factor 2 usage below the 50th percentile (284,344 cells)
- Intermediate: 50th to 75th percentile (142,172)
- SecB: at or above the 75th percentile (142,172)

The ciliated cells (6,678) form a fourth group. See [`../figures/_prep/README.md`](../figures/_prep/README.md).

The final run rendered every atlas panel with `OUTPUT_ROOT` set to a separate atlas output root, not the one used for
the Xenium figures.

## What the stage folders hold

- **`as_run/`**: trimmed copies of the scripts that ran in the atlas project folder `202601_final_atlas` and wrote the
  deposited outputs. The code that wrote the outputs the paper uses is kept as it ran; plots and analyses that are not
  in the paper were removed. Each trimmed copy opens with a note on what was removed; copies with nothing to remove are
  unchanged. The untrimmed scripts stay in `202601_final_atlas`. The copies keep that project's paths: their
  `…/2026_final_atlas/…` and `…/202601_final_atlas/…` paths correspond to `$DATA_ROOT/2026_final_atlas/` here.
- **`01_preprocess_qc/03`–`07`**: the integration scripts that ran on the GPU cluster, with their paths resolved
  through `config/` and their QC plots, metadata JSON and UMAP previews removed. They write to
  `$OUTPUT_ROOT/01_preprocess_qc/integration/`.
- **`03_epithelial_nmf/01_epithelial_nmf.py`**: a re-implementation of
  `03_epithelial_nmf/as_run/11d_epithelial_nmf.py` that runs on the deposited epithelial object and writes to
  `$OUTPUT_ROOT/03_epithelial_nmf/`. `01b_verify_against_deposited.py` compares it with the deposited NMF: the feature
  set is the same, and the loadings and the training-cell usage match to floating-point precision.

The figures read the deposited outputs, not the `$OUTPUT_ROOT` outputs of these scripts.

## Environments

See [`../env/README.md`](../env/README.md).

| Environment | Ran |
|---|---|
| E3 (Python 3.12.3) | data prep, QC and Scrublet |
| E4 (GPU cluster, scvi-tools 1.4.1) | HVGs, CellAssign, scVI and scANVI |
| E2 (Python 3.12.9) | every later Python step |
| E1 (R 4.4.3 / Bioconductor 3.20) | CopyKAT, BayesPrism and the atlas UCell scores |

## Missing code

- **[code missing]** The export of `output/supplemental/T3_epithelial_markeres.csv`, the ciliated-vs-secretory Wilcoxon
  table behind Table S4 (see `02_annotation/README.md` for the part that is traced). The level-1 markers of Table S2
  are written by `02_annotation/as_run/xenium_coverage_report_level1.py`.
- **[code missing]** The original per-cell export `output/supplemental/T1_atlas_metadata.csv` behind Table S1. Its later
  metastatic-site correction is `02_annotation/as_run/fix_metastatic_site_docs.py`.
- `output/spectrum_epithelial_crosstab/geistlinger_cell_annotations.csv` (Figure S9g) is written by
  `figures/supplementary/as_run/geistlinger_cell_annotations.R` from the authors' published SingleCellExperiment objects.
- **[code missing]** Download and preparation of the TCGA-OV inputs `data/cibersort_data_prev/{tcga_ecotyper.txt,
  tcga_hla_clinical.csv}` (429 tumours; the clinical table also holds mutation and HLA-expression columns). They are
  byte-identical copies of files the laboratory has had since August 2024 (an earlier CIBERSORTx project); how those
  were made is not recorded.
- **Downloads (no script).** The 13 public datasets, downloaded from their accessions (Table S1) into
  `pre-processing and integration/raw/<study>/`, and the other public files the figures read
  (`output/spectrum_epithelial_crosstab/GSE180661_GEO_cells.tsv.gz` and the scFEA model files in `tools/scFEA/data/`).

Some intermediate objects are not deposited:
- the level-1 and level-2 objects (`hgsc_atlas_celltype_level1.h5ad`, `hgsc_atlas_celltype_level2.h5ad`) and
  `hgsc_atlas_xenium.h5ad`
- lineage subsets carrying the per-lineage UMAP `X_umap_local`; the deposited `celltype_h5ad/` subsets do not carry it

The steps that read these objects cannot be re-run on the deposit as written. Their results are carried by the
deposited `hgsc_atlas_final.h5ad`, the lineage subsets, `meta.parquet` and the SingleR reference.
