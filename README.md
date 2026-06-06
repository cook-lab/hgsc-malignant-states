# HGSC malignant states — tumour architecture shapes polarized epithelial states

Code to reproduce the analyses and figures in:

> **Tumour architecture shapes polarized epithelial states that predict survival in high-grade serous ovarian cancer.**
> Nersesian\*, Abou-Hamad\*, Durocher, Akiki, Domecq, Southworth, Deng, Meunier, de Ladurantaye, Mes-Masson, Tessier Cloutier, Cook. (2026, preprint)

We integrate scRNA-seq from 13 studies (~2M cells) into an HGSC atlas, resolve a dominant axis
of secretory epithelial polarization — **SecA** (proliferative, progenitor-like) ↔ **SecB**
(quiescent, mucosal-injury/keratin program) — and show with custom Xenium spatial transcriptomics
(8 whole tissues + a 97-patient TMA), patient-derived organoids, and TCGA-OV that this axis is
spatially patterned by a hypoxic gradient, assembles a glycolytic-macrophage / lymphocyte-excluded
niche, and independently predicts overall and progression-free survival.

## Repository layout

```
config/      Central paths, seeds, cohort, polarization schema (config.yml + .py/.R loaders)
shared/      Canonical SecA/SecB gene signatures (noBCAM 7-gene), palettes, helpers
atlas/       scRNA-seq atlas backend pipeline (Python / scanpy)   [01_preprocess_qc … 08_xenium_reference]
spatial/     Xenium/Visium backend pipeline (R / SpatialFeatureExperiment)  [00_setup … 09_external_validation]
figures/     ONE reproduction script per manuscript panel (load intermediate → plot → export)
tables/      Supplemental-data table generators + cleaned outputs
data/        Pointers to deposited objects (not stored here; ~1.4 TB)
docs/        Reproducibility notes, expanded methods, data availability
```

The backend pipelines produce the intermediate objects/caches; each `figures/figureN/` script
loads a cache and renders a publication panel. Panel → script mapping is in `figures/README.md`.

## Reproducing the paper

### Data flow

```
  DEPOSITED DATA (set DATA_ROOT)                         OUTPUT_ROOT (./output, gitignored)
  ├─ entry objects: atlas h5ad, SFE dirs        ──┐
  └─ analysis output caches (per-stage tables,    │   figure scripts ──► output/figures/<figureN>/*.svg|png|pdf
     GAM .rds, NMF/CNV/LIANA/survival outputs) ──┼──►  table scripts  ──► output/tables/*
                                                  │   backend scripts ──► output/<stage>/*  (regeneration / verification)
```

The repo is **deposit-driven**: every figure and table script reads its inputs (entry
objects *and* intermediate caches) from the **deposited data bundle** under `DATA_ROOT`, and
writes only rendered outputs to `OUTPUT_ROOT`. You do **not** need to re-run the heavy backend
to reproduce a figure. The `atlas/` and `spatial/` backend stages are provided to **regenerate
and verify** those caches from the entry objects.

### 1. Set up environments
```bash
conda env create -f environment.yml && conda activate epitype-py     # Python
Rscript -e 'renv::restore(lockfile="renv.lock")'                     # R 4.5.x (410 pinned pkgs)
# (alt R install-from-scratch: Rscript renv_bootstrap.R; versions in docs/r_key_package_versions.csv)
```

### 2. Get the data and point the config at it
Download the deposited bundle (see `data/README.md`) and set the path (no code edits needed):
```bash
export DATA_ROOT=/path/to/deposited/bundle      # contains entry objects + analysis output caches
export OUTPUT_ROOT=./output                      # where renders/regenerated caches go (default)
```

### 3a. Reproduce the FIGURES (primary path — no backend re-run)
```bash
# a few supplementary panels (SF1/SF2) first regenerate small obs caches:
PYTHONPATH=. python figures/_prep/00_extract_atlas_obs.py
PYTHONPATH=. python figures/_prep/00b_extract_integration_umaps.py
# then any panel — each script reproduces its manuscript panel(s):
PYTHONPATH=. python figures/figure2/03_atlas_volcano_secA_secB.py   # Fig 2C  (Python)
Rscript           figures/figure7/01_xenium_forest_cox.R            # Fig 7A/B (R)
```
**To reproduce one specific panel**, look it up in `figures/README.md` (the full Fig 1A…7G /
SF1…SF14 → script map) and run that one script. Each script's header docstring lists its exact
inputs, outputs, and the panel(s) it supports.

### 3b. Reproduce the SUPPLEMENTAL DATA TABLES
```bash
# generators in tables/ (Supp Data 4–7 committed; 1–3 export from the entry objects)
PYTHONPATH=. python tables/<NN>_*.py   ;   Rscript tables/<NN>_*.R
```

### 4. (Optional) Regenerate the analysis from the entry objects
To re-derive the intermediate caches from scratch (e.g. to verify reproducibility), run the
backend stages **in numeric order**, which write to `OUTPUT_ROOT`:
```bash
# Atlas (Python) — the full chain from raw data is included. atlas/01 runs:
#   01_dataprep → 02_aggregate → 03_preprocess_hvg → 04_scvi → 05_cellassign → 06_scanvi
#   → 07_process, with 02_concat_qc_doublets (QC+Scrublet) and 04b_harmony_comparison
#   alongside, then 08_refilter_umap → 09 → 10.
# Steps 04/05/06 (scVI → CellAssign → scANVI) are the ORIGINAL cluster scripts; they need a
# GPU and are OPTIONAL: the integrated atlas object is deposited as an entry object, so you
# can start downstream from 08. The integration code is included so every step is open to
# scrutiny and can be re-executed.
for s in atlas/01_preprocess_qc atlas/02_annotation atlas/03_epithelial_nmf atlas/04_functional \
         atlas/05_cnv atlas/06_cellcomm atlas/07_deconvolution_survival atlas/08_xenium_reference; do
  # run the numbered scripts within each stage in order (see each stage README)
  done
# Spatial (R) — 00_setup is sourced by all; then build → QC → annotation → downstream:
#   spatial/00_setup → 01_build_sfe → 02_qc → 03_annotation_polarization → 04_neighborhood
#   → 05_gradients_gams → 06_spatial_stats → 07_morphometry → 08_clinical_survival → 09_external_validation
```
Then diff `OUTPUT_ROOT/<stage>` against the deposited caches. Each stage dir has a `README.md`
with its ordered scripts and input→output map. **Runtime tiers** are in each script header
(fast / moderate / heavy); the heavy steps (CopyKAT, LIANA, BayesPrism, scFEA, per-sample GAMs,
Lee's L) are flagged. Stochastic steps are seeded from `config.seed` for determinism.

### Integration step (optional GPU re-run)
The repository includes **every step from raw data through final figures, including
multi-study integration** (scVI → CellAssign → scANVI). These are the **original cluster
scripts** (`atlas_00_aggregate` … `atlas_05_process`), migrated as
`atlas/01_preprocess_qc/02_aggregate.py` through `07_process.py` with only the hardcoded
cluster paths replaced by central config — not reconstructions. Re-running the integration is
computationally expensive (originally run on a GPU cluster) and is **not required** to
reproduce downstream results — the integrated atlas object (`hgsc_atlas_scanvi.h5ad`) is
provided as a deposited entry object — but the integration code is included so that every
step is open to scrutiny and can be independently re-executed.

## Integration is included (optional to re-run)

The repository includes the full pipeline from raw data through final figures, **including
the multi-study integration** (scVI → CellAssign → scANVI), as the original cluster scripts
(`atlas/01_preprocess_qc/02_aggregate.py` … `07_process.py`, path-centralised). Re-running it
is computationally expensive (originally run on a GPU cluster) and is **not required** to
reproduce downstream results — the integrated object (`hgsc_atlas_scanvi.h5ad`) is deposited as
an entry object — but the code is included so every step can be independently scrutinized and
re-executed. (See the `atlas/01_preprocess_qc/README.md` for two open reconciliation flags: the
QC/Scrublet placement and the object-naming bridge between the migrated `atlas_05` output and the
deposited integration object.)

## Reproducibility status

This code was independently re-executed and audited (see `docs/REPRODUCIBILITY.md`): the
manuscript is strongly reproducible — the SecA/SecB axis, survival associations, spatial
gradients, and morphometry all regenerate from the deposited inputs. A small number of known
discrepancies (and the seed/cohort hardening applied here) are documented there.

## Citation / contact

David Cook — dacook@ohri.ca · Ottawa Hospital Research Institute.
