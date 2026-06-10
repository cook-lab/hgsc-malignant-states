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

The pipeline has two tiers:

- **Backend** (`atlas/`, `spatial/`) — computationally heavy preprocessing and analysis (raw data
  → integration → annotation → per-analysis results). It produces **standardized preprocessed
  objects** (the integrated atlas, SFE objects, per-stage result tables and model fits).
- **Downstream** (`figures/`, `tables/`) — lightweight scripts that load those objects and render
  the manuscript panels and supplemental tables.

```
  atlas/ + spatial/  ──►  standardized objects  ──►  figures/ + tables/  ──►  output/
   (heavy backend)         (atlas h5ad, SFEs,         (load object →           (panels,
                            result tables/fits)        plot/export)             supp tables)
```

Those backend objects are slow to regenerate (integration needs a GPU; some steps are many
CPU-hours), so they are also **provided directly** (see `data/README.md`): set `DATA_ROOT` to the
data bundle and you can reproduce any figure or table without re-running the backend. All rendered
outputs go to `OUTPUT_ROOT` (`./output`).

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
#   01_dataprep → 02_concat_qc_doublets (QC + Scrublet<0.3) → 03_preprocess_hvg →
#   04_cellassign → {05a_scvi, 05b_scanvi, 05c_mrvi, 05d_sysvi, 05e_harmony} →
#   06_benchmark → 07_finalize (--method scanvi) → 08_refilter_umap → 09 → 10.
# Steps 03–07 are the OFFICIAL cluster scripts (5 integration methods benchmarked; scANVI
# SELECTED). They need a GPU and are OPTIONAL: the integrated atlas object is deposited as an
# entry object, so you can start downstream from 08. The integration code is included so every
# step is open to scrutiny and can be re-executed.
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
multi-study integration**. These are the **official cluster scripts**
(`01_preprocess` … `05_finalize`), migrated as `atlas/01_preprocess_qc/03_preprocess_hvg.py`
through `07_finalize.py` with only the hardcoded cluster paths replaced by central config — not
reconstructions. CellAssign labels (step 04) seed semi-supervised integration; **five methods
are run and benchmarked** (scVI, scANVI, MrVI, SysVI, Harmony; step 06 scib-metrics), with
**scANVI SELECTED** for the manuscript (step 07 `--method scanvi`). Re-running the integration is
computationally expensive (originally run on a GPU cluster) and is **not required** to
reproduce downstream results — the integrated atlas object (`hgsc_atlas_scanvi.h5ad`) is
provided as a deposited entry object — but the integration code is included so that every
step is open to scrutiny and can be independently re-executed (now seeded from `config.seed`;
GPU model training is best-effort-deterministic, **not** bit-identical across hardware, so the
deposited `hgsc_atlas_scanvi.h5ad` is the canonical reference — see docs/REPRODUCIBILITY.md).

## Integration is included (optional to re-run)

The repository includes the full pipeline from raw data through final figures, **including
the multi-study integration**, as the official cluster scripts
(`atlas/01_preprocess_qc/03_preprocess_hvg.py` … `07_finalize.py`, path-centralised). QC and
per-sample Scrublet (< 0.3) are applied in `02_concat_qc_doublets.py`, whose output
`atlas_concatenated_filtered.h5ad` is the raw input to the integration chain; five integration
methods are benchmarked (step 06) and **scANVI is selected** (step 07 → `integrated_scanvi.h5ad`,
obsm `X_scanvi`), feeding `08_refilter_umap.py` which writes the deposited `hgsc_atlas_scanvi.h5ad`.
Re-running it is computationally expensive (originally run on a GPU cluster) and is **not
required** to reproduce downstream results — the integrated object is deposited as an entry
object — but the code is included so every step can be independently scrutinized and re-executed.

## Reproducibility status

This code was independently re-executed and audited (see `docs/REPRODUCIBILITY.md`): the
manuscript is strongly reproducible — the SecA/SecB axis, survival associations, spatial
gradients, and morphometry all regenerate from the deposited inputs. A small number of known
discrepancies (and the seed/cohort hardening applied here) are documented there.

## Citation / contact

David Cook — dacook@ohri.ca · Ottawa Hospital Research Institute.

## License

Released under the **MIT License** — free to use, modify, and redistribute (including
commercially), with no restrictions beyond retaining the copyright/license notice. See
[`LICENSE`](LICENSE).
