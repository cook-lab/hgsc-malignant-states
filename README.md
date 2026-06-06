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

## Reproducing

1. **Environments**
   - Python: `conda env create -f environment.yml && conda activate epitype-py`
   - R (4.5.x): `Rscript -e 'renv::restore(lockfile="renv.lock")'` (pinned lockfile, 410 packages). Alternative install-from-scratch: `Rscript renv_bootstrap.R`. Key versions: `docs/r_key_package_versions.csv`.
2. **Data**: obtain the deposited objects (see `data/README.md`) and set `DATA_ROOT`
   (or edit `config/config.yml`).
3. **Run**: backend stages in numeric order within `atlas/` and `spatial/`, then the
   `figures/` scripts. Each script's header lists its inputs, outputs, and the panel(s) it supports.

## Trust boundary

The scVI/scANVI **integration was run on a compute cluster and is not re-executed here**; its
output (`hgsc_atlas_scanvi.h5ad`) is the entry point for all downstream analysis. Everything
downstream of integration runs locally and is reproducible from the deposited objects.

## Reproducibility status

This code was independently re-executed and audited (see `docs/REPRODUCIBILITY.md`): the
manuscript is strongly reproducible — the SecA/SecB axis, survival associations, spatial
gradients, and morphometry all regenerate from the deposited inputs. A small number of known
discrepancies (and the seed/cohort hardening applied here) are documented there.

## Citation / contact

David Cook — dacook@ohri.ca · Ottawa Hospital Research Institute.
