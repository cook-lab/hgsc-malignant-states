# HGSC malignant states

Code to reproduce the analyses, figures and supplementary tables in:

> **Tumor architecture shapes polarized malignant states that predict outcome and treatment response in high-grade
> serous ovarian cancer.**
> Nersesian\*, Abou-Hamad\*, Durocher, Akiki, Domecq, Southworth, Verey, Deng, Ramos-Valdes, Meunier, de Ladurantaye,
> Shepherd, Provencher, Mes-Masson, Tessier Cloutier, Cook\*\*.
>
> \*Equal contribution. \*\*Corresponding author.

We integrate scRNA-seq from 13 studies (1.98 million cells, 362 samples) into an atlas of high-grade serous ovarian
carcinoma (HGSC). The atlas resolves a dominant axis of malignant secretory polarization, from proliferative,
progenitor-like **SecA** cells to quiescent, injury-like **SecB** cells. We follow this axis with Xenium spatial
transcriptomics (8 whole tumors and tissue microarrays from 378 patients), Visium HD, patient-derived organoids, the
KPCA.B syngeneic mouse model, and bulk transcriptomes from TCGA-OV and the ICON7 trial.

## Repository layout

```
config/      Paths, seed, cohort and polarization thresholds (config.yml + Python / R loaders, study-code lookup)
shared/      SecA / SecB signature genes (signatures.yml), CellAssign marker matrix
atlas/       scRNA-seq atlas backend: QC, integration, annotation, NMF programs, functional scores, CopyKAT,
             TCGA deconvolution, Xenium SingleR reference
spatial/     Xenium backend: SFE build, QC, final TMA cohort, SingleR labels, UCell polarization, pathway scores,
             raw-count export, TMA clinical table
models/      Patient-derived organoids (PDOs, OPTO98) and the KPCA.B mouse model
icon7/       ICON7 trial transcriptomes: cohort build and polarization scores
figures/     One folder per figure: figure1/ … figure7/, supplementary/ (Figures S1–S18), _prep/ (shared loaders)
tables/      Supplementary Tables S1–S12
data/        What the data deposit holds and where the code looks for it (no data are stored here)
env/         The software environment that ran each step, and the scripts that write the lock files
docs/        Versions of the key R packages
```

[`figures/README.md`](figures/README.md) maps every panel, table and reported number to the scripts that produce it.
Figure scripts are named `NN_fig<number><panels>_<what>` and numbered in run order: for example,
`figures/figure4/14_fig4i_local_density.py` computes Figure 4i and `15_fig4i_local_density_plot.py` draws it. Some
output file names do not match the panel letters; each figure README maps outputs to panels.

## Reproducing the paper

### Data availability

No data are stored in this repository. The processed data (atlas objects, Xenium, Visium HD, organoid and mouse
single-cell data) are being deposited, and the accession will be added here. The public inputs are the 13 atlas
studies (accessions in Table S1), TCGA-OV (GDC) and ICON7 (GSE140082). [`data/README.md`](data/README.md) lists what
the deposit holds and the folder layout the code expects.

### Data flow

The pipeline has two tiers:

- **Backend** (`atlas/`, `spatial/`, `models/`, `icon7/`): computationally heavy processing from raw data to
  processed objects and result tables. Integration needs a GPU, and several steps take many CPU-hours.
- **Downstream** (`figures/`, `tables/`): scripts that read those objects and draw each panel or write each table.
  Analyses that belong to a single figure, such as the Figure 6 niche statistics, sit in that figure's folder in run
  order.

```
  atlas/  spatial/   ──►  processed objects,  ──►  figures/  tables/   ──►  $OUTPUT_ROOT
  models/ icon7/          result tables            (load object →          (panels,
  (heavy backend)                                   plot / export)          supplementary tables)
```

The outputs of the atlas, model and ICON7 backends are deposited, and the figures read them from the deposit. The
Xenium figures read tables that the `spatial/` backend writes to `$OUTPUT_ROOT` (SingleR label tables, the final TMA
core set, the raw-count export); each figure README lists the upstream steps it needs.

### 1. Set up the environments

Python 3.12.9 and R 4.4.3 / Bioconductor 3.20 ran everything except the steps listed in
[`env/README.md`](env/README.md).

```bash
conda env create -f environment.yml && conda activate epitype-py   # Python
Rscript -e 'renv::restore(lockfile = "renv.lock")'                 # R, exact versions
# alternative R install without renv: Rscript renv_bootstrap.R; key versions in docs/r_key_package_versions.csv
```

### 2. Point the code at the data

All code outside `as_run/`, `models/` and `icon7/` finds its inputs under `DATA_ROOT` and writes to `OUTPUT_ROOT`
(`config/config.yml`), so no code edits are needed:

```bash
export DATA_ROOT=/path/to/deposit                  # see data/README.md
export OUTPUT_ROOT=/path/to/outputs                # absolute; every step writes here
export ORGANOIDS_ROOT=$DATA_ROOT/202603_organoids  # Figure 3a–e; does not follow DATA_ROOT
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8         # R needs a UTF-8 locale to read config.yml
```

Use one absolute `OUTPUT_ROOT` for every step: later steps read earlier outputs from it, and
`tables/build_all_tables.sh` runs from `tables/`. A few inputs are copied into `$OUTPUT_ROOT` by hand, for example
the probe-QC exclusion list (`spatial/02_qc/`) and the Figure 7a–d tables that Tables S8 and S11 read
(`figures/figure7/`). The folder READMEs list each one.

### 3a. Reproduce the figures

Look the panel up in [`figures/README.md`](figures/README.md), then run that figure folder's scripts in the order of
its README. Run them from the repository root, because R scripts `source()` helpers by repo-relative path:

```bash
python figures/figure4/14_fig4i_local_density.py
python figures/figure4/15_fig4i_local_density_plot.py
Rscript figures/figure7/02_fig7abcd_tma_survival_row.R
```

`run_panel.sh <script>` runs a single script with logging.

### 3b. Reproduce the supplementary tables

```bash
bash tables/build_all_tables.sh                    # Supplementary Tables S1–S12
```

### 4. (Optional) Regenerate the processed objects

To re-derive the processed objects, run the backend stages in the order given in [`atlas/README.md`](atlas/README.md),
[`spatial/README.md`](spatial/README.md), [`models/README.md`](models/README.md) and
[`icon7/README.md`](icon7/README.md). Each stage README lists its scripts in run order with their inputs and outputs.

## Integration is included (optional to re-run)

The repository includes the multi-study integration, as the scripts that ran on the GPU cluster
(`atlas/01_preprocess_qc/03_preprocess_hvg.py` … `07_finalize.py`, with paths resolved through `config/`). Cell QC
and per-sample Scrublet doublet removal run in `atlas/01_preprocess_qc/as_run/02_concat_qc_doublets.ipynb`, whose
output `atlas_concatenated_filtered.h5ad` is the input to integration: HVG selection, CellAssign seed labels, scVI,
then scANVI, written by `07_finalize.py --method scanvi` to `integrated_scanvi.h5ad` (obsm `X_scanvi`).
`as_run/03b_refilter_umap.py` then runs a second doublet pass and the atlas UMAP and writes the deposited
`hgsc_atlas_scanvi.h5ad`.

Re-running integration needs a GPU and is not required to reproduce downstream results, because the integrated object
is deposited. Stochastic steps are seeded from `config.yml` (`seed: 42`). GPU training (scVI / scANVI) is not
bit-identical across hardware, so the deposited integrated object is the reference for everything downstream.

## Scripts as they were run

The `as_run/` folders, `models/` and `icon7/` hold the scripts that wrote the deposited outputs, as they were run.
Plots and analyses not used in the figures or tables were removed; each of these scripts opens with a note on what
was removed, and scripts with nothing to remove are unchanged. They keep the paths of the project folders they ran
in, and the folder READMEs say what to change to run one. Two sets of scripts in `atlas/` resolve their paths through
`config/` instead:

- `atlas/01_preprocess_qc/03`–`07`, the integration scripts described above;
- `atlas/03_epithelial_nmf/01_epithelial_nmf.py`, a re-implementation of the epithelial NMF that runs on the deposited
  epithelial object. `atlas/03_epithelial_nmf/01b_verify_against_deposited.py` checks it against the deposited NMF.

Steps whose code is not in this repository are marked **[code missing]** in the folder READMEs, with the output each
step produced and where that output is.

## Study codes

Samples are named by study code throughout: code, configuration, file names and the deposited data.

| Code | Samples |
|---|---|
| WT01–WT08 | whole-tissue Xenium sections |
| TMA-P### | TMA patients |
| FTE01–FTE15 | fallopian-tube epithelium cores |

No hospital, biobank or pathology identifiers are included. A few scripts that ran on the laboratory's data translate
identifiers through `config/study_codes.{py,R}`. Their key is private and not included; the scripts look for it at
`$STUDY_CODE_KEY`, by default `$OUTPUT_ROOT/private_keys/SuppTable8_code_key_PRIVATE.csv`. Coded names pass through
unchanged, so the key is not needed on the deposited data.

## Contact

David Cook (dacook@ohri.ca), Ottawa Hospital Research Institute.

## License

Released under the **MIT License**: free to use, modify and redistribute (including commercially), with no
restrictions beyond retaining the copyright and license notice. See [`LICENSE`](LICENSE).
