# atlas/05_cnv: CopyKAT copy-number calls and subclones

This stage runs CopyKAT once per atlas sample. Each run uses the sample's epithelial cells plus a diploid reference of
fibroblast, T/NK and endothelial cells from the same sample. Aneuploid epithelial cells are then split into subclones.
The results test whether SecA, Intermediate and SecB cells coexist within genetic clones (Fig 1l,m; Figure S7) and
mark the malignant cells used in Fig 1g.

## Run order

Paths without a root are under `$DATA_ROOT/2026_final_atlas/`. Both steps are the `as_run/` scripts that produced the
deposited calls.

| Step | Script | Reads | Writes (deposited) | Used by |
|---|---|---|---|---|
| 1 | `as_run/19_cnv_extract.py` | `celltype_h5ad/hgsc_atlas_final_{epithelial,fibroblast,t_nk_cell,endothelial}.h5ad`; `output/11d_epithelial_nmf/11d_nmf_usage.csv` | `output/19_cnv/per_sample/<sample>/{counts.mtx.gz, genes.txt, barcodes.csv, ref_barcodes.txt}`; `output/19_cnv/tables/sample_manifest.csv`. A sample is eligible with ≥ 200 epithelial and ≥ 50 reference cells; the reference is capped at 500 cells per sample (seed 42). | step 2; `barcodes.csv` is also read by Fig 1l,m and Figure S7 |
| 2 | `as_run/19_cnv_copykat.R <sample_dir> [n_cores]`, run once per sample folder | one `per_sample/<sample>/` folder | `copykat_prediction.csv` (aneuploid / diploid per cell), `copykat_CNA_results.txt`, `copykat_subclones.csv`, `copykat_subclone_qc.csv`, `DONE.txt` (251 samples completed) | Fig 1g (`figures/_prep/fig1g_polarization_pc_alignment.py`); Fig 1l,m (`figures/figure1/08_fig1lm_clone_composition_and_polarization.py`); Figure S4c (`figures/supplementary/10_figS4c_cnv_aneuploid_by_state.py`); Figure S7 (`figures/supplementary/13_figS7_cnv_clone_validation.py`), then the Discussion number (`figures/_prep/text_state_coexistence.py`) |

**Drivers for step 2.**
- `as_run/19_cnv_run.sh` runs samples in parallel on one machine and resumes after a stop.
- On an HPC cluster, `as_run/19_cnv_bundle.sh` packs the inputs and `as_run/19_cnv_slurm.sh` submits one array task per
  sample.

**Subclones.** Aneuploid cells are clustered on 1 − Pearson distance with Ward.D2. The number of clones k is chosen
in 2–5 by silhouette width. A sample is called monoclonal when the best silhouette is below 0.15. The figure scripts
keep clones with ≥ 20 cells.

**Sample keys.** The `per_sample/` folders are named by the atlas sample id as deposited, before the Nath 2021
relabelling. Figure scripts join on `sample_id_deposited`, which `figures/_prep/nath2021_fix.py` adds.

**Software.** CopyKAT ran in R (E1, with copykat installed from GitHub). Extraction ran in Python.

## `as_run/`

These are the scripts that ran in the atlas project folder `202601_final_atlas`, unchanged: all their code writes or
runs the per-sample inputs and calls, so nothing was trimmed. They keep that project's paths. Comments in
`19_cnv_slurm.sh` and `19_cnv_bundle.sh` name follow-up scripts of that project (status, independence and coexistence
checks) that are not part of the paper or of this repository.

## Missing code

None for this stage.
