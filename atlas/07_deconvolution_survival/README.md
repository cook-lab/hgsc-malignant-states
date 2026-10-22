# atlas/07_deconvolution_survival: TCGA-OV epithelial fraction and SecA / SecB scores

This stage scores TCGA-OV bulk tumours for secretory polarization. Each tumour gets the mean z-score of the 177
SecA-specific and of the 177 SecB-specific NMF genes (`../03_epithelial_nmf`, step 2). BayesPrism, run with an atlas
single-cell reference, estimates each tumour's epithelial fraction, which is used as a covariate in the Cox models of
Fig 7g–j.

## Run order

Paths without a root are under `$DATA_ROOT/2026_final_atlas/`. The three steps are the `as_run/` scripts that wrote
the deposited files.

| Step | Script | Reads | Writes (deposited, in `output/22_tcga_deconvolution/`) | Used by |
|---|---|---|---|---|
| 1 | `as_run/22a_cibersort_reference_v2.py` | `output/17_nmf_labels/celltype_nmf_mapping.csv`; `celltype_h5ad/*.h5ad` (the 13 lineage subsets); `data/cibersort_data_prev/tcga_ecotyper.txt` | `cibersortx_sc_reference_v2.txt` (300 cells per group; genes expressed in ≥ 10% of the cells of at least one group; SecA, Intermediate, SecB and ciliated epithelium plus grouped non-epithelial types), `cibersortx_phenotypes_v2.txt`, `cibersortx_mixture.txt` (a link to the TCGA-OV expression table) | step 2 |
| 2 | `as_run/22b_bayesprism_deconv.R` | step 1 files | `CIBERSORTx_Results.txt` (BayesPrism cell-type fractions per tumour, in CIBERSORTx format) | step 3 |
| 3 | `as_run/22d_signature_survival.py` | `CIBERSORTx_Results.txt`; `output/11e_nmf_characterization/11e_gene_classification.csv`; `data/cibersort_data_prev/{tcga_ecotyper.txt, tcga_hla_clinical.csv}` | `22d_signature_scores.csv`: per tumour, the SecA, SecB and polarization scores (plain and loading-weighted), median groups, the polarization tertile, the epithelial fraction (sum of the epithelial types) and the clinical covariates | `figures/figure7/05_fig7ghij_tcga_km_forest.R` (Fig 7g–j) |

**Polarization sign.** Step 3 writes `polarization` as the SecA score minus the SecB score. Its tertiles
(`polar_tertile`) are labelled SecB-like (lowest), Intermediate and SecA-like (highest). The figure script uses the
tertiles for the Kaplan–Meier curves (SecB-like vs SecA-like) and −`polarization`, i.e. SecB minus SecA, in the Cox
models.

Fig 7g–j reads only the deposited `22d_signature_scores.csv`. The Kaplan–Meier curves and Cox models shown in the
paper are fitted in `figures/figure7/05_fig7ghij_tcga_km_forest.R`.

**Software.** Steps 1 and 3 ran in Python (E2). Step 2 ran in R (E1, BayesPrism; Gibbs sampling, about 10–30 minutes).

## `as_run/`

These are trimmed copies of the scripts that ran in the atlas project folder `202601_final_atlas`. The code that wrote
the files the next step and Fig 7g–j read is kept as it ran. Removed: the gene-check prints and the fallback signature
matrix of `22a` (for the CIBERSORTx web portal; the fractions used come from BayesPrism), the BayesPrism outputs that
nothing reads in `22b` (posterior RDS, fractions CSV, run summary), and the Kaplan–Meier, Cox, plotting and report code
of `22d`. Each copy opens with a note on what was removed. The untrimmed scripts stay in `202601_final_atlas`. The
copies keep that project's paths.

## Missing code

- **[code missing]** Download and preparation of the TCGA-OV inputs `data/cibersort_data_prev/tcga_ecotyper.txt`
  (bulk expression) and `tcga_hla_clinical.csv` (clinical table). Both files are read as they are.
