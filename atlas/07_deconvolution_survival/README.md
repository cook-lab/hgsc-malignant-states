# atlas/07_deconvolution_survival

TCGA-OV bulk deconvolution (BayesPrism), signature-based survival, validation,
and the consensusOV TCGA-subtype mapping.

## Ordered scripts (inputs → outputs)

| # | Script | Inputs | Key outputs |
|---|--------|--------|-------------|
| 01 | `01_cibersort_reference.py` | `03_epithelial_nmf/celltype_nmf_mapping.csv`; per-compartment `celltype_h5ad/*.h5ad`; `tcga_ecotyper.txt` | `cibersortx_sc_reference_v2.txt`, `cibersortx_phenotypes_v2.txt`, `cibersortx_sig_matrix.txt`, `cibersortx_mixture.txt` |
| 02 | `02_bayesprism_deconv.R` | 01 reference + phenotypes + mixture | `CIBERSORTx_Results.txt`, `bayesprism_fractions.csv`, `bayesprism_theta.rds`, summary |
| 03 | `03_bayesprism_merged.R` | 01 reference + `cibersortx_phenotypes_merged.txt` + mixture | `bayesprism_merged_results.txt`, `bayesprism_merged_fractions.csv` (validation run) |
| 04 | `04_survival_analysis.py` | `CIBERSORTx_Results.txt`; `tcga_hla_clinical.csv` | `22c_*` KM/Cox tables + figures; `22c_survival_report.html` |
| 05 | `05_signature_survival.py` | `CIBERSORTx_Results.txt`; `11e_nmf_characterization/11e_gene_classification.csv`; `tcga_ecotyper.txt`; `tcga_hla_clinical.csv` | **`22d_signature_scores.csv`** (KEY); `22d_*` KM/Cox + report |
| 06 | `06_validate_survival.py` | `22d_signature_scores.csv`; `bayesprism_merged_results.txt`; gene class; TCGA expr/clinical | `22f_validation_report.html`; permutation/bootstrap/LOO tables + figures |
| 07 | `07_consensusov_export.py` | `atlas_final` h5ad (X via h5py); `celltype_nmf_mapping.csv` | `consensusov/pseudobulk_{bulk,epi}_counts.tsv.gz`, `pseudobulk_metadata.csv`, `gene_list.tsv` |
| 08 | `08_consensusov_score.R` | 07 pseudobulk | `consensusov/consensusov_calls_{bulk,epi}.csv`, classifier calls, gene id map |
| 09 | `09_consensusov_report.py` | 07 metadata + 08 calls | `consensusov/tables/20c_per_sample_joined.csv` (KEY) + summary/stats; `20_consensusov_report.html` |

## Figures / panels supported

- `22d_signature_scores.csv` → **Fig 7E** (KM OS), **7F** (KM PFS), **7G** (stepwise Cox forest)
- `04_survival_analysis.py` + `06_validate_survival.py` → robustness support for Fig 7E/F/G
- `consensusov/tables/20c_per_sample_joined.csv` + `20c_subtype_composition_summary.csv` → **Fig 3H** (epitype × TCGA subtype)

## Conventions applied

Central config paths; header docstrings; epithelial label **Transitioning →
Intermediate** (BayesPrism `Intermediate_epithelium`; export `pct_intermediate`).
Seeds from config: `random_state=SEED` (reference downsample); `set.seed(CFG$seed)`
for BayesPrism Gibbs sampling (02/03) **and the consensusOV random forest (08)**
— fixing the audit's non-determinism finding; permutation / bootstrap / subgroup
RNGs in 06 all use `SEED`. Superseded `22a_cibersort_reference` (v1) excluded.
