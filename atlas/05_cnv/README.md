# atlas/05_cnv

CopyKAT CNV inference and the within-clone SecA/SecB coexistence / independence
analyses that show the epitype axis is transcriptional, not clonal.

## Ordered scripts (inputs → outputs)

| # | Script | Inputs | Key outputs |
|---|--------|--------|-------------|
| 01 | `01_cnv_extract.py` | per-compartment `celltype_h5ad/*.h5ad` (epi/fibroblast/T-NK/endothelial); `11d_nmf_usage.csv` | `output/05_cnv/per_sample/<id>/{counts.mtx.gz,genes.txt,barcodes.csv,ref_barcodes.txt}`; `tables/sample_manifest.csv` |
| 02 | `02_cnv_copykat.R` (HPC, one process per sample dir) | a `per_sample/<id>/` dir | `<id>/{copykat_prediction.csv,copykat_CNA_results.txt,copykat_subclones.csv,copykat_subclone_qc.csv,DONE.txt}` |
| 03 | `03_cnv_coexistence.py` | `per_sample/*/copykat_subclones.csv` + `barcodes.csv`; `11d_nmf_usage.csv`; epithelial UMAP | `tables/within_clone_coexistence.csv`, `tables/per_sample_verdict.csv`; coexistence figures |
| 04 | `04_cnv_independence.py` | `per_sample/*/copykat_subclones.csv` + `barcodes.csv` | `tables/{chisq_fisher_results,logistic_regression_results,summary_statistics}.csv`; independence figures |

## Figures / panels supported

- `within_clone_coexistence.csv` → **Fig 1J** (alluvial cache), **SF7** (within-clone bars)
- `per_sample_verdict.csv` + independence tables → CNV-independence narrative (Fig 1J / SF4C / SF7)

## Conventions applied

Central config paths; header docstrings; reference subsample / clustering seeded
from config (`np.random.default_rng(SEED)`, `set.seed(CFG$seed)`); epithelial
label **Transitioning → Intermediate** (epitype tidy columns renamed to
`*_intermediate`); CopyKAT and the SLURM/batch shell wrappers preserved as the
documented HPC step (02 is the per-sample worker).
