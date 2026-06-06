# atlas/03_epithelial_nmf

Epithelial NMF — learns the SecA/SecB polarization programs, writes the
canonical `celltype_nmf` epitype schema, and runs the cross-platform UCell
scoring. The NMF `Factor_2` usage and the `celltype_nmf` labels are the
substrate for almost every downstream epitype analysis.

## Ordered scripts (inputs → outputs)

| # | Script | Inputs | Key outputs |
|---|--------|--------|-------------|
| 01 | `01_epithelial_nmf.py` | epithelial h5ad (`celltype_h5ad/hgsc_atlas_final_epithelial.h5ad`); SecA/SecB sigs from `shared/signatures.yml` | `output/03_epithelial_nmf/11d_nmf_usage.csv` (+ `_raw`, `11d_nmf_loadings.csv`, `11d_factor_scores.csv`); factor/polarization figures |
| 02 | `02_prepare_nmf_labels.py` | `atlas_final` h5ad obs + `11d_nmf_usage.csv` | `celltype_nmf` obs column written into `hgsc_atlas_final.h5ad`; `celltype_nmf_mapping.csv`, `celltype_nmf_summary.csv` |
| 03 | `03_ucell_atlas_export.py` | epithelial h5ad + `celltype_nmf_mapping.csv` | `ucell_atlas/atlas_secretory_{counts.mtx.gz,barcodes.tsv,genes.tsv,metadata.csv}` |
| 04 | `04_ucell_atlas_scoring.R` | 03 exports + organoid Seurat object + `shared/signatures.yml` | `ucell_atlas/atlas_ucell_scores.csv`, `shared_gene_list.txt` |
| 05 | `05_ucell_atlas_report.py` | `atlas_ucell_scores.csv` + `atlas_secretory_metadata.csv` | `ucell_atlas/18c_ucell_atlas_report.html`, summary stats + cutoff figures |

## Figures / panels supported

- `11d_nmf_usage.csv` (Factor_2) → **Fig 1G/H, Fig 1I, Fig 2A/C, SF5, SF6, Supp Data 4**
- `celltype_nmf` schema → consumed by stages 05/06/07 (Fig 1J, 3H, 5F, 7E/F/G, Supp Data 7)
- UCell scoring (03–05, OLD 18b — noBCAM-matching, NOT 18b_v2) → cross-platform
  SecB scoring for **Fig 3B / SF11**

## Conventions applied

Central config paths; header docstrings; `np.random.seed(SEED)` / `set.seed(CFG$seed)`
at every NMF/subsample/UCell step; epithelial label **Transitioning → Intermediate**;
SecA/SecB loaded from `shared/signatures.yml`; dead `--figures-only` exploration kept
(it is a documented re-render path), superseded `18b_v2` excluded.
