# atlas/04_functional

Epitype functional characterization: PROGENy pathway activity, Hallmark gene
sets, DoRothEA TF activity, and scFEA metabolic flux for the 4 NMF epitypes
(SecA / Intermediate / SecB / Ciliated).

## Ordered scripts (inputs → outputs)

| # | Script | Inputs | Key outputs |
|---|--------|--------|-------------|
| 01 | `01_epitype_functional_characterization.py` | epithelial h5ad; `output/03_epithelial_nmf/11d_nmf_usage.csv` (Factor_2); `2026_final_atlas/tools/scFEA/data/*` | `output/04_functional/21_{progeny,hallmark,dorothea,flux}_zscored.csv` (+ means/radar/per-cell); heatmaps, violins, radar panels |

## Figures / panels supported

- `21_{progeny,hallmark,dorothea,flux}_zscored.csv` → **Fig 2A** (radars), **SF8** (heatmaps)

## Conventions applied

Central config paths; header docstring; `np.random.seed(SEED)` + `torch.manual_seed(SEED)`
for the scFEA training/downsample; **Transitioning → Intermediate**.

### Documented quirk fixed
- **scFEA flux label-map bug**: the original keyed the readable module label map
  by the integer `Module_id` column (1,2,…) while the flux matrix columns are the
  string module index (`M_1`, `M_2`, …), so labels never matched and heatmaps/radars
  showed raw `M_x`. Fixed by reading the annotation with `index_col=0` (the `M_x`
  index) so labels resolve to `Compound_IN->Compound_OUT`.

### Dependency note
- The flux step trains a **PyTorch** graph neural network (scFEA, Chang et al. 2021);
  `torch` is required (CPU is sufficient). `decoupler` is required for PROGENy/DoRothEA.
