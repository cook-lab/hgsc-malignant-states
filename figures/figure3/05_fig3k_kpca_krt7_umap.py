#!/usr/bin/env python3
"""
Figure 3k -- KPCA.B mouse epithelial Krt7 UMAP
==============================================
PURPOSE
    Krt7 expression on the epithelial UMAP of the KPCA.B mouse tumours, beside the polarization of Fig 3j.

    * Cells: the 12,429 tumour (epithelial) cells with the dedicated epithelial UMAP (epi_umap_1/2), as in Fig 3j.
    * Value: Krt7 in the log-normalised matrix (X). Colour limits run from 0 to the 99th percentile, on the beige ->
      SecB brown map (SECB_CMAP). Cells are drawn in a random order (seed 0).
    * Layout: figures/figure3/_kpca_umap_panel.py (1.8 in square with the colour bar underneath), shared with Fig 3j.

INPUTS (read-only)
    $DATA_ROOT/202608_kpca_secb_mouse/output/33_handoff/kpca_final_annotated.h5ad   ($KPCA_H5AD)

OUTPUTS
    figures_dir/figure3/Fig3_k_mouse_KRT7.{svg,png}

MANUSCRIPT PANEL(S): Fig 3k.
RUNTIME TIER: fast (seconds; reads the 20,598-cell object).
USAGE (repo root): OUTPUT_ROOT=... python figures/figure3/05_fig3k_kpca_krt7_umap.py
"""
import os
import sys

import anndata as ad
import numpy as np
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
from matplotlib.colors import LinearSegmentedColormap, Normalize

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
from config.config import path  # noqa: E402
from _kpca_umap_panel import epi_umap_panel  # noqa: E402

H5AD = os.environ.get("KPCA_H5AD",
                      os.path.join(os.environ.get("DATA_ROOT", "/path/to/deposit"), "202608_kpca_secb_mouse/output/33_handoff/kpca_final_annotated.h5ad"))
FIG = path("figures_dir", "figure3")
os.makedirs(FIG, exist_ok=True)
SEED = 0
SECB_CMAP = LinearSegmentedColormap.from_list(   # beige -> SecB brown
    "beige_secb_brown", [(0.00, "#DDD5CA"), (0.40, "#DDD5CA"), (0.55, "#C0A880"),
                         (0.70, "#9A7D55"), (0.85, "#6E5535"), (1.00, "#3D2B15")], N=256)

a = ad.read_h5ad(H5AD)
tum = a.obs["epi_umap_1"].notna().values
assert tum.sum() == 12_429 and (a.obs.loc[tum, "celltype"].astype(str) == "Tumour").all()
t = a[tum]
ki = t.var_names.get_loc("Krt7")
krt7 = (t.X[:, ki].toarray().ravel() if sp.issparse(t.X) else np.asarray(t.X[:, ki]).ravel()).astype(float)
vmax = float(np.quantile(krt7, 0.99))
ex = t.obs["epi_umap_1"].to_numpy(float); ey = t.obs["epi_umap_2"].to_numpy(float)
eo = np.random.default_rng(SEED).permutation(t.n_obs)
stem = os.path.join(FIG, "Fig3_k_mouse_KRT7")
epi_umap_panel(ex, ey, krt7, eo, SECB_CMAP, Normalize(vmin=0.0, vmax=vmax), "Krt7", "(log-norm)",
               [0.0, vmax], ["0.00", f"{vmax:.2f}"], stem)
print(f"[save] {stem}.svg/.png ({t.n_obs:,} cells; colour limits 0 to {vmax:.3f}, the 99th percentile)")
