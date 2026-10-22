#!/usr/bin/env python3
"""
Fig 2a prep — per-state functional summaries (z-scored) from per-cell outputs
=============================================================================
PURPOSE
    Builds the four z-scored per-state tables read by
    figures/figure2/01_fig2a_radar_characterization.py (and by Figure S8) from the
    PER-CELL activity / flux outputs, joined to the final-atlas epithelial state
    (load_epithelial()['state']): the deposited per-cell PROGENy and DoRothEA
    (decoupler MLM) and Hallmark (scanpy score_genes) scores, and the scFEA flux of
    all 575,366 epithelial cells from fig2a_scfea_all_cells.py.

    Per table: mean per state over the cells present in the per-cell file, then
    z-score each feature across the four states (SecA, Intermediate, SecB,
    Ciliated) with the sample SD (ddof = 1); features with zero variance -> 0.
    All features are kept (PROGENy 14, Hallmark 50, DoRothEA 294, scFEA 168 modules).

    Coverage check: every per-cell barcode must be in the final epithelial set
    (asserted); the number of cells per state used for each table is written out.

INPUTS
    data_root/2026_final_atlas/output/21_epitype_functional_characterization/
        21_{progeny,hallmark,dorothea}_per_cell.parquet
    output_root/figure2/fig2a_scfea_all/21_flux_per_cell.parquet  (fig2a_scfea_all_cells.py)
    figures/_prep/atlas_final_tables.py  (load_epithelial)

OUTPUTS (output_root/figure2/fig2a_functional/)
    21_{progeny,hallmark,dorothea,flux}_zscored.csv   (input of figures/figure2/01_fig2a_radar_characterization.py)
    21_{progeny,hallmark,dorothea,flux}_means.csv     per-state means behind the z-scores
    coverage.csv                                      cells per state per table

MANUSCRIPT PANEL(S): Fig 2a (input).
RUNTIME TIER: fast-moderate (reads a 1.4 GB DoRothEA parquet).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402
from figures._prep.atlas_final_tables import STATE_ORDER, load_epithelial  # noqa: E402

SRC = path("data_root", "2026_final_atlas", "output", "21_epitype_functional_characterization")
# scFEA flux: ALL 575,366 epithelial cells (figures/_prep/fig2a_scfea_all_cells.py); PROGENy, Hallmark and
# DoRothEA per-cell scores come from SRC.
FLUX_SRC = path("output_root", "figure2", "fig2a_scfea_all")
assert os.path.exists(os.path.join(FLUX_SRC, "21_flux_per_cell.parquet")), "run figures/_prep/fig2a_scfea_all_cells.py first"
OUT_DIR = path("output_root", "figure2", "fig2a_functional")
os.makedirs(OUT_DIR, exist_ok=True)
TABLES = ["progeny", "hallmark", "dorothea", "flux"]

state = load_epithelial(columns=["state"])["state"]
cov_rows, diff_rows = [], []
for name in TABLES:
    pc = pd.read_parquet(os.path.join(FLUX_SRC if name == "flux" else SRC, f"21_{name}_per_cell.parquet"))
    pc.index = pc.index.astype(str)
    assert pc.index.is_unique
    missing = int((~pc.index.isin(state.index)).sum())
    assert missing == 0, f"{name}: {missing} per-cell barcodes not in the final epithelial set"
    st = state.loc[pc.index]
    # Native dtype (float32 for scFEA), as in Figure S8: casting to float64 changes the
    # z-score of near-constant flux modules (e.g. M_93) by up to 0.02 after the 4-state scaling.
    means = pc.groupby(st.values, observed=True).mean().reindex(STATE_ORDER)
    means.index.name = "group"
    sd = means.std(ddof=1)
    z = ((means - means.mean()) / sd.replace(0, np.nan)).fillna(0.0)
    means.to_csv(os.path.join(OUT_DIR, f"21_{name}_means.csv"))
    z.to_csv(os.path.join(OUT_DIR, f"21_{name}_zscored.csv"))
    n_by = st.value_counts().reindex(STATE_ORDER)
    cov_rows.append({"table": name, "n_cells": len(pc), "n_features": pc.shape[1],
                     "pct_of_epithelium": round(100 * len(pc) / len(state), 2),
                     **{f"n_{s}": int(n_by[s]) for s in STATE_ORDER}})

pd.DataFrame(cov_rows).to_csv(os.path.join(OUT_DIR, "coverage.csv"), index=False)
print(pd.DataFrame(cov_rows).to_string(index=False))
