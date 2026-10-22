#!/usr/bin/env python3
"""
Figure S4b (data) — canonical-marker dot-plot statistics on the final atlas (chunked)
====================================================================================

PURPOSE
    Chunked (h5py) dot-plot statistics over ALL
    1,980,703 cells of the final atlas. For each level-1 type x marker gene: fraction of cells
    expressing (X > 0), mean expression (adata.X, float64 accumulation), n_cells, and the
    per-gene min-max normalised mean (norm_expression). X is streamed in CSR row chunks and
    only the 27 marker columns are kept; it is never loaded fully.

INPUTS
    obj("atlas_final")  hgsc_atlas_final.h5ad: X (CSR) and var index only.
    figures._prep.atlas_final_tables.load_obs()  celltype_level1 per cell (final-atlas
        metadata; row order asserted identical to the h5ad obs index).

OUTPUTS
    output_root/02_annotation/09f_dotplot_canonical_markers/dotplot_stats.csv  (26 genes)
    output_root/02_annotation/09f_dotplot_canonical_markers/lilra4_stats.csv   (LILRA4 row)
    Both are read by 09_figS4b_canonical_markers_dotplot.py. (The folder name keeps the atlas
    step name 09f_dotplot_canonical_markers.)

MANUSCRIPT PANEL(S)
    Figure S4b (statistics).

RUNTIME TIER
    slow (streams the full CSR X of the 100 GB atlas; ~tens of minutes on the SMB share).
"""
import os
import sys
import time
from pathlib import Path

import h5py
import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import obj, path  # noqa: E402
from figures._prep.atlas_final_tables import load_obs  # noqa: E402

H5 = obj("atlas_final")
OUT = path("output_root", "02_annotation", "09f_dotplot_canonical_markers")
CT = ["Epithelial", "Mesothelial", "Fibroblast", "Smooth muscle", "Pericyte", "Endothelial",
      "T/NK cell", "B cell", "Plasma cell", "Macrophage", "DC", "Neutrophil", "Mast cell"]
GENES = ['EPCAM','PAX8','MSLN','UPK3B','COL1A1','DCN','MYH11','ACTA2','PDGFRB','RGS5','CDH5','VWF',
         'CD3E','NKG7','MS4A1','CD79A','MZB1','XBP1','CD14','C1QA','FCER1A','CD1C','CSF3R','S100A8',
         'KIT','TPSAB1','LILRA4']
t = time.time()
obs = load_obs(columns=["celltype_level1"])
with h5py.File(H5, "r") as f:
    var = ad.io.read_elem(f["var"]).index.astype(str)
    o = f["obs"]
    h5_idx = np.asarray(ad.io.read_elem(o[o.attrs["_index"]])).astype(str)
    assert np.array_equal(h5_idx, obs.index.values.astype(str)), "loader row order != h5ad obs order"
    grp = obs["celltype_level1"].astype(str).values
    gidx = np.array([var.get_loc(g) for g in GENES])
    colmap = np.full(len(var), -1, dtype=np.int64); colmap[gidx] = np.arange(len(GENES))
    X = f["X"]; indptr = X["indptr"][:]; n = len(indptr) - 1
    assert n == len(obs)
    dense = np.zeros((n, len(GENES)), dtype=np.float32)
    CH = 200_000_000; r0 = 0
    while r0 < n:
        s = indptr[r0]
        r1 = min(max(int(np.searchsorted(indptr, s + CH, side="right") - 1), r0 + 1), n)
        e = indptr[r1]
        ind = X["indices"][s:e]; c = colmap[ind]; m = c >= 0
        if m.any():
            dat = X["data"][s:e][m]
            rows = np.repeat(np.arange(r0, r1), np.diff(indptr[r0:r1 + 1]))[m]
            dense[rows, c[m]] = dat
        print(f"rows {r0:,}-{r1:,}/{n:,} ({time.time()-t:.0f}s)", flush=True)
        r0 = r1
recs = []
for ct in CT:
    mk = grp == ct
    sub = dense[mk]
    for j, g in enumerate(GENES):
        v = sub[:, j]
        recs.append(dict(celltype=ct, gene=g, frac_expressing=float(np.mean(v > 0)),
                         mean_expression=float(np.mean(v, dtype=np.float64)), n_cells=int(mk.sum())))
df = pd.DataFrame(recs)
for g in GENES:
    m = df.gene == g; v = df.loc[m, "mean_expression"]
    df.loc[m, "norm_expression"] = (v - v.min()) / (v.max() - v.min()) if v.max() > v.min() else 0.0
df[df.gene != "LILRA4"].to_csv(os.path.join(OUT, "dotplot_stats.csv"), index=False)
df[df.gene == "LILRA4"].to_csv(os.path.join(OUT, "lilra4_stats.csv"), index=False)
print("done", time.time() - t)
