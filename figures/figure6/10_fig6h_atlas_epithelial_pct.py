#!/usr/bin/env python3
# Atlas epithelial object: % of SecA / SecB cells (schema_nmf, the states of the Fig 2g pseudobulk) expressing each Fig 6h
# gene. ALL cells (no subsampling); reads only CSR indices in row blocks. Output: fig6_singler/atlas_epi_pct.csv (dot size of
# the atlas SecB row of Fig 6h, read by 13_fig6_draw_panels.py)
import numpy as np, pandas as pd, h5py
import os
SRC = os.path.join(os.environ.get("DATA_ROOT", "/path/to/deposit"), "2026_final_atlas/celltype_h5ad/hgsc_atlas_final_epithelial.h5ad")
META = os.path.join(os.environ.get("DATA_ROOT", "/path/to/deposit"), "2026_final_atlas/output/fig_secretory_polarization/data/meta.parquet")
GENES = ["VEGFA","BIRC3","NFKB2","CD55","CD47","LGALS9","IL32","LDHA","SLC2A1","NDRG1","SLC16A3","SPP1","CD44","ITGB6","ICAM1","CTSL","MMP7","TGM2","TIMP1"]
f = h5py.File(SRC, "r"); obs_idx = [x.decode() if isinstance(x, bytes) else x for x in f["obs"][f["obs"].attrs["_index"]][:]]
epi = pd.read_parquet(META, columns=["schema_nmf"])["schema_nmf"].reindex(obs_idx).values
g = f["layers/counts"] if "layers" in f and "counts" in f["layers"] else f["X"]; ip = g["indptr"][:]
var = [x.decode() if isinstance(x, bytes) else x for x in f["var"][f["var"].attrs["_index"]][:]]; gi = {v: i for i, v in enumerate(var)}
tg = np.array([gi[x] for x in GENES if x in gi]); names = [x for x in GENES if x in gi]; pos = {t: k for k, t in enumerate(tg)}
cnt = {"SecA": np.zeros(len(tg), np.int64), "SecB": np.zeros(len(tg), np.int64)}; ncell = {"SecA": 0, "SecB": 0}
B = 25000
for s in range(0, len(obs_idx), B):
    e = min(s + B, len(obs_idx)); ind = g["indices"][ip[s]:ip[e]]; rows = np.repeat(np.arange(e - s), np.diff(ip[s:e + 1]))
    m = np.isin(ind, tg); lab = epi[s:e]
    for et in ("SecA", "SecB"):
        sel = lab == et; ncell[et] += int(sel.sum())
        if sel.any():
            mm = m & sel[rows]; k = np.array([pos[t] for t in ind[mm]], int); np.add.at(cnt[et], k, 1)
    if (s // B) % 4 == 0: print(f"{e}/{len(obs_idx)}", flush=True)
out = pd.DataFrame({"gene": names, "pct_SecA": 100 * cnt["SecA"] / ncell["SecA"], "pct_SecB": 100 * cnt["SecB"] / ncell["SecB"]})
out.to_csv(os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "fig6_singler/atlas_epi_pct.csv"), index=False); print(ncell); print(out.round(1).to_string(index=False)); print("EPI_PCT_DONE")
