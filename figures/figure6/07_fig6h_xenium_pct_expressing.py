#!/usr/bin/env python3
# Xenium core-8: % of cells expressing each Fig 6h gene (count > 0) in niche macrophages, SecB epithelium, and SecB cells
# <=30 um of a niche macrophage; per-section percentages averaged across sections. All cells (no subsampling).
# Outputs (<fig6_singler>/): xenium_pct_persample.csv (per-section record), xenium_pct_summary.csv (dot size of the
# Xenium rows of Fig 6h, read by 13_fig6_draw_panels.py).
import gzip, io, numpy as np, pandas as pd
from scipy.io import mmread; from scipy import sparse; from scipy.spatial import cKDTree
import os
R = os.environ.get("OUTPUT_ROOT", "output"); DEG = f"{R}/fig6_viewer/degdata"; SG = f"{R}/singler_recompute"; O = f"{R}/fig6_singler"
WT = ["WT01","WT02","WT03","WT04","WT05","WT06","WT07","WT08"]
GENES = ["VEGFA","BIRC3","NFKB2","CD55","CD47","LGALS9","IL32","LDHA","SLC2A1","NDRG1","SLC16A3","SPP1","CD44","ITGB6","ICAM1","CTSL","MMP7","TGM2","TIMP1"]
ni = pd.read_csv(f"{O}/mac_niche_WT.csv.gz"); rows = []
for s in WT:
    d = pd.read_csv(f"{SG}/{s}_singler_cells.csv.gz", usecols=["barcode", "x", "y", "cell_label_singler"])
    genes = np.array([l.strip() for l in open(f"{DEG}/{s}_genes.txt")]); gi = {g: i for i, g in enumerate(genes)}
    cb = pd.read_csv(f"{DEG}/{s}_cells.csv.gz", usecols=["barcode"]).barcode.values; pos = pd.Series(np.arange(len(cb)), index=cb)
    with gzip.open(f"{DEG}/{s}_counts.mtx.gz", "rb") as f: M = sparse.csr_matrix(mmread(io.BytesIO(f.read()))).T.tocsc()
    d = d[d.barcode.isin(pos.index)]; nn = ni[(ni["sample"] == s) & (ni.group == "niche")].barcode
    niche = d[d.barcode.isin(nn)]; secb = d[d.cell_label_singler == "SecB epithelium"]
    dist, _ = cKDTree(niche[["x", "y"]].values).query(secb[["x", "y"]].values, distance_upper_bound=30) if len(niche) else (np.full(len(secb), np.inf), None)
    groups = {"mac_niche": niche.barcode, "secb": secb.barcode, "secb_contact": secb.barcode[dist <= 30]}
    for gname, bc in groups.items():
        r = pos.loc[bc].values
        if len(r) < 20: continue
        for g in GENES:
            if g in gi: rows.append(dict(sample=s, group=gname, gene=g, n=len(r), pct=100 * (M[r, gi[g]].toarray().ravel() > 0).mean()))
    print(s, {k: len(v) for k, v in groups.items()}, flush=True)
P = pd.DataFrame(rows); P.to_csv(f"{O}/xenium_pct_persample.csv", index=False)
S = P.groupby(["group", "gene"]).pct.mean().reset_index(); S.to_csv(f"{O}/xenium_pct_summary.csv", index=False)
print(S.pivot(index="gene", columns="group", values="pct").round(1).to_string()); print("XPCT_DONE")
