#!/usr/bin/env python3
# Fig 6h atlas macrophage row: per-sample paired contrast, Hypoxic macrophage vs other macrophages, from the full-pass
# pseudobulk accumulators of 08_fig6fh_atlas_macrophage_pseudobulk.py (ALL 309,111 atlas macrophage-lineage cells; no subsampling).
#   per sample (>= 20 cells in each group): effect = mean log1p(CP10k) in Hypoxic macrophages - mean in other macrophages
#   (C1Q tissue-resident, Monocyte-derived, Cycling, Cycling C1Q+ tissue-resident; classical monocytes and cDC2 excluded);
#   summary across samples: n_up, mean effect, two-sided Wilcoxon signed-rank vs 0 (181 paired samples).
# Output: <fig6_singler>/atlas_hypoxic_vs_other_paired.csv (atlas macrophage row of Fig 6h, read by 13_fig6_draw_panels.py).
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
import os
O = os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "fig6_singler")
G = pd.read_csv(f"{O}/atlas_pb/groups.csv"); genes = open(f"{O}/atlas_pb/genes.txt").read().split("\n"); gi = {g: i for i, g in enumerate(genes)}
LS = np.load(f"{O}/atlas_pb/lognorm_sum.npy")
OTH = ["C1Q tissue-resident macrophage", "Monocyte-derived macrophage", "Cycling macrophage", "Cycling C1Q+ tissue-resident macrophage"]
LIST = ["TNF","NFKB1","NFKBIA","BIRC3","BCL2L1","MCL1","VEGFA","CD55","CD47","LGALS9","CD274","C7","SLC16A3","LDHA","SLC2A1","ICAM1","VCAM1",
        "CD44","SPP1","CTSL","MMP7","MMP9","MIF","NDRG1","C1QC","MRC1","CD163","STAB1","APOE","HAVCR2","LAG3"]
D = []
for s, g in G.groupby("sample_id"):
    h = g[g.celltype_level2 == "Hypoxic macrophage"]; o = g[g.celltype_level2.isin(OTH)]
    if h.n_cells.sum() < 20 or o.n_cells.sum() < 20: continue
    D.append(LS[h.index].sum(0) / h.n_cells.sum() - LS[o.index].sum(0) / o.n_cells.sum())
D = np.array(D); rows = []
for g in LIST:
    if g not in gi: continue
    v = D[:, gi[g]]; rows.append(dict(feature=g, n=len(v), n_up=int((v > 0).sum()), mean_effect=v.mean(), p=wilcoxon(v).pvalue if np.any(v != 0) else np.nan))
pd.DataFrame(rows).to_csv(f"{O}/atlas_hypoxic_vs_other_paired.csv", index=False); print("n samples", len(D))
