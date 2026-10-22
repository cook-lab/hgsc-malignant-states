#!/usr/bin/env python3
# PURPOSE: Fig 6f (top block) -- the atlas dot matrix of 15 genes x 6 macrophage / monocyte level-2 types: mean log1p
#   CP10k and % expressing, pooled over every cell of each type from the 08_fig6fh_atlas_macrophage_pseudobulk.py sums.
# Recovered (2026-10-22): this ran inline in the final run's fig6_singler output folder on 2026-09-25 16:34 UTC and was
#   never saved as a file; it is recovered from the analysis session log. Changed: the folder, which now goes through
#   OUTPUT_ROOT; the code is otherwise as it ran.
# OUTPUTS: <fig6_singler>/atlas_subtype_dotmatrix.csv (celltype, gene, mean, pct; read by 13_fig6_draw_panels.py).
# MANUSCRIPT PANEL(S): Fig 6f. RUNTIME TIER: fast.
import os
O=os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "fig6_singler")
import numpy as np, pandas as pd
G=pd.read_csv(f"{O}/atlas_pb/groups.csv"); genes=open(f"{O}/atlas_pb/genes.txt").read().split("\n"); gi={g:i for i,g in enumerate(genes)}
LS=np.load(f"{O}/atlas_pb/lognorm_sum.npy"); NZ=np.load(f"{O}/atlas_pb/nnz.npy")
CT=["Hypoxic macrophage","Cycling macrophage","Monocyte-derived macrophage","C1Q tissue-resident macrophage","Cycling C1Q+ tissue-resident macrophage","Classical monocyte"]
L=["SPP1","VEGFA","SLC2A1","NDRG1","LDHA","CD44","CTSL","TGM2","C1QC","MRC1","A2M","MMP9","MIF","CD163","STAB1"]
M=[dict(celltype=ct,gene=g,mean=LS[(G.celltype_level2==ct).values,gi[g]].sum()/G.n_cells[G.celltype_level2==ct].sum(),
        pct=100*NZ[(G.celltype_level2==ct).values,gi[g]].sum()/G.n_cells[G.celltype_level2==ct].sum()) for ct in CT for g in L if g in gi]
pd.DataFrame(M).to_csv(f"{O}/atlas_subtype_dotmatrix.csv",index=False); print(sorted(set(g for g in L if g not in gi)))
