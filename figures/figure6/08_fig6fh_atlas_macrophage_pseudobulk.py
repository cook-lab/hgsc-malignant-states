#!/usr/bin/env python3
# Full-pass (no subsampling) pseudobulk of the scRNA atlas macrophage object by sample_id x celltype_level2.
# Reads layers/counts (CSR) in row blocks straight from HDF5 (a local copy -- the CookLab exFAT SSD drops out).
# Outputs (atlas_pb/): counts_sum.npy [groups x genes], lognorm_sum.npy (sum of log1p CP10k), nnz.npy, groups.csv, genes.txt
# Nath 2021 cells take their corrected sample (figures/_prep/nath2021_fix.py; run with OUTPUT_ROOT set to the atlas
# output root); block sums use a sparse group-indicator product. Read by 09, 11 and 13 (Fig 6h), and the source of the
# Fig 6f dot matrix (atlas_subtype_dotmatrix.csv, written by 08b_fig6f_atlas_subtype_dotmatrix.py).
import sys, os, numpy as np, pandas as pd, h5py, anndata as ad
from scipy import sparse
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from figures._prep.nath2021_fix import apply_nath2021_fix  # noqa: E402
SRC=sys.argv[1] if len(sys.argv)>1 else os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "atlas_local/hgsc_atlas_final_macrophage.h5ad")
OUT=os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "fig6_singler/atlas_pb"); os.makedirs(OUT,exist_ok=True)
obs=ad.read_h5ad(SRC,backed="r").obs[["study","sample_id","celltype_level2","dataset"]].copy()
obs.index=obs.index.astype(str)
obs=apply_nath2021_fix(obs)
f=h5py.File(SRC,"r"); g=f["layers/counts"]; ip=g["indptr"][:]; nG=int(g.attrs["shape"][1])
genes=[x.decode() if isinstance(x,bytes) else x for x in f["var"][f["var"].attrs["_index"]][:]]
key=obs.sample_id.astype(str)+"||"+obs.celltype_level2.astype(str)
grp,gid=np.unique(key.values,return_inverse=True); nK=len(grp)
CS=np.zeros((nK,nG)); LS=np.zeros((nK,nG)); NZ=np.zeros((nK,nG),np.int64); NC=np.bincount(gid,minlength=nK)
B=20000
for s in range(0,len(obs),B):
    e=min(s+B,len(obs)); a,b=ip[s],ip[e]
    dat=g["data"][a:b].astype(np.float64); ind=g["indices"][a:b]; loc=ip[s:e+1]-a
    rows=np.repeat(np.arange(e-s),np.diff(loc)); lib=np.bincount(rows,weights=dat,minlength=e-s); lib[lib==0]=1
    ln=np.log1p(dat/lib[rows]*1e4)
    X=sparse.csr_matrix((dat,ind,loc),shape=(e-s,nG)); L=sparse.csr_matrix((ln,ind,loc),shape=(e-s,nG))
    Z=sparse.csr_matrix((np.ones_like(dat),ind,loc),shape=(e-s,nG))
    M=sparse.csr_matrix((np.ones(e-s),(gid[s:e],np.arange(e-s))),shape=(nK,e-s))
    CS+=(M@X).toarray(); LS+=(M@L).toarray(); NZ+=(M@Z).toarray().astype(np.int64)
    print(f"{e}/{len(obs)}",flush=True)
G=pd.DataFrame({"group":grp,"n_cells":NC}); G[["sample_id","celltype_level2"]]=G.group.str.split(r"\|\|",expand=True)
G["dataset"]=G.sample_id.map(obs.groupby("sample_id",observed=True).dataset.first().astype(str))
np.save(f"{OUT}/counts_sum.npy",CS); np.save(f"{OUT}/lognorm_sum.npy",LS); np.save(f"{OUT}/nnz.npy",NZ)
G.to_csv(f"{OUT}/groups.csv",index=False); open(f"{OUT}/genes.txt","w").write("\n".join(genes))
print("ATLAS_PB_DONE", nK, "groups")
