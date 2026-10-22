#!/usr/bin/env python3
# Fig 6 cross-platform concordance: each candidate gene x side (macrophage / SecB) called per platform, then a verdict.
# Calls: Xenium (n=8, Wilcoxon on per-section effects): up/down if p<0.05, trend if p<0.1.
#        Visium (WT01 + WT03, descriptive): up/down if both sections |log2FC|>=0.1 same sign with bin p<0.05;
#               mac side "carry-over" if the gene fails the carry-over rule in either section while up.
#        Atlas mac (181 paired samples; hypoxic vs other macrophages): up/down if p<0.05 and |delta log1p|>=0.05.
#        Atlas SecB (Fig 2g pseudobulk DESeq2, 261 tumours, paired ~sample + state): up/down if padj<0.05 and |log2FC|>=0.25.
# Verdict: CONTRADICTED if any two platforms are significant in opposite directions; SUPPORTED if >=2 agree;
#          SINGLE if exactly one platform is significant; else NONE.
# Output: <fig6_singler>/fig6_concordance_table.csv; 13_fig6_draw_panels.py reads its atlas macrophage values for the
# Fig 6h genes that 09 does not list.
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
import os
O=os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "fig6_singler")
X=pd.read_csv(f"{O}/reciprocal_summary.csv"); V=pd.read_csv(f"{O}/visium_candidates_2sections.csv")
D2=pd.read_csv(os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "figure2/atlas_volcano_secA_secB_pseudobulk_deseq2.csv")).set_index("gene")
G=pd.read_csv(f"{O}/atlas_pb/groups.csv"); genes=open(f"{O}/atlas_pb/genes.txt").read().split("\n"); gi={g:i for i,g in enumerate(genes)}
LS=np.load(f"{O}/atlas_pb/lognorm_sum.npy")
OTH=["C1Q tissue-resident macrophage","Monocyte-derived macrophage","Cycling macrophage","Cycling C1Q+ tissue-resident macrophage"]
DA=[]
for s,g in G.groupby("sample_id"):
    h=g[g.celltype_level2=="Hypoxic macrophage"]; o=g[g.celltype_level2.isin(OTH)]
    if h.n_cells.sum()>=20 and o.n_cells.sum()>=20: DA.append(LS[h.index].sum(0)/h.n_cells.sum()-LS[o.index].sum(0)/o.n_cells.sum())
DA=np.array(DA)
CAND=V.gene.unique().tolist()+[g for g in X.feature.unique() if not g.startswith("set:") and g not in V.gene.unique()]
def xcall(cn,g):
    r=X[(X.comparison==cn)&(X.feature==g)]
    if not len(r): return "n/a",""
    r=r.iloc[0]; d="up" if r.mean_effect>0 else "down"; k=r.n_up if r.mean_effect>0 else r.n_sec-r.n_up
    txt=f"{r.mean_effect:+.2f} ({k}/{int(r.n_sec)}, p={r.p:.3f})"
    return (d if r.p<0.05 else ("trend "+d if r.p<0.1 else "ns")), txt
def vcall(side,g):
    r=V[(V.side==side)&(V.gene==g)]
    if len(r)<2: return "n/a",""
    txt=" / ".join(f"{e:+.2f}" for e in r.effect)
    if side=="mac" and ((r.effect>0)&(~r.mac_sourced.astype(bool))).any(): return "carry-over",txt
    sig=(r.p<0.05)&(r.effect.abs()>=0.1)
    if sig.all() and (np.sign(r.effect).nunique()==1): return ("up" if r.effect.iloc[0]>0 else "down"),txt
    if sig.all(): return "mixed",txt
    return "ns",txt
def amac(g):
    if g not in gi: return "n/a",""
    v=DA[:,gi[g]]; p=wilcoxon(v).pvalue if np.any(v!=0) else 1; m=v.mean()
    txt=f"{m:+.2f} ({int((v>0).sum())}/{len(v)})"
    return (("up" if m>0 else "down") if (p<0.05 and abs(m)>=0.05) else "ns"),txt
def asecb(g):
    if g not in D2.index: return "n/a",""
    r=D2.loc[g]; txt=f"{r.log2FoldChange:+.2f} (padj={r.padj:.1g})"
    return (("up" if r.log2FoldChange>0 else "down") if (r.padj<0.05 and abs(r.log2FoldChange)>=0.25) else "ns"),txt
def verdict(calls):
    s=[c for c in calls if c in ("up","down")]
    if "up" in s and "down" in s: return "CONTRADICTED"
    if len(s)>=2: return "SUPPORTED "+s[0]
    if len(s)==1: return "SINGLE "+s[0]
    return "NONE"
rows=[]
for g in dict.fromkeys(CAND):
    xm,xmt=xcall("mac_niche_vs_ctrl",g); vm,vmt=vcall("mac",g); am,amt=amac(g)
    xs,xst=xcall("secb_vs_seca",g); xc,xct=xcall("secb_contact",g); vs,vst=vcall("secb",g); as_,ast=asecb(g)
    rows.append(dict(gene=g,xen_mac=xm,xen_mac_v=xmt,vis_mac=vm,vis_mac_v=vmt,atlas_mac=am,atlas_mac_v=amt,MAC_VERDICT=verdict([xm,vm,am]),
                     xen_secb=xs,xen_secb_v=xst,xen_contact=xc,xen_contact_v=xct,vis_secb=vs,vis_secb_v=vst,atlas_secb=as_,atlas_secb_v=ast,
                     SECB_VERDICT=verdict([xs,vs,as_])))
R=pd.DataFrame(rows); R.to_csv(f"{O}/fig6_concordance_table.csv",index=False)
pd.set_option("display.width",250); pd.set_option("display.max_colwidth",40)
print(R[["gene","xen_mac","vis_mac","atlas_mac","MAC_VERDICT","xen_secb","xen_contact","vis_secb","atlas_secb","SECB_VERDICT"]].to_string(index=False))
