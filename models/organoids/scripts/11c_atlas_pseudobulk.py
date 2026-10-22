#!/usr/bin/env python3
# 11c_atlas_pseudobulk.py
# Pseudobulk the HGSC atlas secretory cells by (patient_id x celltype_nmf) so the
# atlas can be scored as BULK, on the same footing as the bulk organoids. Summed raw
# counts preserve within-sample gene ranking (all UCell needs). Streams the 5.5GB mtx.
import gzip, csv, numpy as np, collections, sys

AD  = "/Volumes/CookLab/Sarah/202601_final_atlas/output/18_ucell_atlas"
OUT = "/Volumes/CookLab/Sarah/202603_organoids/data/11_nextgen"

genes = [l.split("\t")[0].strip() for l in open(f"{AD}/atlas_secretory_genes.tsv")]
barcodes = [l.strip() for l in open(f"{AD}/atlas_secretory_barcodes.tsv")]
ng, nc = len(genes), len(barcodes)
print(f"genes={ng} cells={nc}", flush=True)

meta = {r[list(r.keys())[0]]: r for r in csv.DictReader(open(f"{AD}/atlas_secretory_metadata.csv"))}
# group per barcode (in barcode-file order)
bc_group = []
for bc in barcodes:
    m = meta.get(bc)
    bc_group.append((m['patient_id'], m['celltype_nmf']) if m else None)
cnt = collections.Counter(g for g in bc_group if g)
groups = [g for g, c in cnt.items() if c >= 50]        # >=50 cells
gidx = {g: i for i, g in enumerate(groups)}
ngrp = len(groups)
print(f"pseudobulk groups (>=50 cells): {ngrp}", flush=True)
# per-cell column index (or -1 to drop)
col_of_cell = np.array([gidx.get(g, -1) for g in bc_group], dtype=np.int64)

# accumulate genes x groups via chunked bincount on linear index
acc = np.zeros(ng * ngrp, dtype=np.float64)
import pandas as pd
reader = pd.read_csv(gzip.open(f"{AD}/atlas_secretory_counts.mtx.gz", "rt"),
                     sep=r"\s+", skiprows=3, header=None, names=["g","c","v"],
                     dtype={"g":np.int32,"c":np.int32,"v":np.float32}, chunksize=40_000_000)
tot = 0
for ci, ch in enumerate(reader):
    g = ch["g"].to_numpy() - 1
    c = ch["c"].to_numpy() - 1
    v = ch["v"].to_numpy().astype(np.float64)
    col = col_of_cell[c]
    keep = col >= 0
    lin = g[keep] * ngrp + col[keep]
    acc += np.bincount(lin, weights=v[keep], minlength=ng * ngrp)
    tot += len(ch)
    print(f"  chunk {ci} done; rows~{tot:,}", flush=True)

M = acc.reshape(ng, ngrp)                                # genes x groups
# write matrix (genes x groups) and group metadata
gnames = [f"{p}|{s}" for (p, s) in groups]
import numpy as np
hdr = "gene," + ",".join(gnames)
with open(f"{OUT}/atlas_pseudobulk_counts.csv", "w") as fh:
    fh.write(hdr + "\n")
    for i in range(ng):
        fh.write(genes[i] + "," + ",".join(f"{x:.4g}" for x in M[i]) + "\n")
with open(f"{OUT}/atlas_pseudobulk_meta.csv", "w") as fh:
    fh.write("group,patient_id,celltype_nmf,ncells\n")
    for (p, s) in groups:
        fh.write(f"{p}|{s},{p},{s},{cnt[(p,s)]}\n")
print("WROTE atlas_pseudobulk_counts.csv (%d x %d) + meta" % (ng, ngrp), flush=True)
print("DONE", flush=True)
