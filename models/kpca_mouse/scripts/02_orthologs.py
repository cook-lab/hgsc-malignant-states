#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202608_kpca_secb_mouse/scripts/02_orthologs.py, the script that ran. Kept:
# ortholog_map.csv (Fig 3j, Figure S10, steps 31 and 32) and 02_coverage_summary.csv (checked by the Figure S10
# script). Removed: the per-gene coverage table, which nothing reads. The untrimmed script is in the lab project
# folder.

"""
Step 02 — Human/mouse ortholog map + signature coverage
========================================================
KPCA.B mouse SecA/SecB subanalysis  ·  Cook Lab  ·  2026

Downloads (once, cached) the MGI human-mouse homology report and builds a
human_symbol -> mouse_symbol map retaining ortholog class. Falls back to symbol
case-conversion where MGI has no entry, and reports how often that fires.

The coverage table distinguishes three failure modes, which are NOT the same
thing and get conflated if you only count hits:
    mapped        — ortholog exists and is in the mouse matrix
    off_panel     — ortholog exists but is absent from the Flex probe set
    no_ortholog   — no mouse ortholog at all

Input:   MGI HOM_MouseHumanSequence.rpt (downloaded)
         output/01_kpca_qc/kpca_epithelial.h5ad  (for the observed gene space)
         human 11e gene classification, 11d loadings
Output:  data/external/HOM_MouseHumanSequence.rpt (+ README.md provenance)
         output/02_orthologs/ortholog_map.csv
         output/02_orthologs/02_signature_coverage.csv
         output/02_orthologs/02_coverage_summary.csv

Usage:
    python scripts/02_orthologs.py
    python scripts/02_orthologs.py --force-download
"""

import argparse
import hashlib
import os
import sys
from datetime import datetime, timezone
from importlib import import_module

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
S = import_module("00_setup")

OUT_DIR = os.path.join(S.OUTPUT, "02_orthologs")
EPI_H5AD = os.path.join(S.OUTPUT, "01_kpca_qc", "kpca_epithelial.h5ad")


def download_mgi(force=False):
    """Fetch and cache the MGI homology report, logging provenance."""
    if os.path.exists(S.MGI_RPT) and not force:
        print(f"[1] Using cached MGI report: {S.MGI_RPT}")
        return S.MGI_RPT

    import urllib.request
    os.makedirs(S.EXTERNAL, exist_ok=True)
    print(f"[1] Downloading {S.MGI_URL}")
    urllib.request.urlretrieve(S.MGI_URL, S.MGI_RPT)

    size = os.path.getsize(S.MGI_RPT)
    sha = hashlib.sha256(open(S.MGI_RPT, "rb").read()).hexdigest()
    print(f"    {size/1e6:.1f} MB  sha256={sha[:16]}...")

    with open(os.path.join(S.EXTERNAL, "README.md"), "a") as fh:
        fh.write(
            f"\n## HOM_MouseHumanSequence.rpt\n"
            f"- Source: {S.MGI_URL}\n"
            f"- Downloaded: {datetime.now(timezone.utc).isoformat()}\n"
            f"- Size: {size} bytes\n- SHA256: {sha}\n"
            f"- Used by: scripts/02_orthologs.py\n")
    return S.MGI_RPT


def build_map(rpt):
    """Parse MGI report into human_symbol -> mouse_symbol with class."""
    print("[2] Building ortholog map")
    df = pd.read_csv(rpt, sep="\t", low_memory=False)
    cols = {c.lower().strip(): c for c in df.columns}

    def pick(*cands):
        for c in cands:
            if c in cols:
                return cols[c]
        raise KeyError(f"none of {cands} in {list(df.columns)[:12]}")

    c_org = pick("common organism name", "organism name")
    c_sym = pick("symbol")
    c_key = pick("db class key", "homologene id", "class key")
    print(f"    columns: organism={c_org!r} symbol={c_sym!r} class_key={c_key!r}")

    df = df[[c_org, c_sym, c_key]].dropna()
    df.columns = ["organism", "symbol", "class_key"]
    df["organism"] = df["organism"].str.lower()

    hum = df[df["organism"].str.contains("human")]
    mou = df[df["organism"].str.contains("mouse")]
    print(f"    human rows {len(hum):,} | mouse rows {len(mou):,}")

    # Count members per homology class per species to label 1:1 vs 1:many
    nh = hum.groupby("class_key")["symbol"].nunique()
    nm = mou.groupby("class_key")["symbol"].nunique()

    merged = hum.merge(mou, on="class_key", suffixes=("_human", "_mouse"))
    merged["n_human"] = merged["class_key"].map(nh)
    merged["n_mouse"] = merged["class_key"].map(nm)
    merged["ortholog_class"] = np.where(
        (merged["n_human"] == 1) & (merged["n_mouse"] == 1), "1:1",
        np.where(merged["n_human"] == 1, "1:many",
                 np.where(merged["n_mouse"] == 1, "many:1", "many:many")))

    # Prefer 1:1; among ties keep the alphabetically first mouse symbol so the
    # map is deterministic across runs.
    order = {"1:1": 0, "many:1": 1, "1:many": 2, "many:many": 3}
    merged["_rank"] = merged["ortholog_class"].map(order)
    merged = merged.sort_values(["symbol_human", "_rank", "symbol_mouse"])
    best = merged.drop_duplicates("symbol_human")

    out = best[["symbol_human", "symbol_mouse", "ortholog_class",
                "n_human", "n_mouse", "class_key"]].copy()
    out.columns = ["human_symbol", "mouse_symbol", "ortholog_class",
                   "n_human_in_class", "n_mouse_in_class", "class_key"]
    out["source"] = "mgi"
    print(f"    unique human symbols mapped: {len(out):,}")
    print(out["ortholog_class"].value_counts().to_string())
    return out


def coverage(ortho, mouse_genes):
    """Coverage of every human signature against the mouse gene space."""
    print("[3] Signature coverage")
    gc = pd.read_csv(S.H_GENECLASS)
    loadings = pd.read_csv(S.H_LOADINGS, index_col=0)

    sigs = {
        "SecA_canonical": S.SECA_GENES,
        "SecB_canonical": S.SECB_GENES,
        "SecA_specific_177": gc.loc[gc["class"] == "SecA-specific", "gene"].tolist(),
        "SecB_specific_177": gc.loc[gc["class"] == "SecB-specific", "gene"].tolist(),
        "Shared_124": gc.loc[gc["class"] == "Shared", "gene"].tolist(),
        "NMF_HVG_3007": loadings.columns.tolist(),
    }

    rows, summary = [], []
    for name, genes in sigs.items():
        _, rep = S.map_genes(genes, ortho, restrict_to=mouse_genes)
        rep["signature"] = name
        rows.append(rep)
        vc = rep["status"].value_counts()
        n = len(rep)
        summary.append({
            "signature": name, "n_genes": n,
            "mapped": int(vc.get("mapped", 0)),
            "off_panel": int(vc.get("off_panel", 0)),
            "no_ortholog": int(vc.get("no_ortholog", 0)),
            "pct_mapped": round(100 * vc.get("mapped", 0) / n, 1),
            "pct_via_case_conversion": round(
                100 * ((rep["source"] == "case_conversion") &
                       (rep["status"] == "mapped")).sum() / n, 1),
        })
        print(f"    {name:20s} {vc.get('mapped',0):5d}/{n:5d} mapped "
              f"({100*vc.get('mapped',0)/n:5.1f}%)  "
              f"off_panel={vc.get('off_panel',0):4d}  "
              f"no_ortholog={vc.get('no_ortholog',0):4d}")

    det = pd.concat(rows, ignore_index=True)
    sm = pd.DataFrame(summary)
    return det, sm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force-download", action="store_true")
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    print("=" * 74)
    print("Step 02 — Ortholog map and signature coverage")
    print("=" * 74)

    rpt = download_mgi(args.force_download)
    ortho = build_map(rpt)
    ortho.to_csv(S.ORTHO_MAP, index=False)
    print(f"    wrote {S.ORTHO_MAP}")

    E = sc.read_h5ad(EPI_H5AD)
    mouse_genes = set(E.var_names)
    print(f"    mouse gene space: {len(mouse_genes):,} genes "
          f"({E.shape[0]:,} epithelial cells)")

    det, sm = coverage(ortho, mouse_genes)
    sm.to_csv(os.path.join(OUT_DIR, "02_coverage_summary.csv"), index=False)

    print("\n[4] Canonical marker audit (must be 100% for the analysis to hold)")
    canon = det[det["signature"].isin(["SecA_canonical", "SecB_canonical"])]
    print(canon[["signature", "human_symbol", "mouse_symbol", "source",
                 "status"]].to_string(index=False))
    unmapped = canon[canon["status"] != "mapped"]
    if len(unmapped):
        print(f"\n    WARNING: {len(unmapped)} canonical markers unmapped:")
        print(unmapped[["human_symbol", "mouse_symbol", "status"]].to_string(index=False))
    else:
        print("\n    All canonical SecA/SecB markers mapped and present. OK.")

    print("\nStep 02 complete.")


if __name__ == "__main__":
    main()
