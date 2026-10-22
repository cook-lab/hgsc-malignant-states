#!/usr/bin/env python3
"""
Supplementary Table 7 -- Xenium gene panel composition, panel design rationale and cross-platform probe QC
==========================================================================================================

One row per biological target of the custom Xenium Prime panel (hOvary_477g, design VW8464; 477 genes; the 20
negative-control probes are not listed). Regenerates the table end-to-end from the source files below; nothing is
typed in by hand except the column dictionary.

INPUTS (read-only)
    data/HGSC xenium /<14 runs>/gene_panel.json                         panel definition (gene, Ensembl ID) -- all 14
                                                                        runs are checked to carry the identical 477 genes
    data/HGSC xenium /<run>/cell_feature_matrix.h5                      raw counts (all-cell Xenium detection rates)
    $OUTPUT_ROOT/singler_recompute/<sample>_singler_cells.csv.gz   the analysed (post-QC) cells, per sample
    2026_final_xenium_analysis/output/05_probe_qc/probe_qc_full.csv     cross-platform QC metrics (frozen Feb-2026 run)
    2026_final_xenium_analysis/output/05_probe_qc/flag_qc.csv           the six per-gene QC flags
    2026_final_xenium_analysis/output/05_probe_qc/gene_exclusion_decisions.csv   4-tier decision + reason
    $OUTPUT_ROOT/05_probe_qc/genes_exclude_singler.txt       exclusion list actually read by the SingleR run
    2026_final_atlas/output/xenium_celltype/xenium_celltype_downsampled.h5ad   SingleR reference (var names only)
    $OUTPUT_ROOT/_logs/singler_parallel.log                   SingleR run log (gene counts cross-check)
    2026_final_xenium_analysis/metadata/panel_markers.csv               level-1 / level-2 marker assignments
    2026_final_xenium_analysis/data/Updated Xenium_2024_10.xlsx         panel-design worksheet (sheet "Genelist")
    $OUTPUT_ROOT/9b_scoring/pathway_gene_sets_v2.csv          37 UCell pathway modules
    shared/signatures.yml                                               7-gene SecA / SecB polarization signatures

OUTPUTS
    $OUTPUT_ROOT/final_panels/supplementary/tables/SuppTable7_xenium_gene_panel.xlsx  (README, Gene_panel, Summary)
    + one CSV per data sheet

MANUSCRIPT ITEM(S)
    Table S7; text 477 panel genes, 31 probe-QC failures, 28 excluded from SingleR, 429 SingleR genes, 20 Xenium-only.

RUNTIME: ~3-5 min (streams the sparse raw count matrices of the 14 analysed sections/TMAs from the share).
"""
import json
import os
import re
import sys

import h5py
import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from supp_xlsx import write_table  # noqa: E402

# ------------------------------------------------------------------ paths
SHARE = os.environ.get("DATA_ROOT", "/path/to/deposit")
RAW = f"{SHARE}/202609_final_hgsc_spatial/data/HGSC xenium "
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
QC = f"{SHARE}/2026_final_xenium_analysis/output/05_probe_qc"
RERUN = os.environ.get("OUTPUT_ROOT", "output")
SG = f"{RERUN}/singler_recompute"
EXCL_TXT = f"{RERUN}/05_probe_qc/genes_exclude_singler.txt"
SINGLER_LOG = f"{RERUN}/_logs/singler_parallel.log"
REF_H5AD = f"{SHARE}/2026_final_atlas/output/xenium_celltype/xenium_celltype_downsampled.h5ad"
MARKERS = f"{SHARE}/2026_final_xenium_analysis/metadata/panel_markers.csv"
DESIGN = f"{SHARE}/2026_final_xenium_analysis/data/Updated Xenium_2024_10.xlsx"
PATHWAYS = f"{RERUN}/9b_scoring/pathway_gene_sets_v2.csv"
SIGS = f"{REPO}/shared/signatures.yml"

WT8 = ["WT01", "WT02", "WT03", "WT04", "WT05", "WT06", "WT07", "WT08"]
# analysed sample -> raw Xenium run(s)
RAW_RUN = {s: [s] for s in WT8}
RAW_RUN.update({"WT01": ["WT05_WT01"], "WT05": ["WT05_WT01"],   # dual-tissue slide
                "tma": ["CHUM_TMA_1", "CHUM_TMA_2"],                             # CHUM TMA = TMA A (2 slides)
                "COEUR_TMA_G": ["COEUR_TMA_G"], "COEUR_TMA_H": ["COEUR_TMA_H"]})  # COEUR G+H = TMA B


def log(*a):
    print(*a, flush=True)


# ------------------------------------------------------------------ 1. panel definition (all 14 runs)
runs = sorted(d for d in os.listdir(RAW) if os.path.isdir(os.path.join(RAW, d)) and not d.startswith("."))
panel_sets, panel = {}, None
for r in runs:
    js = [os.path.join(dp, f) for dp, _, fs in os.walk(os.path.join(RAW, r)) for f in fs
          if f == "gene_panel.json"]
    if not js:
        continue
    j = json.load(open(js[0]))
    tg = j["payload"]["targets"]
    genes = [(t["type"]["data"]["name"], t["type"]["data"]["id"], t.get("info", {}).get("gene_coverage"))
             for t in tg if t["type"]["descriptor"] == "gene"]
    n_neg = sum(t["type"]["descriptor"] != "gene" for t in tg)
    panel_sets[r] = frozenset(g for g, _, _ in genes)
    if panel is None:
        panel = pd.DataFrame(genes, columns=["gene", "ensembl_id", "gene_coverage"])
        pinfo = j["payload"]["panel"]
        PANEL_NAME = f'{pinfo["identity"]["name"]} (design {pinfo["identity"]["design_id"]}, v{pinfo["identity"]["version"]})'
        N_NEG = n_neg
assert len(set(panel_sets.values())) == 1, "gene panels differ between runs"
N_PANEL = len(panel)
log(f"panel: {PANEL_NAME}; {N_PANEL} genes + {N_NEG} negative controls; identical across {len(panel_sets)} runs")

# ------------------------------------------------------------------ 2. probe QC (frozen)
full = pd.read_csv(f"{QC}/probe_qc_full.csv")
flags = pd.read_csv(f"{QC}/flag_qc.csv")
dec = pd.read_csv(f"{QC}/gene_exclusion_decisions.csv")
assert set(dec.gene) == set(panel.gene), "QC decisions do not cover the panel"
excl = [g.strip() for g in open(EXCL_TXT) if g.strip()]
dec_excl = set(dec.gene[dec.tier.isin(["EXCLUDE", "HARD_EXCLUDE"])])
assert set(excl) == dec_excl, "SingleR exclusion list != QC decision EXCLUDE/HARD_EXCLUDE"

# ------------------------------------------------------------------ 3. SingleR gene space (reference x panel - excluded)
with h5py.File(REF_H5AD, "r") as f:
    ref_genes = set(f["var"][f["var"].attrs["_index"]][:].astype(str))
in_ref = panel.gene.isin(ref_genes)
singler_genes = set(panel.gene[in_ref]) - set(excl)
N_SHARED, N_SINGLER = int(in_ref.sum()), len(singler_genes)
xen_only = set(panel.gene[~in_ref])
assert xen_only == set(dec.gene[dec.qc_flag == "XENIUM_ONLY"]), "Xenium-only set disagrees with QC"
logtxt = open(SINGLER_LOG).read()
m_sh = re.search(r"Shared genes \(panel ∩ reference\): (\d+)", logtxt)
m_sr = re.search(r"After removing (\d+) excluded genes: (\d+) genes for SingleR", logtxt)
assert m_sh and int(m_sh.group(1)) == N_SHARED and m_sr and int(m_sr.group(2)) == N_SINGLER, "log mismatch"
log(f"SingleR: shared {N_SHARED}, excluded {len(excl)}, used {N_SINGLER} (log agrees)")

# ------------------------------------------------------------------ 4. annotations
mk = pd.read_csv(MARKERS)
mk1 = mk[mk.level == "level1"].groupby("gene").cell_type.apply(lambda s: "; ".join(sorted(set(s))))
mk2 = mk[mk.level == "level2"].groupby("gene").cell_type.apply(lambda s: "; ".join(sorted(set(s))))

des = pd.read_excel(DESIGN, sheet_name="Genelist")
des = des[des["Include in panel"] == True].copy()  # noqa: E712
des["gene"] = des["Genes considered"].astype(str).str.strip()
des = des.set_index("gene")
DES_ONLY = sorted(set(des.index) - set(panel.gene))
PANEL_NO_DESIGN = sorted(set(panel.gene) - set(des.index))

pw = pd.read_csv(PATHWAYS)
pw_by_gene = pw.groupby("gene").pathway.apply(lambda s: "; ".join(sorted(set(s))))
N_MODULES = pw.pathway.nunique()
pw_not_panel = sorted(set(pw.gene) - set(panel.gene))
sigs = yaml.safe_load(open(SIGS))
sig_of = {g: k for k in ("SecA", "SecB") for g in sigs[k]}

# ------------------------------------------------------------------ 5. all-cell Xenium detection rate (analysed cells)
def analysed_barcodes(sample):
    b = pd.read_csv(f"{SG}/{sample}_singler_cells.csv.gz", usecols=["barcode"]).barcode.astype(str)
    return b


gene_idx = None
det = {"WT": np.zeros(N_PANEL, np.int64), "TMA": np.zeros(N_PANEL, np.int64)}
ncell = {"WT": 0, "TMA": 0}
bc_cache = {}
for sample in WT8 + ["tma", "COEUR_TMA_G", "COEUR_TMA_H"]:
    arm = "WT" if sample in WT8 else "TMA"
    want = analysed_barcodes(sample)
    matched = 0
    for run in RAW_RUN[sample]:
        with h5py.File(f"{RAW}/{run}/cell_feature_matrix.h5", "r") as f:
            m = f["matrix"]
            names = m["features/name"][:].astype(str)
            ftype = m["features/feature_type"][:].astype(str)
            gi = {n: i for i, (n, t) in enumerate(zip(names, ftype)) if t == "Gene Expression"}
            order = np.array([gi[g] for g in panel.gene])          # raw feature index of each panel gene
            bcs = pd.Index(m["barcodes"][:].astype(str))
            if sample == "tma":   # CHUM slide-2 barcodes that collide with slide 1 carry a '_TMA2' suffix
                w = want[want.str.endswith("_TMA2")].str.replace("_TMA2", "", regex=False) if run == "CHUM_TMA_2" \
                    else want[~want.str.endswith("_TMA2")]
                if run == "CHUM_TMA_2":
                    w = pd.concat([w, want[~want.str.endswith("_TMA2") & ~want.isin(bc_cache["CHUM_TMA_1"])]])
            else:
                w = want
            sel = bcs.isin(set(w))
            if run == "CHUM_TMA_1":
                bc_cache["CHUM_TMA_1"] = set(bcs[sel])
            matched += int(sel.sum())
            indptr = m["indptr"][:]
            keep = np.repeat(sel, np.diff(indptr))
            idx = m["indices"][:][keep]
            cnt = np.bincount(idx, minlength=len(names))
            det[arm] += cnt[order]
    assert matched == len(want), f"{sample}: matched {matched} of {len(want)} analysed barcodes"
    ncell[arm] += matched
    log(f"  detection: {sample:12s} {matched:>9,} cells")

# ------------------------------------------------------------------ 6. assemble
t = panel[["gene", "ensembl_id"]].copy()
t = t.merge(dec[["gene", "qc_flag", "n_flags", "tier", "decision", "reason", "marker_for"]], on="gene", how="left")
t = t.merge(full[["gene", "rho_corr_of_corr", "top20_jaccard", "sc_detect_rate", "xe_detect_rate",
                  "log2_mean_ratio", "var_ratio", "detect_rate_diff", "mean_nonzero_diff", "overdisp_ratio"]],
            on="gene", how="left")   # the four flag metrics after log2_mean_ratio: added 2026-10-20 (the caption lists them)
FLAG_NAMES = {"flag_corr": "co-expression (rho<0.3 or Jaccard<0.05)", "flag_mean": "mean shift (|log2 ratio|>2.5)",
              "flag_var": "variance ratio (<0.2 or >5)", "flag_detect": "detection shift (|diff|>0.30)",
              "flag_nonzero": "non-zero mean shift (|diff|>0.75)", "flag_overdisp": "overdispersion ratio (<0.5 or >2)"}
fl = flags.set_index("gene")
t["QC flags triggered"] = t.gene.map(lambda g: "; ".join(v for k, v in FLAG_NAMES.items() if g in fl.index and bool(fl.at[g, k])) if g in fl.index else "")
t["xenium_only"] = ~t.gene.isin(ref_genes)
t["in_singler"] = t.gene.isin(singler_genes)
t["excluded_from_singler"] = t.gene.isin(excl)
t["in_ucell_polar_space"] = t.gene.isin(ref_genes)   # UCell SecA/SecB ranked on panel ∩ atlas (457)
STATUS = {"INCLUDE": "Pass", "INCLUDE_MONITOR": "Pass (monitored, WARN)", "KEEP_CAUTION": "Pass (kept with caution)",
          "EXCLUDE": "Excluded from annotation", "HARD_EXCLUDE": "Excluded from annotation (hard)"}
t["qc_status"] = t.tier.map(STATUS)
t.loc[t.xenium_only, "qc_status"] = "Not assessable (absent from atlas reference)"
t["polar_sig"] = t.gene.map(sig_of).fillna("")
t["modules"] = t.gene.map(pw_by_gene).fillna("")
t["n_modules"] = t.gene.map(pw.groupby("gene").size()).fillna(0).astype(int)
t["mk1"] = t.gene.map(mk1).fillna("")
t["mk2"] = t.gene.map(mk2).fillna("")
t["design_reason"] = t.gene.map(des["Inclusion reason"]).fillna("")
t["design_for"] = t.gene.map(des["Considered for"]).fillna("")
t["design_expr"] = t.gene.map(des["Dominant Cell Expressor (from scRNAseq)"]).fillna("")
t["xe_det_wt"] = det["WT"] / ncell["WT"]
t["xe_det_tma"] = det["TMA"] / ncell["TMA"]
t["exclusion_reason"] = np.where(t.excluded_from_singler, t.reason,
                                 np.where(t.xenium_only, "Absent from scRNA-seq reference (cannot enter SingleR)", ""))

out = pd.DataFrame({
    "Gene": t.gene,
    "Ensembl ID": t.ensembl_id,
    "Retained for expression / pathway analysis": "Yes",
    "Used for SingleR annotation": np.where(t.in_singler, "Yes", "No"),
    "Excluded from SingleR by probe QC": np.where(t.excluded_from_singler, "Yes", "No"),
    "Xenium-only (absent from scRNA-seq atlas reference)": np.where(t.xenium_only, "Yes", "No"),
    "QC status": t.qc_status,
    "QC class (PASS/WARN/FAIL)": t.qc_flag.replace({"XENIUM_ONLY": "NA"}),
    "Number of QC flags (of 6)": t.n_flags.where(~t.xenium_only),
    "QC flags triggered": t["QC flags triggered"],
    "Exclusion / caution reason": np.where(t.tier.eq("KEEP_CAUTION"), t.reason, t.exclusion_reason),
    "Correlation of correlations (Spearman rho)": t.rho_corr_of_corr.round(3),
    "Top-20 partner Jaccard": t.top20_jaccard.round(3),
    "Log2 mean ratio (Xenium / scRNA-seq)": t.log2_mean_ratio.round(3),
    "Variance ratio (Xenium / scRNA-seq)": t.var_ratio.round(3),
    "Detection-rate difference (Xenium - scRNA-seq)": t.detect_rate_diff.round(3),
    "Non-zero-mean difference (Xenium - scRNA-seq)": t.mean_nonzero_diff.round(3),
    "Overdispersion ratio (Xenium / scRNA-seq)": t.overdisp_ratio.round(3),
    "scRNA-seq detection rate (reference, 16,000 cells)": t.sc_detect_rate.round(3),
    "Xenium detection rate (QC run, 200k-cell stratified subsample)": t.xe_detect_rate.round(3),
    f"Xenium detection rate, whole tissues (all {ncell['WT']:,} analysed cells)": t.xe_det_wt.round(3),
    f"Xenium detection rate, TMA cores (all {ncell['TMA']:,} analysed cells)": t.xe_det_tma.round(3),
    "Cell-type marker, level 1 (panel design)": t.mk1,
    "Cell-type marker, level 2 (panel design)": t.mk2,
    "Panel design: inclusion reason": t.design_reason,
    "Panel design: considered for": t.design_for,
    "Panel design: dominant expressor (scRNA-seq)": t.design_expr,
    "SecA/SecB polarization signature (7-gene)": t.polar_sig,
    "Number of UCell pathway modules": t.n_modules,
    "UCell pathway module(s)": t.modules,
})
out = out.sort_values(["Used for SingleR annotation", "Gene"], ascending=[False, True]).reset_index(drop=True)

# ------------------------------------------------------------------ 7. summary
vc = dec.tier.value_counts()
summ = pd.DataFrame([
    ("Biological gene targets on panel (gene_panel.json)", N_PANEL, "identical in all 14 raw Xenium runs"),
    ("Negative-control probes (not listed)", N_NEG, ""),
    ("Genes present in scRNA-seq SingleR reference", N_SHARED, "panel ∩ atlas reference"),
    ("Xenium-only genes (absent from reference)", len(xen_only), "cannot be assessed cross-platform; cannot enter SingleR"),
    ("QC PASS (0-1 flags)", int((dec.qc_flag == "PASS").sum()), ""),
    ("QC WARN (2 flags) - kept, monitored", int((dec.qc_flag == "WARN").sum()), ""),
    ("QC FAIL (>=3 flags)", int((dec.qc_flag == "FAIL").sum()), "of which 3 kept with caution (MYL9, TPM2, PF4)"),
    ("  FAIL, excluded (tier EXCLUDE)", int(vc.get("EXCLUDE", 0)), ""),
    ("  FAIL, hard-excluded (rho<0.2)", int(vc.get("HARD_EXCLUDE", 0)), ""),
    ("  FAIL, kept with caution (critical marker)", int(vc.get("KEEP_CAUTION", 0)), ""),
    ("Genes excluded from SingleR (probe QC)", len(excl), "genes_exclude_singler.txt"),
    ("Genes used for SingleR", N_SINGLER, f"{N_SHARED} shared - {len(excl)} excluded (SingleR log agrees)"),
    ("Panel genes minus QC exclusions", N_PANEL - len(excl), "includes the 20 Xenium-only genes, which SingleR cannot use"),
    ("Genes retained for expression / pathway analysis", N_PANEL, "all genes; pathway UCell maxRank = 477"),
    ("Genes in UCell SecA/SecB polarization ranking space", N_SHARED, "panel ∩ atlas"),
    ("UCell pathway modules", N_MODULES, f"{pw.gene.nunique()} distinct genes, {len(pw)} gene-module memberships"),
    ("Panel genes in >=1 UCell module", int((out['Number of UCell pathway modules'] > 0).sum()), ""),
    ("Analysed cells, whole tissues (detection-rate denominator)", ncell["WT"], "8 core whole tissues"),
    ("Analysed cells, TMA (detection-rate denominator)", ncell["TMA"], "CHUM (incl. FTE + off-core) + COEUR G/H, all analysed cells"),
], columns=["Quantity", "Value", "Note"])

readme = [
    f"Panel: custom 10x Genomics Xenium Prime panel {PANEL_NAME}; {N_PANEL} biological gene targets plus {N_NEG} "
    "negative-control probes. The gene list was read from gene_panel.json of every raw run (14 runs; identical).",
    f"Probe QC (cross-platform, scRNA-seq atlas reference vs Xenium): six flags per gene; >=3 flags = FAIL, 2 = WARN, "
    "0-1 = PASS. FAIL genes were excluded from SingleR unless they were a marker for a cell type with <5 reliable "
    "markers (kept with caution: MYL9, TPM2, PF4); FAIL genes with rho<0.2 were hard-excluded. "
    f"Result: {len(excl)} genes excluded from SingleR annotation, all retained for expression/pathway analysis.",
    f"SingleR gene space: panel ∩ reference ({N_SHARED}) minus the {len(excl)} excluded genes = {N_SINGLER} genes "
    f"(confirmed in the SingleR run log). The {len(xen_only)} Xenium-only genes are absent from the reference and "
    f"so never enter SingleR; {N_PANEL}-{len(excl)} = {N_PANEL - len(excl)} is the number of panel genes that are not "
    "QC-excluded, not the SingleR gene count.",
    "QC metrics are those of the frozen probe-QC run that produced the exclusion list (2026_final_xenium_analysis/"
    "output/05_probe_qc; Xenium side = 200,000-cell stratified subsample, scRNA-seq side = the 16,000-cell SingleR "
    "reference). They are reported unchanged because they define the decision. The two 'all analysed cells' "
    "detection-rate columns are recomputed here from the raw count matrices over every analysed cell (no subsampling).",
    "Flag definitions: co-expression = Spearman correlation of a gene's gene-gene correlation profile between "
    "platforms (rho) < 0.3 or top-20 co-expression partner Jaccard < 0.05; mean = |log2 mean ratio| > 2.5; variance "
    "ratio < 0.2 or > 5; detection-rate difference > 0.30; non-zero-mean difference > 0.75; overdispersion ratio < 0.5 "
    "or > 2.",
    "Panel design columns come from the design worksheet 'Updated Xenium_2024_10.xlsx' (sheet Genelist, rows with "
    f"'Include in panel' = TRUE); blank = no design record. Design-sheet genes not on the final panel: {', '.join(DES_ONLY)}; "
    f"panel genes without a design row: {', '.join(PANEL_NO_DESIGN) or 'none'}.",
    "Cell-type marker columns: panel_markers.csv (level 1 = broad types used for marker-coverage risk in QC; "
    "level 2 = subtypes).",
    f"UCell modules: the {N_MODULES} pathway modules of pathway_gene_sets_v2.csv (gene lists in Supplementary Table 9); SecA/SecB = the 7-gene polarization signatures (shared/signatures.yml). "
    + (f"Module genes not on the panel: {', '.join(pw_not_panel)}." if pw_not_panel else "All module genes are on the panel."),
    "Sheets: Gene_panel (one row per gene, sorted SingleR-used first, then alphabetical); Summary (panel-level counts).",
    "Built by tables/07_tableS7_xenium_gene_panel.py.",
]
write_table(7, "xenium_gene_panel",
            title="Supplementary Table 7. Xenium gene panel composition, design rationale and cross-platform probe quality control",
            readme=readme, sheets=[("Gene_panel", out), ("Summary", summ)])
log(summ.to_string(index=False))
