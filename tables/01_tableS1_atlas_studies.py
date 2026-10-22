#!/usr/bin/env python3
"""
Supplementary Table 1 -- per-study / per-sample metadata of the HGSC scRNA-seq atlas
===================================================================================

All counts come from ALL 1,980,703 cells of the final (post-QC, post-exclusion) atlas; no subsampling.

INPUTS (read-only)
    202601_final_atlas/output/supplemental/T1_atlas_metadata.csv        cell-level metadata (post-correction, Jul 30)
    202601_final_atlas/hgsc_atlas_final.h5ad  (obs only, via h5py)        cross-check of totals + library names
    202601_final_atlas/pre-processing and integration/atlas_metadata.xlsx  curated per-library metadata (cross-check)
    202601_final_atlas/pre-processing and integration/raw/<study>/          raw file names -> GSM accessions
    output_root/_prep_caches/nath2021_cell_correction.parquet             Nath 2021 per-cell sample mapping fix
        (figures/_prep/nath2021_fix.py, 2026-10-20): sample, patient, treatment, TP53 and BRCA of each Nath cell
        from its GEO patient and time point. Run with OUTPUT_ROOT set to the atlas output root.

OUTPUTS
    $OUTPUT_ROOT/final_panels/supplementary/tables/SuppTable1_atlas_study_metadata.xlsx
        sheets: README, Studies, Samples, Cells_by_lineage
    $OUTPUT_ROOT/final_panels/supplementary/tables/SuppTable1_atlas_study_metadata_<sheet>.csv (one per sheet)

MANUSCRIPT ITEM(S)
    Table S1; Results / Methods atlas totals (13 studies, 1,980,703 cells, 158 patients, 362 samples).

PUBLICATION FIELDS
    Hard-coded in STUDY_INFO below. Every value was verified on 2026-09-28 against NCBI GEO series/sample records
    (acc.cgi text view), PubMed esummary, PMC full text (Methods / Data availability), the EGA metadata export on disk
    (hornburg_2021/EGAD00001006974-metadata), or the manuscript reference list. The source of each is recorded in the
    "Verification source" column. Fields that could not be verified are left blank and listed in README "To verify".

RUNTIME: ~1 min (CSV read + h5py obs read).
"""
import os
import re
import glob

import sys

import h5py
import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from figures._prep.nath2021_fix import correction as nath2021_correction  # noqa: E402

ATLAS_DIR = os.path.join(os.environ.get("DATA_ROOT", "/path/to/deposit"), "2026_final_atlas")
T1_CSV = f"{ATLAS_DIR}/output/supplemental/T1_atlas_metadata.csv"
H5AD = f"{ATLAS_DIR}/hgsc_atlas_final.h5ad"
XLSX = f"{ATLAS_DIR}/pre-processing and integration/atlas_metadata.xlsx"
RAW = f"{ATLAS_DIR}/pre-processing and integration/raw"
OUT_DIR = os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "final_panels/supplementary/tables")
OUT_XLSX = f"{OUT_DIR}/SuppTable1_atlas_study_metadata.xlsx"

LINEAGES = ["Epithelial", "Mesothelial", "Fibroblast", "Smooth muscle", "Pericyte", "Endothelial",
            "T/NK cell", "B cell", "Plasma cell", "Macrophage", "DC", "Neutrophil", "Mast cell"]
TREAT_NAIVE = {"pre-treatment"}
TREAT_CHEMO = {"post-chemotherapy"}
TREAT_OTHER = {"post-niraparib", "post-chemotherapy_niraparib", "post-chemotherapy_olaparib",
               "post-chemotherapy_pembro", "post-chemotherapy_niraparib_pembro"}

# --------------------------------------------------------------------------------------------------------------
# Verified publication / data metadata (see module docstring). "" = could not be verified -> README "To verify".
# ms_ref = reference number in the manuscript reference list (written as "Manuscript ref #").
# ms_ref_check = an internal citation check; it is not written to the table.
# --------------------------------------------------------------------------------------------------------------
STUDY_INFO = {
    "denisenko_2024": dict(
        label="Denisenko et al., 2024",
        title="Spatial transcriptomics reveals discrete tumour microenvironments and autocrine loops within "
              "ovarian cancer subclones",
        journal="Nat Commun 15, 2860 (2024)", doi="10.1038/s41467-024-47271-y", pmid="38570491",
        accession="GEO GSE211956 (scRNA-seq GSM6506105-GSM6506109)",
        raw_accession="SRA via BioProject PRJNA873051",
        platform="10x Chromium Single Cell 3' v2",
        processing="Whole cell; dissociated ovarian tumour from interval debulking surgery (post-NACT)",
        ms_ref="20", ms_ref_check="Matches data source",
        source="GEO GSE211956 / GSM6506105 (extract protocol: '10x Chromium 3' v2'); PubMed 38570491; "
               "PMC10991508 Methods"),
    "geistlinger_2020": dict(
        label="Geistlinger et al., 2020",
        title="Multiomic Analysis of Subtype Evolution and Heterogeneity in High-Grade Serous Ovarian Carcinoma",
        journal="Cancer Res 80(20), 4335-4345 (2020)", doi="10.1158/0008-5472.CAN-20-0521", pmid="32747365",
        accession="GEO GSE154600 (GSM4675273-GSM4675277)", raw_accession="SRA SRP272287",
        platform="10x Chromium Single Cell 3' v2",
        processing="Whole cell (dissociated omental tumour)",
        ms_ref="6", ms_ref_check="Matches data source",
        source="GEO GSE154600 / GSM4675273 (Chromium Single Cell 3' v2 kit); PubMed 32747365"),
    "hornburg_2021": dict(
        label="Hornburg et al., 2021",
        title="Single-cell dissection of cellular components and interactions shaping the tumor immune "
              "phenotypes in ovarian cancer",
        journal="Cancer Cell 39(7), 928-944.e6 (2021)", doi="10.1016/j.ccell.2021.04.004", pmid="33961783",
        accession="EGA EGAS00001004935 (dataset EGAD00001006974, count matrices)",
        raw_accession="EGA EGAS00001004935 (controlled access)",
        platform="10x Chromium Single Cell, v2 chemistry (3'/5' not stated in EGA record)",
        processing="Whole cell; live cells FACS-sorted into tumour, immune and stromal compartments, each "
                   "sequenced separately (stromal+immune pooled for some desert tumours)",
        ms_ref="9", ms_ref_check="Matches data source",
        source="EGA metadata on disk (studies.csv EGAS00001004935; dataset.csv EGAD00001006974: '10X Genomics "
               "Chromium Single Cell platform (v2 Chemistry)'); PubMed 33961783"),
    "loret_2022": dict(
        label="Loret et al., 2022",
        title="Distinct Transcriptional Programs in Ascitic and Solid Cancer Cells Induce Different Responses to "
              "Chemotherapy in High-Grade Serous Ovarian Cancer",
        journal="Mol Cancer Res 20(10), 1532-1547 (2022)", doi="10.1158/1541-7786.MCR-21-0565", pmid="35749080",
        accession="GEO GSE201047 (GSM6049610-GSM6049631)", raw_accession="SRA via BioProject PRJNA828147",
        platform="10x Genomics, v2 chemistry (3'/5' not stated in GEO)",
        processing="",
        ms_ref="12", ms_ref_check="Matches data source",
        source="GEO GSE201047 / GSM6049610 (characteristic 'platform: 10x Genomics - v2'); PubMed 35749080"),
    "luo_2024": dict(
        label="Luo et al., 2024",
        title="Neoadjuvant PARPi or chemotherapy in ovarian cancer informs targeting effector Treg cells for "
              "homologous-recombination-deficient tumors",
        journal="Cell 187(18), 4905-4925.e24 (2024)", doi="10.1016/j.cell.2024.06.013", pmid="38971151",
        accession="GEO GSE222556 (10x scRNA-seq SubSeries of SuperSeries GSE222557; GSM6925959-GSM6926025)",
        raw_accession="GSA-Human HRA007180 (controlled access)",
        platform="10x Chromium Next GEM Single Cell 3' v3.1 (GEO series design: '10x genomics 3' v3 chemistry')",
        processing="Whole cell, unsorted; gentleMACS tumour dissociation, RBC lysis and dead-cell removal",
        ms_ref="14",
        ms_ref_check="MISMATCH: ms ref 11 is Luo H, Wang K & Li B, Front Immunol 15:1420847 (2024; PMID 38975339), "
                     "which re-analysed GSE211956 + GSE147082. The atlas luo_2024 data (GSE222556) are from "
                     "Luo Y et al., Cell 2024 (PMID 38971151)",
        source="GEO GSE222556 / GSE222557 / GSM6925959 (Chromium Next GEM Single Cell 3' v3.1); PubMed 38971151; "
               "PMC11224428 (Front Immunol paper data source)"),
    "nath_2021": dict(
        label="Nath et al., 2021",
        title="Evolution of core archetypal phenotypes in progressive high grade serous ovarian cancer",
        journal="Nat Commun 12, 3039 (2021)", doi="10.1038/s41467-021-23171-3", pmid="34031395",
        accession="GEO GSE158722", raw_accession="dbGaP phs002294 (controlled access)",
        platform="10x Chromium Single Cell 3' v3 for all atlas-retained samples (the study's iCell8 samples, "
                 "GEO P01-P03, did not pass atlas QC)",
        processing="Malignant ascites/effusion; 10x libraries from whole cells or isolated nuclei (per-sample "
                   "split not stated)",
        ms_ref="13", ms_ref_check="Matches data source",
        source="GEO GSE158722 sample records (platform: 10X / iCell8 per GSM); PMC8144406 Methods "
               "('Chromium Single Cell 3' V3 Kit ... using whole cells or isolated nuclei'); PubMed 34031395"),
    "olalekan_2021": dict(
        label="Olalekan et al., 2021",
        title="Characterizing the tumor microenvironment of metastatic ovarian cancer by single-cell "
              "transcriptomics",
        journal="Cell Rep 35(8), 109165 (2021)", doi="10.1016/j.celrep.2021.109165", pmid="34038734",
        accession="GEO GSE147082 (GSM4416534-GSM4416539)", raw_accession="SRA SRP253041",
        platform="Drop-seq (NextSeq 500)",
        processing="Whole cell; omental metastases dissociated with collagenase IV / DNase I / hyaluronidase",
        ms_ref="16", ms_ref_check="Matches data source",
        source="GEO GSE147082 / GSM4416534 (extract protocol: 'Drop-seq protocol'); PubMed 34038734"),
    "olbrecht_2021": dict(
        label="Olbrecht et al., 2021",
        title="High-grade serous tubo-ovarian cancer refined with single-cell RNA sequencing: specific cell "
              "subtypes influence survival and determine molecular subtype classification",
        journal="Genome Med 13, 111 (2021)", doi="10.1186/s13073-021-00922-x", pmid="34238352",
        accession="Processed counts + metadata: http://blueprint.lambrechtslab.org (files 2093/2095-Olbrecht)",
        raw_accession="EGA EGAS00001004987 (controlled access)",
        platform="10x Chromium Single Cell 3' (kit version not stated in paper; HiSeq 4000)",
        processing="Whole cell (fresh tissue dissociated to single cells); tumour and matched normal "
                   "omentum/peritoneum",
        ms_ref="10", ms_ref_check="Matches data source",
        source="PMC8268616 Methods + Availability of data (EGAS00001004987, blueprint.lambrechtslab.org); "
               "PubMed 34238352"),
    "regner_2021": dict(
        label="Regner et al., 2021",
        title="A multi-omic single-cell landscape of human gynecologic malignancies",
        journal="Mol Cell 81(23), 4924-4941.e10 (2021)", doi="10.1016/j.molcel.2021.10.013", pmid="34739872",
        accession="GEO GSE173682 (HGSOC scRNA-seq GSM5276940, GSM5276943)",
        raw_accession="dbGaP phs002340.v1.p1 (controlled access)",
        platform="10x Chromium Single Cell 3' v3 (matched scATAC-seq not used)",
        processing="Whole cell for scRNA-seq (nuclei were isolated only for scATAC-seq)",
        ms_ref="19", ms_ref_check="Matches data source",
        source="GEO GSE173682 / GSM5276940 (Chromium Single Cell 3' GEM Kit v3); PMC8642316 Methods; "
               "PubMed 34739872"),
    "vazquez_garcia_2022": dict(
        label="Vázquez-García et al., 2022",
        title="Ovarian cancer mutational processes drive site-specific immune evasion",
        journal="Nature 612(7941), 778-786 (2022)", doi="10.1038/s41586-022-05496-1", pmid="36517593",
        accession="GEO GSE180661; Synapse syn25569736 (MSK SPECTRUM)",
        raw_accession="dbGaP phs002857.v1.p1 (controlled access)",
        platform="10x Chromium Single Cell 3' v3",
        processing="Whole cell; fresh tissue dissociated and FACS-sorted into live CD45+ and CD45- fractions "
                   "(fractions from the same site merged into one atlas sample)",
        ms_ref="11", ms_ref_check="Matches data source",
        source="GEO GSE180661 overall design; PMC9771812 Methods ('Chromium Single-Cell 3' Reagent kit v3') + "
               "Data availability; PubMed 36517593"),
    "xu_2022": dict(
        label="Xu et al., 2022",
        title="Single-Cell RNA Sequencing Reveals the Tissue Architecture in Human High-Grade Serous Ovarian "
              "Cancer",
        journal="Clin Cancer Res 28(16), 3590-3602 (2022)", doi="10.1158/1078-0432.CCR-22-0296", pmid="35675036",
        accession="GEO GSE184880 (HGSOC GSM5599225-GSM5599231; normal-ovary samples not used)",
        raw_accession="SRA SRP338994",
        platform="10x Chromium Single Cell 3' v3 (paper) / v3.1 (GEO sample record)",
        processing="Whole cell; enzymatic dissociation (collagenase IV/I, DNase I)",
        ms_ref="17",
        ms_ref_check="MISMATCH: ms ref 12 is Xu J et al., Cancer Lett 593:216928 (2024; PMID 38714290; 34 samples "
                     "/ 18 patients). The atlas xu_2022 data (GSE184880, 7 HGSOC) are from Xu J et al., "
                     "Clin Cancer Res 2022 (PMID 35675036)",
        source="GEO GSE184880 (Series_pubmed_id 35675036) / GSM5599225; PMC9662915 Methods + Data availability"),
    "zhang_2022": dict(
        label="Zhang et al., 2022",
        title="Longitudinal single-cell RNA-seq analysis reveals stress-promoted chemoresistance in metastatic "
              "ovarian cancer",
        journal="Sci Adv 8(8), eabm1831 (2022)", doi="10.1126/sciadv.abm1831", pmid="35196078",
        accession="GEO GSE165897",
        raw_accession="EGA EGAS00001005010 (paper) / EGAD00001006922 (GEO note) (controlled access)",
        platform="10x Chromium Single Cell 3' v2",
        processing="Whole cell; dissociated solid tumour, cryopreserved single-cell suspensions; paired "
                   "treatment-naive and post-NACT samples",
        ms_ref="15", ms_ref_check="Matches data source",
        source="GEO GSE165897; PMC8865800 Methods ('Chromium Single-Cell 3' Reagent Kit v. 2.0') + data "
               "availability; PubMed 35196078"),
    "zheng_2023": dict(
        label="Zheng et al., 2023",
        title="Single-cell analyses implicate ascites in remodeling the ecosystems of primary and metastatic "
              "tumors in ovarian cancer",
        journal="Nat Cancer 4(8), 1138-1156 (2023)", doi="10.1038/s43018-023-00599-8", pmid="37488416",
        accession="Mendeley Data 10.17632/rc47y6m9mp.1 (processed)",
        raw_accession="GSA-Human PRJCA005422",
        platform="10x Chromium Single Cell V(D)J reagent kit (5' gene expression + scTCR-seq)",
        processing="Whole cell; MACS Tumor Dissociation Kit (gentleMACS)",
        ms_ref="18", ms_ref_check="NOT CITED: no Zheng et al. entry in the manuscript reference list",
        source="PMC10447252 Methods + Data availability; PubMed 37488416; Mendeley ID rc47y6m9mp in raw/ file names"),
}

# Per-sample accessions for studies where the atlas sample name maps to a GEO sample file on disk.
# nath_2021: GEO sample record per (patient file, time point), from GSE158722 sample list.
NATH_GSM = {
    "P04_Time1": "GSM4808930", "P04_Time2": "GSM4808931", "P04_Time3": "GSM4808932",
    "P05_Time1": "GSM4808933", "P05_Time2": "GSM4808934", "P05_Time3": "GSM4808935",
    "P06_Time1": "GSM4808936", "P06_Time2": "GSM4808937", "P06_Time3": "GSM4808938",
    "P07_Time1": "GSM4808939", "P07_Time2": "GSM4808940",
    "P08_Time1": "GSM4808941", "P08_Time2": "GSM4808942", "P08_Time3": "GSM4808943",
    "P09_Time1": "GSM4808944", "P09_Time2": "GSM4808945", "P09_Time3": "GSM4808946",
    "P10": "GSM4808947", "P11": "GSM4808948", "P12": "GSM4808949", "P13": "GSM4808950", "P14": "GSM4808951",
    "P15": "GSM4808952", "P16": "GSM4808953", "P17": "GSM4808954", "P18": "GSM4808955", "P19": "GSM4808956",
    "P20": "GSM4808957", "P21": "GSM4808958", "P23": "GSM4808960", "P24": "GSM4808961",
}
# olalekan_2021 / regner_2021: name -> file mapping as coded in 202601_final_atlas/scripts/01_dataprep.ipynb
OLALEKAN_FILE = {"PT-1": "PT-3232", "PT-2": "PT-5150", "PT-3": "PT-6885", "PT-4": "PT-4806", "PT-5": "PT-3401",
                 "PT-6": "PT-2834"}
REGNER_FILE = {"P8": "3BAE2L", "P9": "3E5CFL"}

# Per-sample flags from the GEO cross-check (see README).
OLALEKAN_GEO = {  # GEO characteristics of the file each atlas sample was loaded from
    "PT-1": "GEO PT-3232: age 46, neoadjuvant No, FIGO IIIc",
    "PT-2": "GEO PT-5150: age 62, neoadjuvant Yes, FIGO IVb",
    "PT-3": "GEO PT-6885: age 66, neoadjuvant Yes, FIGO IIIc",
    "PT-4": "GEO PT-4806: age 71, neoadjuvant Yes, FIGO IIIc",
    "PT-5": "GEO PT-3401: age 56, neoadjuvant No, FIGO IIIc",
}


def load_h5ad_obs(path, cols):
    h = h5py.File(path, "r")
    out = {}
    for k in cols:
        g = h["obs"][k]
        if isinstance(g, h5py.Group):
            cats = np.array([c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]], dtype=object)
            out[k] = cats[g["codes"][:]]
        else:
            v = g[:]
            out[k] = np.array([x.decode() for x in v], dtype=object) if v.dtype.kind in "SO" else v
    idx = h["obs"][h["obs"].attrs["_index"]][:]
    h.close()
    df = pd.DataFrame(out)
    df.index = [x.decode() if isinstance(x, bytes) else x for x in idx]
    return df


def clean(v):
    """Missing / placeholder values -> '' (not reported)."""
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return ""
    s = str(v)
    return "" if s in ("nan", "NA", "None", "False", "") else s


def fmt_counts(series):
    vc = series.value_counts()
    return "; ".join(f"{k} {v}" for k, v in vc.items())


def raw_gsm_map():
    """Map (study, atlas sample name) -> GSM from the raw file names on disk."""
    m = {}
    pats = {
        "denisenko_2024": r"^(GSM\d+)_counts_(.+)\.txt\.gz$",
        "geistlinger_2020": r"^(GSM\d+)_(.+)_barcodes\.tsv\.gz$",
        "loret_2022": r"^(GSM\d+)_(.+)_filtered_gene_bc_matrices_h5\.h5$",
        "luo_2024": r"^(GSM\d+)_(.+)_barcodes\.tsv\.gz$",
        "xu_2022": r"^(GSM\d+)_(.+)\.barcodes\.tsv\.gz$",
    }
    for st, p in pats.items():
        for f in os.listdir(f"{RAW}/{st}"):
            r = re.match(p, f)
            if r and not f.startswith("._"):
                m[(st, r.group(2))] = r.group(1)
    for f in os.listdir(f"{RAW}/olalekan_2021"):
        r = re.match(r"^(GSM\d+)_(PT-\d+)\.csv\.gz$", f)
        if r:
            for name, fid in OLALEKAN_FILE.items():
                if fid == r.group(2):
                    m[("olalekan_2021", name)] = r.group(1)
    for f in os.listdir(f"{RAW}/regner_2021"):
        r = re.match(r"^(GSM\d+)_barcodes-(\w+)\.tsv\.gz$", f)
        if r:
            for name, fid in REGNER_FILE.items():
                if fid == r.group(2):
                    m[("regner_2021", name)] = r.group(1)
    return m


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # ---------------- load cell-level metadata (all cells) ----------------
    t1 = pd.read_csv(T1_CSV, low_memory=False)
    assert len(t1) == 1980703, len(t1)
    t1["sample_id"] = t1["sample_id"].astype(int).astype(str)
    t1["patient_id"] = t1["patient_id"].astype(int).astype(str)

    # ---------------- cross-check against final atlas obs ----------------
    obs = load_h5ad_obs(H5AD, ["study", "sample_id", "patient_id", "sample_num", "celltype_level1",
                               "total_counts", "n_genes_by_counts"])
    assert (obs.index.values == t1["cell_id"].values).all(), "cell order differs between T1 and h5ad"
    norm = lambda s: pd.Series(s).astype(str).str.replace(r"\.0$", "", regex=True).values
    for c in ["study", "celltype_level1"]:
        assert (norm(obs[c]) == norm(t1[c])).all(), c
    assert (norm(obs["sample_id"]) == t1["sample_id"].values).all()
    assert (norm(obs["patient_id"]) == t1["patient_id"].values).all()
    assert (obs["total_counts"].values == t1["total_counts"].values).all()
    lib = pd.Series(obs["sample_num"].values).astype(str).str.replace(r"\.txt$", "", regex=True).values

    # Nath 2021 sample mapping fix (2026-10-20): each cell to its GEO patient and time point
    fix = nath2021_correction()
    is_n = (t1["study"] == "nath_2021").values
    ids = t1.loc[is_n, "cell_id"].astype(str).values
    assert len(fix) == is_n.sum() and pd.Index(ids).isin(fix.index).all(), "Nath correction does not cover the atlas"
    for c in ["sample_id", "patient_id", "treatment_status", "TP53_status", "BRCA_status"]:
        t1[c] = t1[c].astype(object)
        t1.loc[is_n, c] = fix.loc[ids, c].values
    lib[is_n] = fix.loc[ids, "sample_num"].values
    t1["library"] = lib

    totals = dict(
        cells=len(t1), studies=t1.study.nunique(),
        patients=t1.groupby(["study", "patient_id"]).ngroups,
        samples=t1.groupby(["study", "sample_id"]).ngroups,
        libraries=t1.groupby(["study", "library"]).ngroups,
        sites=t1.anatomic_site.nunique())
    print("Totals:", totals)

    # curated sheet (cross-check of per-sample fields)
    xl = pd.read_excel(XLSX)
    xl_by_num = xl.groupby("sample_num")

    gsm = raw_gsm_map()

    # ---------------- Samples sheet ----------------
    meta_cols = ["patient_id", "anatomic_site", "metastatic_site", "treatment_status", "stage", "age",
                 "treatment_response", "BRCA_status", "HRD_status", "TP53_status"]
    g = t1.groupby(["study", "sample_id"], sort=False)
    rows = []
    for (st, sid), d in g:
        first = d.iloc[0]
        for c in meta_cols:
            assert d[c].astype(str).nunique() == 1, (st, sid, c)
        libs = list(dict.fromkeys(d["library"]))
        x = xl_by_num.get_group(int(sid))
        assert (x.study == st).all() and (x.patient_id.astype(str) == first.patient_id).all()
        note = ""
        # source accession(s)
        if st == "nath_2021":
            pref = d["barcode"].str.extract(r"^[^:]+:(P\d+(?:_Time\d)?)_")[0]
            parts = sorted(dict.fromkeys(pref))
            acc = ", ".join(NATH_GSM[p] for p in parts)
            src_file = "GSE158722_" + parts[0].split("_")[0] + ".counts.txt.gz"
            assert len(parts) == 1, (sid, parts)   # one GEO sample per atlas sample since the fix
            note = (f"GEO {parts[0].replace('_', ' ')} (file {src_file}); labels from the GEO cell ids "
                    f"(Nath 2021 mapping fix, see README)")
        else:
            acc = ", ".join(dict.fromkeys(gsm.get((st, l), "") for l in libs)).strip(", ")
            if st == "olalekan_2021":
                note = OLALEKAN_GEO.get(libs[0], "")
                note = (note + " -- atlas treatment/age/stage differ from GEO record of the loaded file "
                        "(see README)") if libs[0] in ("PT-1", "PT-2", "PT-4", "PT-5") else note
            if st == "geistlinger_2020" and libs[0] == "T76":
                note = ("Loaded from file GSM4675274_T76_*, whose GEO title is 'Tumor T61' (IIIc, refractory); "
                        "atlas metadata follow GEO 'Tumor T76' (GSM4675275, IVb, resistant)")
        age = clean(first.age)
        rows.append({
            "Study": STUDY_INFO[st]["label"], "Study key": st,
            "Atlas sample ID": int(sid),
            "Original sample name(s)": " + ".join(libs), "n libraries": len(libs),
            "Source sample accession(s)": acc,
            "Atlas patient ID": int(first.patient_id),
            "Anatomic site": clean(first.anatomic_site), "Metastatic site class": clean(first.metastatic_site),
            "Treatment status": clean(first.treatment_status), "Stage": clean(first.stage),
            "Age (years)": int(float(age)) if age else "",
            "Treatment response": clean(first.treatment_response),
            "BRCA status": clean(first.BRCA_status), "HRD status": clean(first.HRD_status),
            "TP53 status": clean(first.TP53_status),
            "n cells": len(d),
            "Median UMI per cell": float(np.median(d.total_counts)),
            "Median genes per cell": float(np.median(d.n_genes_by_counts)),
            "Note": note,
        })
    samples = pd.DataFrame(rows).sort_values(["Study key", "Atlas patient ID", "Atlas sample ID"])
    assert samples["n cells"].sum() == totals["cells"]

    # ---------------- Studies sheet ----------------
    srows = []
    for st, d in t1.groupby("study", sort=True):
        info = STUDY_INFO[st]
        s = samples[samples["Study key"] == st]
        ts = s["Treatment status"]
        srows.append({
            "Study": info["label"], "Study key": st, "Title": info["title"], "Journal": info["journal"],
            "DOI": info["doi"], "PMID": info["pmid"],
            "Data accession(s) used": info["accession"], "Raw data accession": info["raw_accession"],
            "Sequencing platform / chemistry": info["platform"],
            "Tissue processing": info["processing"],
            "n patients": d.groupby("patient_id").ngroups,
            "n samples": len(s), "n libraries": int(s["n libraries"].sum()),
            "n cells": len(d),
            "Median UMI per cell": float(np.median(d.total_counts)),
            "Median genes per cell": float(np.median(d.n_genes_by_counts)),
            "Anatomic sites (n samples)": fmt_counts(s["Anatomic site"]),
            "Metastatic site class (n samples)": fmt_counts(s["Metastatic site class"]),
            "n samples treatment-naive": int(ts.isin(TREAT_NAIVE).sum()),
            "n samples post-chemotherapy": int(ts.isin(TREAT_CHEMO).sum()),
            "n samples other treated (PARPi/ICB +/- chemo)": int(ts.isin(TREAT_OTHER).sum()),
            "n samples treatment unknown": int((ts == "").sum()),
            "Treatment status detail (n samples)": fmt_counts(ts.replace("", "not reported")),
            "Manuscript ref #": info["ms_ref"],
            "Verification source": info["source"],
        })
    studies = pd.DataFrame(srows)
    ts_all = samples["Treatment status"]
    total = {c: "" for c in studies.columns}
    total.update({
        "Study": "Total (13 studies)",
        "n patients": totals["patients"], "n samples": totals["samples"], "n libraries": totals["libraries"],
        "n cells": totals["cells"],
        "Median UMI per cell": float(np.median(t1.total_counts)),
        "Median genes per cell": float(np.median(t1.n_genes_by_counts)),
        "Anatomic sites (n samples)": f"{totals['sites']} sites: " + fmt_counts(samples["Anatomic site"]),
        "Metastatic site class (n samples)": fmt_counts(samples["Metastatic site class"]),
        "n samples treatment-naive": int(ts_all.isin(TREAT_NAIVE).sum()),
        "n samples post-chemotherapy": int(ts_all.isin(TREAT_CHEMO).sum()),
        "n samples other treated (PARPi/ICB +/- chemo)": int(ts_all.isin(TREAT_OTHER).sum()),
        "n samples treatment unknown": int((ts_all == "").sum()),
        "Treatment status detail (n samples)": fmt_counts(ts_all.replace("", "not reported")),
    })
    studies = pd.concat([studies, pd.DataFrame([total])], ignore_index=True)
    assert studies["n cells"].iloc[:-1].sum() == totals["cells"]
    assert studies["n patients"].iloc[:-1].sum() == totals["patients"]

    # ---------------- Cells_by_lineage ----------------
    ct = pd.crosstab(t1.study, t1.celltype_level1)
    assert set(ct.columns) == set(LINEAGES)
    ct = ct[LINEAGES]
    ct.insert(0, "Study", [STUDY_INFO[s]["label"] for s in ct.index])
    ct.index.name = "Study key"
    ct = ct.reset_index()[["Study", "Study key"] + LINEAGES]
    ct["Total"] = ct[LINEAGES].sum(1)
    tot = {"Study": "Total", "Study key": ""}
    tot.update({c: int(ct[c].sum()) for c in LINEAGES + ["Total"]})
    ct = pd.concat([ct, pd.DataFrame([tot])], ignore_index=True)
    assert ct["Total"].iloc[-1] == totals["cells"]
    # obs cross-check of lineage totals
    assert (pd.Series(obs["celltype_level1"]).value_counts()[LINEAGES].values == ct[LINEAGES].iloc[-1].values).all()

    # ---------------- README ----------------
    wt = {st: (int(((samples["Study key"] == st) & (samples["TP53 status"] == "wildtype")).sum()),
               int(samples.loc[(samples["Study key"] == st) & (samples["TP53 status"] == "wildtype"),
                               "Atlas patient ID"].nunique())) for st in ["luo_2024", "nath_2021"]}
    readme = [
        ("Supplementary Table 1", "Metadata of the 13 published scRNA-seq studies integrated into the HGSC atlas."),
        ("Counts", f"All numbers are computed from ALL cells of the final post-QC atlas (hgsc_atlas_final.h5ad; "
                   f"{totals['cells']:,} cells after QC filtering [>=500 UMI, >=300 genes, >=500 cells/sample], "
                   f"doublet removal and exclusion of doublet/artifact clusters). No subsampling."),
        ("Totals", f"{totals['studies']} studies; {totals['patients']} patients; {totals['samples']} tissue samples "
                   f"({totals['libraries']} sequencing libraries); {totals['sites']} anatomic sites; "
                   f"{totals['cells']:,} cells."),
        ("Sample definition", "An atlas 'sample' is one tissue specimen (patient x site x timepoint) = atlas "
                              "sample_id. Libraries from sorted fractions of the same specimen (Vázquez-García "
                              "CD45+/CD45-, Hornburg tumour/immune/stroma) are merged into one sample; 'n "
                              "libraries' counts the underlying sequencing libraries."),
        ("Patient / sample IDs", "Atlas patient_id and sample_id are integer codes unique across the atlas; "
                                 "'Original sample name(s)' are the library names from the source data."),
        ("Blank cells", "Blank = not reported in the source metadata (or not verifiable for publication fields). "
                        "Placeholders stored in the atlas metadata ('False' for Olbrecht HRD/TP53 status, 'nan' for "
                        "Hornburg) are shown as blank."),
        ("TP53 status", "As reported by each source study. The TP53 wild-type calls for Luo et al. 2024 "
                        f"({wt['luo_2024'][0]} samples, {wt['luo_2024'][1]} patients) and Nath et al. 2021 "
                        f"({wt['nath_2021'][0]} samples, {wt['nath_2021'][1]} patients) could not be verified against "
                        "the source publications."),
        ("", ""),
        ("Sheet: Studies", "One row per study plus a Total row."),
        ("  Data accession(s) used", "Repository/accession of the processed count data used to build the atlas."),
        ("  Raw data accession", "Repository of raw reads (often controlled access)."),
        ("  Sequencing platform / chemistry", "Single-cell platform and chemistry as stated in the GEO/EGA record "
                                              "or the paper's Methods (see 'Verification source')."),
        ("  Tissue processing", "Whole cell vs nuclei and any sorting, as stated in the source."),
        ("  n patients / n samples / n libraries / n cells", "Counted in the final atlas (unique (study, patient_id)"
                                                             " and (study, sample_id))."),
        ("  Median UMI / genes per cell", "Median of total_counts / n_genes_by_counts (atlas obs) over all cells."),
        ("  Treatment columns", "naive = pre-treatment; post-chemotherapy = chemotherapy only; other treated = "
                                "post-niraparib (neoadjuvant), post-chemotherapy_niraparib / _olaparib / _pembro / _niraparib_pembro; "
                                "unknown = not reported. Counts are samples."),
        ("  Manuscript ref #", "Number of the data-generating paper in the manuscript reference list."),
        ("Sheet: Samples", "One row per final-atlas sample. Clinical fields as curated in the atlas "
                           "(T1_atlas_metadata.csv, post-correction). 'Source sample accession(s)' = GEO GSM(s) "
                           "where derivable from raw file names on disk (nath_2021: from GEO barcode prefixes)."),
        ("Sheet: Cells_by_lineage", "Cells per study x celltype_level1 (13 lineages)."),
        ("", ""),
        ("FLAGS / TO VERIFY", ""),
        ("nath_2021 sample mapping (corrected)", "GEO GSE158722 count files GSE158722_Pxx are per PATIENT (P01-P09 "
                                     "each contain 2-3 longitudinal time points; P10-P24 one sample each). The atlas "
                                     "dataprep had assigned the 23 files, in order, to the curated names P1-1 ... P8-2, "
                                     "P23, P24, so sample names, patient grouping and treatment did not match GEO. "
                                     "Corrected (figures/_prep/nath2021_fix.py): every Nath cell is assigned "
                                     "from its GEO cell id to its GEO patient and time point, and takes that sample's "
                                     "curated name, patient, treatment, TP53 and BRCA status (the curated design "
                                     "matches GEO). Each Nath sample is now one GEO sample (GSM in 'Source sample "
                                     "accession(s)')."),
        ("olalekan_2021 metadata", "For 4 of 5 atlas samples (PT-1, PT-2, PT-4, PT-5) the curated age / "
                                   "neoadjuvant status / stage differ from the GEO record of the file loaded "
                                   "(GSM4416534-8); PT-3 matches. Likely an ordering mismatch between the curated "
                                   "table and the GEO files; verify against the paper's patient table."),
        ("geistlinger_2020 T76", "GEO file names and titles disagree (GSM4675274 file 'T76' = title 'Tumor T61'; "
                                 "GSM4675275 file 'T77' = title 'Tumor T76'). The atlas loaded file 'T76' but uses "
                                 "the stage/response of GEO 'Tumor T76'."),
        ("Unverified fields", "Loret 2022: 3'/5' chemistry and tissue processing; Hornburg 2021: 3'/5' "
                              "chemistry; Olbrecht 2021: 3' kit version; Zheng 2023: 5' kit version."),
    ]
    readme = pd.DataFrame(readme, columns=["Item", "Description"])

    # ---------------- write ----------------
    sheets = {"README": readme, "Studies": studies, "Samples": samples, "Cells_by_lineage": ct}
    for name, df in sheets.items():
        df.to_csv(f"{OUT_DIR}/SuppTable1_atlas_study_metadata_{name}.csv", index=False)
    write_xlsx(sheets, OUT_XLSX)
    print("Wrote", OUT_XLSX)
    print(studies[["Study", "n patients", "n samples", "n libraries", "n cells"]].to_string(index=False))


def write_xlsx(sheets, path):
    wb = Workbook()
    wb.remove(wb.active)
    hdr_font = Font(bold=True)
    hdr_fill = PatternFill("solid", fgColor="DDDDDD")
    wide = {"Title": 60, "Journal": 30, "Data accession(s) used": 45, "Raw data accession": 35,
            "Sequencing platform / chemistry": 45, "Tissue processing": 50, "Anatomic sites (n samples)": 45,
            "Metastatic site class (n samples)": 35, "Treatment status detail (n samples)": 45,
            "Manuscript citation check": 60, "Verification source": 70, "Note": 70, "Description": 120,
            "Original sample name(s)": 45, "Item": 34}
    for name, df in sheets.items():
        ws = wb.create_sheet(name)
        ws.append(list(df.columns))
        for row in df.itertuples(index=False):
            ws.append([None if (isinstance(v, str) and v == "") else v for v in row])
        for c in ws[1]:
            c.font = hdr_font
            c.fill = hdr_fill
            c.alignment = Alignment(wrap_text=True, vertical="top")
        ws.freeze_panes = "B2" if name != "README" else "A2"
        for j, col in enumerate(df.columns, 1):
            if col in wide:
                w = wide[col]
            else:
                lens = [len(str(col))] + [len(str(v)) for v in df[col].head(500)]
                w = min(max(10, int(np.percentile(lens, 95)) + 2), 40)
            ws.column_dimensions[get_column_letter(j)].width = w
            wrap = col in wide
            for cell in ws.iter_rows(min_row=2, min_col=j, max_col=j):
                cell[0].alignment = Alignment(wrap_text=wrap, vertical="top")
                if isinstance(cell[0].value, (int, float)) and col.startswith(("n ", "Total")) \
                        or col in LINEAGES:
                    cell[0].number_format = "#,##0"
                if col.startswith("Median"):
                    cell[0].number_format = "#,##0"
        if name in ("Studies", "Cells_by_lineage"):
            for c in ws[ws.max_row]:
                c.font = hdr_font
        if name == "README":
            for r in ws.iter_rows(min_row=2):
                if r[0].value and (r[0].value.startswith("Sheet") or r[0].value == "FLAGS / TO VERIFY"):
                    r[0].font = hdr_font
    wb.save(path)


if __name__ == "__main__":
    main()
