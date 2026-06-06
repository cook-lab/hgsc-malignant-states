#!/usr/bin/env python3
"""
Atlas 01 — Step 01: per-study raw data preparation

PURPOSE
    Load each study's raw matrices (10x mtx / per-study h5ad / GEO supplements),
    map gene IDs to symbols (mygene), attach harmonised sample/clinical metadata
    from the atlas metadata sheet, and write one harmonised per-study h5ad. These
    feed step 02 (concat + QC + Scrublet), whose output
    atlas_concatenated_filtered.h5ad is the raw input to the official integration
    chain (steps 03–07: preprocess → CellAssign → 5-method integration →
    benchmark → finalize). See this stage's README.

INPUTS
    DATA_ROOT/2026_final_atlas/raw/<study>/...        (raw per-study matrices; not deposited)
    DATA_ROOT/2026_final_atlas/atlas_metadata.xlsx    (harmonised sample metadata)

OUTPUTS
    DATA_ROOT/2026_final_atlas/processed/<study>.h5ad (one harmonised object per study)

MANUSCRIPT PANEL(S)
    Pre-integration provenance; no panel rendered directly.

RUNTIME TIER
    heavy (raw ingestion of 13 studies; mygene queries require network).

NOTE
    Migrated from 01_dataprep.ipynb. Per-study logic and gene-mapping preserved
    exactly; only hardcoded /Volumes/CookLab/Sarah paths were replaced with
    config-rooted RAW_ROOT / PROCESSED_ROOT / META_XLSX.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path, SEED  # noqa: E402

np.random.seed(SEED)

RAW_ROOT       = path("data_root", "2026_final_atlas", "raw")
PROCESSED_ROOT = path("data_root", "2026_final_atlas", "processed")
META_XLSX      = path("data_root", "2026_final_atlas", "atlas_metadata.xlsx")

import scanpy as sc
import pandas as pd
import anndata as ad
from pathlib import Path
import scipy.sparse as sp
import mygene
from scipy import sparse
import gzip

meta_all = pd.read_excel(f"{META_XLSX}")

#geistlinger
import mygene
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy
mg = mygene.MyGeneInfo()

#load samples
samples = {
    "T59": {
        "barcodes": f"{RAW_ROOT}/geistlinger_2020/GSM4675273_T59_barcodes.tsv.gz",
        "genes": f"{RAW_ROOT}/geistlinger_2020/GSM4675273_T59_genes.tsv.gz",
        "matrix": f"{RAW_ROOT}/geistlinger_2020/GSM4675273_T59_matrix.mtx.gz"
    },
    "T76": {
        "barcodes": f"{RAW_ROOT}/geistlinger_2020/GSM4675274_T76_barcodes.tsv.gz",
        "genes": f"{RAW_ROOT}/geistlinger_2020/GSM4675274_T76_genes.tsv.gz",
        "matrix": f"{RAW_ROOT}/geistlinger_2020/GSM4675274_T76_matrix.mtx.gz"
    },
    "T89": {
        "barcodes": f"{RAW_ROOT}/geistlinger_2020/GSM4675276_T89_barcodes.tsv.gz",
        "genes": f"{RAW_ROOT}/geistlinger_2020/GSM4675276_T89_genes.tsv.gz",
        "matrix": f"{RAW_ROOT}/geistlinger_2020/GSM4675276_T89_matrix.mtx.gz"
    },
    "T90": {
        "barcodes": f"{RAW_ROOT}/geistlinger_2020/GSM4675277_T90_barcodes.tsv.gz",
        "genes": f"{RAW_ROOT}/geistlinger_2020/GSM4675277_T90_genes.tsv.gz",
        "matrix": f"{RAW_ROOT}/geistlinger_2020/GSM4675277_T90_matrix.mtx.gz"
    }
}

adata_list = []
for sample_id, files in samples.items():
    print(f"\nLoading sample {sample_id}...")
    X = sc.read_mtx(files["matrix"]).X.T
    barcodes = pd.read_csv(files["barcodes"], header=None, sep='\t')[0].values
    genes = pd.read_csv(files["genes"], header=None, sep='\t')[0].values
    #convert the ensembl to gene symbols
    query_result = mg.querymany(genes, scopes='ensembl.gene', fields='symbol', species='human')
    ensembl_to_symbol = {entry['query']: entry.get('symbol', entry['query']) for entry in query_result}
    gene_symbols = [ensembl_to_symbol.get(e, e) for e in genes]
    # create AnnData
    adata = ad.AnnData(X)
    adata.obs['barcode'] = barcodes
    adata.obs['sample_id'] = sample_id
    adata.var_names = gene_symbols
    adata.var_names_make_unique()
    # store raw counts
    adata.raw = adata
    adata_list.append(adata)
for adata in adata_list:
    adata.obs_names_make_unique()

#concatenate samples
adata_all = ad.concat(
    adata_list,
    join='outer',
    label=None,  # no batch column
    keys=None,
    index_unique=None
)

#adding metadata
metadata_dict = {
    'T59': {'study': 'geistlinger_2020', 'patient_id': 6, 'sample_num': 6,
            'treatment_status': None, 'histological_subtype': 'Serous', 'stage': 'IV',
            'anatomic_site': 'omentum', 'metastatic_site': 'metastasis',
            'age': None, 'treatment_response': 'resistant', 'BRCA_status': None,
            'HRD_status': None, 'TP53_status': None, 'ref': None},
    'T76': {'study': 'geistlinger_2020', 'patient_id': 8, 'sample_num': 8,
            'treatment_status': None, 'histological_subtype': 'Serous', 'stage': 'IV',
            'anatomic_site': 'omentum', 'metastatic_site': 'metastasis',
            'age': None, 'treatment_response': 'resistant', 'BRCA_status': None,
            'HRD_status': None, 'TP53_status': None, 'ref': None},
    'T89': {'study': 'geistlinger_2020', 'patient_id': 9, 'sample_num': 9,
            'treatment_status': None, 'histological_subtype': 'Serous', 'stage': 'IV',
            'anatomic_site': 'omentum', 'metastatic_site': 'metastasis',
            'age': None, 'treatment_response': 'sensitive', 'BRCA_status': None,
            'HRD_status': None, 'TP53_status': None, 'ref': None},
    'T90': {'study': 'geistlinger_2020', 'patient_id': 10, 'sample_num': 10,
            'treatment_status': None, 'histological_subtype': 'Serous', 'stage': 'IV',
            'anatomic_site': 'omentum', 'metastatic_site': 'metastasis',
            'age': None, 'treatment_response': 'sensitive', 'BRCA_status': None,
            'HRD_status': None, 'TP53_status': None, 'ref': None}
}

for col in metadata_dict['T59'].keys():
    adata_all.obs[col] = adata_all.obs['sample_id'].map(lambda s: metadata_dict[s][col])

#convert all metadata columns except 'barcode' and 'sample_id' to strings
for col in adata_all.obs.columns:
    if col not in ['barcode', 'sample_id']:
        adata_all.obs[col] = adata_all.obs[col].astype(str)

#make obs_names unique
adata_all.obs_names_make_unique()

#save data
adata_all.write(f"{PROCESSED_ROOT}/geistlinger_2020.h5ad")
print("Saved.")

#denisenko
samples = {
    "Y2": f"{RAW_ROOT}/denisenko_2024/GSM6506105_counts_Y2.txt.gz",
    "Y3": f"{RAW_ROOT}/denisenko_2024/GSM6506106_counts_Y3.txt.gz",
    "Y5": f"{RAW_ROOT}/denisenko_2024/GSM6506107_counts_Y5.txt.gz",
    "MJ10": f"{RAW_ROOT}/denisenko_2024/GSM6506108_counts_MJ10.txt.gz",
    "MJ11": f"{RAW_ROOT}/denisenko_2024/GSM6506109_counts_MJ11.txt.gz"
}

# Full metadata table
metadata_table = pd.DataFrame([
    ["denisenko_2024","Y2",1,1,"post-chemotherapy","Serous","III/IV","adnexa","primary","","","","","","GRCh38"],
    ["denisenko_2024","Y3",2,2,"post-chemotherapy","Serous","III/IV","adnexa","primary","","","","","","GRCh38"],
    ["denisenko_2024","Y5",3,3,"post-chemotherapy","Serous","III/IV","adnexa","primary","","","","","","GRCh38"],
    ["denisenko_2024","MJ10",4,4,"post-chemotherapy","Serous","III/IV","adnexa","primary","","","","","","GRCh38"],
    ["denisenko_2024","MJ11",5,5,"post-chemotherapy","Serous","III/IV","adnexa","primary","","","","","","GRCh38"]
], columns=["study","sample_id","patient_id","sample_num","treatment_status","histological_subtype","stage",
            "anatomic_site","metastatic_site","age","treatment_response","BRCA_status","HRD_status","TP53_status","ref"])

# Initialize MyGeneInfo
mg = mygene.MyGeneInfo()
adatas = []

for sample_id, file_path in samples.items():
    print(f"Loading sample {sample_id} from {file_path}...")
    
    # Load count matrix
    counts = pd.read_csv(file_path, index_col=0, sep="\t")
    
    # Transpose to cells x genes
    X = sparse.csr_matrix(counts.T.values)
    
    # Create AnnData
    adata = ad.AnnData(X)
    adata.obs['barcode'] = counts.columns.tolist()
    adata.obs['sample_id'] = sample_id
    adata.var_names = counts.index.tolist()
    adata.raw = adata  # keep raw counts
    
    # Convert gene symbols via MyGene
    genes = adata.var_names.tolist()
    query_result = mg.querymany(
        genes,
        scopes=['ensembl.gene','alias','symbol'],
        fields='symbol',
        species='human'
    )
    ensembl_to_symbol = {entry['query']: entry.get('symbol', entry['query']) for entry in query_result}
    adata.var_names = [ensembl_to_symbol.get(g, g) for g in genes]
    
    # Make gene symbols unique to avoid concat errors
    adata.var_names_make_unique()
    
    adatas.append(adata)

# Concatenate all samples
adata_all = ad.concat(adatas, join='outer', label='sample_id', keys=[s for s in samples.keys()])

# Merge metadata table
adata_all.obs = adata_all.obs.merge(metadata_table, on='sample_id', how='left')

# Convert metadata columns to string
for col in adata_all.obs.columns:
    if col not in ['barcode','sample_id']:
        adata_all.obs[col] = adata_all.obs[col].astype(str)

adata_all.obs_names_make_unique()

# Save AnnData
adata_all.write(f"{PROCESSED_ROOT}/denisenko_2024.h5ad")
print("Saved.")

#loret
# Initialize MyGeneInfo
mg = mygene.MyGeneInfo()

# Sample files
samples = {
    "1_N_OT_PT1": f"{RAW_ROOT}/loret_2022/GSM6049610_1_N_OT_PT1_filtered_gene_bc_matrices_h5.h5",
    "2_N_A_PT1": f"{RAW_ROOT}/loret_2022/GSM6049611_2_N_A_PT1_filtered_gene_bc_matrices_h5.h5",
    "3_N_PER_PT1": f"{RAW_ROOT}/loret_2022/GSM6049612_3_N_PER_PT1_filtered_gene_bc_matrices_h5.h5",
    "4_N_OM_PT1": f"{RAW_ROOT}/loret_2022/GSM6049613_4_N_OM_PT1_filtered_gene_bc_matrices_h5.h5",
    "5_N_BL_PT1": f"{RAW_ROOT}/loret_2022/GSM6049614_5_N_BL_PT1_filtered_gene_bc_matrices_h5.h5",
    "6_T_OT_PT1": f"{RAW_ROOT}/loret_2022/GSM6049615_6_T_OT_PT1_filtered_gene_bc_matrices_h5.h5",
    "7_T_PER_PT1": f"{RAW_ROOT}/loret_2022/GSM6049616_7_T_PER_PT1_filtered_gene_bc_matrices_h5.h5",
    "8_T_OM_PT1": f"{RAW_ROOT}/loret_2022/GSM6049617_8_T_OM_PT1_filtered_gene_bc_matrices_h5.h5",
    "9_T_A_PT1": f"{RAW_ROOT}/loret_2022/GSM6049618_9_T_A_PT1_filtered_gene_bc_matrices_h5.h5",
    "10_N_OT_PT2": f"{RAW_ROOT}/loret_2022/GSM6049619_10_N_OT_PT2_filtered_gene_bc_matrices_h5.h5",
    "11_N_A_PT2": f"{RAW_ROOT}/loret_2022/GSM6049620_11_N_A_PT2_filtered_gene_bc_matrices_h5.h5",
    "12_N_OM_PT2": f"{RAW_ROOT}/loret_2022/GSM6049621_12_N_OM_PT2_filtered_gene_bc_matrices_h5.h5",
    "13_T_OT_PT2": f"{RAW_ROOT}/loret_2022/GSM6049622_13_T_OT_PT2_filtered_gene_bc_matrices_h5.h5",
    "14_T_OM_PT2": f"{RAW_ROOT}/loret_2022/GSM6049623_14_T_OM_PT2_filtered_gene_bc_matrices_h5.h5",
    "15_T_PER_PT2": f"{RAW_ROOT}/loret_2022/GSM6049624_15_T_PER_PT2_filtered_gene_bc_matrices_h5.h5",
    "16_N_OT_PT3": f"{RAW_ROOT}/loret_2022/GSM6049625_16_N_OT_PT3_filtered_gene_bc_matrices_h5.h5",
    "17_N_A_PT3": f"{RAW_ROOT}/loret_2022/GSM6049626_17_N_A_PT3_filtered_gene_bc_matrices_h5.h5",
    "18_N_PER_PT3": f"{RAW_ROOT}/loret_2022/GSM6049627_18_N_PER_PT3_filtered_gene_bc_matrices_h5.h5",
    "19_T_OT_PT3": f"{RAW_ROOT}/loret_2022/GSM6049628_19_T_OT_PT3_filtered_gene_bc_matrices_h5.h5",
    "20_T_A_PT3": f"{RAW_ROOT}/loret_2022/GSM6049629_20_T_A_PT3_filtered_gene_bc_matrices_h5.h5",
    "21_T_OM_PT3": f"{RAW_ROOT}/loret_2022/GSM6049630_21_T_OM_PT3_filtered_gene_bc_matrices_h5.h5",
    "22_T_PER_PT3": f"{RAW_ROOT}/loret_2022/GSM6049631_22_T_PER_PT3_filtered_gene_bc_matrices_h5.h5"
}

# Metadata table
metadata_table = pd.DataFrame([
    ["loret_2022","1_N_OT_PT1",26,26,"pre-treatment","Serous","III","adnexa","primary",58,"resistant","wildtype","","","Hg19"],
    ["loret_2022","2_N_A_PT1",26,27,"pre-treatment","Serous","III","ascites","ascites",58,"resistant","wildtype","","","Hg19"],
    ["loret_2022","3_N_PER_PT1",26,28,"pre-treatment","Serous","III","peritoneum","metastasis",58,"resistant","wildtype","","","Hg19"],
    ["loret_2022","4_N_OM_PT1",26,29,"pre-treatment","Serous","III","omentum","metastasis",58,"resistant","wildtype","","","Hg19"],
    ["loret_2022","5_N_BL_PT1",26,30,"pre-treatment","Serous","III","bladder","metastasis",58,"resistant","wildtype","","","Hg19"],
    ["loret_2022","6_T_OT_PT1",26,31,"post-chemotherapy","Serous","III","adnexa","primary",58,"resistant","wildtype","","","Hg19"],
    ["loret_2022","7_T_PER_PT1",26,32,"post-chemotherapy","Serous","III","peritoneum","metastasis",58,"resistant","wildtype","","","Hg19"],
    ["loret_2022","8_T_OM_PT1",26,33,"post-chemotherapy","Serous","III","omentum","metastasis",58,"resistant","wildtype","","","Hg19"],
    ["loret_2022","9_T_A_PT1",26,34,"post-chemotherapy","Serous","III","ascites","ascites",58,"resistant","wildtype","","","Hg19"],
    ["loret_2022","10_N_OT_PT2",27,35,"pre-treatment","Serous","III","adnexa","primary",62,"sensitive","wildtype","","","Hg19"],
    ["loret_2022","11_N_A_PT2",27,36,"pre-treatment","Serous","III","ascites","ascites",62,"sensitive","wildtype","","","Hg19"],
    ["loret_2022","12_N_OM_PT2",27,37,"pre-treatment","Serous","III","omentum","metastasis",62,"sensitive","wildtype","","","Hg19"],
    ["loret_2022","13_T_OT_PT2",27,38,"post-chemotherapy","Serous","III","adnexa","primary",62,"sensitive","wildtype","","","Hg19"],
    ["loret_2022","14_T_OM_PT2",27,39,"post-chemotherapy","Serous","III","omentum","metastasis",62,"sensitive","wildtype","","","Hg19"],
    ["loret_2022","15_T_PER_PT2",27,40,"post-chemotherapy","Serous","III","peritoneum","metastasis",62,"sensitive","wildtype","","","Hg19"],
    ["loret_2022","16_N_OT_PT3",28,41,"pre-treatment","Serous","III","adnexa","primary",48,"sensitive","wildtype","","","Hg19"],
    ["loret_2022","17_N_A_PT3",28,42,"pre-treatment","Serous","III","ascites","ascites",48,"sensitive","wildtype","","","Hg19"],
    ["loret_2022","18_N_PER_PT3",28,43,"pre-treatment","Serous","III","peritoneum","metastasis",48,"sensitive","wildtype","","","Hg19"],
    ["loret_2022","19_T_OT_PT3",28,44,"post-chemotherapy","Serous","III","adnexa","primary",48,"sensitive","wildtype","","","Hg19"],
    ["loret_2022","20_T_A_PT3",28,45,"post-chemotherapy","Serous","III","ascites","ascites",48,"sensitive","wildtype","","","Hg19"],
    ["loret_2022","21_T_OM_PT3",28,46,"post-chemotherapy","Serous","III","omentum","metastasis",48,"sensitive","wildtype","","","Hg19"],
    ["loret_2022","22_T_PER_PT3",28,47,"post-chemotherapy","Serous","III","peritoneum","metastasis",48,"sensitive","wildtype","","","Hg19"]
], columns=["study","sample_id","patient_id","sample_num","treatment_status","histological_subtype",
            "stage","anatomic_site","metastatic_site","age","treatment_response","BRCA_status","HRD_status","TP53_status","ref"])

# Load each sample, add metadata, concatenate

adatas = []

for sample_id, file_path in samples.items():
    print(f"Loading sample {sample_id}...")
    
    # Load 10X H5
    adata = sc.read_10x_h5(file_path)
    
    # Make gene names unique
    adata.var_names_make_unique()
    
    # Store raw counts
    adata.raw = adata
    
    # Add sample info
    adata.obs['sample_id'] = sample_id
    adata.obs['barcode'] = adata.obs_names
    
    adatas.append(adata)

# Concatenate all Loret samples
adata_all = ad.concat(adatas, join='outer', label='sample_id', keys=list(samples.keys()))

# Merge metadata table
adata_all.obs = adata_all.obs.merge(metadata_table, on='sample_id', how='left')

# Convert metadata columns to string
for col in adata_all.obs.columns:
    if col not in ['barcode','sample_id']:
        adata_all.obs[col] = adata_all.obs[col].astype(str)

adata_all.obs_names_make_unique()

# Save AnnData
adata_all.write(f"{PROCESSED_ROOT}/loret_2022.h5ad")
print("Saved.")

#luo
# Initialize MyGeneInfo
mg = mygene.MyGeneInfo()

# Load metadata
meta_all = pd.read_excel(f"{META_XLSX}")

# Filter for Luo 2024 study
meta_luo = meta_all[meta_all['study'] == 'luo_2024']

# List of Luo 2024 sample IDs
sample_ids = meta_luo['sample_id'].unique()

# Base folder containing Luo 2024 files
base_folder = Path(f"{RAW_ROOT}/luo_2024/")

# Build samples dictionary dynamically
samples = {}
for sid in sample_ids:
    # Files are in format: GSMXXXXXX_<sample_id>_barcodes/features/matrix
    barcodes_file = list(base_folder.glob(f"*{sid}_barcodes.tsv.gz"))[0]
    features_file = list(base_folder.glob(f"*{sid}_features.tsv.gz"))[0]
    matrix_file = list(base_folder.glob(f"*{sid}_matrix.mtx.gz"))[0]
    samples[sid] = {
        "barcodes": str(barcodes_file),
        "genes": str(features_file),
        "matrix": str(matrix_file)
    }

adata_list = []

# Loop over samples
for sample_id, files in samples.items():
    print(f"\nLoading sample {sample_id}...")

    # Load matrix and transpose to cells x genes
    X = sc.read_mtx(files["matrix"]).X.T

    # Load barcodes and genes
    barcodes = pd.read_csv(files["barcodes"], header=None, sep='\t')[0].values
    genes = pd.read_csv(files["genes"], header=None, sep='\t')[0].values

    # Convert Ensembl -> gene symbols
    query_result = mg.querymany(genes, scopes='ensembl.gene', fields='symbol', species='human')
    ensembl_to_symbol = {entry['query']: entry.get('symbol', entry['query']) for entry in query_result}
    gene_symbols = [ensembl_to_symbol.get(e, e) for e in genes]

    # Create AnnData
    adata = ad.AnnData(X)
    adata.obs['barcode'] = barcodes
    adata.obs['sample_id'] = sample_id
    adata.var_names = gene_symbols
    adata.var_names_make_unique()

    # Keep raw counts
    adata.raw = adata

    # Append to list
    adata_list.append(adata)

# Make obs names unique for each sample
for adata in adata_list:
    adata.obs_names_make_unique()

# Concatenate all samples (union of genes)
adata_all = ad.concat(
    adata_list,
    join='outer',
    label=None,  # no batch column
    keys=None,
    index_unique=None
)

# Build metadata dictionary from Excel
metadata_dict = {}
for _, row in meta_luo.iterrows():
    metadata_dict[row['sample_id']] = {
        'study': row['study'],
        'patient_id': row['patient_id'],
        'sample_num': row['sample_num'],
        'treatment_status': row['treatment_status'],
        'histological_subtype': row['histological_subtype'],
        'stage': row['stage'],
        'anatomic_site': row['anatomic_site'],
        'metastatic_site': row['metastatic_site'],
        'age': row['age'],
        'treatment_response': row['treatment_response'],
        'BRCA_status': row['BRCA_status'],
        'HRD_status': row['HRD_status'],
        'TP53_status': row['TP53_status'],
        'ref': row['ref']
    }

# Add metadata columns
for col in metadata_dict[sample_ids[0]].keys():
    adata_all.obs[col] = adata_all.obs['sample_id'].map(lambda s: metadata_dict[s][col])

# Convert all metadata columns except 'barcode' and 'sample_id' to strings
for col in adata_all.obs.columns:
    if col not in ['barcode', 'sample_id']:
        adata_all.obs[col] = adata_all.obs[col].astype(str)

# Make obs_names unique
adata_all.obs_names_make_unique()

# Save finalized Luo 2024 AnnData
adata_all.write(f"{PROCESSED_ROOT}/luo_2024.h5ad")
print("Saved.")

#check all to this point
# -----------------------------
# 1) Load AnnData objects
# -----------------------------
geist = sc.read_h5ad(f"{PROCESSED_ROOT}/geistlinger_2020.h5ad")
denis = sc.read_h5ad(f"{PROCESSED_ROOT}/denisenko_2024.h5ad")
loret = sc.read_h5ad(f"{PROCESSED_ROOT}/loret_2022.h5ad")
luo = sc.read_h5ad(f"{PROCESSED_ROOT}/luo_2024.h5ad")

datasets = {
    "Geistlinger": geist,
    "Denisenko": denis,
    "Loret": loret,
    "Luo_2024": luo
}

# -----------------------------
# 2) Check X and raw
# -----------------------------
for name, ad in datasets.items():
    print(f"\n{name}:")
    print("X type:", type(ad.X))
    
    # Safe raw check — avoids errors
    raw_exists = ad.raw is not None
    raw_type_matches = raw_exists and isinstance(ad.raw.X, (np.ndarray, type(ad.X)))
    print("raw counts preserved:", raw_exists and raw_type_matches)

    print("shape (cells x genes):", ad.shape)

# -----------------------------
# 3) Gene symbols and overlaps
# -----------------------------
genes_dict = {name: set(ad.var_names) for name, ad in datasets.items()}

all_genes_union = set.union(*genes_dict.values())
all_genes_intersection = set.intersection(*genes_dict.values())

print("\n\n===== Gene Overlap Summary =====")
for name, genes in genes_dict.items():
    print(f"{name}: {len(genes)} genes")

print(f"\nGenes overlapping in ALL FOUR datasets: {len(all_genes_intersection)}")

# Pairwise gene overlaps
pairs = [
    ("Geistlinger", "Denisenko"),
    ("Geistlinger", "Loret"),
    ("Geistlinger", "Luo_2024"),
    ("Denisenko", "Loret"),
    ("Denisenko", "Luo_2024"),
    ("Loret", "Luo_2024"),
]

for a, b in pairs:
    overlap = genes_dict[a] & genes_dict[b]
    print(f"Overlap {a} & {b}: {len(overlap)}")

# -----------------------------
# 4) Metadata columns overlap
# -----------------------------
obs_cols_dict = {name: set(ad.obs.columns) for name, ad in datasets.items()}

all_metadata_union = set.union(*obs_cols_dict.values())
all_metadata_intersection = set.intersection(*obs_cols_dict.values())

print("\n\n===== Metadata Columns Summary =====")
for name, cols in obs_cols_dict.items():
    print(f"{name}: {len(cols)} columns")
    print(cols)

print("\nColumns overlapping in ALL FOUR datasets:")
print(all_metadata_intersection)

# Pairwise overlaps
for a, b in pairs:
    overlap = obs_cols_dict[a] & obs_cols_dict[b]
    only_a = obs_cols_dict[a] - obs_cols_dict[b]
    only_b = obs_cols_dict[b] - obs_cols_dict[a]

    print(f"\n===== {a} vs {b} =====")
    print("Overlap:", overlap)
    print(f"Only in {a}:", only_a)
    print(f"Only in {b}:", only_b)

import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
import numpy as np
import re

# Nath 2021 sample info
samples = {
    "P1-1": f"{RAW_ROOT}/nath_2021/GSE158722_P01.counts.txt.gz",
    "P1-2": f"{RAW_ROOT}/nath_2021/GSE158722_P02.counts.txt.gz",
    "P2-1": f"{RAW_ROOT}/nath_2021/GSE158722_P03.counts.txt.gz",
    "P2-2": f"{RAW_ROOT}/nath_2021/GSE158722_P04.counts.txt.gz",
    "P2-3": f"{RAW_ROOT}/nath_2021/GSE158722_P05.counts.txt.gz",
    "P3-1": f"{RAW_ROOT}/nath_2021/GSE158722_P06.counts.txt.gz",
    "P3-2": f"{RAW_ROOT}/nath_2021/GSE158722_P07.counts.txt.gz",
    "P3-3": f"{RAW_ROOT}/nath_2021/GSE158722_P08.counts.txt.gz",
    "P4-1": f"{RAW_ROOT}/nath_2021/GSE158722_P09.counts.txt.gz",
    "P4-2": f"{RAW_ROOT}/nath_2021/GSE158722_P10.counts.txt.gz",
    "P4-3": f"{RAW_ROOT}/nath_2021/GSE158722_P11.counts.txt.gz",
    "P5-1": f"{RAW_ROOT}/nath_2021/GSE158722_P12.counts.txt.gz",
    "P5-2": f"{RAW_ROOT}/nath_2021/GSE158722_P13.counts.txt.gz",
    "P5-3": f"{RAW_ROOT}/nath_2021/GSE158722_P14.counts.txt.gz",
    "P6-1": f"{RAW_ROOT}/nath_2021/GSE158722_P15.counts.txt.gz",
    "P6-2": f"{RAW_ROOT}/nath_2021/GSE158722_P16.counts.txt.gz",
    "P6-3": f"{RAW_ROOT}/nath_2021/GSE158722_P17.counts.txt.gz",
    "P7-1": f"{RAW_ROOT}/nath_2021/GSE158722_P18.counts.txt.gz",
    "P7-2": f"{RAW_ROOT}/nath_2021/GSE158722_P19.counts.txt.gz",
    "P8-1": f"{RAW_ROOT}/nath_2021/GSE158722_P20.counts.txt.gz",
    "P8-2": f"{RAW_ROOT}/nath_2021/GSE158722_P21.counts.txt.gz",
    "P23":  f"{RAW_ROOT}/nath_2021/GSE158722_P23.counts.txt.gz",
    "P24":  f"{RAW_ROOT}/nath_2021/GSE158722_P24.counts.txt.gz"
}

# Nath 2021 metadata table
metadata_table = pd.DataFrame([
    ["nath_2021","P1-1",63,115,"pre-treatment","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P1-2",63,116,"post-chemotherapy","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P2-1",64,117,"post-chemotherapy","Serous","","ascites","ascites","","","wildtype","","wildtype",""],
    ["nath_2021","P2-2",64,118,"post-chemotherapy","Serous","","ascites","ascites","","","wildtype","","wildtype",""],
    ["nath_2021","P2-3",64,119,"post-chemotherapy","Serous","","ascites","ascites","","","wildtype","","wildtype",""],
    ["nath_2021","P3-1",65,120,"post-chemotherapy","Serous","","ascites","ascites","","","mutated","","wildtype",""],
    ["nath_2021","P3-2",65,121,"post-chemotherapy","Serous","","ascites","ascites","","","mutated","","wildtype",""],
    ["nath_2021","P3-3",65,122,"post-chemotherapy","Serous","","ascites","ascites","","","mutated","","wildtype",""],
    ["nath_2021","P4-1",66,123,"post-chemotherapy","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P4-2",66,124,"post-chemotherapy","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P4-3",66,125,"post-chemotherapy","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P5-1",67,126,"pre-treatment","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P5-2",67,127,"post-chemotherapy","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P5-3",67,128,"post-chemotherapy","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P6-1",68,129,"post-chemotherapy","Serous","","ascites","ascites","","","wildtype","","wildtype",""],
    ["nath_2021","P6-2",68,130,"post-chemotherapy","Serous","","ascites","ascites","","","wildtype","","wildtype",""],
    ["nath_2021","P6-3",68,131,"post-chemotherapy","Serous","","ascites","ascites","","","wildtype","","wildtype",""],
    ["nath_2021","P7-1",69,132,"pre-treatment","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P7-2",69,133,"post-chemotherapy_olaparib","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P8-1",70,134,"pre-treatment","Serous","","ascites","ascites","","","wildtype","","wildtype",""],
    ["nath_2021","P8-2",70,135,"post-chemotherapy","Serous","","ascites","ascites","","","wildtype","","wildtype",""],
    ["nath_2021","P8-3",70,136,"post-chemotherapy","Serous","","ascites","ascites","","","wildtype","","wildtype",""],
    ["nath_2021","P9-1",71,137,"post-chemotherapy_niraparib","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P9-2",71,138,"post-chemotherapy_niraparib","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P9-3",71,139,"post-chemotherapy_niraparib","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P10",72,140,"pre-treatment","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P11",73,141,"pre-treatment","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P12",74,142,"pre-treatment","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P13",75,143,"pre-treatment","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P14",76,144,"pre-treatment","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P15",77,145,"pre-treatment","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P16",78,146,"pre-treatment","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P17",79,147,"pre-treatment","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P18",80,148,"post-chemotherapy","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P19",81,149,"post-chemotherapy_niraparib_pembro","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P20",82,150,"post-chemotherapy_niraparib","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P21",83,151,"post-chemotherapy","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P22",84,152,"post-chemotherapy","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P23",85,153,"post-chemotherapy_niraparib","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P24",86,154,"post-chemotherapy_pembro","Serous","","ascites","ascites","","","","","",""]
], columns=["study","sample_id","patient_id","sample_num","treatment_status","histological_subtype","stage",
            "anatomic_site","metastatic_site","age","treatment_response","BRCA_status","HRD_status","TP53_status","ref"])

# ----------------------------
# Loader
# ----------------------------
adatas = []

for sample_id, file_path in samples.items():
    print(f"\nLoading sample {sample_id}")

    # Load counts
    counts = pd.read_csv(file_path, sep="\t", index_col=0)
    
    # Drop the numeric Gene ID column if present
    if "Gene ID" in counts.columns:
        counts = counts.drop(columns=["Gene ID"])
    
    # Force numeric
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0)
    print(f"  Initial shape (genes × cells): {counts.shape}")

    # Remove all-zero genes
    nonzero_mask = counts.sum(axis=1) > 0
    dropped = (~nonzero_mask).sum()
    counts = counts.loc[nonzero_mask]
    print(f"  After zero-gene filter: {counts.shape} (dropped {dropped})")

    # Clean gene IDs (remove version suffix if any)
    genes_clean = counts.index.to_series().str.replace(r"\.\d+$", "", regex=True)

    # Transpose to cells × genes
    X = sparse.csr_matrix(counts.values.T)

    # Build AnnData object
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=counts.columns),
        var=pd.DataFrame(index=genes_clean)
    )

    # Add metadata
    adata.obs["barcode"] = adata.obs.index
    adata.obs["sample_id"] = sample_id
    adata.var_names_make_unique()
    adata.raw = adata

    adatas.append(adata)

# Concatenate all samples
adata_all = ad.concat(
    adatas,
    join="outer",
    label="sample_id",
    keys=list(samples.keys())
)

# Merge metadata
adata_all.obs = adata_all.obs.merge(metadata_table, on="sample_id", how="left")
for col in adata_all.obs.columns:
    adata_all.obs[col] = adata_all.obs[col].astype(str)
adata_all.obs_names_make_unique()

# Save
output_path = f"{PROCESSED_ROOT}/nath_2021.h5ad"
adata_all.write(output_path)
print(f"\nSaved Nath 2021 dataset to:\n{output_path}")

#olalekan
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
import mygene

# -------------------------------------------------------------------
# Samples for Olalekan 2021
# -------------------------------------------------------------------
samples = {
    "PT-1": f"{RAW_ROOT}/olalekan_2021/GSM4416534_PT-3232.csv.gz",
    "PT-2": f"{RAW_ROOT}/olalekan_2021/GSM4416535_PT-5150.csv.gz",
    "PT-3": f"{RAW_ROOT}/olalekan_2021/GSM4416536_PT-6885.csv.gz",
    "PT-4": f"{RAW_ROOT}/olalekan_2021/GSM4416537_PT-4806.csv.gz",
    "PT-5": f"{RAW_ROOT}/olalekan_2021/GSM4416538_PT-3401.csv.gz",
    "PT-6": f"{RAW_ROOT}/olalekan_2021/GSM4416539_PT-2834.csv.gz"
}

# -------------------------------------------------------------------
# Metadata table for Olalekan 2021
# EXACT column structure as your Nath metadata table
# -------------------------------------------------------------------
metadata_table = pd.DataFrame([
    ["olalekan_2021","PT-1",87,155,"post-chemotherapy","Serous","IV",
     "omentum","metastasis",62,"","","","","GRCh38"],
    ["olalekan_2021","PT-2",88,156,"pre-treatment","Serous","III",
     "omentum","metastasis",56,"","","","","GRCh38"],
    ["olalekan_2021","PT-3",89,157,"post-chemotherapy","Serous","III",
     "omentum","metastasis",66,"","","","","GRCh38"],
    ["olalekan_2021","PT-4",90,158,"pre-treatment","Serous","III",
     "omentum","metastasis",46,"","","","","GRCh38"],
    ["olalekan_2021","PT-5",91,159,"post-chemotherapy","Serous","III",
     "omentum","metastasis",71,"","","","","GRCh38"],
    ["olalekan_2021","PT-6",92,160,"post-chemotherapy","Mixed Mullerian","III",
     "omentum","metastasis",66,"","","","","GRCh38"]
], columns=[
    "study","sample_id","patient_id","sample_num","treatment_status",
    "histological_subtype","stage","anatomic_site","metastatic_site","age",
    "treatment_response","BRCA_status","HRD_status","TP53_status","ref"
])

# -------------------------------------------------------------------
# Initialize MyGeneInfo
# -------------------------------------------------------------------
mg = mygene.MyGeneInfo()

adatas = []

# -------------------------------------------------------------------
# Load each sample
# -------------------------------------------------------------------
for sample_id, file_path in samples.items():
    print(f"Loading sample {sample_id} from {file_path}...")

    # Load the matrix
    counts = pd.read_csv(file_path, sep=",", index_col=0)

    # Force numeric (avoid dtype object)
    counts = counts.apply(pd.to_numeric, errors='coerce').fillna(0)

    # Convert to sparse matrix (transpose to cells × genes)
    X = sparse.csr_matrix(counts.values.T)

    # Create AnnData
    adata = ad.AnnData(X)
    adata.obs['barcode'] = counts.columns.tolist()
    adata.obs['sample_id'] = sample_id
    adata.var_names = counts.index.tolist()

    # Save raw counts
    adata.raw = adata

    # Gene symbol conversion via MyGene (same procedure as Nath)
    genes = adata.var_names.tolist()
    query_result = mg.querymany(
        genes,
        scopes=['ensembl.gene', 'alias', 'symbol'],
        fields='symbol',
        species='human'
    )
    ensembl_to_symbol = {entry['query']: entry.get('symbol', entry['query']) for entry in query_result}

    adata.var_names = [ensembl_to_symbol.get(g, g) for g in genes]
    adata.var_names_make_unique()

    adatas.append(adata)

# -------------------------------------------------------------------
# Concatenate all samples 
# -------------------------------------------------------------------
adata_all = ad.concat(adatas, join='outer', label='sample_id', keys=[s for s in samples.keys()])

# -------------------------------------------------------------------
# Merge metadata into .obs
# -------------------------------------------------------------------
adata_all.obs = adata_all.obs.merge(metadata_table, on='sample_id', how='left')

# -------------------------------------------------------------------
# Convert metadata columns to string (except barcode, sample_id)
# -------------------------------------------------------------------
for col in adata_all.obs.columns:
    if col not in ['barcode','sample_id']:
        adata_all.obs[col] = adata_all.obs[col].astype(str)

adata_all.obs_names_make_unique()

# -------------------------------------------------------------------
# Save final dataset
# -------------------------------------------------------------------
output_path = f"{PROCESSED_ROOT}/olalekan_2021.h5ad"
adata_all.write(output_path)

print("Saved.")

#check nath and olalekan
import scanpy as sc
import numpy as np
import pandas as pd

# -----------------------------
# 1) Load AnnData objects
# -----------------------------
nath = sc.read_h5ad(f"{PROCESSED_ROOT}/nath_2021.h5ad")
olalekan = sc.read_h5ad(f"{PROCESSED_ROOT}/olalekan_2021.h5ad")

datasets = {
    "Nath_2021": nath,
    "Olalekan_2021": olalekan
}

# -----------------------------
# 2) Check X and raw
# -----------------------------
for name, ad in datasets.items():
    print(f"\n{name}:")
    print("X type:", type(ad.X))
    
    raw_exists = ad.raw is not None
    raw_type_matches = raw_exists and isinstance(ad.raw.X, (np.ndarray, type(ad.X)))
    print("raw counts preserved:", raw_exists and raw_type_matches)

    print("shape (cells x genes):", ad.shape)

# -----------------------------
# 3) Gene symbols and overlaps
# -----------------------------
genes_dict = {name: set(ad.var_names) for name, ad in datasets.items()}
all_genes_union = set.union(*genes_dict.values())
all_genes_intersection = set.intersection(*genes_dict.values())

print("\n\n===== Gene Overlap Summary =====")
for name, genes in genes_dict.items():
    print(f"{name}: {len(genes)} genes")

print(f"Genes overlapping in BOTH datasets: {len(all_genes_intersection)}")

# -----------------------------
# 4) Metadata columns overlap
# -----------------------------
obs_cols_dict = {name: set(ad.obs.columns) for name, ad in datasets.items()}
all_metadata_union = set.union(*obs_cols_dict.values())
all_metadata_intersection = set.intersection(*obs_cols_dict.values())

print("\n\n===== Metadata Columns Summary =====")
for name, cols in obs_cols_dict.items():
    print(f"{name}: {len(cols)} columns")
    print(sorted(cols))

print("\nColumns overlapping in BOTH datasets:")
print(sorted(all_metadata_intersection))

#regner
import scanpy as sc
import pandas as pd
import anndata as ad
from scipy import io
import mygene

# -------------------------------------------------------------------
# 1) Initialize MyGeneInfo
# -------------------------------------------------------------------
mg = mygene.MyGeneInfo()

# -------------------------------------------------------------------
# 2) Define Regner 2021 samples and file paths
# -------------------------------------------------------------------
samples = {
    "P8": {
        "barcodes": f"{RAW_ROOT}/regner_2021/GSM5276940_barcodes-3BAE2L.tsv.gz",
        "genes": f"{RAW_ROOT}/regner_2021/GSM5276940_features-3BAE2L.tsv.gz",
        "matrix": f"{RAW_ROOT}/regner_2021/GSM5276940_matrix-3BAE2L.mtx.gz"
    },
    "P9": {
        "barcodes": f"{RAW_ROOT}/regner_2021/GSM5276943_barcodes-3E5CFL.tsv.gz",
        "genes": f"{RAW_ROOT}/regner_2021/GSM5276943_features-3E5CFL.tsv.gz",
        "matrix": f"{RAW_ROOT}/regner_2021/GSM5276943_matrix-3E5CFL.mtx.gz"
    }
}

# -------------------------------------------------------------------
# 3) Define full metadata columns and data
# -------------------------------------------------------------------
full_metadata_cols = [
    "study","sample_id","patient_id","sample_num","treatment_status",
    "histological_subtype","stage","anatomic_site","metastatic_site",
    "age","treatment_response","BRCA_status","HRD_status","TP53_status","ref"
]

metadata_dict = {
    "P8": {
        "study": "regner_2021",
        "sample_id": "P8",
        "patient_id": 101,
        "sample_num": 174,
        "treatment_status": "pre-treatment",
        "histological_subtype": "Serous",
        "stage": "II",
        "anatomic_site": "adnexa",
        "metastatic_site": "primary",
        "age": 61,
        "treatment_response": "",
        "BRCA_status": "",
        "HRD_status": "",
        "TP53_status": "",
        "ref": "GRCh38"
    },
    "P9": {
        "study": "regner_2021",
        "sample_id": "P9",
        "patient_id": 102,
        "sample_num": 175,
        "treatment_status": "pre-treatment",
        "histological_subtype": "Serous",
        "stage": "III",
        "anatomic_site": "adnexa",
        "metastatic_site": "primary",
        "age": 59,
        "treatment_response": "",
        "BRCA_status": "",
        "HRD_status": "",
        "TP53_status": "",
        "ref": "GRCh38"
    }
}

# -------------------------------------------------------------------
# 4) Load each sample
# -------------------------------------------------------------------
adata_list = []

for sample_id, files in samples.items():
    print(f"Loading sample {sample_id}...")

    # Load matrix and transpose to cells × genes
    X = sc.read_mtx(files["matrix"]).X.T

    # Load barcodes and genes
    barcodes = pd.read_csv(files["barcodes"], header=None, sep='\t')[0].values
    genes = pd.read_csv(files["genes"], header=None, sep='\t')[0].values

    # Convert Ensembl -> gene symbols
    query_result = mg.querymany(genes, scopes='ensembl.gene', fields='symbol', species='human')
    ensembl_to_symbol = {entry['query']: entry.get('symbol', entry['query']) for entry in query_result}
    gene_symbols = [ensembl_to_symbol.get(e, e) for e in genes]

    # Create AnnData
    adata = ad.AnnData(X)
    adata.obs['barcode'] = barcodes
    adata.obs['sample_id'] = sample_id
    adata.var_names = gene_symbols
    adata.var_names_make_unique()

    # Save raw counts
    adata.raw = adata

    # Make obs_names unique now to avoid warnings
    adata.obs_names_make_unique()

    adata_list.append(adata)

# -------------------------------------------------------------------
# 5) Concatenate all samples (union of genes)
# -------------------------------------------------------------------
adata_all = ad.concat(adata_list, join='outer', label='sample_id', keys=[s for s in samples.keys()])

# -------------------------------------------------------------------
# 6) Add full metadata (all columns preserved)
# -------------------------------------------------------------------
for col in full_metadata_cols:
    adata_all.obs[col] = adata_all.obs['sample_id'].map(lambda s: metadata_dict[s].get(col, ""))

# Convert metadata columns to string (except 'barcode' and 'sample_id')
for col in adata_all.obs.columns:
    if col not in ['barcode', 'sample_id']:
        adata_all.obs[col] = adata_all.obs[col].astype(str)

# Make obs_names unique again after merging metadata
adata_all.obs_names_make_unique()

# -------------------------------------------------------------------
# 7) Save final AnnData
# -------------------------------------------------------------------
output_path = f"{PROCESSED_ROOT}/regner_2021.h5ad"
adata_all.write(output_path)
print("Saved.")

#xu
import scanpy as sc
import pandas as pd
import anndata as ad
from scipy import io
import mygene

# -------------------------------------------------------------------
# 1) Initialize MyGeneInfo
# -------------------------------------------------------------------
mg = mygene.MyGeneInfo()

# -------------------------------------------------------------------
# 2) Define Xu 2022 samples and file paths
# -------------------------------------------------------------------
samples = {
    "Cancer1": {
        "barcodes": f"{RAW_ROOT}/xu_2022/GSM5599225_Cancer1.barcodes.tsv.gz",
        "genes": f"{RAW_ROOT}/xu_2022/GSM5599225_Cancer1.genes.tsv.gz",
        "matrix": f"{RAW_ROOT}/xu_2022/GSM5599225_Cancer1.matrix.mtx.gz"
    },
    "Cancer2": {
        "barcodes": f"{RAW_ROOT}/xu_2022/GSM5599226_Cancer2.barcodes.tsv.gz",
        "genes": f"{RAW_ROOT}/xu_2022/GSM5599226_Cancer2.genes.tsv.gz",
        "matrix": f"{RAW_ROOT}/xu_2022/GSM5599226_Cancer2.matrix.mtx.gz"
    },
    "Cancer3": {
        "barcodes": f"{RAW_ROOT}/xu_2022/GSM5599227_Cancer3.barcodes.tsv.gz",
        "genes": f"{RAW_ROOT}/xu_2022/GSM5599227_Cancer3.genes.tsv.gz",
        "matrix": f"{RAW_ROOT}/xu_2022/GSM5599227_Cancer3.matrix.mtx.gz"
    },
    "Cancer4": {
        "barcodes": f"{RAW_ROOT}/xu_2022/GSM5599228_Cancer4.barcodes.tsv.gz",
        "genes": f"{RAW_ROOT}/xu_2022/GSM5599228_Cancer4.genes.tsv.gz",
        "matrix": f"{RAW_ROOT}/xu_2022/GSM5599228_Cancer4.matrix.mtx.gz"
    },
    "Cancer5": {
        "barcodes": f"{RAW_ROOT}/xu_2022/GSM5599229_Cancer5.barcodes.tsv.gz",
        "genes": f"{RAW_ROOT}/xu_2022/GSM5599229_Cancer5.genes.tsv.gz",
        "matrix": f"{RAW_ROOT}/xu_2022/GSM5599229_Cancer5.matrix.mtx.gz"
    },
    "Cancer6": {
        "barcodes": f"{RAW_ROOT}/xu_2022/GSM5599230_Cancer6.barcodes.tsv.gz",
        "genes": f"{RAW_ROOT}/xu_2022/GSM5599230_Cancer6.genes.tsv.gz",
        "matrix": f"{RAW_ROOT}/xu_2022/GSM5599230_Cancer6.matrix.mtx.gz"
    },
    "Cancer7": {
        "barcodes": f"{RAW_ROOT}/xu_2022/GSM5599231_Cancer7.barcodes.tsv.gz",
        "genes": f"{RAW_ROOT}/xu_2022/GSM5599231_Cancer7.genes.tsv.gz",
        "matrix": f"{RAW_ROOT}/xu_2022/GSM5599231_Cancer7.matrix.mtx.gz"
    }
}

# -------------------------------------------------------------------
# 3) Define full metadata columns and data
# -------------------------------------------------------------------
full_metadata_cols = [
    "study","sample_id","patient_id","sample_num","treatment_status",
    "histological_subtype","stage","anatomic_site","metastatic_site",
    "age","treatment_response","BRCA_status","HRD_status","TP53_status","ref"
]

metadata_dict = {
    "Cancer1": {"study":"xu_2022","sample_id":"Cancer1","patient_id":103,"sample_num":176,"treatment_status":"pre-treatment","histological_subtype":"Serous","stage":"III","anatomic_site":"adnexa","metastatic_site":"primary","age":50,"treatment_response":"","BRCA_status":"mutated","HRD_status":"","TP53_status":"","ref":"GRCh38"},
    "Cancer2": {"study":"xu_2022","sample_id":"Cancer2","patient_id":104,"sample_num":177,"treatment_status":"pre-treatment","histological_subtype":"Serous","stage":"II","anatomic_site":"adnexa","metastatic_site":"primary","age":51,"treatment_response":"","BRCA_status":"wildtype","HRD_status":"","TP53_status":"","ref":"GRCh38"},
    "Cancer3": {"study":"xu_2022","sample_id":"Cancer3","patient_id":105,"sample_num":178,"treatment_status":"pre-treatment","histological_subtype":"Serous","stage":"I","anatomic_site":"adnexa","metastatic_site":"primary","age":41,"treatment_response":"","BRCA_status":"wildtype","HRD_status":"","TP53_status":"","ref":"GRCh38"},
    "Cancer4": {"study":"xu_2022","sample_id":"Cancer4","patient_id":106,"sample_num":179,"treatment_status":"pre-treatment","histological_subtype":"Serous","stage":"I","anatomic_site":"adnexa","metastatic_site":"primary","age":47,"treatment_response":"","BRCA_status":"wildtype","HRD_status":"","TP53_status":"","ref":"GRCh38"},
    "Cancer5": {"study":"xu_2022","sample_id":"Cancer5","patient_id":107,"sample_num":180,"treatment_status":"pre-treatment","histological_subtype":"Serous","stage":"II","anatomic_site":"adnexa","metastatic_site":"primary","age":57,"treatment_response":"","BRCA_status":"mutated","HRD_status":"","TP53_status":"","ref":"GRCh38"},
    "Cancer6": {"study":"xu_2022","sample_id":"Cancer6","patient_id":108,"sample_num":181,"treatment_status":"pre-treatment","histological_subtype":"Serous","stage":"III","anatomic_site":"adnexa","metastatic_site":"primary","age":48,"treatment_response":"","BRCA_status":"wildtype","HRD_status":"","TP53_status":"","ref":"GRCh38"},
    "Cancer7": {"study":"xu_2022","sample_id":"Cancer7","patient_id":109,"sample_num":182,"treatment_status":"pre-treatment","histological_subtype":"Serous","stage":"I","anatomic_site":"adnexa","metastatic_site":"primary","age":53,"treatment_response":"","BRCA_status":"wildtype","HRD_status":"","TP53_status":"","ref":"GRCh38"}
}

# -------------------------------------------------------------------
# 4) Load each sample
# -------------------------------------------------------------------
adata_list = []

for sample_id, files in samples.items():
    print(f"Loading sample {sample_id}...")

    # Load matrix and transpose to cells × genes
    X = sc.read_mtx(files["matrix"]).X.T

    # Load barcodes and genes
    barcodes = pd.read_csv(files["barcodes"], header=None, sep='\t')[0].values
    genes = pd.read_csv(files["genes"], header=None, sep='\t')[0].values

    # Convert Ensembl -> gene symbols
    query_result = mg.querymany(genes, scopes='ensembl.gene', fields='symbol', species='human')
    ensembl_to_symbol = {entry['query']: entry.get('symbol', entry['query']) for entry in query_result}
    gene_symbols = [ensembl_to_symbol.get(e, e) for e in genes]

    # Create AnnData
    adata = ad.AnnData(X)
    adata.obs['barcode'] = barcodes
    adata.obs['sample_id'] = sample_id
    adata.var_names = gene_symbols
    adata.var_names_make_unique()

    # Save raw counts
    adata.raw = adata

    # Make obs_names unique per sample
    adata.obs_names_make_unique()

    adata_list.append(adata)

# -------------------------------------------------------------------
# 5) Concatenate all samples (union of genes)
# -------------------------------------------------------------------
adata_all = ad.concat(adata_list, join='outer', label='sample_id', keys=[s for s in samples.keys()])

# -------------------------------------------------------------------
# 6) Add full metadata (all columns preserved)
# -------------------------------------------------------------------
for col in full_metadata_cols:
    adata_all.obs[col] = adata_all.obs['sample_id'].map(lambda s: metadata_dict[s].get(col, ""))

# Convert metadata columns to string (except 'barcode', 'sample_id')
for col in adata_all.obs.columns:
    if col not in ['barcode', 'sample_id']:
        adata_all.obs[col] = adata_all.obs[col].astype(str)

# Make obs_names unique again after merging metadata
adata_all.obs_names_make_unique()

# -------------------------------------------------------------------
# 7) Save final AnnData
# -------------------------------------------------------------------
output_path = f"{PROCESSED_ROOT}/xu_2022.h5ad"
adata_all.write(output_path)
print("Saved Xu 2022 concatenated AnnData successfully.")

#check nath olalekan regner and xu 
import scanpy as sc
import numpy as np
import pandas as pd

# -----------------------------
# 1) Load AnnData objects
# -----------------------------
nath = sc.read_h5ad(f"{PROCESSED_ROOT}/nath_2021.h5ad")
olalekan = sc.read_h5ad(f"{PROCESSED_ROOT}/olalekan_2021.h5ad")
regner = sc.read_h5ad(f"{PROCESSED_ROOT}/regner_2021.h5ad")
xu = sc.read_h5ad(f"{PROCESSED_ROOT}/xu_2022.h5ad")

datasets = {
    "Nath_2021": nath,
    "Olalekan_2021": olalekan,
    "Regner_2021": regner,
    "Xu_2022": xu
}

# -----------------------------
# 2) Check X and raw
# -----------------------------
for name, ad in datasets.items():
    print(f"\n{name}:")
    print("X type:", type(ad.X))
    
    raw_exists = ad.raw is not None
    raw_type_matches = raw_exists and isinstance(ad.raw.X, (np.ndarray, type(ad.X)))
    print("raw counts preserved:", raw_exists and raw_type_matches)

    print("shape (cells x genes):", ad.shape)

# -----------------------------
# 3) Gene symbols and overlaps
# -----------------------------
genes_dict = {name: set(ad.var_names) for name, ad in datasets.items()}
all_genes_union = set.union(*genes_dict.values())
all_genes_intersection = set.intersection(*genes_dict.values())

print("\n\n===== Gene Overlap Summary =====")
for name, genes in genes_dict.items():
    print(f"{name}: {len(genes)} genes")
print(f"Genes overlapping in ALL FOUR datasets: {len(all_genes_intersection)}")

# -----------------------------
# 4) Metadata columns overlap
# -----------------------------
obs_cols_dict = {name: set(ad.obs.columns) for name, ad in datasets.items()}
all_metadata_union = set.union(*obs_cols_dict.values())
all_metadata_intersection = set.intersection(*obs_cols_dict.values())

print("\n\n===== Metadata Columns Summary =====")
for name, cols in obs_cols_dict.items():
    print(f"{name}: {len(cols)} columns")
    print(sorted(cols))

print("\nColumns overlapping in ALL FOUR datasets:")
print(sorted(all_metadata_intersection))

#olbrecht
import scanpy as sc

# Load Olbrecht dataset
adata = sc.read_h5ad(f"{RAW_ROOT}/olbrecht_2021/olbrecht_2021_raw_counts.h5ad")

# -----------------------------
# Remove extra metadata columns
# -----------------------------
columns_to_drop = ['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'label_short', 'ident']
adata.obs = adata.obs.drop(columns=[c for c in columns_to_drop if c in adata.obs.columns])

# -----------------------------
# Keep raw counts
# -----------------------------
adata.raw = adata

# -----------------------------
# Make obs_names unique
# -----------------------------
adata.obs_names_make_unique()

# -----------------------------
# Ensure all metadata columns are strings (except barcode, sample_id)
# -----------------------------
for col in adata.obs.columns:
    if col not in ['barcode', 'sample_id']:
        adata.obs[col] = adata.obs[col].astype(str)

# -----------------------------
# Save final AnnData
# -----------------------------
output_path = f"{PROCESSED_ROOT}/olbrecht_2021.h5ad"
adata.write(output_path)

print("Saved.")

#zhang
# -----------------------------
# Load full counts matrix
# -----------------------------
umi_file = f"{RAW_ROOT}/zhang_2022/GSE165897_UMIcounts_HGSOC.tsv.gz"
df = pd.read_csv(umi_file, sep="\t", index_col=0)

# Extract sample identifiers from column names
# Assumes format: CELLBARCODE-SAMPLEID
sample_ids = [c.split('-', 1)[1] for c in df.columns]

# Add sample IDs as a row for easier splitting
df.columns = pd.MultiIndex.from_arrays([sample_ids, df.columns], names=['sample_id','barcode'])

# Unique sample IDs
unique_samples = sorted(set(sample_ids))

# -----------------------------
# Manually define metadata table
# -----------------------------
metadata_table = pd.DataFrame([
    ["zhang_2022","EOC1005_pPer",110,183,"pre-treatment","Serous","IV","peritoneum","metastasis",73,"","","","", "GRCh38"],
    ["zhang_2022","EOC1005_iTum2",110,184,"post-treatment","Serous","IV","adnexa","primary",73,"","","","", "GRCh38"],
    ["zhang_2022","EOC136_pMes1",111,185,"pre-treatment","Serous","IV","omentum","metastasis",64,"","","","", "GRCh38"],
    ["zhang_2022","EOC136_iOme",111,186,"post-treatment","Serous","IV","omentum","metastasis",64,"","","","", "GRCh38"],
    ["zhang_2022","EOC153_pOme",112,187,"pre-treatment","Serous","IV","omentum","metastasis",78,"","","","", "GRCh38"],
    ["zhang_2022","EOC153_iOme1",112,188,"post-treatment","Serous","IV","omentum","metastasis",78,"","","","", "GRCh38"],
    ["zhang_2022","EOC227_pOme1",113,189,"pre-treatment","Serous","IV","omentum","metastasis",74,"","","","", "GRCh38"],
    ["zhang_2022","EOC227_iOme1",113,190,"post-treatment","Serous","IV","omentum","metastasis",74,"","","","", "GRCh38"],
    ["zhang_2022","EOC3_pPer1",114,191,"pre-treatment","Serous","IV","peritoneum","metastasis",67,"","","","", "GRCh38"],
    ["zhang_2022","EOC3_iOme2",114,192,"post-treatment","Serous","IV","omentum","metastasis",67,"","","","", "GRCh38"],
    ["zhang_2022","EOC349_pPer2",115,193,"pre-treatment","Serous","IV","peritoneum","metastasis",67,"","","","", "GRCh38"],
    ["zhang_2022","EOC349_iOme1",115,194,"post-treatment","Serous","IV","omentum","metastasis",67,"","","","", "GRCh38"],
    ["zhang_2022","EOC372_pPer",116,195,"pre-treatment","Serous","III","peritoneum","metastasis",68,"","","","", "GRCh38"],
    ["zhang_2022","EOC372_iPer",116,196,"post-treatment","Serous","III","peritoneum","metastasis",68,"","","","", "GRCh38"],
    ["zhang_2022","EOC443_pOme",117,197,"pre-treatment","Serous","IV","omentum","metastasis",54,"","","","", "GRCh38"],
    ["zhang_2022","EOC443_iOme1",117,198,"post-treatment","Serous","IV","omentum","metastasis",54,"","","","", "GRCh38"],
    ["zhang_2022","EOC540_p2Ome",118,199,"pre-treatment","Serous","IV","omentum","metastasis",62,"","","","", "GRCh38"],
    ["zhang_2022","EOC540_iOme2",118,200,"post-treatment","Serous","IV","omentum","metastasis",62,"","","","", "GRCh38"],
    ["zhang_2022","EOC733_pPer",119,201,"pre-treatment","Serous","IV","peritoneum","metastasis",72,"","","","", "GRCh38"],
    ["zhang_2022","EOC733_iOme",119,202,"post-treatment","Serous","IV","omentum","metastasis",72,"","","","", "GRCh38"],
    ["zhang_2022","EOC87_pPer1_2",120,203,"pre-treatment","Serous","III","peritoneum","metastasis",62,"","","","", "GRCh38"],
    ["zhang_2022","EOC87_iOme1",120,204,"post-treatment","Serous","III","omentum","metastasis",62,"","","","", "GRCh38"]
], columns=[
    "study","sample_id","patient_id","sample_num","treatment_status","histological_subtype",
    "stage","anatomic_site","metastatic_site","age","treatment_response","BRCA_status","HRD_status",
    "TP53_status","ref"
])


# -----------------------------
# Create per-sample AnnData objects
# -----------------------------
adatas = []

for sample_id in unique_samples:
    print(f"Processing sample {sample_id}...")

    # Select columns for this sample
    sample_cols = df.loc[:, sample_id]
    
    # Convert to sparse matrix
    X = sparse.csr_matrix(sample_cols.values)
    
    # Create AnnData
    adata = ad.AnnData(X.T)  # transpose to cells x genes
    adata.obs['barcode'] = sample_cols.columns
    adata.obs['sample_id'] = sample_id
    adata.var_names = sample_cols.index.tolist()
    
    # Save raw counts
    adata.raw = adata
    
    adatas.append(adata)

# -----------------------------
# Concatenate all samples
# -----------------------------
adata_all = ad.concat(adatas, join='outer', label='sample_id', keys=unique_samples)

# -----------------------------
# Merge metadata
# -----------------------------
adata_all.obs = adata_all.obs.merge(metadata_table, on='sample_id', how='left')

# Ensure metadata columns are strings except barcode and sample_id
for col in adata_all.obs.columns:
    if col not in ['barcode', 'sample_id']:
        adata_all.obs[col] = adata_all.obs[col].astype(str)

# Make obs_names unique
adata_all.obs_names_make_unique()

# -----------------------------
# Save final AnnData
# -----------------------------
output_path = f"{PROCESSED_ROOT}/zhang_2022.h5ad"
adata_all.write(output_path)

print("Saved.")

#zheng
import scanpy as sc
import pandas as pd

adata = sc.read_h5ad(f"{RAW_ROOT}/zheng_2023/zheng_2023_raw_counts.h5ad")

metadata_table = pd.DataFrame([
    ["zheng_2023","HGSOC1_PT",121,205,"pre-treatment","Serous","III","adnexa","primary",57,"sensitive","","","", "GRCh38"],
    ["zheng_2023","HGSOC1_MT",121,206,"pre-treatment","Serous","III","omentum","metastasis",57,"sensitive","","","", "GRCh38"],
    ["zheng_2023","HGSOC1_AS",121,207,"pre-treatment","Serous","III","ascites","ascites",57,"sensitive","","","", "GRCh38"],
    ["zheng_2023","HGSOC2_PT",122,208,"pre-treatment","Serous","II","adnexa","primary",43,"sensitive","","","", "GRCh38"],
    ["zheng_2023","HGSOC2_AS",122,209,"pre-treatment","Serous","II","ascites","ascites",43,"sensitive","","","", "GRCh38"],
    ["zheng_2023","HGSOC3_PT",123,210,"pre-treatment","Serous","III","adnexa","primary",68,"resistant","","","", "GRCh38"],
    ["zheng_2023","HGSOC3_MT",123,211,"pre-treatment","Serous","III","omentum","metastasis",68,"resistant","","","", "GRCh38"],
    ["zheng_2023","HGSOC3_AS",123,212,"pre-treatment","Serous","III","ascites","ascites",68,"resistant","","","", "GRCh38"],
    ["zheng_2023","HGSOC4_PT",124,213,"pre-treatment","Serous","III","adnexa","primary",71,"sensitive","","","", "GRCh38"],
    ["zheng_2023","HGSOC4_MT",124,214,"pre-treatment","Serous","III","omentum","metastasis",71,"sensitive","","","", "GRCh38"],
    ["zheng_2023","HGSOC5_PT",125,215,"pre-treatment","Serous","III","adnexa","primary",50,"sensitive","","","", "GRCh38"],
    ["zheng_2023","HGSOC5_AS",125,216,"pre-treatment","Serous","III","ascites","ascites",50,"sensitive","","","", "GRCh38"],
    ["zheng_2023","HGSOC6_PT",126,217,"pre-treatment","Serous","III","adnexa","primary",67,"resistant","","","", "GRCh38"],
    ["zheng_2023","HGSOC6_MT",126,218,"pre-treatment","Serous","III","omentum","metastasis",67,"resistant","","","", "GRCh38"],
    ["zheng_2023","HGSOC6_AS",126,219,"pre-treatment","Serous","III","ascites","ascites",67,"resistant","","","", "GRCh38"],
    ["zheng_2023","HGSOC7_PT",127,220,"pre-treatment","Serous","III","adnexa","primary",47,"resistant","","","", "GRCh38"],
    ["zheng_2023","HGSOC8_PT",128,221,"pre-treatment","Serous","III","adnexa","primary",63,"sensitive","","","", "GRCh38"],
    ["zheng_2023","HGSOC8_AS",128,222,"pre-treatment","Serous","III","ascites","ascites",63,"sensitive","","","", "GRCh38"],
    ["zheng_2023","HGSOC9_PT",129,223,"pre-treatment","Serous","III","adnexa","primary",48,"sensitive","","","", "GRCh38"],
    ["zheng_2023","HGSOC9_AS",129,224,"pre-treatment","Serous","III","ascites","ascites",48,"sensitive","","","", "GRCh38"],
    ["zheng_2023","HGSOC10_PT",130,225,"pre-treatment","Serous","III","adnexa","primary",66,"sensitive","","","", "GRCh38"],
    ["zheng_2023","HGSOC10_AS",130,226,"pre-treatment","Serous","III","ascites","ascites",66,"sensitive","","","", "GRCh38"]
], columns=["study","sample_id","patient_id","sample_num","treatment_status","histological_subtype",
            "stage","anatomic_site","metastatic_site","age","treatment_response","BRCA_status","HRD_status","TP53_status","ref"])

# -----------------------------
# 3. Subset AnnData to keep only metadata-matching cells
# -----------------------------
adata = adata[adata.obs['Samples'].isin(metadata_table['sample_id']), :].copy()

# -----------------------------
# 4. Merge metadata into adata.obs
# -----------------------------
adata.obs = adata.obs.merge(metadata_table, left_on='Samples', right_on='sample_id', how='left')

# -----------------------------
# 5. Remove unnecessary columns for compatibility
# -----------------------------
columns_to_keep = [
    "barcode", "study", "sample_id", "patient_id", "sample_num",
    "treatment_status", "histological_subtype", "stage", "anatomic_site",
    "metastatic_site", "age", "treatment_response", "BRCA_status",
    "HRD_status", "TP53_status", "ref"
]


# If 'barcode' does not exist yet, create it from obs_names
if 'barcode' not in adata.obs.columns:
    adata.obs['barcode'] = adata.obs_names

# Keep only compatible columns
adata.obs = adata.obs[columns_to_keep]

# -----------------------------
# 6. Ensure raw counts
# -----------------------------
adata.raw = adata

# -----------------------------
# 7. Make unique obs and var names
# -----------------------------
adata.obs_names_make_unique()
adata.var_names_make_unique()

# -----------------------------
# 8. Save cleaned AnnData
# -----------------------------
output_path = f"{PROCESSED_ROOT}/zheng_2023.h5ad"
adata.write(output_path)
print("Saved.")

#check all after nath
import scanpy as sc
import numpy as np
import pandas as pd

# -----------------------------
# 1) Load AnnData objects
# -----------------------------
nath = sc.read_h5ad(f"{PROCESSED_ROOT}/nath_2021.h5ad")
olalekan = sc.read_h5ad(f"{PROCESSED_ROOT}/olalekan_2021.h5ad")
regner = sc.read_h5ad(f"{PROCESSED_ROOT}/regner_2021.h5ad")
xu = sc.read_h5ad(f"{PROCESSED_ROOT}/xu_2022.h5ad")
olbrecht = sc.read_h5ad(f"{PROCESSED_ROOT}/olbrecht_2021.h5ad")
zheng = sc.read_h5ad(f"{PROCESSED_ROOT}/zheng_2023.h5ad")
zhang = sc.read_h5ad(f"{PROCESSED_ROOT}/zhang_2022.h5ad")

datasets = {
    "Nath_2021": nath,
    "Olalekan_2021": olalekan,
    "Regner_2021": regner,
    "Xu_2022": xu,
    "Olbrecht_2021": olbrecht,
    "Zheng_2023": zheng,
    "Zhang_2022": zhang
}

# -----------------------------
# 2) Check X and raw
# -----------------------------
print("\n===== Raw Counts & X Check =====")
for name, ad in datasets.items():
    print(f"\n{name}:")
    print("X type:", type(ad.X))
    
    raw_exists = ad.raw is not None
    raw_type_matches = raw_exists and isinstance(ad.raw.X, (np.ndarray, type(ad.X)))
    print("raw counts preserved:", raw_exists and raw_type_matches)

    print("shape (cells x genes):", ad.shape)

# -----------------------------
# 3) Gene symbols and overlaps
# -----------------------------
genes_dict = {name: set(ad.var_names) for name, ad in datasets.items()}
all_genes_union = set.union(*genes_dict.values())
all_genes_intersection = set.intersection(*genes_dict.values())

print("\n\n===== Gene Overlap Summary =====")
for name, genes in genes_dict.items():
    print(f"{name}: {len(genes)} genes")
print(f"Genes overlapping in ALL SEVEN datasets: {len(all_genes_intersection)}")

# -----------------------------
# 4) Metadata columns overlap
# -----------------------------
obs_cols_dict = {name: set(ad.obs.columns) for name, ad in datasets.items()}
all_metadata_union = set.union(*obs_cols_dict.values())
all_metadata_intersection = set.intersection(*obs_cols_dict.values())

print("\n\n===== Metadata Columns Summary =====")
for name, cols in obs_cols_dict.items():
    print(f"{name}: {len(cols)} columns")
    print(sorted(cols))

print("\nColumns overlapping in ALL SEVEN datasets:")
print(sorted(all_metadata_intersection))

#vazquez_garcia inspect
import h5py

file_path = f"{RAW_ROOT}/vazquez_garcia_2022/GSE180661_matrix.h5"

with h5py.File(file_path, 'r') as f:
    print("Keys at root:", list(f.keys()))
    
    # Inspect the 'obs' dataset
    print("obs dtype:", f['obs'].dtype)
    print("obs shape:", f['obs'].shape)
    
    # Inspect the 'var' dataset
    print("var dtype:", f['var'].dtype)
    print("var shape:", f['var'].shape)
    
    # Inspect first row of obs
    print("First obs entry:", f['obs'][0])
    
    # Inspect first row of var
    print("First var entry:", f['var'][0])

#vazquez_garcia inspect
import scanpy as sc

file_path = f"{RAW_ROOT}/vazquez_garcia_2022/GSE180661_matrix.h5"

# Load the file
adata = sc.read(file_path)

# Inspect the data
print(adata)  # summary: # cells x # genes, layers, etc.

# View first few rows of obs (cell metadata)
print(adata.obs.head())

# View first few rows of var (gene metadata)
print(adata.var.head())

# Inspect first few entries of the matrix
print(adata.X[:5, :5])

#load in metadata
import pandas as pd

adata.obs['sample_id'] = adata.obs.index.str.rsplit('_', n=1).str[0]

# Load your metadata
meta_all = pd.read_excel(f"{META_XLSX}")

# Extract sample_id from barcode in adata
adata.obs['sample_id'] = adata.obs.index.str.rsplit('_', n=1).str[0]

# Filter metadata to only include the samples present in adata
meta_subset = meta_all[meta_all['sample_id'].isin(adata.obs['sample_id'])]

# Merge metadata into adata.obs
adata.obs = adata.obs.merge(meta_subset, on='sample_id', how='left')

# Optional: reorder columns or fill missing values if needed
adata.obs.head()

#integrate metadata to vazquez_garcia
barcode_to_sample = {bc: bc.rsplit('_', 1)[0] for bc in adata.obs_names}
meta_obs = pd.DataFrame(index=adata.obs_names)

for col in ['study', 'sample_id', 'patient_id', 'sample_num', 'treatment_status', 
            'histological_subtype', 'stage', 'anatomic_site', 'metastatic_site', 
            'age', 'treatment_response', 'BRCA_status', 'HRD_status', 'TP53_status', 'ref']:
    # Pull metadata for each barcode, convert to string, fill missing with empty string
    meta_obs[col] = [str(meta_subset.loc[meta_subset['sample_id'] == barcode_to_sample[bc], col].values[0])
                     if not meta_subset.loc[meta_subset['sample_id'] == barcode_to_sample[bc], col].empty 
                     else '' 
                     for bc in adata.obs_names]

# Replace the AnnData obs
adata.obs = meta_obs

# Ensure var names are unique
adata.var_names_make_unique()

# Now write safely
save_path = f"{PROCESSED_ROOT}/vazquez_garcia_2022.h5ad"
adata.write(save_path)

#vazquez_garcia
import scanpy as sc
import pandas as pd

# -----------------------------
# Paths
# -----------------------------
file_path = f"{RAW_ROOT}/vazquez_garcia_2022/GSE180661_matrix.h5"
metadata_path = f"{META_XLSX}"
save_path = f"{PROCESSED_ROOT}/vazquez_garcia_2022.h5ad"

# -----------------------------
# Load the data
# -----------------------------
adata = sc.read(file_path)

# Quick inspection
print(adata)
print(adata.obs.head())
print(adata.var.head())
print(adata.X[:5, :5])

# -----------------------------
# Extract sample_id from barcode
# -----------------------------
adata.obs['sample_id'] = adata.obs.index.str.rsplit('_', n=1).str[0]

# -----------------------------
# Load metadata
# -----------------------------
meta_all = pd.read_excel(metadata_path)

# Subset to only the samples present in adata
meta_subset = meta_all[meta_all['sample_id'].isin(adata.obs['sample_id'])]

# Map each barcode to its sample_id
barcode_to_sample = {bc: bc.rsplit('_', 1)[0] for bc in adata.obs_names}

# -----------------------------
# Build obs with all 15 metadata columns
# -----------------------------
metadata_columns = [
    'study', 'sample_id', 'patient_id', 'sample_num', 'treatment_status', 
    'histological_subtype', 'stage', 'anatomic_site', 'metastatic_site', 
    'age', 'treatment_response', 'BRCA_status', 'HRD_status', 'TP53_status', 'ref'
]

meta_obs = pd.DataFrame(index=adata.obs_names)

for col in metadata_columns:
    values = []
    for bc in adata.obs_names:
        sample_id = barcode_to_sample[bc]
        subset = meta_subset.loc[meta_subset['sample_id'] == sample_id, col]
        if not subset.empty:
            values.append(str(subset.values[0]))
        else:
            values.append('')
    meta_obs[col] = values

# Replace adata.obs with the new metadata
adata.obs = meta_obs

# -----------------------------
# Preserve raw counts
# -----------------------------
adata.raw = adata.copy()

# -----------------------------
# Ensure var names are unique
# -----------------------------
adata.var_names_make_unique()

# -----------------------------
# Save processed dataset
# -----------------------------
adata.write(save_path)
print(f"Processed dataset saved to {save_path}")

# -----------------------------
# Optional: Verify obs metadata
# -----------------------------
print("Metadata summary per column:")
for col in metadata_columns:
    print(f"{col}: {adata.obs[col].nunique()} unique values, {adata.obs[col].isna().sum()} missing")

#check all
import scanpy as sc
import numpy as np
import pandas as pd

# -----------------------------
# 1) Dataset paths
# -----------------------------
paths = {
    "Denisenko_2024": f"{PROCESSED_ROOT}/denisenko_2024.h5ad",
    "Geistlinger_2020": f"{PROCESSED_ROOT}/geistlinger_2020.h5ad",
    "Loret_2022": f"{PROCESSED_ROOT}/loret_2022.h5ad",
    "Luo_2024": f"{PROCESSED_ROOT}/luo_2024.h5ad",
    "Nath_2021": f"{PROCESSED_ROOT}/nath_2021.h5ad",
    "Olalekan_2021": f"{PROCESSED_ROOT}/olalekan_2021.h5ad",
    "Olbrecht_2021": f"{PROCESSED_ROOT}/olbrecht_2021.h5ad",
    "Regner_2021": f"{PROCESSED_ROOT}/regner_2021.h5ad",
    "Vazquez_Garcia_2022": f"{PROCESSED_ROOT}/vazquez_garcia_2022.h5ad",
    "Xu_2022": f"{PROCESSED_ROOT}/xu_2022.h5ad",
    "Zhang_2022": f"{PROCESSED_ROOT}/zhang_2022.h5ad",
    "Zheng_2023": f"{PROCESSED_ROOT}/zheng_2023.h5ad"
}

datasets = {}
zero_count_cells = {}

# -----------------------------
# 2) Load datasets & check X/raw
# -----------------------------
for name, path in paths.items():
    ad = sc.read_h5ad(path)
    datasets[name] = ad
    
    # X type
    print(f"\n{name}:")
    print("Shape (cells x genes):", ad.shape)
    print("X type:", type(ad.X))
    
    # Raw counts preserved
    raw_exists = ad.raw is not None
    raw_type_matches = raw_exists and isinstance(ad.raw.X, (np.ndarray, type(ad.X)))
    print("Raw counts preserved:", raw_exists and raw_type_matches)
    
    # Zero-count cells
    zero_cells = np.sum(ad.X.sum(axis=1) == 0)
    zero_count_cells[name] = zero_cells
    print(f"Zero-count cells: {zero_cells} / {ad.n_obs}")

# -----------------------------
# 3) Gene symbols and overlaps
# -----------------------------
genes_dict = {name: set(ad.var_names) for name, ad in datasets.items()}
all_genes_union = set.union(*genes_dict.values())
all_genes_intersection = set.intersection(*genes_dict.values())

print("\n\n===== Gene Overlap Summary =====")
for name, genes in genes_dict.items():
    print(f"{name}: {len(genes)} genes")

print(f"\nGenes overlapping in ALL datasets: {len(all_genes_intersection)}")

# Pairwise overlaps
from itertools import combinations
for a, b in combinations(datasets.keys(), 2):
    overlap = genes_dict[a] & genes_dict[b]
    print(f"Overlap {a} & {b}: {len(overlap)}")

# -----------------------------
# 4) Metadata columns overlap
# -----------------------------
obs_cols_dict = {name: set(ad.obs.columns) for name, ad in datasets.items()}
all_metadata_union = set.union(*obs_cols_dict.values())
all_metadata_intersection = set.intersection(*obs_cols_dict.values())

print("\n\n===== Metadata Columns Summary =====")
for name, cols in obs_cols_dict.items():
    print(f"{name}: {len(cols)} columns")
    print(cols)

print("\nColumns overlapping in ALL datasets:")
print(all_metadata_intersection)

# Pairwise metadata overlaps
for a, b in combinations(datasets.keys(), 2):
    overlap = obs_cols_dict[a] & obs_cols_dict[b]
    only_a = obs_cols_dict[a] - obs_cols_dict[b]
    only_b = obs_cols_dict[b] - obs_cols_dict[a]

    print(f"\n===== {a} vs {b} =====")
    print("Overlap:", overlap)
    print(f"Only in {a}:", only_a)
    print(f"Only in {b}:", only_b)

# -----------------------------
# 5) Zero-count cells summary
# -----------------------------
print("\n\n===== Zero-count Cells Summary =====")
for name, zero_cells in zero_count_cells.items():
    print(f"{name}: {zero_cells} / {datasets[name].n_obs}")

#new nath 
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import numpy as np

# Nath 2021 sample info
samples = {
    "P1-1": f"{RAW_ROOT}/nath_2021/GSE158722_P01.counts.txt.gz",
    "P1-2": f"{RAW_ROOT}/nath_2021/GSE158722_P02.counts.txt.gz",
    "P2-1": f"{RAW_ROOT}/nath_2021/GSE158722_P03.counts.txt.gz",
    "P2-2": f"{RAW_ROOT}/nath_2021/GSE158722_P04.counts.txt.gz",
    "P2-3": f"{RAW_ROOT}/nath_2021/GSE158722_P05.counts.txt.gz",
    "P3-1": f"{RAW_ROOT}/nath_2021/GSE158722_P06.counts.txt.gz",
    "P3-2": f"{RAW_ROOT}/nath_2021/GSE158722_P07.counts.txt.gz",
    "P3-3": f"{RAW_ROOT}/nath_2021/GSE158722_P08.counts.txt.gz",
    "P4-1": f"{RAW_ROOT}/nath_2021/GSE158722_P09.counts.txt.gz",
    "P4-2": f"{RAW_ROOT}/nath_2021/GSE158722_P10.counts.txt.gz",
    "P4-3": f"{RAW_ROOT}/nath_2021/GSE158722_P11.counts.txt.gz",
    "P5-1": f"{RAW_ROOT}/nath_2021/GSE158722_P12.counts.txt.gz",
    "P5-2": f"{RAW_ROOT}/nath_2021/GSE158722_P13.counts.txt.gz",
    "P5-3": f"{RAW_ROOT}/nath_2021/GSE158722_P14.counts.txt.gz",
    "P6-1": f"{RAW_ROOT}/nath_2021/GSE158722_P15.counts.txt.gz",
    "P6-2": f"{RAW_ROOT}/nath_2021/GSE158722_P16.counts.txt.gz",
    "P6-3": f"{RAW_ROOT}/nath_2021/GSE158722_P17.counts.txt.gz",
    "P7-1": f"{RAW_ROOT}/nath_2021/GSE158722_P18.counts.txt.gz",
    "P7-2": f"{RAW_ROOT}/nath_2021/GSE158722_P19.counts.txt.gz",
    "P8-1": f"{RAW_ROOT}/nath_2021/GSE158722_P20.counts.txt.gz",
    "P8-2": f"{RAW_ROOT}/nath_2021/GSE158722_P21.counts.txt.gz",
    "P23":  f"{RAW_ROOT}/nath_2021/GSE158722_P23.counts.txt.gz",
    "P24":  f"{RAW_ROOT}/nath_2021/GSE158722_P24.counts.txt.gz"
}

# Nath 2021 metadata table
metadata_table = pd.DataFrame([
    ["nath_2021","P1-1",63,115,"pre-treatment","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P1-2",63,116,"post-chemotherapy","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P2-1",64,117,"post-chemotherapy","Serous","","ascites","ascites","","","wildtype","","wildtype",""],
    ["nath_2021","P2-2",64,118,"post-chemotherapy","Serous","","ascites","ascites","","","wildtype","","wildtype",""],
    ["nath_2021","P2-3",64,119,"post-chemotherapy","Serous","","ascites","ascites","","","wildtype","","wildtype",""],
    ["nath_2021","P3-1",65,120,"post-chemotherapy","Serous","","ascites","ascites","","","mutated","","wildtype",""],
    ["nath_2021","P3-2",65,121,"post-chemotherapy","Serous","","ascites","ascites","","","mutated","","wildtype",""],
    ["nath_2021","P3-3",65,122,"post-chemotherapy","Serous","","ascites","ascites","","","mutated","","wildtype",""],
    ["nath_2021","P4-1",66,123,"post-chemotherapy","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P4-2",66,124,"post-chemotherapy","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P4-3",66,125,"post-chemotherapy","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P5-1",67,126,"pre-treatment","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P5-2",67,127,"post-chemotherapy","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P5-3",67,128,"post-chemotherapy","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P6-1",68,129,"post-chemotherapy","Serous","","ascites","ascites","","","wildtype","","wildtype",""],
    ["nath_2021","P6-2",68,130,"post-chemotherapy","Serous","","ascites","ascites","","","wildtype","","wildtype",""],
    ["nath_2021","P6-3",68,131,"post-chemotherapy","Serous","","ascites","ascites","","","wildtype","","wildtype",""],
    ["nath_2021","P7-1",69,132,"pre-treatment","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P7-2",69,133,"post-chemotherapy_olaparib","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P8-1",70,134,"pre-treatment","Serous","","ascites","ascites","","","wildtype","","wildtype",""],
    ["nath_2021","P8-2",70,135,"post-chemotherapy","Serous","","ascites","ascites","","","wildtype","","wildtype",""],
    ["nath_2021","P8-3",70,136,"post-chemotherapy","Serous","","ascites","ascites","","","wildtype","","wildtype",""],
    ["nath_2021","P9-1",71,137,"post-chemotherapy_niraparib","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P9-2",71,138,"post-chemotherapy_niraparib","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P9-3",71,139,"post-chemotherapy_niraparib","Serous","","ascites","ascites","","","mutated","","mutated",""],
    ["nath_2021","P10",72,140,"pre-treatment","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P11",73,141,"pre-treatment","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P12",74,142,"pre-treatment","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P13",75,143,"pre-treatment","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P14",76,144,"pre-treatment","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P15",77,145,"pre-treatment","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P16",78,146,"pre-treatment","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P17",79,147,"pre-treatment","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P18",80,148,"post-chemotherapy","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P19",81,149,"post-chemotherapy_niraparib_pembro","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P20",82,150,"post-chemotherapy_niraparib","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P21",83,151,"post-chemotherapy","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P22",84,152,"post-chemotherapy","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P23",85,153,"post-chemotherapy_niraparib","Serous","","ascites","ascites","","","","","",""],
    ["nath_2021","P24",86,154,"post-chemotherapy_pembro","Serous","","ascites","ascites","","","","","",""]
], columns=["study","sample_id","patient_id","sample_num","treatment_status","histological_subtype","stage",
            "anatomic_site","metastatic_site","age","treatment_response","BRCA_status","HRD_status","TP53_status","ref"])



# -----------------------------
# Loader for one Nath sample
# -----------------------------
def load_nath_sample(sample_id: str, file_path: str) -> ad.AnnData:
    print(f"\nLoading {sample_id}: {file_path}", flush=True)

    # Load raw table (do not set index_col initially)
    df = pd.read_csv(file_path, sep="\t")
    if "Gene Symbol" not in df.columns:
        raise ValueError(
            f"[{sample_id}] Expected 'Gene Symbol' column missing. "
            f"Columns: {df.columns.tolist()[:20]}"
        )

    # Set gene symbol as index
    df = df.set_index("Gene Symbol")

    # Drop NON-expression annotation columns if present
    drop_cols = [c for c in df.columns if c.lower() in {"gene id", "gene_id", "ensembl id", "ensembl_id"}]
    if drop_cols:
        df = df.drop(columns=drop_cols)

    # Force numeric expression only
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)

    # Remove all-zero genes
    df = df.loc[df.sum(axis=1) > 0]

    # Build sparse AnnData (cells × genes)
    X = sp.csr_matrix(df.values.T)

    # Use raw barcodes as initial obs index
    obs = pd.DataFrame(index=df.columns.astype(str))
    var = pd.DataFrame(index=df.index.astype(str))

    adata = ad.AnnData(X=X, obs=obs, var=var)

    # Minimal obs/var hygiene
    adata.obs["sample_id"] = str(sample_id)
    adata.obs["barcode"] = adata.obs_names.astype(str)

    # Make obs_names globally unique across ALL samples BEFORE concatenation
    adata.obs_names = (adata.obs["sample_id"] + ":" + adata.obs["barcode"]).values
    adata.obs_names_make_unique()

    # Make genes unique if needed
    adata.var_names_make_unique()

    # Sanity checks
    print(f"  -> genes × cells (input): {df.shape[0]} × {df.shape[1]}")
    print(
        f"  -> AnnData: {adata.n_obs} cells × {adata.n_vars} genes | "
        f"X dtype={adata.X.dtype} | sparse={sp.issparse(adata.X)}"
    )

    for g in ["ACTB", "GAPDH", "B2M"]:
        if g in adata.var_names:
            j = adata.var_names.get_loc(g)
            nnz = adata.X[:, j].nnz
            print(f"     {g}: {nnz} / {adata.n_obs} cells nonzero")

    return adata

# -----------------------------
# Load all samples
# -----------------------------
adatas = []
failed = []

for sid, fpath in samples.items():
    try:
        adatas.append(load_nath_sample(sid, fpath))
    except Exception as e:
        failed.append((sid, str(e)))
        print(f"!! FAILED {sid}: {e}", flush=True)

print(f"\nLoaded {len(adatas)} / {len(samples)} samples.", flush=True)
if failed:
    print("\nFailures:", flush=True)
    for sid, err in failed:
        print(f" - {sid}: {err}", flush=True)


# -----------------------------
# CONCATENATE ALL SAMPLES
# -----------------------------
# Keep sample_id as a column from concat as well (label=...)
# Use outer join to keep union of genes across samples
adata_all = ad.concat(
    adatas,
    join="outer",
    label="sample_id",
    keys=[a.obs["sample_id"].iloc[0] for a in adatas],
)

print("\nConcatenation complete.")
print(adata_all)


# -----------------------------
# MERGE METADATA
# -----------------------------

obs = adata_all.obs.copy()

# Guard: if any duplicate column names slipped in, drop duplicates (keeps first)
obs = obs.loc[:, ~obs.columns.duplicated()].copy()

# Merge metadata on sample_id (safe)
obs = obs.merge(metadata_table, on="sample_id", how="left", suffixes=("", "_meta"))

# Ensure obs index matches AnnData obs_names (keep your existing unique IDs)
obs.index = adata_all.obs_names.astype(str)

# Convert metadata to string-only values (h5ad-safe)
for col in obs.columns:
    obs[col] = obs[col].astype(str).fillna("NA")

adata_all.obs = obs
adata_all.obs_names_make_unique()

# -----------------------------
# FINAL SANITY + SAVE
# -----------------------------
print("\nFinal AnnData summary:")
print(adata_all)
print("X dtype:", adata_all.X.dtype)
print("Sparse:", sp.issparse(adata_all.X))

save_path = f"{PROCESSED_ROOT}/nath_2021.h5ad"
adata_all.write(save_path)
print(f"\nSaved processed Nath 2021 dataset (raw counts in X) to: {save_path}", flush=True)


obs = adata_all.obs.copy()

# Guard: if any duplicate column names slipped in, drop duplicates (keeps first)
obs = obs.loc[:, ~obs.columns.duplicated()].copy()

# Merge metadata on sample_id (safe)
obs = obs.merge(metadata_table, on="sample_id", how="left", suffixes=("", "_meta"))

# Ensure obs index matches AnnData obs_names (keep your existing unique IDs)
obs.index = adata_all.obs_names.astype(str)

# Convert metadata to string-only values (h5ad-safe)
for col in obs.columns:
    obs[col] = obs[col].astype(str).fillna("NA")

adata_all.obs = obs
adata_all.obs_names_make_unique()

# -----------------------------
# FINAL SANITY + SAVE
# -----------------------------
print("\nFinal AnnData summary:")
print(adata_all)
print("X dtype:", adata_all.X.dtype)
print("Sparse:", sp.issparse(adata_all.X))

save_path = f"{PROCESSED_ROOT}/nath_2021.h5ad"
adata_all.write(save_path)
print(f"\nSaved processed Nath 2021 dataset (raw counts in X) to: {save_path}", flush=True)

#hornburg
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from pathlib import Path

# 1) Load metadata
meta_all = pd.read_excel(f"{META_XLSX}")
meta_h = meta_all[meta_all["study"] == "hornburg_2021"].copy()
print(f"Metadata rows: {len(meta_h)}, patients: {meta_h['patient_id'].nunique()}")

# 2) Build file-path mapping — scan EGAF dirs, deduplicate by filename, keep only metadata matches
base = Path(path("data_root", "2026_final_atlas", "hornburg"))
all_files = sorted(base.glob("EGAF*/RNA-*.txt"))
all_files = [f for f in all_files if not f.name.endswith(".md5")]

# Deduplicate by filename (keep first occurrence)
seen = {}
for f in all_files:
    if f.name not in seen:
        seen[f.name] = f

meta_ids = set(meta_h["sample_id"].tolist())
samples = {name: path for name, path in seen.items() if name in meta_ids}
print(f"Files matched to metadata: {len(samples)}")

# 3) Load each count matrix
adatas = []
for sample_id, file_path in sorted(samples.items()):
    print(f"Loading {sample_id} ...")
    counts = pd.read_csv(file_path, sep="\t", index_col=0)
    # genes x cells -> transpose to cells x genes
    X = sp.csr_matrix(counts.T.values)
    adata = ad.AnnData(X)
    adata.obs_names = counts.columns.str.strip('"').tolist()
    adata.var_names = counts.index.str.strip('"').tolist()
    adata.var_names_make_unique()
    adata.obs["barcode"] = adata.obs_names.tolist()
    adata.obs["sample_id"] = sample_id
    adata.raw = adata
    adatas.append(adata)
    print(f"  {adata.shape[0]} cells x {adata.shape[1]} genes")

# 4) Concatenate
adata_all = ad.concat(adatas, join="outer")
adata_all.obs_names_make_unique()
print(f"\nConcatenated: {adata_all.shape}")

# 5) Merge metadata
adata_all.obs = adata_all.obs.reset_index(drop=True)
adata_all.obs = adata_all.obs.merge(meta_h, on="sample_id", how="left")
for col in adata_all.obs.columns:
    if col not in ["barcode", "sample_id"]:
        adata_all.obs[col] = adata_all.obs[col].astype(str)
adata_all.obs_names_make_unique()

# 6) Save
out_path = f"{PROCESSED_ROOT}/hornburg_2021.h5ad"
adata_all.write(out_path)
print(f"\nSaved to {out_path}")

# 7) Verify
print(f"Shape: {adata_all.shape}")
print(f"X type: {type(adata_all.X)}")
print(f"Raw: {adata_all.raw is not None}")
print(f"NaN in study: {adata_all.obs['study'].isna().sum()}")
print(f"Zero-count cells: {(adata_all.X.sum(axis=1) == 0).sum()}")
print(f"\nSample ID value counts:\n{adata_all.obs['sample_id'].value_counts()}")
print(f"\nStudy values: {adata_all.obs['study'].unique()}")
