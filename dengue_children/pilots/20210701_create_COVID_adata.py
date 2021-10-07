
import h5py
import numpy as np
import pandas as pd
import anndata
import scipy.sparse as sp
import scanpy as sc
import random
from collections import defaultdict

############### create the dataset by randomly selecting certain quantity cells from each cell subtype
inna = anndata.read_h5ad('/home/yike/phd/dengue/data/dataset_from_google/COVID PBMC/innate.h5ad')

a = inna.obs['cell_type'].value_counts().tolist()
b = [int(i/10) for i in a]
cells_n = {ct: n for ct, n in zip(inna.obs['cell_type'].value_counts().index.tolist(), b)}

cells_n['platelet'] = 1000
cells_n['conventional dendritic cell'] = 900
cells_n['plasmacytoid dendritic cell'] = 700
cells_n['platelet'] = 1000
cells_n['granulocyte'] = 500

cells = []
for key in cells_n.keys():
    cells.append(random.sample(inna[inna.obs['cell_type'] == key].obs_names.tolist(), int(cells_n[key])))

cell_list = []
for cell in cells:
    cell_list = cell_list + cell

in_adata = inna[inna.obs_names.isin(cell_list)]
in_adata.obs = in_adata.obs[['severity', 'timepoint', 'outcome', 'days_since_onset', 'sex', 'age', 'donar', 'cell_type']]
in_adata.obs.columns = ['severity', 'timepoint', 'outcome', 'days_since_onset', 'sex', 'age', 'donar', 'cell_subtype']

in_adata.obs['cell_type'] = in_adata.obs['cell_subtype'].replace({
    'non-classical monocyte': 'Monocyte',
     'classical monocyte': 'Monocyte',
     'NK_CD16hi': 'NK_cells',
     'plasmacytoid dendritic cell': 'pDCs',
     'NK_CD56loCD16lo': 'NK_cells',
     'platelet': 'Platelets',
     'NK_CD56hiCD16lo': 'NK_cells',
     'conventional dendritic cell': 'cDCs',
     'granulocyte': 'Granulocytes',
     'intermediate monocyte': 'Monocytes'
})

in_adata.write('/home/yike/phd/dengue/data/dataset_from_google/COVID PBMC/innate_scaled_filter.h5ad')

################################################################
adaptive = anndata.read_h5ad('/home/yike/phd/dengue/data/dataset_from_google/COVID PBMC/adaptive.h5ad')

a = adaptive.obs['cell_type'].value_counts().tolist()
b = [int(i/10) for i in a]
cells_n = {ct: n for ct, n in zip(adaptive.obs['cell_type'].value_counts().index.tolist(), b)}

cells_n['plasmablast'] = 1000

cells = []
for key in cells_n.keys():
    cells.append(random.sample(adaptive[adaptive.obs['cell_type'] == key].obs_names.tolist(), int(cells_n[key])))

cell_list = []
for cell in cells:
    cell_list = cell_list + cell    

ad_adata = adaptive[adaptive.obs_names.isin(cell_list)]

ad_adata.obs = ad_adata.obs[['severity', 'timepoint', 'outcome', 'days_since_onset', 'sex', 'age', 'cell_type']]
ad_adata.obs.columns = ['severity', 'timepoint', 'outcome', 'days_since_onset', 'sex', 'age', 'cell_subtype']

ad_adata.obs['cell_type'] = ad_adata.obs['cell_subtype'].replace({
    'naive B cell': 'B_cells',
     'plasmablast': 'Plasmablasts',
     'CD4-positive, alpha-beta memory T cell': 'T_cells',
     'memory B cell': 'B_cells',
     'regulatory T cell': 'T_cells',
     'gamma-delta T cell': 'T_cells',
     'naive CD8+ T cell': 'T_cells',
     'TissueResMemT': 'TissueResMemT',
     'naive CD4+ T cell': 'T_cells',
     'CD8-positive, alpha-beta memory T cell': 'T_cells',
     'mucosal invariant T cell (MAIT)': 'T_cells',
     'double-positive T cell (DPT)': 'T_cells',
     'double negative T cell (DNT)': 'T_cells',
     'TCRVbeta13.1pos': 'TCRVbeta13.1pos'
})

ad_adata.write('/home/yike/phd/dengue/data/dataset_from_google/COVID PBMC/adaptive_scaled_filter.h5ad')
######################################################

adata_C = anndata.concat([in_adata, ad_adata])
adata_C.write('/home/yike/phd/dengue/data/dataset_from_google/COVID PBMC/adata_scaled_filter.h5ad')

############### create the dataset by randomly selecting 50 cells from each cell subtype of each patient, if the number is less than 50, then select all cells
inna = anndata.read_h5ad('/home/yike/phd/dengue/data/dataset_from_google/COVID PBMC/innate.h5ad')

cells = []
for ID in inna.obs['donor'].unique().tolist():
    inna_ID = inna[inna.obs['donor'] == ID]
    for ct in inna.obs['cell_type'].unique().tolist():
        n_cells = len(inna_ID[inna_ID.obs['cell_type'] == ct].obs_names.tolist())
        if n_cells >= 50:
            cells.append(random.sample(inna_ID[inna_ID.obs['cell_type'] == ct].obs_names.tolist(), 50))
        else:
            cells.append(inna_ID[inna_ID.obs['cell_type'] == ct].obs_names.tolist())

cell_list = []
for cell in cells:
    cell_list = cell_list + cell

in_adata = inna[inna.obs_names.isin(cell_list)]
in_adata.obs = in_adata.obs.rename(columns={'cell_type':'cell_subtype'})

in_adata.obs['cell_type'] = in_adata.obs['cell_subtype'].replace({
    'non-classical monocyte': 'Monocytes',
     'classical monocyte': 'Monocytes',
     'NK_CD16hi': 'NK_cells',
     'plasmacytoid dendritic cell': 'pDCs',
     'NK_CD56loCD16lo': 'NK_cells',
     'platelet': 'Platelets',
     'NK_CD56hiCD16lo': 'NK_cells',
     'conventional dendritic cell': 'cDCs',
     'granulocyte': 'Granulocytes',
     'intermediate monocyte': 'Monocytes'
})

in_adata.write('/home/yike/phd/dengue/data/dataset_from_google/COVID PBMC/50__subtype__ID/innate_filter.h5ad')

#########################################
adaptive = anndata.read_h5ad('/home/yike/phd/dengue/data/dataset_from_google/COVID PBMC/adaptive.h5ad')

cells = []
for ID in adaptive.obs['donor'].unique().tolist():
    adaptive_ID = adaptive[adaptive.obs['donor'] == ID]
    for ct in adaptive.obs['cell_type'].unique().tolist():
        n_cells = len(adaptive_ID[adaptive_ID.obs['cell_type'] == ct].obs_names.tolist())
        if n_cells >= 50:
            cells.append(random.sample(adaptive_ID[adaptive_ID.obs['cell_type'] == ct].obs_names.tolist(), 50))
        else:
            cells.append(adaptive_ID[adaptive_ID.obs['cell_type'] == ct].obs_names.tolist())

cell_list = []
for cell in cells:
    cell_list = cell_list + cell

ad_adata = adaptive[adaptive.obs_names.isin(cell_list)]
ad_adata.obs = ad_adata.obs.rename(columns={'cell_type':'cell_subtype'})

ad_adata.obs['cell_type'] = ad_adata.obs['cell_subtype'].replace({
    'naive B cell': 'B_cells',
     'plasmablast': 'Plasmablasts',
     'CD4-positive, alpha-beta memory T cell': 'T_cells',
     'memory B cell': 'B_cells',
     'regulatory T cell': 'T_cells',
     'gamma-delta T cell': 'T_cells',
     'naive CD8+ T cell': 'T_cells',
     'TissueResMemT': 'T_cells',
     'naive CD4+ T cell': 'T_cells',
     'CD8-positive, alpha-beta memory T cell': 'T_cells',
     'mucosal invariant T cell (MAIT)': 'T_cells',
     'double-positive T cell (DPT)': 'T_cells',
     'double negative T cell (DNT)': 'T_cells',
     'TCRVbeta13.1pos': 'TCRVbeta13.1pos'
})

ad_adata.write('/home/yike/phd/dengue/data/dataset_from_google/COVID PBMC/50__subtype__ID/adaptive_filter.h5ad')

#################################
adata_C = anndata.concat([in_adata, ad_adata])
adata_C.write('/home/yike/phd/dengue/data/dataset_from_google/COVID PBMC/50__subtype__ID/adata_filter.h5ad')

################################ How to transfer h5 file to adata
fdn = '/home/yike/phd/dengue/data/dataset_from_google/COVID PBMC/GSE161918_RAW/GSM4929081_B1_10xlane1_RNA_filtered_feature_bc_matrix.h5'
f = h5py.File(fdn,'r') 

for key in ['barcodes', 'data', 'features', 'indices', 'indptr', 'shape']:
    if key == 'features':
        for k in f['matrix'][key].keys():
            print(f['matrix'][key][k].shape)
    else:
        print(f['matrix'][key].shape)
        
data = f['matrix']['data'][:]
indices = f['matrix']['indices'][:]
indptr = f['matrix']['indptr'][:]

X = sp.csr_matrix((data, indices, indptr))

obs = pd.DataFrame([], index = f['matrix']['barcodes'][:])
var = pd.DataFrame([], index = f['matrix']['features']['name'][:])

obs.index = obs.index.str.decode('utf8')
var.index = obs.index.str.decode('utf8')

adata = anndata.AnnData(X=X, obs=obs, var=var)
adata.var_names_make_unique