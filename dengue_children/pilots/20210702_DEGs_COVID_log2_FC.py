
import h5py
import numpy as np
import pandas as pd
import anndata
import scipy.sparse as sp
import scanpy as sc
import random
import scipy.sparse
from scipy.stats import ks_2samp

print('select data of sick dengue kids')
fn_h5ad = '/home/yike/phd/dengue/data/mergedata_20210519.h5ad'
adata_D = anndata.read_h5ad(fn_h5ad)
adata_D = (adata_D[adata_D.obs['dataset'] == 'child'])[(adata_D[adata_D.obs['dataset'] == 'child']).obs['sick'] =='sick']

print('select data of sick COVID patients')
adata_C = sc.read_h5ad('/home/yike/phd/dengue/data/dataset_from_google/COVID PBMC/50__subtype__ID/adata_filter.h5ad')
adata_C = adata_C[adata_C.obs['sick'] == 'sick']

sc.pp.normalize_total(adata_D, target_sum=1e6) #normalize data to CPM (counts per million)
sc.pp.normalize_total(adata_C, target_sum=1e6) 

genes_D = adata_D.var_names
genes_C = adata_C.var_names
genes = [gene for gene in genes_D if gene in genes_C] 
############### get the average expression and log2_fold_change between COVID and dengue
cell_types = ['Monocytes', 'Plasmablasts', 'B_cells', 'T_cells', 'NK_cells', 'cDCs', 'pDCs']

m = len(genes)
ress = pd.DataFrame([])
for cell_type in cell_types:
    adatag_D = adata_D[adata_D.obs['cell_type'] == cell_type][:, genes]
    adatag_C = adata_C[adata_C.obs['cell_type'] == cell_type][:, genes]
    X_D = adatag_D.X
    X_C = adatag_C.X
    avg_D = X_D.mean(axis=0)
    avg_C = X_C.mean(axis=0)
    if scipy.sparse.issparse(X_D):
        avg_D = np.asarray(avg_D).reshape(-1)
    if scipy.sparse.issparse(X_C):
        avg_C = np.asarray(avg_C).reshape(-1)
    res = pd.DataFrame([], index=adatag_D.var_names)
    
    res['avg_D'] = avg_D
    res['avg_C'] = avg_C
    
    # Compute log2 fold changes
    log2_fc = np.log2(avg_C + 0.1) - np.log2(avg_D + 0.1)
    res['log2_fold_change'] = log2_fc
    
    res['cell_type'] = cell_type
    
    ress = pd.concat([ress, res])
    
ress.to_csv('/home/yike/phd/dengue/data/tables/DEGs_vs_COVID.tsv')

############### get p value from ks test between COVID and dengue
#ress = pd.read_csv('/home/yike/phd/dengue/data/tables/DEGs_vs_COVID.tsv', index_col=[0, 'cell_type'])

for cell_type in cell_types:
    adatag_D = adata_D[adata_D.obs['cell_type'] == cell_type][:, genes]
    adatag_C = adata_C[adata_C.obs['cell_type'] == cell_type][:, genes]#    X_D = adatag_D.X
    X_C = adatag_C.X
    for i, gene in enumerate(genes):
        data_D = X_D[:, i].toarray()[:, 0]
        data_C = X_C[:, i].toarray()[:, 0]
        res = ks_2samp(data_D, data_C, alternative='two-sided', mode='auto')
        ress.loc[gene, cell_type]['statistic'] = res[0]
        ress.loc[gene, cell_type]['pvalue'] = res[1]
ress.to_csv('/home/yike/phd/dengue/data/tables/DEGs_vs_COVID.tsv')