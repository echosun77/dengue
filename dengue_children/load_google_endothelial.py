################################### common dataset

%config Completer.use_jedi = False

import numpy as np
import pandas as pd

import anndata
import scanpy as sc

data = pd.read_csv('/home/yike/phd/dengue/data/dataset_from_google/human_lung/mouse/mouse_droplet_TMS_UMIs.csv', index_col=0)
dataT = data.T
matrix = dataT.values

obsname = pd.read_csv('/home/yike/phd/dengue/data/dataset_from_google//human_lung/mouse/mouse_droplet_TMS_metadata.csv', index_col=0)
obsname.index.name = None

varname = pd.DataFrame([], index=dataT.columns)
varname.index.name = None

adata_SS2 = anndata.AnnData(X=matrix, obs=obsname, var=varname)

save_file = '/home/yike/phd/dengue/data/dataset_from_google/human_lung/mouse/mouse_droplet.h5ad'
adata_SS2.write(save_file)

adata_droplet = sc.read_h5ad('/home/yike/phd/dengue/data/dataset_from_google/human_lung/mouse/mouse_droplet.h5ad')
adata_facs = sc.read_h5ad('/home/yike/phd/dengue/data/dataset_from_google/human_lung/mouse/mouse_facs.h5ad')
adata = anndata.concat([adata_droplet, adata_facs], join='outer')
adata.write('/home/yike/phd/dengue/data/dataset_from_google/human_lung/mouse/mouse_lung.h5ad')






######################################## eg., load human lung, pick some cells
%config Completer.use_jedi = False

import numpy as np
import pandas as pd

import anndata
import scanpy as sc

obsname = pd.read_csv('/home/yike/phd/dengue/data/dataset_from_google//human_lung/human/krasnow_hlca_10x_metadata.csv', index_col=0)
obsname.index.name = None

obs_immune = obsname[obsname['compartment'] == 'Immune']
obs_immune['free_annotation'].value_counts()

number = {
    'Macrophage': 500,
    'Natural Killer': 500,
    'CD4+ Memory/Effector T': 500,
    'Classical Monocyte': 500,
    'CD8+ Naive T': 500,
    'CD8+ Memory/Effector T': 300,
    'CD4+ Naive T': 300,
    'Nonclassical Monocyte': 200,
    'B': 200,
    'Basophil/Mast 2': 200,
    'Natural Killer T': 200,
    'IGSF21+ Dendritic': 100,
    'Myeloid Dendritic Type 2': 100, 
    'Proliferating Macrophage': 100,
    'OLR1+ Classical Monocyte': 100,
    'Intermediate Monocyte': 100,
    'Plasma': 100,
    'TREM2+ Dendritic': 75,
    'EREG+ Dendritic': 75,
    'Plasmacytoid Dendritic': 75,
    'Myeloid Dendritic Type 1': 75,
    'Proliferating NK/T': 60,
    'Platelet/Megakaryocyte': 40, 
    }

dic = {}
for key in obs_immune['free_annotation'].value_counts().keys():
    dic[key] = []
 
cell_list = []
for i, (cell_name, row) in enumerate(obsname['free_annotation'].items()):
    if row in number.keys():
        n = number[row]
        if len(dic[row]) < n:
           dic[row].append(i)
           cell_list.append(i)
    else:
        pass

from collections import defaultdict
data = defaultdict(list)

for key in dic.keys():
    cell_ID = dic[key]
    index = []
    with open('/home/yike/phd/dengue/data/dataset_from_google/human_lung/human/krasnow_hlca_10x_UMIs.csv', 'r') as f:
        columns = np.array(f.readline().split(','))[cell_ID].tolist()
        for line in f:
            fill = line.split(',')
            index.append(fill[0])
            cell_ID = dic[key]
            data[key].append(np.array(fill)[1:][cell_ID])
        data[key] = pd.DataFrame(data[key], columns=columns, index=index).astype(np.float32)

data_immune = pd.DataFrame([])
for key in data.keys():
    data_immune = pd.concat([data_immune, data[key]], axis=1)

data_immune.to_csv('/home/yike/phd/dengue/data/dataset_from_google/human_lung/human/krasnow_hlca_10x_UMIs_immune.csv')

dataT = data_immune.T
matrix = dataT.values

obs_immunename = obsname.iloc[cell_list]
obs_immunename.index.name = None

varname = pd.DataFrame([], index=dataT.columns)
varname.index.name = None

adata_10X_immune = anndata.AnnData(X=matrix, obs=obs_immunename, var=varname)

save_file = '/home/yike/phd/dengue/data/dataset_from_google/human_lung/human/human_lung_10X_immune.h5ad'
adata_10X_immune.write(save_file)

adata_10X_endo = sc.read_h5ad('/home/yike/phd/dengue/data/dataset_from_google/human_lung/human/human_lung_10X_endo.h5ad')
adata_10X_immune = sc.read_h5ad('/home/yike/phd/dengue/data/dataset_from_google/human_lung/human/human_lung_10X_immune.h5ad')
adata_10X_rest = sc.read_h5ad('/home/yike/phd/dengue/data/dataset_from_google/human_lung/human/human_lung_10X_rest.h5ad')

adata = anndata.concat([adata_10X_endo, adata_10X_endo, adata_10X_rest], join='outer')
adata.write('/home/yike/phd/dengue/data/dataset_from_google/human_lung/human/human_lung_10X.h5ad')


########################################### load human cell atlas, pick specific cells
import os
import gzip

fdn = "../data/dataset_from_google/human_cell_landscape/"
fdn_counts = fdn+"dge_rmbatch_data/"
fdn_metadata = fdn+"annotation_rmbatch_data_revised417/"
fns = os.listdir(fdn_metadata)
data = []
fns = [x for x in fns if 'Adult' in x]
fns = fns[41:]
for ifn, fn in enumerate(fns):
    print('Sample', ifn+1, 'out of', len(fns))
    fn_meta = fdn_metadata+fn
    fn_full = fdn_counts+fn.replace('-', '').split('_')[0]+'.rmbatchdge.txt.gz'
    if fn.replace('-', '').split('_')[0]+'.rmbatchdge.txt.gz' not in os.listdir(fdn_counts):
        continue
    print('Reading metadata')
    meta = pd.read_csv(fn_meta, sep=',', index_col=3) # should be index_col=3
    print('Finding endothelial cells')
    #idx = meta['Celltype'].isin(cell_types_whitelist)
    #idx = meta.index[idx]
    #if len(idx) < 50:
    #    continue
    with gzip.open(fn_full, 'rt') as f:
        print('Getting cell names')
        cellnames = f.readline().rstrip('\n').split()####### first line  # .rstrip() 删右边， .lstrip() 删左边 # first line == cellnames
        #cellnames = np.array([x.strip('"') for x in cellnames])[idx] # .strip('"') 去掉字符串中的"
        cellnames = np.array([x.strip('"') for x in cellnames])
        counts = []
        genes = []
        print('Getting counts')
        for i, line in enumerate(f):
            if (i % 100) == 0:
                print(f'Gene {i}', end='\r')
            fields = line.rstrip('\n').split()
            gene = fields[0].strip('"')
            #counts_gene = np.array(fields[1:]).astype(np.int64)[idx]
            counts_gene = np.array(fields[1:]).astype(np.int64)
            genes.append(gene)
            counts.append(counts_gene)

    counts = np.array(counts)
    counts = pd.DataFrame(counts, index=genes, columns=cellnames).T
    
    ###############################
    counts = counts.loc[meta.index]
    ###############################
    #idx = counts.index[counts['PECAM1'] > 1] 
    if 'PECAM1' not in counts.columns:
        continue
    idx = counts['PECAM1'] > 1
    idx &= (counts > 0).sum(axis=1) >= 250 
     
    meta = meta.loc[idx]
    counts = counts.loc[idx]
    data.append({
        'sample': fn.split('.')[0],
        'counts': counts,
        'metadata': meta,
    })PE
    #if ifn == 10:
        #break
        
###########################################################################
meta_all = pd.concat([x['metadata'] for x in data], join='outer')
#meta_all.to_csv(fdn+'meta.tsv')
data_all = pd.concat([x['counts'] for x in data], join='outer').fillna(0)
meta_all = meta_all.loc[data_all.index]
adata = anndata.AnnData(
    X=data_all.values, # perhaps .T
    obs=meta_all,
    var=pd.DataFrame([], index=data_all.columns) # perhaps data_all.T.columns
)
adata.write(fdn+'adata_0520_3.h5ad')
#adata1 = sc.read_h5ad(fdn+'adata_0520_1.h5ad')
#adata2 = sc.read_h5ad(fdn+'adata_0520_2.h5ad')
#adata3 = sc.read_h5ad(fdn+'adata_0520_3.h5ad')
#adata = anndata.concat([adata1, adata2, adata3], join='outer')
####################################################
fdn = "../data/dataset_from_google/human_cell_landscape/"
adata = sc.read_h5ad(fdn+'adata_YK.h5ad')

sc.pp.filter_cells(adata, min_genes=250) 
sc.pp.filter_genes(adata, min_cells=3) 
# normalize adata
sc.pp.normalize_total(adata, target_sum=1e6) # 1e4
# logarithmize the data
sc.pp.log1p(adata)
print('Number of genes: {:d}'.format(adata.n_vars))

mt_gene = adata.var_names.str.startswith('MT-')
hla_gene = adata.var_names.str.startswith('HLA-')
rp_gene = adata.var_names.str.startswith('RPS') | adata.var_names.str.startswith('RPL') 

bad_genes = mt_gene | hla_gene| rp_gene
good_genes = ~bad_genes

adata = adata[:, good_genes]
print('Number of genes after cell filter: {:d}'.format(adata.n_vars))

sc.pp.highly_variable_genes(adata, min_mean=0.5, max_mean=5, min_disp=1)
sc.pl.highly_variable_genes(adata)
print('get {:d} highly variable cells'.format(adata.var.highly_variable.sum()))

sc.tl.pca(adata, use_highly_variable=True, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['PECAM1', 'CDH5', 'GJA5', 'VWF'])

import leidenalg
#sc.tl.leiden(adata, resolution=0.00015, partition_type=leidenalg.CPMVertexPartition, key_added='leiden_0.00015')
#sc.tl.leiden(adata, resolution=0.0002, partition_type=leidenalg.CPMVertexPartition, key_added='leiden_0.0002')
#sc.tl.leiden(adata, resolution=0.0001, partition_type=leidenalg.CPMVertexPartition, key_added='leiden_0.0001')
sc.tl.leiden(adata, resolution=0.00005, partition_type=leidenalg.CPMVertexPartition, key_added='leiden_0.00005')

sc.pl.umap(adata, color=['leiden_0.0001', ])

sc.tl.rank_genes_groups(adata, 'leiden_0.0001', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

marker_genes = ['CLDN5', 'PECAM1', 'CDH5', 'ERG', # pan-endothelial cells
    'GJA5', # arterial ECS
    'ACKR1', 'VWF', # Venous ECs 
    'PDPN', 'LYVE1', # Lymphatic ECs
    'ICAM1', 'VCAM1', 'SELE', 'ACKR3', 'ACKR4',
    'PTPRC', 'CD3D', 'MS4A1', 'CD14', 'FCGR3A', 
                 'EPCAM', 'TFF3'
                ]
fig,ax = plt.subplots(figsize=[4, 2.5], dpi=150)
sc.pl.dotplot(adata, marker_genes, groupby='leiden_0.0001', use_raw=False, ax=ax)
############################################################

