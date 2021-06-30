%config Completer.use_jedi = False
import os

import numpy as np
import pandas as pd

import anndata
import scanpy as sc

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib import gridspec
import matplotlib as mpl
import seaborn as sns

from collections import defaultdict
############################################################
print('Load high-quality cells only')
fn_h5ad = '/home/yike/phd/dengue/data/mergedata_20210519.h5ad'
adata = anndata.read_h5ad(fn_h5ad)
adata = adata[adata.obs['cell_quality'] == 'high']
adata = adata[adata.obs['cell_type'].isin(['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs'])]
adata.obs['dataset'] = adata.obs['platform'].replace({
    '10X': 'child',
    'plate': 'adult'   
})
adata.obs['sick'] = adata.obs['Condition'].replace({
    'S_dengue': 'sick',
    'dengue': 'sick',
    'DWS': 'sick'
})
sc.pp.normalize_total(adata, target_sum=1e6) #normalize data to CPM (counts per million)

###########################################################################
adata_children = adata[adata.obs['dataset'] == 'child']
adata_adults = adata[adata.obs['dataset'] == 'adult']

conditions = list(adata.obs['Condition'].astype('category').cat.categories)
datasets = list(adata.obs['dataset'].astype('category').cat.categories)
sicks = list(adata.obs['sick'].astype('category').cat.categories)

adata_dic = {}
n_ID = {}
l_ID = {}

for condition in conditions:
    for dataset in datasets:
        adata_dic[(condition, dataset)] = adata[adata.obs['Condition'] == condition][adata[adata.obs['Condition'] == condition].obs['dataset'] == dataset]
        n_ID[(condition, dataset)] = len(adata_dic[(condition, dataset)].obs['ID'].astype('category').cat.categories)
        l_ID[(condition, dataset)] = list(adata_dic[(condition, dataset)].obs['ID'].astype('category').cat.categories)
        
for sick in sicks:
    for dataset in datasets:
        adata_dic[(sick, dataset)] = adata[adata.obs['sick'] == sick][adata[adata.obs['sick'] == sick].obs['dataset'] == dataset]
        n_ID[(sick, dataset)] = len(adata_dic[(sick, dataset)].obs['ID'].astype('category').cat.categories)
        l_ID[(sick, dataset)] = list(adata_dic[(sick, dataset)].obs['ID'].astype('category').cat.categories)

###############################################################################
del_if = []
for cd in ['Healthy', 'dengue', 'DWS', 'S_dengue']:
    adata_n = adata_dic[cd, 'child']
    n_cell = adata_n.obs[['cell_type', 'ID']].groupby(['cell_type', 'ID']).size().unstack()
    idx = np.argwhere(n_cell.values < 5)
    del_i= [[cd, n_cell.index[i], n_cell.columns[j]] for i, j in zip(idx[:, 0], idx[:, 1])]
    for i in range(len(del_i)): 
        del_if.append(del_i[i]) 
adata_new = adata[~adata.obs['ID'].isin([i[2] for i in del_if])]
adata_cDCs = adata_new[adata_new.obs['cell_type'] == 'cDCs']
adata_cDCs.write('/home/yike/phd/dengue/data/mergedata_20210519_cDCs.h5ad')

#######################################
# delete information
#[['Healthy', 'cDCs', '3_074_01'],
# ['DWS', 'cDCs', '5_089_01'],
# ['S_dengue', 'cDCs', '1_140_01'],
# ['S_dengue', 'cDCs', '5_193_01']]
 ######################################

fn_int = '/home/yike/phd/dengue/data/interaction_source_file/interactions_DB.tsv'
interactions = pd.read_csv(fn_int, sep=',')[['gene_name_a', 'gene_name_b']]
genes = np.unique(interactions)
#genes = [i for i in genes if i in adata.var_names]
genes = [i for i in genes if i not in ['CCL3L3', 'CCL4L1', 'CCN6', 'KIR3DS1', 'YARS1']]

#####################################################################################
from log2_FC_functions import log2_FC_all

log2_fc_ave = log2_FC_all(genes, adata_cDCs, 'S_dengue', 'dengue', 'child')[1]
log2_fc_ave['cell_type']=['cDCs']*log2_fc_all.shape[0]

log2_fc_ave['gene'] = log2_fc_ave.index.values
log2_fc_ave = log2_fc_ave.set_index('cell_type')
log2_fc_ave.to_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/log2_fc_ave_kids_cDCs.tsv')

#####################################################################################
adatag = adata_cDCs[:, genes]
from numpy import * # 调用numpy所有函数
adatag_children = adatag[adatag.obs['dataset'] == 'child']

IDs = list(adatag_children.obs['ID'].astype('category').cat.categories)

#cell_types = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
conditions = ['S_dengue', 'dengue', 'DWS', 'Healthy']
from collections import defaultdict
adatag_ch_ct_cd_ID = {}
for condition in conditions:
    for ID in IDs:
        adatag_ch_ct_cd_ID[(condition, ID)] = adatag_children[adatag_children.obs['Condition'] == condition][adatag_children[adatag_children.obs['Condition'] == condition].obs['ID'] == ID]

gene_exp_ave = {}
ave_exp_ave = {}
for condition in conditions:
    exp_fra = []
    ave_exp = []
    for ID in list(adatag_children[adatag_children.obs['Condition'] == condition].obs['ID'].astype('category').cat.categories):
        exp_fra.append((adatag_ch_ct_cd_ID[(condition, ID)].X > 0).toarray().mean(axis=0))
        ave_exp.append(adatag_ch_ct_cd_ID[(condition, ID)].X.toarray().mean(axis=0))
        
    exp_fra = np.array(exp_fra).mean(axis=0)
    ave_exp = np.array(ave_exp).mean(axis=0)
    
    gene_exp_ave[condition] = pd.DataFrame(exp_fra, index=adatag_children[adatag_children.obs['Condition'] == condition].var.index, columns=['gene_expre'])
    ave_exp_ave[condition] = pd.DataFrame(ave_exp, index=adatag_children[adatag_children.obs['Condition'] == condition].var.index, columns=['ave_exp'])
    
exp_fra_ave = pd.DataFrame([])
ave_exp = pd.DataFrame([])
for key in gene_exp_ave.keys():
    gene_exp_ave[key]['cell_type'] = ['cDCs']*(gene_exp_ave[key].shape[0])
    gene_exp_ave[key]['condition'] = [key]*(gene_exp_ave[key].shape[0])
    exp_fra_ave = pd.concat([exp_fra_ave, gene_exp_ave[key]])
    exp_fra_ave['gene'] = exp_fra_ave.index
    exp_fra_ave.to_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/exp_fra_ave_cDCs.tsv')
    
    ave_exp_ave[key]['cell_type'] = ['cDCs']*(ave_exp_ave[key].shape[0])
    ave_exp_ave[key]['condition'] = [key]*(ave_exp_ave[key].shape[0])
    ave_exp = pd.concat([ave_exp, ave_exp_ave[key]])
    ave_exp['gene'] = ave_exp.index
    ave_exp.to_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/ave_exp_ave_cDCs.tsv')

##################################################################################
from scipy import stats as ss
from collections import defaultdict

df = adatag_children.obs[['cell_type', 'Condition']].copy()
stas = []
p_value = []
for gene in log2_fc_cDCs.loc['cDCs'].index.tolist():
    df[gene] = adatag_children[:, gene].X.toarray()[:, 0]
    SD = df[(df['Condition'] == 'S_dengue')][gene]
    D = df[(df['Condition'] == 'dengue')][gene]
    res = ss.ks_2samp(SD, D, alternative='two-sided')
    stas.append(res[0])
    p_value.append(res[1])
log2_fc_cDCs['statistic'] = stas
log2_fc_cDCs['p_value'] = p_value
log2_fc_cDCs.to_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/log2_fc_ave_kids_cDCs.tsv')   
###############################################################################
log2_fc_cDCs = pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/log2_fc_ave_kids_cDCs.tsv', index_col=['cell_type', 'gene'])
gene_exp_cDCs = pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/exp_fra_ave_cDCs.tsv', index_col=['cell_type', 'gene', 'condition'])
ave_exp_cDCs = pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/ave_exp_ave_cDCs.tsv', index_col=['cell_type', 'gene', 'condition'])
gene_exp_cDCs.drop(['Unnamed: 0'], axis=1, inplace=True)
ave_exp_cDCs.drop(['Unnamed: 0'], axis=1, inplace=True)

idx_SD = [(i[0], i[1], 'S_dengue') for i in log2_fc_cDCs.index]
idx_D = [(i[0], i[1], 'dengue') for i in log2_fc_cDCs.index]

gene_exp_cDCs_SD = gene_exp_cDCs.loc[idx_SD]
gene_exp_cDCs_D = gene_exp_cDCs.loc[idx_D]
ave_exp_cDCs_SD = ave_exp_cDCs.loc[idx_SD]
ave_exp_cDCs_D = ave_exp_cDCs.loc[idx_D]

log2_fc_cDCs['SD_exp_frac'] = gene_exp_cDCs_SD['gene_expre'].tolist()
log2_fc_cDCs['D_exp_frac'] = gene_exp_cDCs_D['gene_expre'].tolist()
log2_fc_cDCs['SD_ave_exp'] = ave_exp_cDCs_SD['ave_exp'].tolist()
log2_fc_cDCs['D_ave_exp'] = ave_exp_cDCs_D['ave_exp'].tolist()

log2_fc_cDCs.to_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/log2_fc_ave_kids_cDCs.tsv')  
##############################################################################
log2_fc = pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/log2_fc_ave_kids.tsv', index_col=['cell_type', 'gene'])
gene_exp = pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/exp_fra_ave.tsv', index_col=['cell_type', 'gene', 'condition'])
ave_exp = pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/ave_exp_ave.tsv', index_col=['cell_type', 'gene', 'condition'])

log2_fc_new = log2_fc.loc[['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'pDCs']]
log2_fc = pd.concat([log2_fc_new, log2_fc_cDCs])
log2_fc.to_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/log2_fc_ave_kids.tsv')

gene_exp_new = gene_exp.loc[['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'pDCs']]
gene_exp = pd.concat([gene_exp_new, gene_exp_cDCs])
gene_exp.to_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/exp_fra_ave.tsv')

ave_exp_new = ave_exp.loc[['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'pDCs']]
ave_exp = pd.concat([ave_exp_new, ave_exp_cDCs])
ave_exp.to_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/ave_exp_ave.tsv')
#########################################################################
def update_cDCs(fdn): 
    
    inters_im_all = pd.read_csv(fdn+'.tsv', index_col=0)
    
    inters_im_cDCs = inters_im_all[(inters_im_all['cta'] == 'cDCs') | (inters_im_all['ctb'] == 'cDCs')]
    inters_im_other = inters_im_all[~ ((inters_im_all['cta'] == 'cDCs') | (inters_im_all['ctb'] == 'cDCs'))]
    log2_fc =pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/log2_fc_ave_kids.tsv', index_col=['gene','cell_type', ])

    inters_im_cDCs_cta = inters_im_cDCs[(inters_im_cDCs['cta'] == 'cDCs') & (inters_im_cDCs['ctb'] != 'cDCs')].copy()

    inters_im_cDCs_cta = inters_im_cDCs_cta.set_index(['ga', 'cta'])
    inters_im_cDCs_cta['ga_log2FC'] = log2_fc.loc[inters_im_cDCs_cta.index]['fold_2_change']
    inters_im_cDCs_cta['ga_comfrac'] = log2_fc.loc[inters_im_cDCs_cta.index]['comp_frac']
    inters_im_cDCs_cta['ga_D_exp_frac'] = log2_fc.loc[inters_im_cDCs_cta.index]['D_exp_frac']
    inters_im_cDCs_cta['ga_SD_exp_frac'] = log2_fc.loc[inters_im_cDCs_cta.index]['SD_exp_frac']

    inters_im_cDCs_cta = inters_im_cDCs_cta.reset_index()

    inters_im_cDCs_ctb = inters_im_cDCs[(inters_im_cDCs['ctb'] == 'cDCs') & (inters_im_cDCs['cta'] != 'cDCs')].copy()

    inters_im_cDCs_ctb = inters_im_cDCs_ctb.set_index(['gb', 'ctb'])
    inters_im_cDCs_ctb['gb_log2FC'] = log2_fc.loc[inters_im_cDCs_ctb.index]['fold_2_change']
    inters_im_cDCs_ctb['gb_comfrac'] = log2_fc.loc[inters_im_cDCs_ctb.index]['comp_frac']
    inters_im_cDCs_ctb['gb_D_exp_frac'] = log2_fc.loc[inters_im_cDCs_ctb.index]['D_exp_frac']
    inters_im_cDCs_ctb['gb_SD_exp_frac'] = log2_fc.loc[inters_im_cDCs_ctb.index]['SD_exp_frac']

    inters_im_cDCs_ctb = inters_im_cDCs_ctb.reset_index()

    inters_im_cDCs_ctb = inters_im_cDCs_ctb[inters_im_cDCs_cta.columns]

    inters_im_cDCs_ctab = inters_im_cDCs[(inters_im_cDCs['ctb'] == 'cDCs') & (inters_im_cDCs['cta'] == 'cDCs')].copy()

    inters_im_cDCs_ctab = inters_im_cDCs_ctab.set_index(['ga', 'cta'])
    inters_im_cDCs_ctab['ga_log2FC'] = log2_fc.loc[inters_im_cDCs_ctab.index]['fold_2_change']
    inters_im_cDCs_ctab['ga_comfrac'] = log2_fc.loc[inters_im_cDCs_ctab.index]['comp_frac']
    inters_im_cDCs_ctab['ga_D_exp_frac'] = log2_fc.loc[inters_im_cDCs_ctab.index]['D_exp_frac']
    inters_im_cDCs_ctab['ga_SD_exp_frac'] = log2_fc.loc[inters_im_cDCs_ctab.index]['SD_exp_frac']

    inters_im_cDCs_ctab = inters_im_cDCs_ctab.reset_index()

    inters_im_cDCs_ctab = inters_im_cDCs_ctab.set_index(['gb', 'ctb'])
    inters_im_cDCs_ctab['gb_log2FC'] = log2_fc.loc[inters_im_cDCs_ctab.index]['fold_2_change']
    inters_im_cDCs_ctab['gb_comfrac'] = log2_fc.loc[inters_im_cDCs_ctab.index]['comp_frac']
    inters_im_cDCs_ctab['gb_D_exp_frac'] = log2_fc.loc[inters_im_cDCs_ctab.index]['D_exp_frac']
    inters_im_cDCs_ctab['gb_SD_exp_frac'] = log2_fc.loc[inters_im_cDCs_ctab.index]['SD_exp_frac']

    inters_im_cDCs_ctab = inters_im_cDCs_ctab.reset_index()

    inters_im_cDCs_ctab = inters_im_cDCs_ctab[inters_im_cDCs_cta.columns]

    inters_im_cDCs_new = pd.concat([inters_im_cDCs_cta, inters_im_cDCs_ctb, inters_im_cDCs_ctab])

    inters_im_all_new = pd.concat([inters_im_cDCs_new, inters_im_other])

    inters_im_all_new = inters_im_all_new[inters_im_all.columns]
    inters_im_all_new.to_csv(fdn+'_new.tsv')
    return inters_im_all_new

update_cDCs('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/inters/inters_mix_all_002ave')

