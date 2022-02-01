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
import random
import itertools
from numpy import *  


import sys
sys.path.append('/home/yike/phd/dengue/dengue_children') #enter the YK_util file directory
sys.path.append('/home/yike/phd/dengue/dengue_children/paper_figures') #enter the YK_util file directory
import YK_util, optimized_pair_comparison_new, functions_CCC_figure_fra_pair
from YK_util import *
from optimized_pair_comparison_new import *
from functions_CCC_figure_fra_pair import *

path = '/home/yike/phd/dengue/data/mergedata_20211001.h5ad'
adata = sc.read_h5ad(path)
adata_kid = subsetdata(adata)
adata_kid = normalizedata(adata_kid, log1p=True) # log1p = 2
adata_kid = removegenes(adata_kid)
# batch 1:
# patient 6_029_01: 
# only contains cells in ['Plasmablasts', 'doublets', 'B_cells', 'T_cells', 'NK_cells']

print('Load interaction') 
# DB_fn_int = '/home/yike/phd/dengue/data/interaction_source_file/interactions_DB.tsv'
# DB_interactions = pd.read_csv(DB_fn_int, sep=',')[['gene_name_a', 'gene_name_b']]

fn_int = '/home/yike/phd/dengue/data/interaction_source_file/omni_DB_inters.tsv'
interactions = pd.read_csv(fn_int, sep='\t')[['genesymbol_intercell_source', 'genesymbol_intercell_target']]

genes = np.unique(interactions)
genes = [gene for gene in genes if gene in adata_kid.var_names]

save_tabels = '/home/yike/phd/dengue/data/tables/dataset_20211001/'
save_figures = '/home/yike/phd/dengue/figures/paper_figure/dataset_20211001/'

#pair comparison generation for cell types of kids 
cts = adata_kid.obs['cell_type_new'].unique().tolist()
cts.remove('doublets')
cts.remove('unknown')

pair = pair_comparison(adata_kid, cts, gene_cut_off=False, log1p=2)

pair['pair_res'].to_csv(save_tabels + 'ct_pair_gene_cut_0.tsv', sep='\t', index=False)
pair['log_FCs'].to_csv(save_tabels + 'ct_log2FC_gene_cut_0.tsv', sep='\t')

fra_av = fra_avg(adata_kid, cts, log1p=2)
fra_av['fra'].to_csv(save_tabels + 'ct_fra_gene_cut_0.tsv', sep='\t', index=False)
fra_av['avg'].to_csv(save_tabels + 'ct_avg_gene_cut_0.tsv', sep='\t', index=False)

data = combination(pair['pair_res'], fra_av['fra'], fra_av['avg'])
data.to_csv(save_tabels + 'ct_data_pair_gene_cut_0.tsv', sep='\t')

#pair comparison generation for cell subtypes of kids
csts = adata_kid.obs['cell_subtype_new'].unique().tolist()
csts.remove('doublets')
csts.remove('megakaryocytes')
csts.remove('plasmacytoid DCs')
csts.remove('conventional DCs')
csts.remove('unknown')

pair = pair_comparison(adata_kid, csts, gene_cut_off=False, log1p=2)

pair['pair_res'].to_csv(save_tabels + 'cst_pair_gene_cut_0.tsv', sep='\t', index=False)
pair['log_FCs'].to_csv(save_tabels + 'cst_log2FC_gene_cut_0.tsv', sep='\t')

fra_av = fra_avg(adata_kid, csts, log1p=2)
fra_av['fra'].to_csv(save_tabels + 'cst_fra_gene_cut_0.tsv', sep='\t', index=False)
fra_av['avg'].to_csv(save_tabels + 'cst_avg_gene_cut_0.tsv', sep='\t', index=False)

data = combination(pair['pair_res'], fra_av['fra'], fra_av['avg'])
data.to_csv(save_tabels + 'cst_data_pair_gene_cut_0.tsv', sep='\t')

ct_fra = pd.read_csv(save_tabels + 'ct_fra_gene_cut_0.tsv', sep='\t', index_col=['cell_type_new', 'gene', 'condition'])
ct_avg = pd.read_csv(save_tabels + 'ct_avg_gene_cut_0.tsv', sep='\t', index_col=['cell_type_new', 'gene', 'condition'])

cst_fra = pd.read_csv(save_tabels + 'cst_fra_gene_cut_0.tsv', sep='\t', index_col=['cell_subtype_new', 'gene', 'condition'])
cst_avg = pd.read_csv(save_tabels + 'cst_avg_gene_cut_0.tsv', sep='\t', index_col=['cell_subtype_new', 'gene', 'condition'])

ct_pair = pd.read_csv(save_tabels + 'ct_data_pair_gene_cut_0.tsv', sep='\t', index_col=['cell_subtype', 'gene'])
ct_log2FC = pd.read_csv(save_tabels + 'ct_log2FC_gene_cut_0.tsv', sep='\t', index_col=0)

cst_pair = pd.read_csv(save_tabels + 'cst_data_pair_gene_cut_0.tsv', sep='\t', index_col=['cell_subtype', 'gene'])
cst_log2FC = pd.read_csv(save_tabels + 'cst_log2FC_gene_cut_0.tsv', sep='\t', index_col=0)

cst_pair = pd.concat([cst_pair, ct_pair.loc[['megakaryocytes', 'plasmacytoid DCs', 'conventional DCs']]])
cst_log2FC = pd.concat([cst_log2FC, ct_log2FC.loc[['megakaryocytes', 'plasmacytoid DCs', 'conventional DCs']]])
cst_fra = pd.concat([cst_fra, ct_fra.loc[['megakaryocytes', 'plasmacytoid DCs', 'conventional DCs']]])
cst_avg = pd.concat([cst_avg, ct_avg.loc[['megakaryocytes', 'plasmacytoid DCs', 'conventional DCs']]])

cst_pair.to_csv(save_tabels + 'cst_data_pair_gene_cut_0.tsv', sep='\t')
cst_log2FC.to_csv(save_tabels + 'cst_log2FC_gene_cut_0.tsv', sep='\t')
cst_fra.to_csv(save_tabels + 'cst_fra_gene_cut_0.tsv', sep='\t')
cst_avg.to_csv(save_tabels + 'cst_avg_gene_cut_0.tsv', sep='\t')

ct_pair['log2_fold_change'] = np.log2(ct_pair['S_avg'] + 1) - np.log2(ct_pair['NS_avg'] + 1)
cst_pair['log2_fold_change'] = np.log2(cst_pair['S_avg'] + 1) - np.log2(cst_pair['NS_avg'] + 1)
ct_pair.to_csv(save_tabels + 'ct_data_pair_gene_cut_0.tsv', sep='\t')
cst_pair.to_csv(save_tabels + 'cst_data_pair_gene_cut_0.tsv', sep='\t')

import anndataks

cts = adata_kid.obs['cell_type_new'].unique().tolist()
cts.remove('doublets')
cts.remove('unknown')

csts = adata_kid.obs['cell_subtype_new'].unique().tolist()
csts.remove('doublets')
csts.remove('megakaryocytes')
csts.remove('plasmacytoid DCs')
csts.remove('conventional DCs')
csts.remove('unknown')

subcts = cts + csts
conditions = ['S_dengue', 'dengue']

results = {}
for subct in subcts:
    if subct in adata_kid.obs['cell_type_new'].astype('category').cat.categories:
        adata_ct = adata_kid[adata_kid.obs['cell_type_new'] == subct]
    elif subct in adata_kid.obs['cell_subtype_new'].astype('category').cat.categories:
        adata_ct = adata_kid[adata_kid.obs['cell_subtype_new'] == subct]
        
    adata_SD = adata_ct[adata_ct.obs['Condition'] == 'S_dengue']
    adata_D = adata_ct[adata_ct.obs['Condition'] == 'dengue']
    # while calculating ks test pvalue, the adata is log1ped, so the argument log1p=2
    results[subct] = anndataks.compare(adata_D, adata_SD, log1p=2, mode='asymp') # log2_fold_change: adata_Sd vs adata_D

ks_res = pd.DataFrame([])
for subct in subcts:
    results[subct]['cell_subtype'] = [subct] * results[subct].shape[0]
    ks_res = pd.concat([ks_res, results[subct]])

ks_res.set_index(['cell_subtype', ks_res.index], inplace=True)
ks_res.to_csv(save_tabels + 'ks_pvalue.tsv', sep='\t')

ks_res = pd.read_csv(save_tabels + 'ks_pvalue.tsv', sep='\t', index_col=['cell_subtype', 'Unnamed: 1'])
ct_pair = pd.concat([ct_pair, ks_res.loc[ct_pair.index][['statistic', 'pvalue']]], axis=1)
ct_pair.to_csv(save_tabels + 'ct_data_pair_gene_cut_0.tsv', sep='\t')

cst_pair = pd.concat([cst_pair, ks_res.loc[cst_pair.index][['statistic', 'pvalue']]], axis=1)
cst_pair.to_csv(save_tabels + 'cst_data_pair_gene_cut_0.tsv', sep='\t')