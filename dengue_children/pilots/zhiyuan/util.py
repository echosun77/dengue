import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
import scanpy as sc #for scanpy >= 1.3.7
import anndata as ann
import numpy as np
import scipy as sp
import pandas as pd
import logging
from copy import copy
import glob
import os
import anndataks

def getdata(path):
    """Load data from file path
    
    Parameters:
    ----------
    path: str
        absolute file path (.loom).

    Returns:
    -------
    adata: scanpy.adata
        scanpy adata object
    """
    adata = sc.read_loom(path, var_names = 'var_names', obs_names = 'obs_names')
    return(adata)

def subsetdata(adata, quality='high', platform ='10X', doublets='no'):
    """subset data by cell quality, platform and non_doublets
    
    Parameters:
    ----------
    adata: scanpy.adata
        scanpy adata object
    quality: str
        'high' or 'low'
    platform: str
        '10X' for children, 'plate' for adult
    doublets: str
        'yes' or 'no'

    Returns:
    -------
    temp: scanny.adata
        scanpy adata object
    """
    temp = adata[(adata.obs.cell_quality == quality) & 
                 (adata.obs.platform == platform) & 
                 (adata.obs.doublets == doublets),].copy()
    return(temp)

def normalizedata(adata):
    """normalize the dataset
    
    Parameters:
    ----------
    adata: scanpy.adata
        scanpy adata object

    Returns:
    -------
    adata: scanpy.adata
        scanpy adata object
    """
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata, base=2)
    return(adata)

def removegenes(adata):
    """remove human HLA genes from the dataset 
    Parameters:
    ----------
    adata: scanpy.adata
        scanpy adata object

    Returns:
    -------
    temp: scanpy.adata
        scanpy adata object
    """
    IGKV = [x for x in adata.var_names if x.startswith('IGKV')]
    IGHV = [x for x in adata.var_names if x.startswith('IGHV')]
    IGLV = [x for x in adata.var_names if x.startswith('IGLV')]
    IGLC = [x for x in adata.var_names if x.startswith('IGLC')]
    IGKC = [x for x in adata.var_names if x.startswith('IGKC')]
    IGHC = [x for x in adata.var_names if x.startswith('IGHC')]
    IGHC = [x for x in adata.var_names if x.startswith('IGHC')]
    TRAV = [x for x in adata.var_names if x.startswith('TRAV')]
    TRBV = [x for x in adata.var_names if x.startswith('TRBV')]
    
    #try removing IGHG genes and MZB1 and JCHAIN
#     IGHG = [x for x in adata.var_names if x.startswith('IGHG')]
    exclude = IGKV + IGHV + IGLV + IGLC + IGKC + IGHC + TRAV + TRBV 
    gene = [x for x in adata.var_names if x not in exclude]
    temp = adata[:,gene].copy()
    return(temp)



def cluster (adata):
    sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=2000)
    sc.pp.pca(adata, n_comps=40, use_highly_variable=True, svd_solver='arpack')
    sc.pp.neighbors(adata, n_pcs = 15)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.3, key_added = 'leiden_r0.3')
    sc.tl.leiden(adata, resolution=0.5, key_added = 'leiden_r0.5')
    sc.tl.rank_genes_groups(adata, groupby='leiden_r0.3', key_added='rank_genes_r0.3')
    sc.tl.rank_genes_groups(adata, groupby='leiden_r0.5', key_added='rank_genes_r0.5')
#     sc.tl.louvain(adata, resolution=0.2, key_added = 'louvain_r0.2')
#     sc.tl.louvain(adata, resolution=0.3, key_added = 'louvain_r0.3')
#     sc.tl.rank_genes_groups(adata, groupby='louvain_r0.2', key_added='rank_genes_r0.2')
#     sc.tl.rank_genes_groups(adata, groupby='louvain_r0.3', key_added='rank_genes_r0.3')


def KS(adata):
    adata_1 = adata[adata.obs.loc[:,'DENV_reads'] == 0,]
    adata_2 = adata[adata.obs.loc[:,'DENV_reads'] != 0,]
    results = anndataks.compare(adata_1, adata_2, log1p=2)
    results_sort = results.sort_values(by = 'statistic', ascending=False)
#     results_sort_2 = results_sort[(results_sort.pvalue < 0.05)]
    return(results_sort)

def virus_nor(adata):
    adata.obs['virus'] = 'no'
    adata.obs.loc[adata.obs.DENV_reads != 0, 'virus'] = 'yes'
    temp = adata[adata.obs.DENV_reads != 0,]
    temp.obs['DENV_reads_cpm'] = (temp.obs['DENV_reads'] / temp.obs['n_counts']) * 1e6
    temp.obs['DENV_reads_nor'] = np.log2(temp.obs['DENV_reads_cpm'] + 1)
    adata.obs['DENV_reads_nor'] = 0
    ind = temp.obs.index
    adata.obs.loc[ind, 'DENV_reads_nor'] = temp.obs.loc[ind, 'DENV_reads_nor']

from scipy import stats
def spearman(adata):
    df = pd.DataFrame(index=adata.var.index, columns=['spearman_cor', 'p_value'])
    DENV = list(adata.obs.loc[:, 'DENV_reads_nor'])
    for gene in df.index:
        gene_exp = list(adata[:,gene].X.todense().flatten().A1)
#         print(gene_exp)
#         print(len(gene_exp))
#         print(DENV)
#         print(len(DENV))
        temp = stats.spearmanr(gene_exp, DENV)
        df.loc[gene, 'spearman_cor'] = temp[0]
        df.loc[gene, 'p_value'] = temp[1]
    return(df)

def scatter(adata, df, gene, path):
    name = path+'/'+gene+'.png'
    rcParams['figure.figsize']=(7,7)
    rcParams['font.size']= 18
    df_temp = pd.DataFrame(index = adata.obs.index, columns = ['DENV', gene])
    df_temp.loc[:,'DENV'] = adata.obs.loc[:, 'DENV_reads_nor']
    df_temp.loc[:,gene] = adata[:,gene].X.todense().flatten().A1
    ax = df_temp.plot.scatter(x = 'DENV', y = gene)
    text = 'rho='+ str(round(df.loc[gene, 'spearman_cor'],2))
    print(text)
    ax.text(max(df_temp.DENV) - 1,min(df_temp[gene]) + 0.5, text, fontsize = 14,
           bbox=dict(facecolor='none', edgecolor='black', boxstyle='round'))
    ax.set_xlabel('DENV_reads (log2(cpm+1))')
    ax.set_ylabel(gene+' (log2(cpm+1))')
    fig = ax.get_figure()
    fig.savefig(name, bbox_inches = 'tight')

def histo(df):
    rcParams['figure.figsize']=(7,7)
    rcParams['font.size']= 18
    ax = df.spearman_cor.plot.hist(bins=50, alpha=0.5)
    ax.set_xlabel('spearman_rho')
    ax.set_ylabel('gene_number')
    fig = ax.get_figure()
    
def pcf_plotting(adata,gene,path, filename = ''):
    if filename == '':
        name = gene + '.png'
    else:
        name = filename + '.png'
    rcParams['font.size']=18
    n_bins = 100
    fig, ax = plt.subplots(figsize=(8, 8))
    x1 = adata[adata.obs.loc[:,'virus'] == 'yes',gene].X.todense().flatten().A1
    x2 = adata[adata.obs.loc[:,'virus'] == 'no',gene].X.todense().flatten().A1

    ax.hist(x1, n_bins, density=True, histtype='step',
                           cumulative=True, label='virus_harboring_cells')
    ax.hist(x2, n_bins, density=True, histtype='step',
                           cumulative=True, label='bystanders')

    ax.grid(True)
    ax.legend(loc='lower center')
    ax.set_title(gene)
    ax.set_xlabel('log2(cpm+1)')
    ax.set_ylabel('Fraction of cells < x')
    fig.savefig(os.path.join(path, name), bbox_inches = 'tight')


def scatter_2d(df_merge, title, fold_upper, fold_lower, rho_upper, rho_lower, path, filename):

    ind_1 = df_merge[((df_merge.log2_fold_change > fold_upper) | (df_merge.log2_fold_change < fold_lower)) & 
              ((df_merge.spearman_rho > rho_upper) | (df_merge.spearman_rho < rho_lower))].index
    ind_2 = df_merge[((df_merge.log2_fold_change >fold_upper) | (df_merge.log2_fold_change < fold_lower)) & 
              ((df_merge.spearman_rho <= rho_upper) & (df_merge.spearman_rho >= rho_lower))].index
    ind_3 = df_merge[((df_merge.log2_fold_change <= fold_upper) & (df_merge.log2_fold_change >= fold_lower)) & 
              ((df_merge.spearman_rho > rho_upper) | (df_merge.spearman_rho < rho_lower))].index

    df_merge.loc[:, 'outlier'] = 'no'
    df_merge.loc[ind_1, 'outlier'] = '1'
    df_merge.loc[ind_2, 'outlier'] = '2'
    df_merge.loc[ind_3, 'outlier'] = '3'
    
    
    rcParams['figure.figsize']=(14,14)
    fig, ax = plt.subplots()
    sb.scatterplot(data = df_merge, x = 'log2_fold_change', y = 'spearman_rho', hue='outlier', 
               hue_order = ['no', '1', '2', '3'], ax = ax)
    ax.set_title(title)
    ax.get_legend().remove()
    ax.set_xlabel('log2 fold change (Virus_harboring_cells vs bystanders)')
    for index in (list(ind_1) + list(ind_2) + list(ind_3)):
        ax.text(df_merge.loc[index, 'log2_fold_change']+0.1, 
        df_merge.loc[index, 'spearman_rho'],
        index, color = 'black',
        fontsize = 14)
        
    fig.savefig(os.path.join(path, filename), bbox_inches = 'tight')
    return(ind_1, ind_2, ind_3)