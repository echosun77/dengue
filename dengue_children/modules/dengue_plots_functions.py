import pandas as pd 
import numpy as np

import anndata
import scanpy as sc

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

from collections import defaultdict

################### define a function to show gene average expression in different cell types at several conditions
def bar_ave_update(adata, gene, ax2_t, ax1_b, ax1_t):

    df = adata.obs[['cell_type', 'Condition']].copy()
    df[gene] = adata[:, gene].X.toarray()[:, 0]
    df['is_expressed'] = df[gene] > 0
    gby = df.groupby(['cell_type', 'Condition']).mean()
    # mean_expression = gby.loc['cell_type', 'CXCL10'] # cell_type
    # frac_expressed =  gby.loc[cell_type, 'is_expressed']

    ############################### 
    a1 = 1; a2 = 2
    gs = gridspec.GridSpec(2,1, height_ratios=[a1, a2], hspace=0.05)

    fig, axs = plt.subplots(dpi=300, figsize=(3.6, 2.4))
    ax1 = plt.subplot(gs[0, 0:])
    ax2 = plt.subplot(gs[1, 0:], sharex = ax1)

    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.axes.get_xaxis().set_visible(False)
    
    ax1.tick_params(labeltop='off')  # don't put tick labels at the top
    ax2.xaxis.tick_bottom() # move ticks and ticklabels (if present) to the bottom of the axes

    ax1.set_ylim(ax1_b, ax1_t)  # outliers only
    ax2.set_ylim(0, ax2_t)  # most of the data

    d = .015
    oa1,oa2=(a1+a2)/a1,(a1+a2)/a2
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d, +d), (-d*oa1, +d*oa1), **kwargs, label='_nolegend_')        # top-left diagonal
    #ax2.plot((1 - d, 1 + d), (-d*oa1, +d*oa1), **kwargs)  # top-right diagonal

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d*oa2, 1 + d*oa2), **kwargs)  # bottom-left diagonal
    #ax2.plot((1 - d, 1 + d), (1 - d*oa2, 1 + d*oa2), **kwargs)  # bottom-right diagonal
    
    ############################### 
    gby.unstack(1)[gene][['Healthy', 'dengue', 'S_dengue']].plot.bar(ax=ax1, width=0.8)
    gby.unstack(1)[gene][['Healthy', 'dengue', 'S_dengue']].plot.bar(ax=ax2, legend=None, width=0.8)
    
    ########################################## expression fraction
    # gby.unstack(1)['is_expressed'][['Healthy', 'dengue', 'S_dengue']].plot.bar(ax=ax, width=0.8) 
    ##########################################
    
    ax1.legend(labels=['Healthy', 'Dengue', 'Severe dengue'], fontsize=8)
    ax2.set_xlabel(None)
    #ax2.set_ylabel('Fraction of cells \n expressing ' + gene, fontsize=10)
    ax2.set_ylabel('gene expression[cpm]', fontsize=10)
    ax2.set_xticklabels(['B cells', 'Monocytes', 'NK cells', 'Plasmablasts', 'T cells', 'cDCs', 'pDCs'])
    plt.xticks(fontsize=10)
    ax2.yaxis.set_label_coords(-0.15, 0.8)
    ax1.set_title(gene)

    plt.show()
    
    return {'fig': fig, 'ax': axs}

################### define a function to show gene average expression in different cell types at several conditions

def bar_ave(adata, gene, ylim):

    df = adata.obs[['cell_type', 'Condition', 'ID']].copy()
    df[gene] = adata[:, gene].X.toarray()[:, 0]
    df['is_expressed'] = df[gene] > 0
    gby = df.groupby(['cell_type', 'Condition']).mean()
    gby_new = df.groupby(['cell_type', 'Condition', 'ID']).mean()

    ave_exp = []
    exp_frac = []
    for idx in list(gby.index):
        ave_exp.append(gby_new.loc[idx][gene].mean())
        exp_frac.append(gby_new.loc[idx]['is_expressed'].mean())

    gby['ave_exp'] = ave_exp
    gby['exp_frac'] = exp_frac

    fig, ax = plt.subplots(dpi=300, figsize=(3.6, 2.4))
    gby.unstack(1)['ave_exp'][['Healthy', 'dengue', 'S_dengue']].plot.bar(ax=ax, width=0.8)
    
    ########################################## expression fraction
    # gby.unstack(1)['is_expressed'][['Healthy', 'dengue', 'S_dengue']].plot.bar(ax=ax, width=0.8) 
    ##########################################
    
    ax.set_xlabel(None)
    ax.set_xticklabels(['B cells', 'Monocytes', 'NK cells', 'Plasmablasts', 'T cells', 'cDCs', 'pDCs'])
    
    ax.set_ylim(0, ylim)
    ax.set_ylabel('gene expression[cpm]', fontsize=10)
    ax.legend(labels=['Healthy', 'Dengue', 'Severe dengue'], fontsize=8)
    plt.xticks(fontsize=10)
    plt.title(gene)
    plt.legend(loc='upper right')
    plt.show()
    
    return {'fig': fig, 'ax': ax}

def bar_exp(adata, gene, ylim):

    df = adata.obs[['cell_type', 'Condition', 'ID']].copy()
    df[gene] = adata[:, gene].X.toarray()[:, 0]
    df['is_expressed'] = df[gene] > 0
    gby = df.groupby(['cell_type', 'Condition']).mean()
    gby_new = df.groupby(['cell_type', 'Condition', 'ID']).mean()

    ave_exp = []
    exp_frac = []
    for idx in list(gby.index):
        ave_exp.append(gby_new.loc[idx][gene].mean())
        exp_frac.append(gby_new.loc[idx]['is_expressed'].mean())

    gby['ave_exp'] = ave_exp
    gby['exp_frac'] = exp_frac

    fig, ax = plt.subplots(dpi=300, figsize=(3.6, 2.4))
    gby.unstack(1)['exp_frac'][['Healthy', 'dengue', 'S_dengue']].plot.bar(ax=ax, width=0.8)

    ax.set_xlabel(None)
    ax.set_xticklabels(['B cells', 'Monocytes', 'NK cells', 'Plasmablasts', 'T cells', 'cDCs', 'pDCs'])
    
    ax.set_ylim(0, ylim)
    ax.set_ylabel('Fraction of cells\n expressing '+gene, fontsize=10)
    ax.legend(labels=['Healthy', 'Dengue', 'Severe dengue'], fontsize=8)
    plt.xticks(fontsize=10)
    plt.title(gene)
    plt.legend(loc='upper right')
    plt.show()
    
    return {'fig': fig, 'ax': ax}

# average gene expression in one cell type at different conditions
def bar_ct(adata, gene, cell_type):

    df = adata.obs[['cell_type', 'Condition', 'ID']].copy()
    df[gene] = adata[:, gene].X.toarray()[:, 0]
    df['is_expressed'] = df[gene] > 0
    gby = df.groupby(['cell_type', 'Condition']).mean()
    gby_new = df.groupby(['cell_type', 'Condition', 'ID']).mean()

    ave_exp = []
    exp_frac = []
    for idx in list(gby.index):
        ave_exp.append(gby_new.loc[idx][gene].mean())
        exp_frac.append(gby_new.loc[idx]['is_expressed'].mean())

    gby['ave_exp'] = ave_exp
    gby['exp_frac'] = exp_frac

    fig, ax = plt.subplots(dpi=300, figsize=(3.6, 2.4))
    gby.unstack(1)['ave_exp'][['Healthy', 'dengue', 'S_dengue']].loc[cell_type].plot.bar(ax=ax, width=0.5, label='_nolegend_', color='orange')
    ax.set_ylabel('gene expression[cpm]', fontsize=10)
    ax.set_xticklabels(['Healthy', 'Dengue', 'Severe dengue'], rotation=90)
    ax.set_xlabel(None)
    plt.xticks(fontsize=10)
    plt.title(gene + ' in ' + cell_type)
    plt.show()
    
    return {'fig': fig, 'ax': ax}

def heatmappy(adata, gene, columns=['cell_type', 'Condition'], **kwargs):
    dfi = pd.Series(adata[:, gene].X.toarray()[:, 0], index=adata.obs_names, name=gene).to_frame()
    for col in columns:
        dfi[col] = adata.obs[col]
    dfi['is_expressed'] = dfi[gene] > 0 # fraction; if without > 0, it will be average expression
    tmp = dfi.groupby(columns).mean()['is_expressed'] * 100
    ax = sns.heatmap(tmp.unstack(1)[['Healthy', 'dengue', 'DWS', 'S_dengue']], **kwargs)
    return ax

def KS_test(adata, cell_type, gene):
    from scipy import stats as ss
    from collections import defaultdict

    #adata_ch = adata_children[adata_children.obs['Condition'].isin(['Healthy', 'dengue', 'S_dengue'])]
    #adata_ch = adata_ch[adata_ch.obs['cell_type'].isin(['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs'])]
    
    df = adata.obs[['cell_type', 'Condition']].copy()
    df[gene] = adata[:, gene].X.toarray()[:, 0]
    SD = df[(df['Condition'] == 'S_dengue') & (df['cell_type'] == cell_type)][gene]
    D = df[(df['Condition'] == 'dengue') & (df['cell_type'] == cell_type)][gene]

    res = ss.ks_2samp(SD, D, alternative='two-sided')
    
    return {'statistic': res[0], 'pvalue': res[1]}

# cmap: bwr, viridis_r, Blues, RdBu_r, viridis, RdBu_r
def vs_dot(adata, cond1, cond2, gb, inters, vmax_n):
    conditions = list(adata.obs['Condition'].astype('category').cat.categories)
    datasets = list(adata.obs['dataset'].astype('category').cat.categories)
    sicks = list(adata.obs['sick'].astype('category').cat.categories)
    cell_types = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']

    from collections import defaultdict
    adata_dic = {}

    for dataset in datasets:                                                                                                                          
        for condition in conditions:
            adata_dic[(condition, dataset)] = adata[adata.obs['Condition'] == condition][adata[adata.obs['Condition'] == condition].obs['dataset'] == dataset]
        for sick in sicks:
            adata_dic[(sick, dataset)] = adata[adata.obs['sick'] == sick][adata[adata.obs['sick'] == sick].obs['dataset'] == dataset]
        for cell_type in cell_types:
            adata_dic[(cell_type, dataset)] = adata[adata.obs['cell_type'] == cell_type][adata[adata.obs['cell_type'] == cell_type].obs['dataset'] == dataset]
    
    fig,axs = plt.subplots(2, 2, figsize=(10,8), gridspec_kw={'wspace':-0.4}, sharex=True, facecolor='white', dpi=300, sharey=True)
    for age, ax_row in zip(['child', 'adult'], axs):
        for cond, ax in zip([cond1, cond2], ax_row):
            sc.pl.dotplot(adata_dic[(cond, age)], inters, groupby=gb, ax=ax, show=False, cmap='plasma', vmin=0, vmax=vmax_n + 0.5 *(age=='adult'), dot_min=0, dot_max=1)
            # vmax=1 + 0.5 *(age=='adult')
            ax.set_title(cond + ' ' + age, fontsize=15, y=0.8, loc='left')
            
    axes = fig.get_axes()
    #del axes[6:7], axes[14:15]
    axes[6].remove()
    axes[7].remove()
    axes[14].remove()
    axes[15].remove()
    #axes[8].set_yticks([])
    #axes[16].set_yticks([])
    fig.tight_layout()
    return {'figure':fig, 'subplots':axs}

# violin plot
def vs_violin(adata, cond1, cond2, gb, gene):
    conditions = list(adata.obs['Condition'].astype('category').cat.categories)
    datasets = list(adata.obs['dataset'].astype('category').cat.categories)
    sicks = list(adata.obs['sick'].astype('category').cat.categories)
    cell_types = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']

    from collections import defaultdict
    adata_dic = {}

    for dataset in datasets:                                                                                                                          
        for condition in conditions:
            adata_dic[(condition, dataset)] = adata[adata.obs['Condition'] == condition][adata[adata.obs['Condition'] == condition].obs['dataset'] == dataset]
        for sick in sicks:
            adata_dic[(sick, dataset)] = adata[adata.obs['sick'] == sick][adata[adata.obs['sick'] == sick].obs['dataset'] == dataset]
        for cell_type in cell_types:
            adata_dic[(cell_type, dataset)] = adata[adata.obs['cell_type'] == cell_type][adata[adata.obs['cell_type'] == cell_type].obs['dataset'] == dataset]
    
    fig,axs = plt.subplots(2, 2, figsize=(10,6), gridspec_kw={'wspace':0.2, 'hspace':0.75}, facecolor='white', dpi=300, sharey=True)
    for age, ax_row in zip(['child', 'adult'], axs):
        for cond, ax in zip([cond1, cond2], ax_row):
            sc.pl.violin(adata_dic[(cond, age)], gene, groupby=gb, rotation=45, ax=ax, show=False)
            #sc.pl.dotplot(adata_dic[(cond, age)], inters, groupby=gb, ax=ax, show=False, cmap='RdBu_r', vmin=0, vmax=vmax_n + 0.5 *(age=='adult'), dot_min=0, dot_max=1)
            ax.set_title(cond + ' ' + age, fontsize=15, y=1, loc='center')
    fig.tight_layout()
    return {'figure':fig, 'subplots':axs}

# stacked violin plot
def vs_stacked_violin(adata, cond1, cond2, gb, genes, vmax_n):
    conditions = list(adata.obs['Condition'].astype('category').cat.categories)
    datasets = list(adata.obs['dataset'].astype('category').cat.categories)
    sicks = list(adata.obs['sick'].astype('category').cat.categories)
    cell_types = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']

    from collections import defaultdict
    adata_dic = {}

    for dataset in datasets:                                                                                                                          
        for condition in conditions:
            adata_dic[(condition, dataset)] = adata[adata.obs['Condition'] == condition][adata[adata.obs['Condition'] == condition].obs['dataset'] == dataset]
        for sick in sicks:
            adata_dic[(sick, dataset)] = adata[adata.obs['sick'] == sick][adata[adata.obs['sick'] == sick].obs['dataset'] == dataset]
        for cell_type in cell_types:
            adata_dic[(cell_type, dataset)] = adata[adata.obs['cell_type'] == cell_type][adata[adata.obs['cell_type'] == cell_type].obs['dataset'] == dataset]
    
    fig,axs = plt.subplots(2, 2, figsize=(10,5), gridspec_kw={'wspace':-0.2, 'hspace':0.3}, facecolor='white', dpi=300, sharex=True, sharey=True)
    for age, ax_row in zip(['child', 'adult'], axs):
        for cond, ax in zip([cond1, cond2], ax_row):
            sc.pl.stacked_violin(adata_dic[(cond, age)], genes, groupby=gb, var_group_rotation=45, ax=ax, swap_axes=True, show=False, cmap='RdBu_r', vmin=0, vmax=vmax_n)
            ax.set_title(cond + ' ' + age, fontsize=15, y=0.6, loc='left')
    
    axes = fig.get_axes()  
    axes[4].set_xticks([])
    axes[7+len(genes)].set_xticks([])
    axes[6+len(genes)].remove()
    axes[12+3*len(genes)].remove()
    
    
    fig.tight_layout()
    return {'figure':fig, 'subplots':axs}
