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

def get_inters_fra_pair(data_ct, cts, fra_pair_cut_off, fra_exp_cut_off):    
    fn_int = '/home/yike/phd/dengue/data/interaction_source_file/inters_YK_20220324.tsv'
    interactions = pd.read_csv(fn_int, sep='\t')[['genesymbol_intercell_source', 'genesymbol_intercell_target']]
    res = []
    for _, row in interactions.iterrows():
        ga = row['genesymbol_intercell_source']
        gb = row['genesymbol_intercell_target']
        if (ga not in genes) | (gb not in genes):
            continue
        for cst in cts:
            cst_med = {gene: data_ct.loc[cst, gene]['med_pair'] for gene in [ga, gb]} 
            cst_fra_pair = {gene: data_ct.loc[cst, gene]['fra_pair'] for gene in [ga, gb]}
            cst_neg_fra_pair = {gene: data_ct.loc[cst, gene]['neg_fra_pair'] for gene in [ga, gb]}
            cst_fra = {(gene, cd): data_ct.loc[cst, gene][cd+'_fra'] for gene in [ga, gb] for cd in ['S', 'NS']}
            for ct in cts:
                ct_med = {gene: data_ct.loc[ct, gene]['med_pair'] for gene in [ga, gb]}
                ct_fra_pair = {gene: data_ct.loc[ct, gene]['fra_pair'] for gene in [ga, gb]}
                ct_neg_fra_pair = {gene: data_ct.loc[ct, gene]['neg_fra_pair'] for gene in [ga, gb]}
                ct_fra = {(gene, cd): data_ct.loc[ct, gene][cd+'_fra'] for gene in [ga, gb] for cd in ['S', 'NS']}

                ##################### upregulated
                if (cst_fra_pair[ga] >= fra_pair_cut_off) & (ct_fra_pair[gb] >= fra_pair_cut_off) & (cst_fra[ga, 'S'] > fra_exp_cut_off) & (ct_fra[gb, 'S'] > fra_exp_cut_off):
                    res.append([
                        ga, cst, cst_med[ga], cst_fra_pair[ga], cst_neg_fra_pair[ga], cst_fra[ga, 'S'], cst_fra[ga, 'NS'],
                        gb, ct, ct_med[gb], ct_fra_pair[gb], ct_neg_fra_pair[gb], ct_fra[gb, 'S'], ct_fra[gb, 'NS'], 'up'
                    ])
                ##################### downregulated
                if (cst_neg_fra_pair[ga] <= 1 - fra_pair_cut_off) & (ct_neg_fra_pair[gb] <= 1 - fra_pair_cut_off) & (cst_fra[ga, 'NS'] > fra_exp_cut_off) & (ct_fra[gb, 'NS'] > fra_exp_cut_off):
                    res.append([
                        ga, cst, cst_med[ga], cst_fra_pair[ga], cst_neg_fra_pair[ga], cst_fra[ga, 'S'], cst_fra[ga, 'NS'],
                        gb, ct, ct_med[gb], ct_fra_pair[gb], ct_neg_fra_pair[gb], ct_fra[gb, 'S'], ct_fra[gb, 'NS'], 'down'
                    ])
                ##################### mixregulated
                if (cst_fra_pair[ga] >= fra_pair_cut_off) & (ct_neg_fra_pair[gb] <= 1 - fra_pair_cut_off) & (cst_fra[ga, 'S'] > fra_exp_cut_off) & (ct_fra[gb, 'NS'] > fra_exp_cut_off):
                    res.append([
                        ga, cst, cst_med[ga], cst_fra_pair[ga], cst_neg_fra_pair[ga], cst_fra[ga, 'S'], cst_fra[ga, 'NS'],
                        gb, ct, ct_med[gb], ct_fra_pair[gb], ct_neg_fra_pair[gb], ct_fra[gb, 'S'], ct_fra[gb, 'NS'], 'mix'
                    ])

                if (cst_neg_fra_pair[ga] <= 1 - fra_pair_cut_off) & (ct_fra_pair[gb] >= fra_pair_cut_off) & (cst_fra[ga, 'NS'] > fra_exp_cut_off) & (ct_fra[gb, 'S'] > fra_exp_cut_off):
                    res.append([
                        gb, ct, ct_med[gb], ct_fra_pair[gb], ct_neg_fra_pair[gb], ct_fra[gb, 'S'], ct_fra[gb, 'NS'], 
                        ga, cst, cst_med[ga], cst_fra_pair[ga], cst_neg_fra_pair[ga], cst_fra[ga, 'S'], cst_fra[ga, 'NS'], 'mix'
                    ])

    res = pd.DataFrame(res, columns=[
        'ga', 'csta', 'ga_med_pair', 'ga_fra_pair','ga_neg_fra_pair','ga_SD_fra', 'ga_D_fra', 
        'gb', 'cstb', 'gb_med_pair', 'gb_fra_pair','gb_neg_fra_pair','gb_SD_fra', 'gb_D_fra',
        'inter_type'])
    
    return res

def get_inters_med_pair(data_ct, cts, med_pair_cut_off, fra_exp_cut_off):    
    fn_int = '/home/yike/phd/dengue/data/interaction_source_file/inters_YK_20220324.tsv'
    interactions = pd.read_csv(fn_int, sep='\t')[['genesymbol_intercell_source', 'genesymbol_intercell_target']]
    res = []
    for _, row in interactions.iterrows():
        ga = row['genesymbol_intercell_source']
        gb = row['genesymbol_intercell_target']
        if (ga not in genes) | (gb not in genes):
            continue
        for cst in cts:
            cst_med = {gene: data_ct.loc[cst, gene]['med_pair'] for gene in [ga, gb]} 
            cst_fra_pair = {gene: data_ct.loc[cst, gene]['fra_pair'] for gene in [ga, gb]} 
            cst_neg_fra_pair = {gene: data_ct.loc[cst, gene]['neg_fra_pair'] for gene in [ga, gb]}
            cst_fra = {(gene, cd): data_ct.loc[cst, gene][cd+'_fra'] for gene in [ga, gb] for cd in ['S', 'NS']}
            for ct in cts:
                ct_med = {gene: data_ct.loc[ct, gene]['med_pair'] for gene in [ga, gb]}
                ct_fra_pair = {gene: data_ct.loc[ct, gene]['fra_pair'] for gene in [ga, gb]}
                ct_neg_fra_pair = {gene: data_ct.loc[ct, gene]['neg_fra_pair'] for gene in [ga, gb]}
                ct_fra = {(gene, cd): data_ct.loc[ct, gene][cd+'_fra'] for gene in [ga, gb] for cd in ['S', 'NS']}
                
                #####################
                if (cst_med[ga] > med_pair_cut_off) & (ct_med[gb] > med_pair_cut_off) & (cst_fra[ga, 'S'] > fra_exp_cut_off) & (ct_fra[gb, 'S'] > fra_exp_cut_off):
                    res.append([
                        ga, cst, cst_med[ga], cst_fra_pair[ga], cst_neg_fra_pair[ga], cst_fra[ga, 'S'], cst_fra[ga, 'NS'],
                        gb, ct, ct_med[gb], ct_fra_pair[gb], ct_neg_fra_pair[gb], ct_fra[gb, 'S'], ct_fra[gb, 'NS'], 'up'
                    ])
                #####################
                if (cst_med[ga] < -med_pair_cut_off) & (ct_med[gb] < -med_pair_cut_off) & (cst_fra[ga, 'NS'] > fra_exp_cut_off) & (ct_fra[gb, 'NS'] > fra_exp_cut_off):
                    res.append([
                        ga, cst, cst_med[ga], cst_fra_pair[ga], cst_neg_fra_pair[ga], cst_fra[ga, 'S'], cst_fra[ga, 'NS'],
                        gb, ct, ct_med[gb], ct_fra_pair[gb], ct_neg_fra_pair[gb], ct_fra[gb, 'S'], ct_fra[gb, 'NS'], 'down'
                    ])
                #####################
                if (cst_med[ga] > med_pair_cut_off) & (ct_med[gb] < -med_pair_cut_off) & (cst_fra[ga, 'S'] > fra_exp_cut_off) & (ct_fra[gb, 'NS'] > fra_exp_cut_off):
                    res.append([
                        ga, cst, cst_med[ga], cst_fra_pair[ga], cst_neg_fra_pair[ga], cst_fra[ga, 'S'], cst_fra[ga, 'NS'],
                        gb, ct, ct_med[gb], ct_fra_pair[gb], ct_neg_fra_pair[gb], ct_fra[gb, 'S'], ct_fra[gb, 'NS'], 'mix'
                    ])

                if (cst_med[ga] < -med_pair_cut_off) & (ct_med[gb] > med_pair_cut_off) & (cst_fra[ga, 'NS'] > fra_exp_cut_off) & (ct_fra[gb, 'S'] > fra_exp_cut_off):
                    res.append([
                        gb, ct, ct_med[gb], ct_fra_pair[gb], ct_neg_fra_pair[gb], ct_fra[gb, 'S'], ct_fra[gb, 'NS'], 
                        ga, cst, cst_med[ga], cst_fra_pair[ga], cst_neg_fra_pair[ga], cst_fra[ga, 'S'], cst_fra[ga, 'NS'], 'mix'
                    ])
    res = pd.DataFrame(res, columns=[
        'ga', 'csta', 'ga_med_pair', 'ga_fra_pair','ga_neg_fra_pair','ga_SD_fra', 'ga_D_fra', 
        'gb', 'cstb', 'gb_med_pair', 'gb_fra_pair','gb_neg_fra_pair','gb_SD_fra', 'gb_D_fra',
        'inter_type'])
    
    return res

def randomization(genes, adata_kid):
    log2fc = defaultdict(list)
    r = defaultdict(list)
    pvalue = {}
    
    adata_children = adata_kid[adata_kid.obs['Condition'].isin(['S_dengue', 'dengue'])].copy()
    
    for i, inter in enumerate(genes):
        if i % 100 == 0:
            print(i)
        ga = inter['ga']
        cta = inter['cta']
        gb = inter['gb']
        ctb = inter['ctb']

        csts = adata_kid.obs['cell_subtype_new'].unique().tolist()
        csts.remove('doublets')
        csts.remove('megakaryocytes')
        csts.remove('plasmacytoid DCs')
        csts.remove('conventional DCs')
        csts.remove('unknown')
        
        ct_obs = {x: ['cell_subtype_new', 'cell_type_new'][x not in csts] for x in [cta, ctb]}
        
        adata_g = {}
        for gene, ct in zip([ga, gb], [cta, ctb]):
            adata_ct = adata_children[adata_children.obs[ct_obs[ct]] == ct]
            adata_g[gene] = adata_ct[:, gene]  
            
        avg = {(gene, cd): adata_g[gene][adata_g[gene].obs['Condition'] == cd].X.toarray().mean() for gene in [ga, gb] for cd in ['S_dengue', 'dengue']}
        lfc = {gene: np.log2(avg[gene, 'S_dengue'] + 0.1) - np.log2(avg[gene, 'dengue'] + 0.1) for gene in [ga, gb]}
        log2fc[(ga, cta, gb, ctb)].append([lfc[ga], lfc[gb]])
        r0 = (float(lfc[ga])**2 + float(lfc[gb])**2)**0.5
        r[(ga, cta, gb, ctb)].append(r0)
        
        p = 0
        for i in range(1000):
            adata_i = adata_g
            raw = {gene: adata_i[gene].obs['Condition'].tolist() for gene in [ga, gb]}
            for gene in [ga, gb]:
                random.shuffle(raw[gene])
                adata_i[gene].obs['Condition'] = raw[gene]
            avg_i = {(gene, cd): adata_i[gene][adata_i[gene].obs['Condition'] == cd].X.toarray().mean() for gene in [ga, gb] for cd in ['S_dengue', 'dengue']}
            log2fc_i = {gene: np.log2(avg_i[gene, 'S_dengue'] + 0.1) - np.log2(avg_i[gene, 'dengue'] + 0.1) for gene in [ga, gb]}
            log2fc[(ga, cta, gb, ctb)].append([log2fc_i[ga], log2fc_i[gb]])
            ri = (float(log2fc_i[ga])**2 + float(log2fc_i[gb])**2)**0.5
            r[(ga, cta, gb, ctb)].append(ri)
            if ri >= r0:
                p += 1
        pvalue[(ga, cta, gb, ctb)] = p * 0.001

    res = pd.DataFrame([])
    for key in log2fc.keys():
        log2fc[key] = pd.DataFrame(log2fc[key], columns = ['log2fc_ga', 'log2fc_gb'])
        log2fc[key]['r'] = r[key]
        log2fc[key]['pvalue'] = pvalue[key]
        for j, s in enumerate(['ga', 'cta', 'gb', 'ctb']):
            log2fc[key][s] = key[j]
        res = pd.concat([res, log2fc[key]])
    res = res.set_index(['ga', 'cta', 'gb', 'ctb'])
        
    return res

def ran_filter(ran_res, genes):
    sig_res = pd.DataFrame([], columns = ['log2fc_ga', 'log2fc_gb', 'r', 'pvalue', 'ga', 'csta', 'gb', 'cstb'])
    i = 0
    
    for inter in genes:
        ga = inter['ga']
        csta = inter['cta']
        gb = inter['gb']
        cstb = inter['ctb']
        if ran_res.loc[ga, csta, gb, cstb][:1]['pvalue'][0] == 0:
            i += 1
            loc = ran_res.loc[ga, csta, gb, cstb][:1].loc[ga, csta, gb, cstb].tolist()
            for j in [ga, csta, gb, cstb]:
                loc.append(j)
            sig_res.loc[i] = loc
    return sig_res

def randomization_plot(fdn, ran_filter, ran_res):
    ran_genes = [{'ga': idx[0], 'cta': idx[1], 'gb': idx[2], 'ctb': idx[3]} for idx in ran_filter.index]
    for inter in ran_genes:
        fig, ax = plt.subplots(figsize=[3, 3], dpi=300)
        ga = inter['ga']
        cta = inter['cta']
        gb = inter['gb']
        ctb = inter['ctb']

        log2fc = ran_res.loc[ga, cta, gb, ctb]

        x0 = log2fc['log2fc_ga'].tolist()[0]
        y0 = log2fc['log2fc_gb'].tolist()[0]    

        x = log2fc['log2fc_ga'].tolist()[1:]
        y = log2fc['log2fc_gb'].tolist()[1:]

        ax.scatter(x0, y0, c='r', s=10, label='original data')
        ax.scatter(x, y, c='gray', s=10, label='randomized data', alpha=0.5)
        ax.legend(loc='lower left')
        
        ax.axvline(0, c='gray', zorder=-3, lw=0.5)
        ax.axhline(0, c='gray', zorder=-3, lw=0.5)
        ax.set_ylim(-4.5, 4.5)
        ax.set_xlim(-4.5, 4.5)
        ax.set_xlabel('Log2fc of ' + ga + ' in ' + cta.replace('_', ' '))
        ax.set_ylabel('Log2fc of ' + gb + ' in ' + ctb.replace('_', ' '))
        p = ran_res.loc[ga, cta, gb, ctb]['pvalue'][0]
        if p == 0:
            ax.set_title('p value < 0.001')
        else:
            ax.set_title('p value: ' + str(p))
        plt.savefig(fdn + ga + '_in_' + cta.replace('/', '_') + '&' + gb + '_in_' + ctb.replace('/', '_') + '.png', bbox_inches = 'tight')
#######################################################################

def violin(gene, cell_type, fdn=None):
    from scipy.stats import gaussian_kde as gs_kde

    #adata_v = adata_kid[adata_kid.obs['Condition'].isin(['S_dengue', 'dengue'])]
    adata_ct = adata_kid[adata_kid.obs['cell_type_new'] == cell_type]

    SD_IDs = list(adata_ct[adata_ct.obs['Condition'] == 'S_dengue'].obs['ID'].astype('category').cat.categories)
    D_IDs = list(adata_ct[adata_ct.obs['Condition'] == 'dengue'].obs['ID'].astype('category').cat.categories)

    df = {}
    for ID in D_IDs + SD_IDs: # from dengue to severe dengue
        if ID in D_IDs:
            ID_info = adata_ct[(adata_ct.obs['ID'] == ID) & (adata_ct.obs['Condition'] == 'dengue')][:, gene].X.toarray()[:, 0]
        else:
            ID_info = adata_ct[(adata_ct.obs['ID'] == ID) & (adata_ct.obs['Condition'] == 'S_dengue')][:, gene].X.toarray()[:, 0]
        df[ID] = (ID_info + 1).tolist()

    ##########################################################################
    def violin_show(df, n):
        fig, ax = plt.subplots(figsize=(5, 1), dpi=300)
#         hline = np.log10(mean(list(ave.values())))
        hline = np.log10(adata_ct[adata_ct.obs['Condition'].isin(['S_dengue', 'dengue'])][:, gene].X.toarray()[:, 0].mean() + 1)
        ax.axhline(hline, color='gray', zorder=0.5, lw=0.5, ls='--') # average expression line 
        for i, (key, col) in enumerate(df.items()):
            xmin = 0
            xmax = 6
            x = np.logspace(xmin, xmax, 1000)
            x = np.log10(x)

            num = np.log10(col)
            if num.max() == num.min():
                y = np.zeros(len(x))
                y[0] = 1
            else:
                kde = gs_kde(num, bw_method=1, weights=None) 
                #weights: f None (default), the samples are assumed to be equally weighted
                y = kde(x)

            if key in D_IDs:
                if key == D_IDs[0]:
                    ax.fill_betweenx(x, -y+5.5*i, y+5.5*i, facecolor='c', edgecolor='black', lw=0.4, label='D')
                else:
                    ax.fill_betweenx(x, -y+5.5*i, y+5.5*i, facecolor='c', edgecolor='black', lw=0.4)
            else:
                if key == SD_IDs[0]:
                    ax.fill_betweenx(x, -y+5.5*i, y+5.5*i, facecolor='pink', edgecolor='black', lw=0.4, label='SD')
                else:  
                    ax.fill_betweenx(x, -y+5.5*i, y+5.5*i, facecolor='pink', edgecolor='black', lw=0.4)

        ax.legend(bbox_to_anchor=(1.21, 1.05))
        ax.set_yticks([0, 1, 2])
        ax.set_yticklabels(['0.1', '$10$', '$10^2$'])
        ax.set_ylim([-0.1, 2.1])
        #ax.set_xticks([a * 5.5 for a in range(len(D_IDs + SD_IDs))])
        #ax.set_xticklabels(df.keys(), rotation=90)
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.set_ylabel('Gene exp\n[cpm+1]')

        ax.set_title('%s in %s'%(gene, cell_type.replace('_', ' ')))

        if fdn is not None:
            fig.savefig(fdn + '%s_in_%s_%s.png'%(gene, cell_type, n), bbox_inches = 'tight')

    ########################################################################################
    ave = {ID: mean(val) for ID, val in df.items()}
    ave = {k: v for k, v in sorted(ave.items(), key=lambda item: item[1])}
    df = {ID: df[ID] for ID in ave.keys()}
    ######################################################
    new_keys = [ID for ID in ave.keys() if ID in D_IDs] + [ID for ID in ave.keys() if ID in SD_IDs]
    df_new = {ID: df[ID] for ID in new_keys}

    violin_show(df, '1')
    violin_show(df_new, '2')

def inter_number(adata, inters_df, vmax, trend, cmap):

    it_n = inters_df.groupby(['csta', 'cstb']).size().unstack(fill_value=0)
        
    cell_types = adata.obs['cell_type_new'].unique().tolist()
    cell_types.remove('doublets')
    cell_types.remove('unknown')
    
    #cell_types = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
    idx_ap = list(set(cell_types) - set(it_n.index.tolist()))
    col_ap = list(set(cell_types) - set(it_n.columns.tolist()))
    for idx in idx_ap:
        it_n.loc[idx] = [0] * it_n.shape[1]
    for col in col_ap:
        it_n[col] = [0] * it_n.shape[0]
    
    it_n = it_n[it_n.index]

    pairs = []
    for i, cell_type in enumerate(cell_types):
        cts = cell_types[i+1:]
        for ct in cts:
            pairs.append([cell_type, ct])

    for [cta, ctb] in pairs:
        it_n.loc[ctb][cta] = it_n.loc[cta][ctb]

    from scipy.spatial.distance import pdist
    from scipy.cluster.hierarchy import linkage, leaves_list
    import matplotlib.patches as mpatches
    lkg_idx = linkage(pdist(it_n.values), optimal_ordering=True)
    best_idx = leaves_list(lkg_idx)
    best_idx = it_n.index[best_idx].tolist()

    if it_n.loc[best_idx[-1]][best_idx[-1]] > it_n.loc[best_idx[0]][best_idx[0]]:
        best_idx = best_idx
    else:
        best_idx.reverse()

    it_n = it_n.loc[best_idx]
    it_n = it_n[best_idx]

    fig, ax = plt.subplots(figsize=[3, 2], dpi=300)
    sns.heatmap(it_n.T, ax=ax, cmap=cmap, linecolor='w', linewidths=1, vmin=0, vmax=vmax) # cmap='plasma', 'magma'

    for x in range(len(best_idx)):
        for y in range(len(best_idx)):
            if y < x:
                dots = [[x, y],
                        [x, y+1],
                        [x+1, y+1],
                        [x+1, y],
                ]
                e = mpatches.Polygon(np.array(dots), color='w')
                ax.add_patch(e)

    ax.axvline(0, c='black')
    ax.axhline(0, c='black')

    ax.axvline(len(cell_types), c='black')
    ax.axhline(len(cell_types), c='black')
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    xlabels = [i.get_text() for i in ax.get_xticklabels()]
    for i, ct in enumerate(xlabels):
        if ct in ['B_cells', 'T_cells', 'NK_cells']:
            xlabels[i] = ct.replace('_', ' ')
    ax.set_xticklabels(xlabels)
    #xlabels.reverse()
    ax.set_yticklabels(xlabels)
    ax.text(len(cell_types) + 2.05, 3.8, 'Number of interactions', verticalalignment='center', rotation=90)
    ax.set_title('%sregulated'%trend)
    return {'figure': fig, 'ax': ax}

def inter_mix_number(adata, inters_df, vmax, cmap):

    cell_types = adata.obs['cell_type_new'].unique().tolist()
    cell_types.remove('doublets')
    cell_types.remove('unknown')
    
    it_n = inters_df.groupby(['csta', 'cstb']).size().unstack(fill_value=0)

    idx_ap = list(set(cell_types) - set(it_n.index.tolist()))
    col_ap = list(set(cell_types) - set(it_n.columns.tolist()))
    for idx in idx_ap:
        it_n.loc[idx] = [0] * it_n.shape[1]
    for col in col_ap:
        it_n[col] = [0] * it_n.shape[0]

    from scipy.spatial.distance import pdist
    from scipy.cluster.hierarchy import linkage, leaves_list
    import matplotlib.patches as mpatches
    lkg_idx = linkage(pdist(it_n.values), optimal_ordering=True)
    best_idx = leaves_list(lkg_idx)
    best_idx = it_n.index[best_idx].tolist()

    if it_n.loc[best_idx[-1]].sum() > it_n.loc[best_idx[0]].sum():
        best_idx = best_idx
    else:
        best_idx.reverse()
    
    it_n = it_n.loc[best_idx]
    lkg_col = linkage(pdist(it_n.T.values), optimal_ordering=True)
    best_col = leaves_list(lkg_col)
    best_col = it_n.T.index[best_col].tolist()
    
    if it_n[best_col[-1]].sum() > it_n[best_col[0]].sum():
        best_col = best_col
    else:
        best_col.reverse()

    it_n = it_n[best_col]

    fig, ax = plt.subplots(figsize=[3, 2], dpi=300)
    sns.heatmap(it_n.T, ax=ax, cmap=cmap, linecolor='w', linewidths=1, vmin=0, vmax=vmax) # cmap='plasma', 'magma'
    ax.set_xlabel(None)
    ax.set_ylabel(None)

    ax.axvline(0, c='black')
    ax.axhline(0, c='black')

    ax.axvline(len(cell_types), c='black')
    ax.axhline(len(cell_types), c='black')

    xlabels = [i.get_text() for i in ax.get_xticklabels()]
    ylabels = [i.get_text() for i in ax.get_yticklabels()]
    for i, ct in enumerate(xlabels):
        if ct in ['B_cells', 'T_cells', 'NK_cells']:
            xlabels[i] = ct.replace('_', ' ')
            
    for i, ct in enumerate(ylabels):
        if ct in ['B_cells', 'T_cells', 'NK_cells']:
            ylabels[i] = ct.replace('_', ' ')        
    
    ax.text(len(cell_types) + 2.05, 3.8, 'Number of interactions', verticalalignment='center', rotation=90)
    ax.set_title('%sregulated (Down vs Up)'%'Mix')

    ax.set_xticklabels(xlabels, c='k')
    ax.set_yticklabels(ylabels, c='k')

    return {'figure': fig, 'ax': ax}

def cst_number(Mon_cst, color):
    cst_its = pd.read_csv('/home/yike/phd/dengue/data/tables/dataset_20211001/' + 'cst_inters_in_cts_strict.tsv', sep='\t', index_col=['ga', 'csta', 'gb', 'cstb'])
    cst_its = cst_its[~ cst_its.duplicated()]
    cst_inters = cst_its.index

    Mon_n = {cst: 0 for cst in Mon_cst}

    for inter in cst_inters:
        for cst in Mon_cst:
            if inter[1] == cst:
                Mon_n[cst] += 1
            elif inter[3] == cst:
                Mon_n[cst] += 1

    Mon_it_n = pd.DataFrame(Mon_n.values(), columns=['Number of interactions'])
    Mon_it_n['Cell subtype'] = [cst.replace('_', ' ') for cst in Mon_n.keys()]
    Mon_it_n.sort_values('Number of interactions', ascending=False, inplace=True)
    
    fig, ax = plt.subplots(figsize=[0.9, 1], dpi=300)
    sns.barplot(data=Mon_it_n, x='Cell subtype', y='Number of interactions', ax=ax, color=color)
    ax.set_xticklabels([ct.replace('_', ' ') for ct in Mon_cst], rotation=90)
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    ylim = (Mon_it_n['Number of interactions'].max() // 5 + 1) * 5 + 1
    ax.set_ylim(0, ylim)
    ax.set_yticks(range(0, ylim, 5))
    ax.set_yticklabels(range(0, ylim, 5))
    
    return {'figure': fig, 'ax': ax}

def s_mushrooms(genes, vmax=3):
    '''
    genes = [{'ITGAX': ['B_cells', 'NK_cells'],
          'ITGB2': ['cDCs'],
          'ICAM1': ['Plasmablasts']},
         {'CCL4L2': ['Monocytes'], 'VSIR': ['pDCs']}]
    '''
    from matplotlib.patches import Wedge
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import math
    import numpy as np
    import pandas as pd
    import itertools

    #save_tabels = '/home/yike/phd/dengue/data/tables/dataset_20211001/'

    conditions = ['S_dengue', 'dengue']
    cmap = plt.cm.get_cmap('viridis')
    vmin = -1
    threshold = 0.05
    frac_ct = pd.read_csv('/home/yike/phd/dengue/data/tables/dataset_20211001/' + 'ct_fra_gene_cut_0.tsv', 
                         sep='\t', index_col=['cell_type_new', 'condition', 'gene'], squeeze=True)
    avg_ct = pd.read_csv('/home/yike/phd/dengue/data/tables/dataset_20211001/' + 'ct_avg_gene_cut_0.tsv', 
                         sep='\t', index_col=['cell_type_new', 'condition', 'gene'], squeeze=True)

    frac_cst = pd.read_csv('/home/yike/phd/dengue/data/tables/dataset_20211001/' + 'cst_fra_gene_cut_0.tsv', 
                          sep='\t', index_col=['cell_subtype_new', 'condition', 'gene'], squeeze=True)
    avg_cst = pd.read_csv('/home/yike/phd/dengue/data/tables/dataset_20211001/' + 'cst_avg_gene_cut_0.tsv', 
                          sep='\t', index_col=['cell_subtype_new', 'condition', 'gene'], squeeze=True)

    yl = sum([len(list(itertools.chain.from_iterable(genesi.values()))) for genesi in genes])
    fig = plt.figure(figsize=((1 + 0.8 * 2) * 0.6, (1 + yl)* 0.6), dpi=300)

    grid = plt.GridSpec(yl , 2, wspace=0.1, hspace=0.1)
    
    cell_types = ['B_cells',
                 'T_cells',
                 'NK_cells',
                 'Monocytes',
                 'Plasmablasts',
                 'plasmacytoid DCs',
                 'conventional DCs',
                 'megakaryocytes']
    
    cell_subtypes = ['naive B cells',
                     'activated B cells',
                     'memory B cells',
                     'XCL_low NK cells',
                     'XCL_high NK cells',
                     'CD8+ effector T cells',
                     'CD8+ naive/memory T cells',
                     'CD4+ T cells',
                     'macrophages',
                     'non_classical monocytes',
                     'classical monocytes',
                     'non_cycling Plasmablasts',
                     'cycling Plasmablasts',
                     ]

    axs = []
    for i in range(len(genes)):
         axs.append(plt.subplot(grid[sum(len(list(itertools.chain.from_iterable(genesi.values()))) for genesi in genes[: i]): sum(len(list(itertools.chain.from_iterable(genesi.values()))) for genesi in genes[: i+1]), 0: 1]))
    size_bar = plt.subplot(grid[0: 5, 1: 2])

    datap = []
    for genesi, ax in zip(genes, axs):
        cts = list(genesi.values())
        gs = list(genesi.keys())
        yticklabels = []
        for i, (csts, gene) in enumerate(zip(cts, gs)):
            for cst in csts:
                avgs = []
                yticklabels.append(gene + ' in\n' + cst.replace('_', ' '))
                for k, cond in enumerate(conditions):
                    if cst in cell_types:
                        fr = frac_ct.loc[(cst, cond, gene)]
                        av = np.log10(avg_ct.loc[(cst, cond, gene)] + 0.1)
                    elif cst in cell_subtypes:
                        fr = frac_cst.loc[(cst, cond, gene)]
                        av = np.log10(avg_cst.loc[(cst, cond, gene)] + 0.1)
                    avgs.append(av)
                    r = 0.5 * fr**0.3
                    color = cmap((min(vmax, av) - vmin) / (vmax - vmin))
                    theta0, theta1 = 180 * (k > 0), 180 + 180 * (k > 0)
                    datap.append({
                        'r': r,
                        'facecolor': color,
                        'center': (0, len(yticklabels)-1),
                        'theta': (theta0, theta1),
                        'ax': ax,
                    })
                if avgs[0] - avgs[1] > threshold:
                    datap[-2]['edgecolor'] = 'red'
                    datap[-1]['edgecolor'] = 'none'
                elif avgs[0] - avgs[1] < -threshold:
                    datap[-1]['edgecolor'] = 'red'
                    datap[-2]['edgecolor'] = 'none'
                else:
                    datap[-1]['edgecolor'] = 'none'
                    datap[-2]['edgecolor'] = 'none'   


        ax.set_yticks(np.arange(len(list(itertools.chain.from_iterable(genesi.values())))))
        ax.set_yticklabels(yticklabels)
        ax.set_ylim(-0.6, len(list(itertools.chain.from_iterable(genesi.values()))) - 0.4)        
        ax.set_xticks([])
        ax.set_xlim(-0.6, 1 - 0.4)

    for datum in datap:
        ax = datum['ax']
        r = datum['r']
        color = datum['facecolor']
        center = datum['center']
        theta0, theta1 = datum['theta']
        ec = datum['edgecolor']

        h = Wedge(
            center, r, theta0, theta1, facecolor=color, edgecolor=ec
        )
        ax.add_artist(h)
        ax.set_aspect(1)

    size_bar.set_ylim(-0.6, 5 - 0.4)        
    c = [(0.5, i) for i in range(5)]
    radius = [0.5 * fr**0.3 for fr in [0.05, 0.1, 0.2, 0.4, 0.8]]
    for c, r in zip(c, radius):
        e = Wedge(c, r, 0, 180, facecolor='gray',)
        size_bar.add_artist(e)
    size_bar.set_aspect(1)
    size_bar.set_yticks([])
    size_bar.set_yticks(range(5))
    size_bar.set_yticklabels(['5', '10', '20', '40', '80'])
    size_bar.yaxis.tick_right()
    size_bar.yaxis.set_label_position('right')
    size_bar.set_ylabel('Gene exp frac')
    size_bar.set_xticks([])
    size_bar.spines['bottom'].set_visible(False)
    size_bar.spines['top'].set_visible(False)
    size_bar.spines['right'].set_visible(False)
    size_bar.spines['left'].set_visible(False)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax) 
    cmap = plt.cm.get_cmap('viridis')
    position = fig.add_axes([0.7, 0.2, 0.05, 2/yl])
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=position, ax=axs[-1], label='Gene exp \n(log10[cpm+0.1])')

    fig.tight_layout()
    
    return {'figure': fig, 'axs': ax}

def randomization_gene(adata, genes, ):
    '''
    genes = {
    ct1: [],
    ct2: [],
    ct3: [],
    }
    '''
    import random
    log2fc = defaultdict(list)
    r = defaultdict(list)
    pvalue = {}

    csts = adata.obs['cell_subtype_new'].unique().tolist()
    csts.remove('doublets')
    csts.remove('megakaryocytes')
    csts.remove('plasmacytoid DCs')
    csts.remove('conventional DCs')
    csts.remove('unknown')

    for cell_type in genes.keys():
        ct_obs = 'cell_subtype_new'
        if cell_type not in csts:
            ct_obs = 'cell_type_new'

        adata_T = {cd: adata[(adata.obs[ct_obs] == cell_type) & (adata.obs['Condition'] == cd)] for cd in ['S_dengue', 'dengue']} # cell type tested for significance
        adata_C = {cd: adata[(adata.obs[ct_obs] != cell_type) & (adata.obs['Condition'] == cd)] for cd in ['S_dengue', 'dengue']} # control cell types

        for gene in genes[cell_type]:
            X_T = {cd: adata_T[cd][:, gene].X.toarray()[:, 0] for cd in ['S_dengue', 'dengue']}
            n_T = {cd: adata_T[cd].shape[0] for cd in ['S_dengue', 'dengue']}

            X_C = {cd: adata_C[cd][:, gene].X.toarray()[:, 0] for cd in ['S_dengue', 'dengue']}
            n_C = {cd: adata_C[cd].shape[0] for cd in ['S_dengue', 'dengue']}

            avg_T = {cd: X_T[cd].mean() for cd in ['S_dengue', 'dengue']}
            avg_C = {cd: X_C[cd].mean() for cd in ['S_dengue', 'dengue']}
            lfc0 = {t: np.log2(avg['S_dengue'] + 0.1) - np.log2(avg['dengue'] + 0.1) for avg, t in
                    zip([avg_T, avg_C], ['T', 'C'])}
            log2fc[(cell_type, gene)].append([lfc0['T'], lfc0['C']])
            r0 = (float(lfc0['T']) ** 2 + float(lfc0['C']) ** 2) ** 0.5
            r[(cell_type, gene)].append(r0)

            p = 0
            adata_Ti = adata[(adata.obs[ct_obs] == cell_type) & (adata.obs['Condition'].isin(['S_dengue', 'dengue']))]
            adata_Ci = adata[(adata.obs[ct_obs] != cell_type) & (adata.obs['Condition'].isin(['S_dengue', 'dengue']))]

            for i in range(1000):
                lfci = {}
                for name, n, adatai in zip(['T', 'C'], [n_T, n_C], [adata_Ti, adata_Ci]):
                    Xi = adatai[:, gene].X.toarray().tolist() 
                    idx_SDi = random.sample(range(len(Xi)), n['S_dengue'])
                    X_SDi = [Xi[i] for i in idx_SDi]
                    remain_idx = tuple(set(range(len(Xi))) - set(idx_SDi))
                    X_Di = [Xi[i] for i in remain_idx]
                    lfci[name] = np.log2(np.mean(X_SDi) + 0.1) - np.log2(np.mean(X_Di) + 0.1)

                log2fc[(cell_type, gene)].append([lfci['T'], lfci['C']])
                ri = (float(lfci['T']) ** 2 + float(lfci['C']) ** 2) ** 0.5
                r[(cell_type, gene)].append(ri)

                if ri >= r0:
                    p += 1
                if i % 100 == 0:
                    print(i)
            pvalue[(cell_type, gene)] = p * 0.001

    res = pd.DataFrame([])
    for key in log2fc.keys():
        log2fc[key] = pd.DataFrame(log2fc[key], columns = ['log2fc_ct', 'log2fc_unct'])
        log2fc[key]['r'] = r[key]
        log2fc[key]['pvalue'] = pvalue[key]
        for i, s in enumerate(['cell_type', 'gene']):
            log2fc[key][s] = key[i]
        res = pd.concat([res, log2fc[key]])
    res = res.set_index(['cell_type', 'gene'])

    return res

def randomization_gene_plot(fdn, ran_genes, ran_res):
    for inter in ran_genes:
        fig, ax = plt.subplots(figsize=[3, 3], dpi=300)
        ct = inter['cell_type']
        gene = inter['gene']

        log2fc = ran_res.loc[ct, gene]

        x0 = log2fc['log2fc_ct'].tolist()[0]
        y0 = log2fc['log2fc_unct'].tolist()[0]    

        x = log2fc['log2fc_ct'].tolist()[1:]
        y = log2fc['log2fc_unct'].tolist()[1:]

        ax.scatter(x0, y0, c='r', s=10, label='original data')
        ax.scatter(x, y, c='gray', s=10, label='randomized data', alpha=0.5)
        ax.legend(loc='lower left')
        
        ax.axvline(0, c='gray', zorder=-3, lw=0.5)
        ax.axhline(0, c='gray', zorder=-3, lw=0.5)
        ax.set_ylim(-4.5, 4.5)
        ax.set_xlim(-4.5, 4.5)
        ax.set_xlabel('Log2fc of ' + gene + ' in ' + ct.replace('_', ' '))
        ax.set_ylabel('Log2fc of ' + gene + ' in ' + 'bystanders')
        p = ran_res.loc[ct, gene]['pvalue'][0]
        if p == 0:
            ax.set_title('p value < 0.001')
        else:
            ax.set_title('p value: ' + str(p))
        plt.savefig(fdn + gene + '_in_' + ct.replace('/', '_') + '.png', bbox_inches = 'tight')
#######################################################################

def cst_mushrooms(genes, vmax=3):
    '''
    genes = [{'ITGAX': ['B_cells', 'NK_cells'],
          'ITGB2': ['cDCs'],
          'ICAM1': ['Plasmablasts']},
         {'CCL4L2': ['Monocytes'], 'VSIR': ['pDCs']}]
    '''
    from matplotlib.patches import Wedge
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import math
    import numpy as np
    import pandas as pd
    import itertools

    conditions = ['S_dengue', 'dengue']
    cmap = plt.cm.get_cmap('viridis')
    vmin = -1
    threshold = 0.05

    #save_tabels = '/home/yike/phd/dengue/data/tables/dataset_20211001/'
    frac_ct = pd.read_csv('/home/yike/phd/dengue/data/tables/dataset_20211001/' + 'ct_fra_gene_cut_0.tsv', 
                         sep='\t', index_col=['cell_type_new', 'condition', 'gene'], squeeze=True)
    avg_ct = pd.read_csv('/home/yike/phd/dengue/data/tables/dataset_20211001/' + 'ct_avg_gene_cut_0.tsv', 
                         sep='\t', index_col=['cell_type_new', 'condition', 'gene'], squeeze=True)

    frac_cst = pd.read_csv('/home/yike/phd/dengue/data/tables/dataset_20211001/' + 'cst_fra_gene_cut_0.tsv', 
                          sep='\t', index_col=['cell_subtype_new', 'condition', 'gene'], squeeze=True)
    avg_cst = pd.read_csv('/home/yike/phd/dengue/data/tables/dataset_20211001/' + 'cst_avg_gene_cut_0.tsv', 
                          sep='\t', index_col=['cell_subtype_new', 'condition', 'gene'], squeeze=True)

    yl = sum([len(list(itertools.chain.from_iterable(genesi.values()))) for genesi in genes])
    fig = plt.figure(figsize=((1 + 0.8 * 2) * 0.6, (1 + yl)* 0.6), dpi=300)

    grid = plt.GridSpec(yl , 2, wspace=0.1, hspace=0.1)
    
    cell_types = ['B_cells',
                 'T_cells',
                 'NK_cells',
                 'Monocytes',
                 'Plasmablasts',
                 'plasmacytoid DCs',
                 'conventional DCs',
                 'megakaryocytes']
    
    cell_subtypes = ['naive B cells',
                     'activated B cells',
                     'memory B cells',
                     'XCL_low NK cells',
                     'XCL_high NK cells',
                     'CD8+ effector T cells',
                     'CD8+ naive/memory T cells',
                     'CD4+ T cells',
                     'macrophages',
                     'non_classical monocytes',
                     'classical monocytes',
                     'non_cycling Plasmablasts',
                     'cycling Plasmablasts',
                     ]

    axs = []
    for i in range(len(genes)):
         axs.append(plt.subplot(grid[sum(len(list(itertools.chain.from_iterable(genesi.values()))) for genesi in genes[: i]): sum(len(list(itertools.chain.from_iterable(genesi.values()))) for genesi in genes[: i+1]), 0: 1]))
    size_bar = plt.subplot(grid[0: 5, 1: 2])

    datap = []
    for genesi, ax in zip(genes, axs):
        cts = list(genesi.values())
        gs = list(genesi.keys())
        yticklabels = []
        for i, (csts, gene) in enumerate(zip(cts, gs)):
            for cst in csts:
                avgs = []
                yticklabels.append(cst.replace('_', ' '))
                for k, cond in enumerate(conditions):
                    if cst in cell_types:
                        fr = frac_ct.loc[(cst, cond, gene)]
                        av = np.log10(avg_ct.loc[(cst, cond, gene)] + 0.1)
                    elif cst in cell_subtypes:
                        fr = frac_cst.loc[(cst, cond, gene)]
                        av = np.log10(avg_cst.loc[(cst, cond, gene)] + 0.1)
                    avgs.append(av)
                    r = 0.5 * fr**0.3
                    color = cmap((min(vmax, av) - vmin) / (vmax - vmin))
                    theta0, theta1 = 180 * (k > 0), 180 + 180 * (k > 0)
                    datap.append({
                        'r': r,
                        'facecolor': color,
                        'center': (0, len(yticklabels)-1),
                        'theta': (theta0, theta1),
                        'ax': ax,
                    })
                if avgs[0] - avgs[1] > threshold:
                    datap[-2]['edgecolor'] = 'red'
                    datap[-1]['edgecolor'] = 'none'
                elif avgs[0] - avgs[1] < -threshold:
                    datap[-1]['edgecolor'] = 'red'
                    datap[-2]['edgecolor'] = 'none'
                else:
                    datap[-1]['edgecolor'] = 'none'
                    datap[-2]['edgecolor'] = 'none'   

        axs[0].set_title(gene)
        ax.set_yticks(np.arange(len(list(itertools.chain.from_iterable(genesi.values())))))
        ax.set_yticklabels(yticklabels)
        ax.set_ylim(-0.6, len(list(itertools.chain.from_iterable(genesi.values()))) - 0.4)        
        ax.set_xticks([])
        ax.set_xlim(-0.6, 1 - 0.4)

    for datum in datap:
        ax = datum['ax']
        r = datum['r']
        color = datum['facecolor']
        center = datum['center']
        theta0, theta1 = datum['theta']
        ec = datum['edgecolor']

        h = Wedge(
            center, r, theta0, theta1, facecolor=color, edgecolor=ec
        )
        ax.add_artist(h)
        ax.set_aspect(1)

    size_bar.set_ylim(-0.6, 5 - 0.4)        
    c = [(0.5, i) for i in range(5)]
    radius = [0.5 * fr**0.3 for fr in [0.05, 0.1, 0.2, 0.4, 0.8]]
    for c, r in zip(c, radius):
        e = Wedge(c, r, 0, 180, facecolor='gray',)
        size_bar.add_artist(e)
    size_bar.set_aspect(1)
    size_bar.set_yticks([])
    size_bar.set_yticks(range(5))
    size_bar.set_yticklabels(['5', '10', '20', '40', '80'])
    size_bar.yaxis.tick_right()
    size_bar.yaxis.set_label_position('right')
    size_bar.set_ylabel('Gene exp frac')
    size_bar.set_xticks([])
    size_bar.spines['bottom'].set_visible(False)
    size_bar.spines['top'].set_visible(False)
    size_bar.spines['right'].set_visible(False)
    size_bar.spines['left'].set_visible(False)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax) 
    cmap = plt.cm.get_cmap('viridis')
    position = fig.add_axes([0.7, 0.2, 0.05, 2/yl])
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=position, ax=axs[-1], label='Gene exp \n(log10[cpm+0.1])')
    fig.tight_layout()
    
    return {'figure': fig, 'axs': ax}

def com_mushrooms(genes, cst_plots, vmax=3):

    from matplotlib.patches import Wedge
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import math
    import numpy as np
    import pandas as pd
    import itertools

    conditions = ['S_dengue', 'dengue']
    cmap = plt.cm.get_cmap('viridis')
    vmin = -1
    # vmax = 3
    threshold = 0.05

    #save_tabels = '/home/yike/phd/dengue/data/tables/dataset_20211001/'
    frac_ct = pd.read_csv('/home/yike/phd/dengue/data/tables/dataset_20211001/' + 'ct_fra_gene_cut_0.tsv', 
                         sep='\t', index_col=['cell_type_new', 'condition', 'gene'], squeeze=True)
    avg_ct = pd.read_csv('/home/yike/phd/dengue/data/tables/dataset_20211001/' + 'ct_avg_gene_cut_0.tsv', 
                         sep='\t', index_col=['cell_type_new', 'condition', 'gene'], squeeze=True)

    frac_cst = pd.read_csv('/home/yike/phd/dengue/data/tables/dataset_20211001/' + 'cst_fra_gene_cut_0.tsv', 
                          sep='\t', index_col=['cell_subtype_new', 'condition', 'gene'], squeeze=True)
    avg_cst = pd.read_csv('/home/yike/phd/dengue/data/tables/dataset_20211001/' + 'cst_avg_gene_cut_0.tsv', 
                          sep='\t', index_col=['cell_subtype_new', 'condition', 'gene'], squeeze=True)

    if len(cst_plots) < 5:
        length = 5
    else:
        length = len(cst_plots)

    fig = plt.figure(figsize=(2 + 0.5 * len(genes), 1 + length * 0.6), dpi=300)
    grid = plt.GridSpec(1 , len(genes) + 1, wspace=0.1, hspace=0.1)

    # fig, axs = plt.subplots(1, len(genes) + 1, 
    #                         figsize=((1 + 0.48 * len(genes)), (1 + length)* 0.6), dpi=300)

    cell_types = ['B_cells',
                 'T_cells',
                 'NK_cells',
                 'Monocytes',
                 'Plasmablasts',
                 'plasmacytoid DCs',
                 'conventional DCs',
                 'megakaryocytes']

    cell_subtypes = ['naive B cells',
                     'activated B cells',
                     'memory B cells',
                     'XCL_low NK cells',
                     'XCL_high NK cells',
                     'CD8+ effector T cells',
                     'CD8+ naive/memory T cells',
                     'CD4+ T cells',
                     'macrophages',
                     'non_classical monocytes',
                     'classical monocytes',
                     'non_cycling Plasmablasts',
                     'cycling Plasmablasts',
                     ]
    axs = []
    for i in range(len(genes)):
         axs.append(plt.subplot(grid[0: len(cst_plots), i: i+1]))

    axs.append(plt.subplot(grid[0: 5, len(genes): len(genes) + 1])) # size_bar
    #axs.append(plt.subplot(grid[0: 2, len(genes) + 1: len(genes) + 2])) # color_bar

    for gene, ax in zip(genes, axs[:-1]):
        datap = []
        yticklabels = []
        for cst in cst_plots:
            avgs = []
            yticklabels.append(cst.replace('_', ' '))
            for k, cond in enumerate(conditions):
                if cst in cell_types:
                    fr = frac_ct.loc[(cst, cond, gene)]
                    av = np.log10(avg_ct.loc[(cst, cond, gene)] + 0.1)
                elif cst in cell_subtypes:
                    fr = frac_cst.loc[(cst, cond, gene)]
                    av = np.log10(avg_cst.loc[(cst, cond, gene)] + 0.1)
                avgs.append(av)
                r = 0.5 * fr**0.3
                color = cmap((min(vmax, av) - vmin) / (vmax - vmin))
                theta0, theta1 = 180 * (k > 0), 180 + 180 * (k > 0)
                datap.append({
                    'r': r,
                    'facecolor': color,
                    'center': (0, len(yticklabels)-1),
                    'theta': (theta0, theta1),
                    'ax': ax,
                })
            if avgs[0] - avgs[1] > threshold:
                datap[-2]['edgecolor'] = 'red'
                datap[-1]['edgecolor'] = 'none'
            elif avgs[0] - avgs[1] < -threshold:
                datap[-1]['edgecolor'] = 'red'
                datap[-2]['edgecolor'] = 'none'
            else:
                datap[-1]['edgecolor'] = 'none'
                datap[-2]['edgecolor'] = 'none'   

        ax.set_title(gene)

        if ax == axs[0]:
            ax.set_yticks(np.arange(len(cst_plots)))
            ax.set_yticklabels(yticklabels)
        else:
            ax.set_yticks([])
            ax.set_yticklabels([])

        ax.set_ylim(-0.6, len(cst_plots) - 0.4)        
        ax.set_xticks([])
        ax.set_xlim(-0.6, 1 - 0.4)

        for datum in datap:
            ax = datum['ax']
            r = datum['r']
            color = datum['facecolor']
            center = datum['center']
            theta0, theta1 = datum['theta']
            ec = datum['edgecolor']

            h = Wedge(
                center, r, theta0, theta1, facecolor=color, edgecolor=ec
            )
            ax.add_artist(h)
            ax.set_aspect(1)

    axs[-1].set_ylim(-0.6, 5 - 0.4)        
    c = [(0.5, i) for i in range(5)]
    radius = [0.5 * fr**0.3 for fr in [0.05, 0.1, 0.2, 0.4, 0.8]]
    for c, r in zip(c, radius):
        e = Wedge(c, r, 0, 180, facecolor='gray',)
        axs[-1].add_artist(e)
    axs[-1].set_aspect(1)
    axs[-1].set_yticks([])
    axs[-1].set_yticks(range(5))
    axs[-1].set_yticklabels(['5', '10', '20', '40', '80'])
    axs[-1].yaxis.tick_right()
    axs[-1].yaxis.set_label_position('right')
    axs[-1].set_ylabel('Gene exp frac')
    axs[-1].set_xticks([])
    axs[-1].spines['bottom'].set_visible(False)
    axs[-1].spines['top'].set_visible(False)
    axs[-1].spines['right'].set_visible(False)
    axs[-1].spines['left'].set_visible(False)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax) 
    cmap = plt.cm.get_cmap('viridis')
    position = fig.add_axes([0.98, 0.4, 0.01, 0.3]) 
    # The dimensions [left, bottom, width, height] of the new Axes.
    # All quantities are in fractions of figure width and height.
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs[-1], cax=position, 
                        label='Gene exp \n(log10[cpm+0.1])')
    cbar.set_ticks(range(-1, 4))
    cbar.set_ticklabels(range(-1, 4))

    fig.tight_layout()
    
    return {'figure': fig, 'axs': axs}

def endo_flowers(genes):
    '''
    genes = [{'ITGAX': ['B_cells', 'NK_cells'],
          'ITGB2': ['Lymphatic ECs', 'Arterial ECs'],},
         {'CCL4L2': ['Monocytes'], 'VSIR': ['Endothelial cells']}]
    '''
    from matplotlib.patches import Wedge
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import math
    import numpy as np
    import pandas as pd
    import itertools

    conditions = ['S_dengue', 'dengue']
    cmap = plt.cm.get_cmap('viridis')
    vmin, vmax = -1, 3
    threshold = 0.05
    frac_cst = pd.read_csv('/home/yike/phd/dengue/data/tables/cell_subtype/original/fra.tsv', 
                      sep='\t', index_col=['cell_subtype_2', 'condition', 'gene'], squeeze=True)
    avg_cst = pd.read_csv('/home/yike/phd/dengue/data/tables/cell_subtype/original/avg.tsv', sep='\t',
                      index_col=['cell_subtype_2', 'condition', 'gene'], squeeze=True)
    
    frac_ct = pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/20210625_figure_4_code/fra.tsv', 
                          index_col=['cell_type', 'condition', 'gene'], squeeze=True)
    avg_ct = pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/20210625_figure_4_code/exp.tsv', 
                         index_col=['cell_type', 'condition', 'gene'], squeeze=True)
    
    #fra_avg_endo = pd.read_csv('/home/yike/phd/dengue/data/tables/endos/fra_avg.tsv', sep='\t', index_col='gene')
    fra_avg_endo_cst = pd.read_csv('/home/yike/phd/dengue/data/tables/endos/fra_avg_cst.tsv',
                                  sep='\t', index_col=['cell_subtype', 'gene'])
    

    yl = sum([len(list(itertools.chain.from_iterable(genesi.values()))) for genesi in genes])
    fig = plt.figure(figsize=((1 + 0.8 * 2) * 0.6, (1 + yl)* 0.6), dpi=300)

    grid = plt.GridSpec(yl , 2, wspace=0.1, hspace=0.1)
    
    cell_types = ['B_cells', 'T_cells', 'NK_cells', 'cDCs', 'pDCs', 'Monocytes', 'Plasmablasts']
    cell_subtypes = ['Naive_B_cells', 'NK', 'CD4_T_cells', 'NKT', 'Memory_B_cells', 'CD8_T_cells',
                     'Macrophages', 'non_classical_monocytes', 'IgA', 'IgG1_proliferate', 'pDCs',
                     'cDC_IFN', 'Classical_monocytes', 'IgG1_IgG2', 'cDC2', 'IgG1', 'IgM', 'cDC1']

    axs = []
    for i in range(len(genes)):
         axs.append(plt.subplot(grid[sum(len(list(itertools.chain.from_iterable(genesi.values()))) for genesi in genes[: i]): sum(len(list(itertools.chain.from_iterable(genesi.values()))) for genesi in genes[: i+1]), 0: 1]))
    size_bar = plt.subplot(grid[0: 5, 1: 2])

    datap = []
    for genesi, ax in zip(genes, axs): # {'ITGAX': ['B_cells', 'NK_cells'], 'ITGB2': ['Lymphatic ECs', 'Arterial ECs'],}
        cts = list(genesi.values()) # [['B_cells', 'NK_cells'], ['Lymphatic ECs', 'Arterial ECs']]
        gs = list(genesi.keys()) # ['ITGAX', 'ITGB2']
        yticklabels = []
        
        for i, (csts, gene) in enumerate(zip(cts, gs)): # ['B_cells', 'NK_cells'], 'ITGAX'
            for cst in csts: # 'B_cells'
                avgs = []
                yticklabels.append(gene + ' in\n' + cst.replace('_', ' '))
                
                if cst in cell_types + cell_subtypes:
                    for j, cond in enumerate(conditions):
                        if cst in cell_types:
                            fr = frac_ct.loc[(cst, cond, gene)]
                            av = np.log10(avg_ct.loc[(cst, cond, gene)] + 0.1)
                        elif cst in cell_subtypes:
                            fr = frac_cst.loc[(cst, cond, gene)]
                            av = np.log10(avg_cst.loc[(cst, cond, gene)] + 0.1)

                        avgs.append(av)

                        r = 0.5 * fr**0.3
                        color = cmap((min(vmax, av) - vmin) / (vmax - vmin))
                        theta0, theta1 = 180 * (j > 0), 180 + 180 * (j > 0)
                        datap.append({
                            'r': r,
                            'facecolor': color,
                            'center': (0, len(yticklabels)-1),
                            'theta': (theta0, theta1),
                            'ax': ax,
                        })
                                       
                    if avgs[0] - avgs[1] > threshold:
                        datap[-2]['edgecolor'] = 'red'
                        datap[-1]['edgecolor'] = 'none'
                    elif avgs[0] - avgs[1] < -threshold:
                        datap[-1]['edgecolor'] = 'red'
                        datap[-2]['edgecolor'] = 'none'
                    else:
                        datap[-1]['edgecolor'] = color
                        datap[-2]['edgecolor'] = color
                    
                    
                elif cst == 'Endothelial cells':

                    for k, cst in enumerate(['Lymphatic ECs', 'Venous ECs', 'Capillary ECs', 'Arterial ECs']):
                        fr = fra_avg_endo_cst.loc[cst, gene]['fra']
                        av = np.log10(fra_avg_endo_cst.loc[cst, gene]['avg'] + 0.1)

                        r = 0.5 * fr**0.3
                        color = cmap((min(vmax, av) - vmin) / (vmax - vmin))
                        theta0, theta1 = 90 * k, 90 * (k + 1)
                        datap.append({
                            'r': r,
                            'facecolor': color,
                            'center': (0, len(yticklabels)-1),
                            'theta': (theta0, theta1),
                            'ax': ax,
                            'edgecolor': 'white'
                        })
                    
        ax.set_yticks(np.arange(len(list(itertools.chain.from_iterable(genesi.values())))))
        ax.set_yticklabels(yticklabels)
        ax.set_ylim(-0.6, len(list(itertools.chain.from_iterable(genesi.values()))) - 0.4)        
        ax.set_xticks([])
        ax.set_xlim(-0.6, 1 - 0.4)

    for datum in datap:
        ax = datum['ax']
        r = datum['r']
        color = datum['facecolor']
        center = datum['center']
        theta0, theta1 = datum['theta']
        ec = datum['edgecolor']

        h = Wedge(
            center, r, theta0, theta1, facecolor=color, edgecolor=ec, linewidth=0.4
        )
        ax.add_artist(h)
        ax.set_aspect(1)

    size_bar.set_ylim(-0.6, 5 - 0.4)        
    c = [(0.5, i) for i in range(5)]
    radius = [0.5 * fr**0.3 for fr in [0.05, 0.1, 0.2, 0.4, 0.8]]
    for c, r in zip(c, radius):
        e = Wedge(c, r, 0, 90, facecolor='gray',)
        size_bar.add_artist(e)
    size_bar.set_aspect(1)
    size_bar.set_yticks([])
    size_bar.set_yticks(range(5))
    size_bar.set_yticklabels(['5', '10', '20', '40', '80'])
    size_bar.yaxis.tick_right()
    size_bar.yaxis.set_label_position('right')
    size_bar.set_ylabel('Gene exp frac')
    size_bar.set_xticks([])
    #size_bar.axvline(0.5, c='white')
    size_bar.spines['bottom'].set_visible(False)
    size_bar.spines['top'].set_visible(False)
    size_bar.spines['right'].set_visible(False)
    size_bar.spines['left'].set_visible(False)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax) 
    cmap = plt.cm.get_cmap('viridis')
    position = fig.add_axes([0.7, 0.2, 0.05, 2/yl])
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=position, ax=axs[-1], label='Gene exp \n(log10[cpm+0.1])')

    fig.tight_layout()
    return {'fig': fig, 'ax': axs}