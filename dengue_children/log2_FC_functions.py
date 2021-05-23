import pandas as pd 
import numpy as np

import anndata
import scanpy as sc

import matplotlib.pyplot as plt
import seaborn as sns

def log2_FC(adata, cond1, cond2, datas):
    '''the log2 fold_change of genes by comparision between each cond1 and cond2, eg., 'sick' vs 'Healthy', 'DWS' vs 'Healthy'
    Args:
        adata: the adata used for analysis, eg., adata, adata_B_cells
        cond1, cond2: comparision between cond1 and cond2, eg., 'sick' vs 'Healthy', 'DWS' vs 'Healthy'
                i.e., log2_FC = np.log2(ffoldchange[cond1, datas] + 0.1) - np.log2(foldchange[cond2, datas] + 0.1)
        datas: dataset used for analysis, eg, 'child' or 'adult'
    Returns:
        dict containing DataFrame with log2 fold_change of genes by comparision between cond1 and cond2, keys are IDs [cond1_ID, cond2_ID]
    '''
    
    conditions = list(adata.obs['Condition'].astype('category').cat.categories)
    datasets = list(adata.obs['dataset'].astype('category').cat.categories)
    sicks = list(adata.obs['sick'].astype('category').cat.categories)

    from collections import defaultdict
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
            
    from collections import defaultdict
    
    fc = {}
    IDs = l_ID[(cond1, datas)] + l_ID[(cond2, datas)]
    for ID in IDs:
        if ID in l_ID[(cond1, datas)]:
            sick = cond1
        else:
            sick = cond2
        adata_ID = adata_dic[(sick, datas)][adata_dic[(sick, datas)].obs['ID'] == ID]
        fc[ID] = (adata_ID.X.toarray()).mean(axis=0)

    log2_fc = {}
    for ID1 in l_ID[(cond1, datas)]:
        for ID2 in l_ID[(cond2, datas)]:
            log2_fc[(ID1, ID2)] = np.log2(fc[ID1] + 0.1) - np.log2(fc[ID2] + 0.1)
            log2_fc[(ID1, ID2)] = pd.DataFrame(log2_fc[(ID1, ID2)], columns=['fold_2_change'], index = adata_dic[(cond1, datas)].var_names) 
    
    return log2_fc

def log2_FC_all(genes, adata, cond1, cond2, datas, save=None):
    '''the log2 fold_change of genes by comparision between each cond1 and cond2, eg., 'sick' vs 'Healthy', 'DWS' vs 'Healthy'
    Args:
        genes: list containing genes to calculate log2 fold_change
        adata, cond1, cond2, datas: parameters for function 'log2_FC'
            log_FC(adata, cond1, cond2, datas) returns the dict containing DataFrame obtained from function'log2_FC'
    Returns:
        list containing two tuple, 
        tuple 1: dict containing list with log2 fold_change of genes from log2_fc, keys are gene_names
        tuple 2: DataFrame containing the average of log2 fold_change by comparision beteen each sick and each healthy, index are gene_names
    '''
    
    from collections import defaultdict

    list_log2_fc = defaultdict(list)
    log2_fc_all = []

    log2_fc = log2_FC(adata, cond1, cond2, datas)
    
    for gene in genes:
        n_pos_fc = 0
        for key in log2_fc.keys():
            list_log2_fc[gene].append(log2_fc[key].loc[gene][0])
            list_log2_fc[gene].sort()
            if log2_fc[key].loc[gene][0] > 0:
                n_pos_fc += 1
            comp_frac = n_pos_fc/len(log2_fc.keys())
        log2_fc_all.append([np.mean(list_log2_fc[gene]), comp_frac])

    log2_fc_all = pd.DataFrame(log2_fc_all, columns=['fold_2_change', 'comp_frac'], index=pd.Index(genes)).sort_values('fold_2_change', ascending=False)
    
    df_log2_fc = pd.DataFrame([])
    for gene in genes:
        df_log2_fc[gene]=list_log2_fc[gene]
    df_log2_fc = df_log2_fc.set_index(pd.Index(log2_fc.keys()))
    df_log2_fc = df_log2_fc.T

    if save is not None:
        df_log2_fc.to_excel(save + 'log2_fc_by_genes.xls')
        log2_fc_all.to_excel(save + 'log2_fc_all.xls')

    return [list_log2_fc, log2_fc_all, df_log2_fc]

def top_inters(interactions,genes, adata, cond1, cond2, datas, top_n, save=None):
    '''the interactions with one gene with the top_n largest log2_FC by comparision between each cond1 and cond2, eg., 'sick' vs 'Healthy', 'DWS' vs 'Healthy'
    Args:
        interactions: list containing pairs of ligand and receptor genes for cell-cell communication, usually obtained from the 'interaction_unpacked_mouse.tsv'
        genes, adata, cond1, cond2, datas： parameters for function 'log2_FC_all'
            log2_FC_all(genes, adata, cond1, cond2, datas) returns list containing two tuple, 
                tuple 1: dict containing list with log2 fold_change of genes from log2_fc, keys are gene_names
                tuple 2: DataFrame containing the average of log2 fold_change by comparision beteen each sick and each healthy, index are gene_names
        top_n: int that the top number of genes to be analyzed
    Returns:
        DataFrame containing interactions with at least one gene in the top_n log2 fold_change list, and their average log2 fold_change
    '''
    
    from collections import defaultdict
    
    interactions_new = interactions.drop_duplicates()
    
    log2_fc_all = log2_FC_all(genes, adata, cond1, cond2, datas)[1]
    top_log2_fc = list(log2_fc_all[:top_n].index)

    top_inters = []
    for _,row in interactions_new.iterrows():
        gene_a = row['gene_name_a']
        gene_b = row['gene_name_b']
        if gene_a in top_log2_fc:
            top_inters.append([gene_a, gene_b, log2_fc_all.loc[gene_a][0], log2_fc_all.loc[gene_b][0]])
        elif gene_b in top_log2_fc:
            top_inters.append([gene_b, gene_a, log2_fc_all.loc[gene_b][0], log2_fc_all.loc[gene_a][0]])
    top_inter_list = pd.DataFrame(top_inters, columns = ['gene_name_a', 'gene_name_b', 'log2_fc_ga', 'log2_fc_gb']).drop_duplicates()
    top_inter_list['sum_log2_fc'] = top_inter_list['log2_fc_ga'] + top_inter_list['log2_fc_gb']
    top_inter_list = top_inter_list.sort_values('sum_log2_fc', ascending = False)
    
    if save is not None:
        top_inter_list.to_excel(save + cond1 + ' in ' + datas + '_top_inters.xls')
    
    return top_inter_list

def last_inters(interactions,genes, adata, cond1, cond2, datas, last_n, save=None):
    '''the interactions with one gene with the last_n smallest log2_FC by comparision between each cond1 and cond2, eg., 'sick' vs 'Healthy', 'DWS' vs 'Healthy'
    Args:
        interactions: list containing pairs of ligand and receptor genes for cell-cell communication, usually obtained from the 'interaction_unpacked_mouse.tsv'
        genes, adata, cond1, cond2, datas： parameters for function 'log2_FC_all'
            log2_FC_all(genes, adata, cond1, cond2, datas) returns list containing two tuple, 
                tuple 1: dict containing list with log2 fold_change of genes from log2_fc, keys are gene_names
                tuple 2: DataFrame containing the average of log2 fold_change by comparision beteen each sick and each healthy, index are gene_names
        last_n: int that the last number of genes to be analyzed
    Returns:
        DataFrame containing interactions with at least one gene in the last_n log2 fold_change list, and their average log2 fold_change
    '''
    
    from collections import defaultdict
    
    interactions_new = interactions.drop_duplicates()
    
    log2_fc_all = log2_FC_all(genes, adata, cond1, cond2, datas)[1]
    last_log2_fc = list(log2_fc_all[-last_n:].index)

    last_inters = []
    for _,row in interactions_new.iterrows():
        gene_a = row['gene_name_a']
        gene_b = row['gene_name_b']
        if gene_a in last_log2_fc:
            last_inters.append([gene_a, gene_b, log2_fc_all.loc[gene_a][0], log2_fc_all.loc[gene_b][0]])
        elif gene_b in last_log2_fc:
            last_inters.append([gene_b, gene_a, log2_fc_all.loc[gene_b][0], log2_fc_all.loc[gene_a][0]])
    last_inter_list = pd.DataFrame(last_inters, columns = ['gene_name_a', 'gene_name_b', 'log2_fc_ga', 'log2_fc_gb']).drop_duplicates()
    last_inter_list['sum_log2_fc'] = last_inter_list['log2_fc_ga'] + last_inter_list['log2_fc_gb']
    last_inter_list = last_inter_list.sort_values('sum_log2_fc', ascending = False)
    
    if save is not None:
        last_inter_list.to_csv(save + cond1 + ' in ' + datas + '_last_inters.tsv')
    
    return last_inter_list

def top_inters_gbs(interactions,genes, adata, cond1, cond2, datas, top_n, groupby, save=None):
    '''the interactions with one gene with the top_n largest log2_FC in each groupby by comparision between each cond1 and cond2, eg., 'sick' vs 'Healthy', 'DWS' vs 'Healthy' 
    Args:
        interactions: list containing pairs of ligand and receptor genes for cell-cell communication, usually obtained from the 'interaction_unpacked_mouse.tsv'
        genes, adata, cond1, cond2, datas： parameters for function 'log2_FC_all'
            log2_FC_all(genes, adata, cond1, cond2, datas) returns list containing two tuple, 
                tuple 1: dict containing list with log2 fold_change of genes from log2_fc, keys are gene_names
                tuple 2: DataFrame containing the average of log2 fold_change by comparision beteen each sick and each healthy, index are gene_names
        top_n: int that the top number of genes to be analyzed
        groupby: name of group used to analyze the data, eg., 'cell_type', 'cell_subtype' 
    Returns:
        DataFrame containing interactions with the gene & its log2_FC in the top_n log2 fold_change list, and the log2_FC of the other gene in each groupby
        List containing Datafram interactions with the gene & its log2_FC in the top_n log2 fold_change list, and the log2_FC of the other gene, with groupby as key
    '''
    
    from collections import defaultdict
    
    if 'doublets' in list(adata.obs[groupby].astype('category').cat.categories):
        gbs = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
    else:
        gbs = list(adata.obs[groupby].astype('category').cat.categories)

    adata_gb = {}
    log2_fc_all = {}
    top_log2_fc = {}
    inters = {}
    for gb in gbs:
        adata_gb[gb] = adata[adata.obs[groupby] == gb]
        log2_fc_all[gb] = log2_FC_all(genes, adata_gb[gb], cond1, cond2, datas)[1]
        top_log2_fc[gb] = list(log2_fc_all[gb][:top_n].index)
        inters[gb] = top_inters(interactions, genes, adata_gb[gb], cond1, cond2, datas, top_n)
        
    top_its = defaultdict(list)
    for gb1 in gbs:
        for _,row in inters[gb1].iterrows():
            gene_a = row['gene_name_a']
            gene_b = row['gene_name_b']
            
            top_inter = [gene_a, gb1, log2_fc_all[gb1].loc[gene_a][0], gene_b]
            global cls
            cls = ['gene_name_a',
                   'ga_ct',
                   'log2_fc_ga',
                   'gene_name_b']
            for gb2 in gbs:
                top_inter.extend([log2_fc_all[gb2].loc[gene_b][0]])
                cls.extend(['log2_fc_gb_' + gb2])
            top_its[gb1].append(top_inter)

    top_inter_list = pd.DataFrame([])
    for gb in gbs:
        top_its[gb] = pd.DataFrame(top_its[gb], columns=cls, index=pd.Index([gb]*len(top_its[gb]))).sort_values('log2_fc_ga', ascending=False)
        top_inter_list = pd.concat([top_inter_list, top_its[gb]])

    if save is not None:
        top_inter_list.to_csv(save + cond1 + '_in_' + datas + '_by_' + groupby + '_top_inters.tsv')

    return [top_inter_list, top_its]

def last_inters_gbs(interactions,genes, adata, cond1, cond2, datas, last_n, groupby, save=None):
    '''the interactions with one gene with the last_n largest log2_FC in each groupby by comparision between each cond1 and cond2, eg., 'sick' vs 'Healthy', 'DWS' vs 'Healthy' 
    Args:
        interactions: list containing pairs of ligand and receptor genes for cell-cell communication, usually obtained from the 'interaction_unpacked_mouse.tsv'
        genes, adata, cond1, cond2, datas： parameters for function 'log2_FC_all'
            log2_FC_all(genes, adata, cond1, cond2, datas) returns list containing two tuple, 
                tuple 1: dict containing list with log2 fold_change of genes from log2_fc, keys are gene_names
                tuple 2: DataFrame containing the average of log2 fold_change by comparision beteen each sick and each healthy, index are gene_names
        last_n: int that the last number of genes to be analyzed
        groupby: name of group used to analyze the data, eg., 'cell_type', 'cell_subtype' 
    Returns:
        DataFrame containing interactions with the gene & its log2_FC in the last_n log2 fold_change list, and the log2_FC of the other gene in each groupby
        List containing Datafram interactions with the gene & its log2_FC in the last_n log2 fold_change list, and the log2_FC of the other gene, with groupby as key
    '''
    
    from collections import defaultdict
    
    if 'doublets' in list(adata.obs[groupby].astype('category').cat.categories):
        gbs = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
    else:
        gbs = list(adata.obs[groupby].astype('category').cat.categories)

    adata_gb = {}
    log2_fc_all = {}
    last_log2_fc = {}
    inters = {}
    for gb in gbs:
        adata_gb[gb] = adata[adata.obs[groupby] == gb]
        log2_fc_all[gb] = log2_FC_all(genes, adata_gb[gb], cond1, cond2, datas)[1]
        last_log2_fc[gb] = list(log2_fc_all[gb][-last_n:].index)
        inters[gb] = last_inters(interactions, genes, adata_gb[gb], cond1, cond2, datas, last_n)
        
    last_its = defaultdict(list)
    for gb1 in gbs:
        for _,row in inters[gb1].iterrows():
            gene_a = row['gene_name_a']
            gene_b = row['gene_name_b']
            
            last_inter = [gene_a, gb1, log2_fc_all[gb1].loc[gene_a][0], gene_b]
            global cls
            cls = ['gene_name_a',
                   'ga_ct',
                   'log2_fc_ga',
                   'gene_name_b']
            for gb2 in gbs:
                last_inter.extend([log2_fc_all[gb2].loc[gene_b][0]])
                cls.extend(['log2_fc_gb_' + gb2])
            last_its[gb1].append(last_inter)

    last_inter_list = pd.DataFrame([])
    for gb in gbs:
        last_its[gb] = pd.DataFrame(last_its[gb], columns=cls, index=pd.Index([gb]*len(last_its[gb]))).sort_values('log2_fc_ga', ascending=False)
        last_inter_list = pd.concat([last_inter_list, last_its[gb]])

    if save is not None:
        last_inter_list.to_csv(save + cond1 + '_in_' + datas + '_by_' + groupby + '_last_inters.tsv')

    return [last_inter_list, last_its]

def cul_plot_log2_FC(genes, adata, cond1, cond2, datas, top_n, save=None):
    '''cumulative plots of top_n genes with largest log2 fold_change by comparision between cond1 and cond2 in datas
    Args:
        genes, adata, cond1, cond2, datas： parameters for function 'log2_FC_all'
            log2_FC_all(genes, adata, cond1, cond2, datas) returns list containing two tuple, 
                tuple 1: dict containing list with log2 fold_change of genes from log2_fc, keys are gene_names
                tuple 2: DataFrame containing the average of log2 fold_change by comparision beteen each sick and each healthy, index are gene_names
        top_n: int that the top number of genes to be plotted
    Returns: cumulative plots of top_n genes shows 'Fraction of comparison between sick and healthy with log2FC > x' vs 'log2_FC'
    '''
    from collections import defaultdict

    conditions = list(adata.obs['Condition'].astype('category').cat.categories)
    datasets = list(adata.obs['dataset'].astype('category').cat.categories)
    sicks = list(adata.obs['sick'].astype('category').cat.categories)

    n_ID = {}
    adata_dic = {}
    for condition in conditions:
        for dataset in datasets:
            adata_dic[(condition, dataset)] = adata[adata.obs['Condition'] == condition][adata[adata.obs['Condition'] == condition].obs['dataset'] == dataset]
            n_ID[(condition, dataset)] = len(adata_dic[(condition, dataset)].obs['ID'].astype('category').cat.categories)

    for sick in sicks:
        for dataset in datasets:
            adata_dic[(sick, dataset)] = adata[adata.obs['sick'] == sick][adata[adata.obs['sick'] == sick].obs['dataset'] == dataset]
            n_ID[(sick, dataset)] = len(adata_dic[(sick, dataset)].obs['ID'].astype('category').cat.categories)
    
    list_log2_fc = log2_FC_all(genes, adata, cond1, cond2, datas)[0]
    log2_fc_all = log2_FC_all(genes, adata, cond1, cond2, datas)[1]
    top_log2_fc = list(log2_fc_all[:top_n].index)
    
    n = n_ID[(cond1, datas)] * n_ID[(cond2, datas)]
    y = 1.0 - np.linspace(0, 1, n)
    for gene in top_log2_fc:
        fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=300, facecolor='white')
        ax.plot(list_log2_fc[gene], y, lw=3, color='red', zorder=6)
        ax.axvline(0, lw=2, color='grey', ls='--')
        ax.axvline(np.mean(list_log2_fc[gene]), lw=2, color='gold')
        ax.set_title(gene, fontsize=20)
        ax.grid(True)
        ax.set_xlim(-10, 15)
        ax.set_ylim(-0.01, 1.01)
        ax.set_ylabel('Fraction of comparison between\n'  + cond1 +' ' + 'and ' + cond2 + ' with log2FC > x', fontsize=15)
        ax.set_xlabel('log2FC', fontsize=15)
        if save is not None:
            plt.savefig(save + gene +'.png')

def cul_plot_log2_FC_gb(genes, adata, cond1, cond2, datas, top_n, groupby, save=None):
    '''cumulative plots of top_n genes with largest log2 fold_change by comparision between cond1 and cond2 in datas
    Args:
        genes, adata, cond1, cond2, datas： parameters for function 'log2_FC_all'
            log2_FC_all(genes, adata, cond1, cond2, datas) returns list containing two tuple, 
                tuple 1: dict containing list with log2 fold_change of genes from log2_fc, keys are gene_names
                tuple 2: DataFrame containing the average of log2 fold_change by comparision beteen each sick and each healthy, index are gene_names
        top_n: int that the top number of genes to be plotted
        groupby: name of group used to analyze the data, eg., 'cell_type', 'cell_subtype'
    Returns: set groupby as 'cell_type' as the example, 
        returns cumulative plots of top_n genes in each cell_type shows 'Fraction of comparison between sick and healthy with log2FC > x' vs 'log2_FC', 
        each figure only includes culmulative plots of the gene in the cell_type
    '''
    from collections import defaultdict
    
    if 'doublets' in list(adata.obs[groupby].astype('category').cat.categories):
        gbs = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
    else:
        gbs = list(adata.obs[groupby].astype('category').cat.categories)

    n_gbs = len(gbs)
    
    adata_gb = {}
    list_log2_fc = {}
    log2_fc_all = {}
    top_log2_fc = {}
    n_ID = {}
    adata_dic = {}
    
    for gb in gbs:
        adata_gb[gb] = adata[adata.obs[groupby] == gb]
        list_log2_fc[gb] = log2_FC_all(genes, adata_gb[gb], cond1, cond2, datas)[0]
        log2_fc_all[gb] = log2_FC_all(genes, adata_gb[gb], cond1, cond2, datas)[1]
        top_log2_fc[gb] = list(log2_fc_all[gb][:top_n].index)
        
        conditions = list(adata_gb[gb].obs['Condition'].astype('category').cat.categories)
        datasets = list(adata_gb[gb].obs['dataset'].astype('category').cat.categories)
        sicks = list(adata_gb[gb].obs['sick'].astype('category').cat.categories)

        for condition in conditions:
            adata_f = adata_gb[gb][adata_gb[gb].obs['Condition'] == condition]
            for dataset in datasets:
                adata_dic[(gb, condition, dataset)] = adata_f[adata_f.obs['dataset'] == dataset]
                n_ID[(gb, condition, dataset)] = len(adata_dic[(gb, condition, dataset)].obs['ID'].astype('category').cat.categories)

        for sick in sicks:
            adata_f = adata_gb[gb][adata_gb[gb].obs['sick'] == sick]
            for dataset in datasets:
                adata_dic[(gb, sick, dataset)] = adata_f[adata_f.obs['dataset'] == dataset]
                n_ID[(gb, sick, dataset)] = len(adata_dic[(gb, sick, dataset)].obs['ID'].astype('category').cat.categories)
    
    for gb in gbs:
        for gene in top_log2_fc[gb]:
            fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=300, facecolor='white')
            n = n_ID[(gb, cond1, datas)] * n_ID[(gb, cond2, datas)]
            y = 1.0 - np.linspace(0, 1, n)
            ax.plot(list_log2_fc[gb][gene], y, lw=3, zorder=6)
            ax.axvline(0, lw=2, color='grey', ls='--')
            ax.axvline(np.mean(list_log2_fc[gb][gene]), lw=2, color='gold')
            ax.set_title(gene + ' in ' + gb, fontsize=20)
            ax.grid(True)
            ax.set_xlim(-10, 15)
            ax.set_ylim(-0.01, 1.01)
            ax.set_ylabel('Fraction of comparison between\n'  + cond1 +' ' + 'and ' + cond2 + ' with log2FC > x', fontsize=15)
            ax.set_xlabel('log2FC', fontsize=15) 

            if save is not None:
                plt.savefig(save + gene + ' in ' + gb +'.png')

def cul_plot_log2_FC_gbs(genes, adata, cond1, cond2, datas, top_n, groupby, save=None):
    '''cumulative plots of top_n genes with largest log2 fold_change by comparision between cond1 and cond2 in datas
    Args:
        genes, adata, cond1, cond2, datas： parameters for function 'log2_FC_all'
            log2_FC_all(genes, adata, cond1, cond2, datas) returns list containing two tuple, 
                tuple 1: dict containing list with log2 fold_change of genes from log2_fc, keys are gene_names
                tuple 2: DataFrame containing the average of log2 fold_change by comparision beteen each sick and each healthy, index are gene_names
        top_n: int that the top number of genes to be plotted
        groupby: name of group used to analyze the data, eg., 'cell_type', 'cell_subtype'
    Returns: set groupby as 'cell_type' as the example, 
        returns cumulative plots of top_n genes in each cell_type shows 'Fraction of comparison between sick and healthy with log2FC > x' vs 'log2_FC', 
        each figure includes culmulative plots of the gene in all cell_types
    '''
    from collections import defaultdict
    
    if 'doublets' in list(adata.obs[groupby].astype('category').cat.categories):
        gbs = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
    else:
        gbs = list(adata.obs[groupby].astype('category').cat.categories)

    n_gbs = len(gbs)
    colors = sns.color_palette('hls', n_gbs)
    
    adata_gb = {}
    list_log2_fc = {}
    log2_fc_all = {}
    top_log2_fc = {}
    n_ID = {}
    adata_dic = {}
    
    for gb in gbs:
        adata_gb[gb] = adata[adata.obs[groupby] == gb]
        list_log2_fc[gb] = log2_FC_all(genes, adata_gb[gb], cond1, cond2, datas)[0]
        log2_fc_all[gb] = log2_FC_all(genes, adata_gb[gb], cond1, cond2, datas)[1]
        top_log2_fc[gb] = list(log2_fc_all[gb][:top_n].index)
        
        conditions = list(adata_gb[gb].obs['Condition'].astype('category').cat.categories)
        datasets = list(adata_gb[gb].obs['dataset'].astype('category').cat.categories)
        sicks = list(adata_gb[gb].obs['sick'].astype('category').cat.categories)

        for condition in conditions:
            adata_f = adata_gb[gb][adata_gb[gb].obs['Condition'] == condition]
            for dataset in datasets:
                adata_dic[(gb, condition, dataset)] = adata_f[adata_f.obs['dataset'] == dataset]
                n_ID[(gb, condition, dataset)] = len(adata_dic[(gb, condition, dataset)].obs['ID'].astype('category').cat.categories)

        for sick in sicks:
            adata_f = adata_gb[gb][adata_gb[gb].obs['sick'] == sick]
            for dataset in datasets:
                adata_dic[(gb, sick, dataset)] = adata_f[adata_f.obs['dataset'] == dataset]
                n_ID[(gb, sick, dataset)] = len(adata_dic[(gb, sick, dataset)].obs['ID'].astype('category').cat.categories)
    
    for gb1 in gbs:
        for gene in top_log2_fc[gb1]:
            fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=300, facecolor='white')
            for gb2, cl in zip(gbs, colors):
                n = n_ID[(gb2, cond1, datas)] * n_ID[(gb2, cond2, datas)]
                y = 1.0 - np.linspace(0, 1, n)
                ax.plot(list_log2_fc[gb2][gene], y, lw=3, zorder=6, label=gb2, color=cl)
                ax.axvline(0, lw=2, color='grey', ls='--')
                ax.axvline(np.mean(list_log2_fc[gb2][gene]), lw=2, color=cl, ls='--')
            ax.legend(loc='upper right')
            ax.set_title(gene + ' in ' + gb1, fontsize=20)
            ax.grid(True)
            ax.set_xlim(-10, 15)
            ax.set_ylim(-0.01, 1.01)
            ax.set_ylabel('Fraction of comparison between\n'  + cond1 +' ' + 'and ' + cond2 + ' with log2FC > x', fontsize=15)
            ax.set_xlabel('log2FC', fontsize=15) 

            if save is not None:
                plt.savefig(save + gene + ' in ' + gb1 +'.png')

def cul_plot_log2_FC_inters(interactions, genes, adata, cond1, cond2, datas, top_n, groupby, save=None):
    '''cumulative plots of top_n genes with largest log2 fold_change by comparision between cond1 and cond2 in datas
    Args:
        genes, adata, cond1, cond2, datas： parameters for function 'log2_FC_all'
            log2_FC_all(genes, adata, cond1, cond2, datas) returns list containing two tuple, 
                tuple 1: dict containing list with log2 fold_change of genes from log2_fc, keys are gene_names
                tuple 2: DataFrame containing the average of log2 fold_change by comparision beteen each sick and each healthy, index are gene_names
        top_n: int that the top number of genes to be plotted
        groupby: name of group used to analyze the data, eg., 'cell_type', 'cell_subtype'
    Returns: set groupby as 'cell_type' as the example, 
        returns cumulative plots of top_n genes in each cell_type shows 'Fraction of comparison between sick\nand healthy with log2FC > x' vs 'log2_FC', 
        each figure includes culmulative plots of the gene in all cell_types
    '''
    from collections import defaultdict
    
    if 'doublets' in list(adata.obs[groupby].astype('category').cat.categories):
        gbs = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
    else:
        gbs = list(adata.obs[groupby].astype('category').cat.categories)

    n_gbs = len(gbs)
    colors = sns.color_palette('hls', n_gbs)
    
    adata_gb = {}
    list_log2_fc = {}
    n_ID = {}
    adata_dic = {}
    
    for gb in gbs:
        adata_gb[gb] = adata[adata.obs[groupby] == gb]
        list_log2_fc[gb] = log2_FC_all(genes, adata_gb[gb], cond1, cond2, datas)[0]
        
        conditions = list(adata_gb[gb].obs['Condition'].astype('category').cat.categories)
        datasets = list(adata_gb[gb].obs['dataset'].astype('category').cat.categories)
        sicks = list(adata_gb[gb].obs['sick'].astype('category').cat.categories)

        for condition in conditions:
            adata_f = adata_gb[gb][adata_gb[gb].obs['Condition'] == condition]
            for dataset in datasets:
                adata_dic[(gb, condition, dataset)] = adata_f[adata_f.obs['dataset'] == dataset]
                n_ID[(gb, condition, dataset)] = len(adata_dic[(gb, condition, dataset)].obs['ID'].astype('category').cat.categories)

        for sick in sicks:
            adata_f = adata_gb[gb][adata_gb[gb].obs['sick'] == sick]
            for dataset in datasets:
                adata_dic[(gb, sick, dataset)] = adata_f[adata_f.obs['dataset'] == dataset]
                n_ID[(gb, sick, dataset)] = len(adata_dic[(gb, sick, dataset)].obs['ID'].astype('category').cat.categories)
    
    
    top_inters = top_inters_gbs(interactions,genes, adata, cond1, cond2, datas, top_n, groupby)[1]
    
    for gb1 in gbs:
        for _,row in top_inters[gb1].iterrows():
            fig, axs = plt.subplots(1, 2, figsize=(15, 5), dpi=300, facecolor='white', sharex=True, sharey=True)
            gene_a = row['gene_name_a']
            gene_b = row['gene_name_b']

            n1 = n_ID[(gb1, cond1, datas)] * n_ID[(gb1, cond2, datas)]
            y1 = 1.0 - np.linspace(0, 1, n1)

            axs[0].plot(list_log2_fc[gb1][gene_a], y1, lw=3, zorder=6, label=gb1)
            axs[0].axvline(np.mean(list_log2_fc[gb1][gene_a]), lw=2, ls='--', color='gold')
            axs[0].set_title(gene_a + ' in ' + gb1, fontsize=20)

            for gb2, cl2 in zip(gbs, colors):
                n2 = n_ID[(gb2, cond1, datas)] * n_ID[(gb2, cond2, datas)]
                y2 = 1.0 - np.linspace(0, 1, n2)
                axs[1].plot(list_log2_fc[gb2][gene_b], y2, lw=3, zorder=6, label=gb2, color=cl2)
                axs[1].axvline(np.mean(list_log2_fc[gb2][gene_b]), lw=2, color=cl2, ls='--')
                axs[1].set_title(gene_b, fontsize=20)

            for i in [0, 1]:
                axs[i].axvline(0, lw=2, color='grey') 
                axs[i].grid(True)
                axs[i].set_xlim(-10, 15)
                axs[i].set_ylim(-0.01, 1.01)
                axs[i].legend(loc='upper right')

            fig.text(0.5, 0, 'log2FC', ha='center', fontsize=15)
            fig.text(0, 0.5,'Fraction of comparison between\n'  + cond1 +' ' + 'and ' + cond2 + ' with log2FC > x', va='center', rotation='vertical', fontsize=15)

            fig.tight_layout(rect=[0.04,0.04,1.02,1], pad=0, h_pad=0, w_pad=0)

            if save is not None:
                plt.savefig(save + gene_a + '_in_' + gb1 + '_&_' + gene_b + '.png', bbox_inches='tight')

def dotplot_log2_FC(adata, cond1, cond2, datas, top_genes, inters, genes_increase, genes_decrease, save=None):
    sc.pp.log1p(adata)

    from collections import defaultdict

    conditions = list(adata.obs['Condition'].astype('category').cat.categories)
    datasets = list(adata.obs['dataset'].astype('category').cat.categories)
    sicks = list(adata.obs['sick'].astype('category').cat.categories)

    n_ID = {}
    adata_dic = {}
    for condition in conditions:
        for dataset in datasets:
            adata_dic[(condition, dataset)] = adata[adata.obs['Condition'] == condition][adata[adata.obs['Condition'] == condition].obs['dataset'] == dataset]
            n_ID[(condition, dataset)] = len(adata_dic[(condition, dataset)].obs['ID'].astype('category').cat.categories)

    for sick in sicks:
        for dataset in datasets:
            adata_dic[(sick, dataset)] = adata[adata.obs['sick'] == sick][adata[adata.obs['sick'] == sick].obs['dataset'] == dataset]
            n_ID[(sick, dataset)] = len(adata_dic[(sick, dataset)].obs['ID'].astype('category').cat.categories)
    
    fig,axs = plt.subplots(2, 1, figsize=(25,8), gridspec_kw={'wspace':0.3}, dpi=300, facecolor='white')
    mainplot = {}
    ax=axs[0]
    mainplot[0] = sc.pl.dotplot(adata_dic[(cond1, datas)], inters, groupby='cell_type', ax=ax, show=False, cmap='RdBu_r')
    ax.set_title(cond1 + ' ' + datas, fontsize=12, y=0.7, loc='right')

    ax=axs[1]
    mainplot[1] = sc.pl.dotplot(adata_dic[(cond2, datas)], inters, groupby='cell_type', ax=ax, show=False, cmap='RdBu_r')
    ax.set_title(cond2 + ' ' + datas, fontsize=12, y=0.7, loc='right')

    for i in [0,1]: 
        for xlabel in mainplot[i]['mainplot_ax'].get_xticklabels():
            gene = xlabel.get_text()
            if gene in top_genes:
                mainplot[i]['mainplot_ax'].axvline(xlabel.get_position()[0], color='grey', alpha=0.3, zorder = 0)

            for ylabel in mainplot[i]['mainplot_ax'].get_yticklabels():
                cell_type = ylabel.get_text()
                if [gene, cell_type] in genes_increase:
                    x = [xlabel.get_position()[0]]
                    y = [ylabel.get_position()[1]]
                    mainplot[i]['mainplot_ax'].scatter(x, y, alpha=0.3, zorder=0, s=200, marker='s', c='', edgecolors='r', linewidths=1.5)
                elif [gene, cell_type] in genes_decrease:
                    x = [xlabel.get_position()[0]]
                    y = [ylabel.get_position()[1]]
                    mainplot[i]['mainplot_ax'].scatter(x, y, alpha=0.3, zorder=0, s=200, marker='s', c='', edgecolors='b', linewidths=1.5)
    if save is not None:
        plt.savefig(save)

def cul_plot_genes(genes, adata, cond1, cond2, datas, groupby, save=None):

    from collections import defaultdict
    
    if 'doublets' in list(adata.obs[groupby].astype('category').cat.categories):
        gbs = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
    else:
        gbs = list(adata.obs[groupby].astype('category').cat.categories)

    n_gbs = len(gbs)
    colors = sns.color_palette('hls', n_gbs)

    adata_gb = {}
    list_log2_fc = {}
    log2_fc_all = {}
    top_log2_fc = {}
    n_ID = {}
    adata_dic = {}
    
    for gb in gbs:
        adata_gb[gb] = adata[adata.obs[groupby] == gb]
        list_log2_fc[gb] = log2_FC_all(genes, adata_gb[gb], cond1, cond2, datas)[0]
        
        conditions = list(adata_gb[gb].obs['Condition'].astype('category').cat.categories)
        datasets = list(adata_gb[gb].obs['dataset'].astype('category').cat.categories)
        sicks = list(adata_gb[gb].obs['sick'].astype('category').cat.categories)

        for condition in conditions:
            adata_f = adata_gb[gb][adata_gb[gb].obs['Condition'] == condition]
            for dataset in datasets:
                adata_dic[(gb, condition, dataset)] = adata_f[adata_f.obs['dataset'] == dataset]
                n_ID[(gb, condition, dataset)] = len(adata_dic[(gb, condition, dataset)].obs['ID'].astype('category').cat.categories)

        for sick in sicks:
            adata_f = adata_gb[gb][adata_gb[gb].obs['sick'] == sick]
            for dataset in datasets:
                adata_dic[(gb, sick, dataset)] = adata_f[adata_f.obs['dataset'] == dataset]
                n_ID[(gb, sick, dataset)] = len(adata_dic[(gb, sick, dataset)].obs['ID'].astype('category').cat.categories)
    
    
    for gene in genes:
        fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=300)
        for gb, cl in zip(gbs, colors):
            n = n_ID[(gb, cond1, datas)] * n_ID[(gb, cond2, datas)]
            y = 1.0 - np.linspace(0, 1, n)
            ax.plot(list_log2_fc[gb][gene], y, lw=3, zorder=6, color=cl, label=gb)
            ax.legend(loc='upper right', fontsize=15)
            ax.axvline(0, lw=2, color='grey', ls='--')
            ax.axvline(np.mean(list_log2_fc[gb][gene]), lw=2, color=cl)
            ax.set_title(gene, fontsize=25)
            ax.grid(True)
            ax.set_xlim(-10, 15)
            ax.set_ylim(-0.01, 1.01)
            ax.set_xticklabels([-10, -5, 0, 5, 10, 15], fontsize=15)
            ax.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=15)
            ax.set_ylabel('Fraction of comparison between\n'  + cond1 +' ' + 'and ' + cond2 + ' with log2FC > x', fontsize=15)
            ax.set_xlabel('log2FC', fontsize=25) 

            if save is not None:
                plt.savefig(save + gene +'.png')


def cul_plot_inters(interactions, adata, cond1, cond2, datas, groupby, save=None):

    from collections import defaultdict
    
    if 'doublets' in list(adata.obs[groupby].astype('category').cat.categories):
        gbs = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
    else:
        gbs = list(adata.obs[groupby].astype('category').cat.categories)

    n_gbs = len(gbs)
    colors = sns.color_palette('hls', n_gbs)
    
    adata_gb = {}
    list_log2_fc = {}
    n_ID = {}
    adata_dic = {}
    genes = np.unique(interactions)

    for gb in gbs:
        adata_gb[gb] = adata[adata.obs[groupby] == gb]
        list_log2_fc[gb] = log2_FC_all(genes, adata_gb[gb], cond1, cond2, datas)[0]
        
        conditions = list(adata_gb[gb].obs['Condition'].astype('category').cat.categories)
        datasets = list(adata_gb[gb].obs['dataset'].astype('category').cat.categories)
        sicks = list(adata_gb[gb].obs['sick'].astype('category').cat.categories)

        for condition in conditions:
            adata_f = adata_gb[gb][adata_gb[gb].obs['Condition'] == condition]
            for dataset in datasets:
                adata_dic[(gb, condition, dataset)] = adata_f[adata_f.obs['dataset'] == dataset]
                n_ID[(gb, condition, dataset)] = len(adata_dic[(gb, condition, dataset)].obs['ID'].astype('category').cat.categories)

        for sick in sicks:
            adata_f = adata_gb[gb][adata_gb[gb].obs['sick'] == sick]
            for dataset in datasets:
                adata_dic[(gb, sick, dataset)] = adata_f[adata_f.obs['dataset'] == dataset]
                n_ID[(gb, sick, dataset)] = len(adata_dic[(gb, sick, dataset)].obs['ID'].astype('category').cat.categories)
    
    for inter in interactions:
        fig, axs = plt.subplots(1, 2, figsize=(15, 5), dpi=300, facecolor='white', sharex=True, sharey=True)
        for gb, cl in zip(gbs, colors):
            n = n_ID[(gb, cond1, datas)] * n_ID[(gb, cond2, datas)]
            y = 1.0 - np.linspace(0, 1, n)
            
            for i in [0, 1]:
                axs[i].plot(list_log2_fc[gb][inter[i]], y, lw=3, zorder=6, label=gb, color=cl)
                axs[i].axvline(np.mean(list_log2_fc[gb][inter[i]]), lw=2, color=cl, ls='--')
                axs[i].set_title(inter[i], fontsize=20)
                
                axs[i].axvline(0, lw=2, color='grey') 
                axs[i].grid(True)
                axs[i].set_xlim(-10, 15)
                axs[i].set_ylim(-0.01, 1.01)
                axs[i].legend(loc='upper right')

        fig.text(0.5, 0, 'log2FC', ha='center', fontsize=15)
        fig.text(0, 0.5,'Fraction of comparison between\n'  + cond1 +' ' + 'and ' + cond2 + ' with log2FC > x', va='center', rotation='vertical', fontsize=15)

        fig.tight_layout(rect=[0.04,0.04,1.02,1], pad=0, h_pad=0, w_pad=0)

        if save is not None:
            plt.savefig(save + inter[0] + '_&_' + inter[1] + '.png', bbox_inches='tight')

def cul_plot_gene(genes, adata, cond1, cond2, datas, cell_type, save=None):

    from collections import defaultdict
    
    if 'doublets' in list(adata.obs['cell_type'].astype('category').cat.categories):
        gbs = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
    else:
        gbs = list(adata.obs['cell_type'].astype('category').cat.categories)
    
    adata_gb = {}
    list_log2_fc = {}
    n_ID = {}
    adata_dic = {}
    
    for gb in gbs:
        adata_gb[gb] = adata[adata.obs['cell_type'] == gb]
        list_log2_fc[gb] = log2_FC_all(genes, adata_gb[gb], cond1, cond2, datas)[0]
        
        conditions = list(adata_gb[gb].obs['Condition'].astype('category').cat.categories)
        datasets = list(adata_gb[gb].obs['dataset'].astype('category').cat.categories)
        sicks = list(adata_gb[gb].obs['sick'].astype('category').cat.categories)

        for condition in conditions:
            adata_f = adata_gb[gb][adata_gb[gb].obs['Condition'] == condition]
            for dataset in datasets:
                adata_dic[(gb, condition, dataset)] = adata_f[adata_f.obs['dataset'] == dataset]
                n_ID[(gb, condition, dataset)] = len(adata_dic[(gb, condition, dataset)].obs['ID'].astype('category').cat.categories)

        for sick in sicks:
            adata_f = adata_gb[gb][adata_gb[gb].obs['sick'] == sick]
            for dataset in datasets:
                adata_dic[(gb, sick, dataset)] = adata_f[adata_f.obs['dataset'] == dataset]
                n_ID[(gb, sick, dataset)] = len(adata_dic[(gb, sick, dataset)].obs['ID'].astype('category').cat.categories)
    
    n = n_ID[(cell_type, cond1, datas)] * n_ID[(cell_type, cond2, datas)]
    y = 1.0 - np.linspace(0, 1, n)
    for gene in genes:
        fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=300, facecolor='white')
        ax.plot(list_log2_fc[cell_type][gene], y, lw=3, zorder=6)
        ax.axvline(0, lw=2, color='grey')
        ax.axvline(np.mean(list_log2_fc[cell_type][gene]), lw=2, color='gold', ls='--')
        ax.set_title(gene + ' in ' + cell_type, fontsize=20)
        ax.grid(True)
        ax.set_xlim(-10, 15)
        ax.set_ylim(-0.01, 1.01)
        ax.set_ylabel('Fraction of comparison between\n'  + cond1 +' ' + 'and ' + cond2 + ' with log2FC > x', fontsize=15)

        ax.set_xlabel('log2FC', fontsize=15)
        if save is not None:
            plt.savefig(save + gene +'.png')

def cul_plot_gene_cts(genes, adata, cond1, cond2, datas, cell_types, save=None):

    from collections import defaultdict
    
    if 'doublets' in list(adata.obs['cell_type'].astype('category').cat.categories):
        gbs = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
    else:
        gbs = list(adata.obs['cell_type'].astype('category').cat.categories)
    
    adata_gb = {}
    list_log2_fc = {}
    n_ID = {}
    adata_dic = {}
    
    for gb in gbs:
        adata_gb[gb] = adata[adata.obs['cell_type'] == gb]
        list_log2_fc[gb] = log2_FC_all(genes, adata_gb[gb], cond1, cond2, datas)[0]
        
        conditions = list(adata_gb[gb].obs['Condition'].astype('category').cat.categories)
        datasets = list(adata_gb[gb].obs['dataset'].astype('category').cat.categories)
        sicks = list(adata_gb[gb].obs['sick'].astype('category').cat.categories)

        for condition in conditions:
            adata_f = adata_gb[gb][adata_gb[gb].obs['Condition'] == condition]
            for dataset in datasets:
                adata_dic[(gb, condition, dataset)] = adata_f[adata_f.obs['dataset'] == dataset]
                n_ID[(gb, condition, dataset)] = len(adata_dic[(gb, condition, dataset)].obs['ID'].astype('category').cat.categories)

        for sick in sicks:
            adata_f = adata_gb[gb][adata_gb[gb].obs['sick'] == sick]
            for dataset in datasets:
                adata_dic[(gb, sick, dataset)] = adata_f[adata_f.obs['dataset'] == dataset]
                n_ID[(gb, sick, dataset)] = len(adata_dic[(gb, sick, dataset)].obs['ID'].astype('category').cat.categories)
    
    colors = sns.color_palette('plasma', len(cell_types))
    for gene in genes:
        fig, ax = plt.subplots(1, 1, figsize=(3.6, 2.4), dpi=300)
        for cell_type, cl in zip(cell_types, colors):
            n = n_ID[(cell_type, cond1, datas)] * n_ID[(cell_type, cond2, datas)]
            y = 1.0 - np.linspace(0, 1, n)
            
            if cell_type == 'B_cells':
                ct = 'B cells'
            elif cell_type == 'T_cells':
                ct = 'T cells'
            elif cell_type == 'NK_cells':
                ct = 'NK cells'
            else:
                ct = cell_type
                
            ax.plot(list_log2_fc[cell_type][gene], y, lw=3, zorder=6, label=ct, color=cl)
            ax.axvline(np.mean(list_log2_fc[cell_type][gene]), lw=2, color=cl, ls='--')
        ax.legend(loc='upper right', fontsize=8)
        ax.axvline(0, lw=2, color='grey')
        ax.set_title(gene)
        ax.grid(True)
        ax.set_xlim(-10, 15)
        ax.set_ylim(-0.01, 1.01)
        ax.set_ylabel('Fraction of comparison between\n'  + cond1 +' ' + 'and ' + cond2 + ' with log2FC > x', fontsize=10)

        ax.set_xlabel('log2FC', fontsize=10)
        if save is not None:
            plt.savefig(save + gene +'.png')
            
        