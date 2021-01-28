import pandas as pd 
import numpy as np 
# import adata_utilis as au  # au.split

def interaction_analysis_same(datas1, datas2, cond1, cond2, ct, res):
    '''analysis of the interactions between ct and other cell_types,
    Args:
        datas1, cond1 --> requiremets for interaction_a
        datas2, cond2 --> requiremets for interaction_b
        ct --> analysis of ct and other cell_types between interaction_a and interaction_b
        res --> the dict of dataframe containing columns we want to analyze, eg. 
                res = average_fraction(adatag, ['dataset', 'Condition', 'cell_type'], 
                                       axis='obs', greater_than=0, log=False)
    Returns:
        dict containing DataFrame with interactions expressed both in interaction_a and interaction_b
    '''
    
    from collections import defaultdict
    
    cell_types = ['B_cells','NK_cells' , 'T_cells', 'Monocytes', 'Plasmablasts', 'pDCs', 'cDCs']
    same_inters = defaultdict(list)
    
    inter_list = res.keys()
    for cell_type in cell_types:
        if (datas1, cond1, cell_type, ct) in inter_list:
            interaction_a = res[(datas1, cond1, cell_type, ct)]
            #pd.read_excel(fdn + datas1 + '_' + cond1 + '_' + cell_type + '_' + ct + '.xls') 
            interaction_b = res[(datas2, cond2, cell_type, ct)]
        else:
            interaction_a = res[(datas1, cond1, ct, cell_type)]
            interaction_b = res[(datas2, cond2, ct, cell_type)]
        
        for _,row1 in interaction_a.iterrows():
            gene_a1 = row1['gene_name_a']
            gene_a2 = row1['gene_name_b']
            for _,row2 in interaction_b.iterrows():
                gene_b1 = row2['gene_name_a']
                gene_b2 = row2['gene_name_b']
                if gene_a1 == gene_b1 and gene_a2 == gene_b2:
                    same_inters[cell_type].append(row1.tolist() + row2.tolist())

    for cell_type in cell_types:
        same_inters[cell_type]=pd.DataFrame(
            same_inters[cell_type],
            columns=['dataset', 'Condition',
                     'cell_type1','cell_type2', 
                     'gene_name_a', 'gene_name_b',
                     'frac1', 'frac2',
                     'avg1','avg2', 
                     'fra_sum', 'av_sum',
                     'dataset', 'Condition',
                     'cell_type1','cell_type2', 
                     'gene_name_a', 'gene_name_b',
                     'frac1', 'frac2',
                     'avg1','avg2', 
                     'fra_sum', 'av_sum'])

    return same_inters

def interaction_analysis_only_a(datas1, datas2, cond1, cond2, ct, res):
    '''analysis of the interactions between ct and other cell_types,
    Args:
        datas1, cond1 --> requiremets for interaction_a
        datas2, cond2 --> requiremets for interaction_b
        ct --> analysis of ct and other cell_types between interaction_a and interaction_b
        res --> the dict of dataframe containing columns we want to analyze, eg. 
                res = average_fraction(adatag, ['dataset', 'Condition', 'cell_type'], 
                                       axis='obs', greater_than=0, log=False)
    Returns:
    dict containing DataFrame with interactions expressed only in interaction_a
    '''
    from collections import defaultdict
    
    cell_types = ['B_cells','NK_cells' , 'T_cells', 'Monocytes', 'Plasmablasts', 'pDCs', 'cDCs']
    only_a_inters = defaultdict(list)
    
    inter_list = res.keys()
    for cell_type in cell_types:
        if (datas1, cond1, cell_type, ct) in inter_list:
            interaction_a = res[(datas1, cond1, cell_type, ct)]
            interaction_b = res[(datas2, cond2, cell_type, ct)]
        else:
            interaction_a = res[(datas1, cond1, ct, cell_type)]
            interaction_b = res[(datas2, cond2, ct, cell_type)]
        
        gene_bs = []
        for _,row in interaction_b.iterrows(): 
            gene_bs.append([row['gene_name_a'],row['gene_name_b']]) # 'gene_name_a', 'gene_name_b'

        for _,row in interaction_a.iterrows(): 
            gene_a1 = row['gene_name_a']
            gene_a2 = row['gene_name_b']           
            if [gene_a1,gene_a2] not in gene_bs:
                only_a_inters[cell_type].append(row.tolist()) 

    for cell_type in cell_types:
        only_a_inters[cell_type]=pd.DataFrame(
            only_a_inters[cell_type],
            columns=['dataset', 'Condition',
                     'cell_type1','cell_type2', 
                     'gene_name_a', 'gene_name_b',
                     'frac1', 'frac2',
                     'avg1','avg2', 
                     'fra_sum', 'av_sum'])
    
    return only_a_inters

def interaction_analysis_only_b(datas1, datas2, cond1, cond2, ct, res):
    '''analysis of the interactions between ct and other cell_types,
    Args:
        datas1, cond1 --> requiremets for interaction_a
        datas2, cond2 --> requiremets for interaction_b
        ct --> analysis of ct and other cell_types between interaction_a and interaction_b
        res --> the dict of dataframe containing columns we want to analyze, eg. 
                res = average_fraction(adatag, ['dataset', 'Condition', 'cell_type'], 
                                       axis='obs', greater_than=0, log=False)
    Returns:
        dict containing DataFrame with interactions expressed only in interaction_b
    '''

    from collections import defaultdict
    
    cell_types = ['B_cells','NK_cells' , 'T_cells', 'Monocytes', 'Plasmablasts', 'pDCs', 'cDCs']
    only_b_inters = defaultdict(list)
    
    inter_list = res.keys()
    for cell_type in cell_types:
        if (datas1, cond1, cell_type, ct) in inter_list:
            interaction_a = res[(datas1, cond1, cell_type, ct)]
            interaction_b = res[(datas2, cond2, cell_type, ct)]
        else:
            interaction_a = res[(datas1, cond1, ct, cell_type)]
            interaction_b = res[(datas2, cond2, ct, cell_type)]
        
        gene_as = []
        for _,row in interaction_a.iterrows(): 
            gene_as.append([row['gene_name_a'],row['gene_name_b']]) 
       
        for _,row in interaction_b.iterrows(): 
            gene_b1 = row['gene_name_a']
            gene_b2 = row['gene_name_b']
            if [gene_b1,gene_b2] not in gene_as:
                    only_b_inters[cell_type].append(row.tolist())

    for cell_type in cell_types:
        only_b_inters[cell_type]=pd.DataFrame(
            only_b_inters[cell_type],
            columns=['dataset', 'Condition',
                     'cell_type1','cell_type2', 
                     'gene_name_a', 'gene_name_b',
                     'frac1', 'frac2',
                     'avg1','avg2', 
                     'fra_sum', 'av_sum'])
    
    return only_b_inters

def interaction_analysis(datas1, datas2, cond1, cond2, ct, res):
    '''analysis of the interactions between ct and other cell_types,
    Args:
        datas1, cond1 --> requiremets for interaction_a
        datas2, cond2 --> requiremets for interaction_b
        ct --> analysis of ct and other cell_types between interaction_a and interaction_b
        res --> the dict of dataframe containing columns we want to analyze, eg. 
                res = average_fraction(adatag, ['dataset', 'Condition', 'cell_type'], 
                                       axis='obs', greater_than=0, log=False)
    Returns:
        dict containing DataFrame with interactions expressed only in interaction_a
        and dict containing DataFrame with interactions expressed only in interaction_b
        and dict containing DataFrame with interactions expressed both in interaction_a and interaction_b
    '''
    
    from collections import defaultdict
    
    cell_types = ['B_cells','NK_cells' , 'T_cells', 'Monocytes', 'Plasmablasts', 'pDCs', 'cDCs']
    same_inters = defaultdict(list)
    only_a_inters = defaultdict(list)
    only_b_inters = defaultdict(list)
    
    inter_list = res.keys()
    for cell_type in cell_types:
        if (datas1, cond1, cell_type, ct) in inter_list:
            interaction_a = res[(datas1, cond1, cell_type, ct)]
            interaction_b = res[(datas2, cond2, cell_type, ct)]
        else:
            interaction_a = res[(datas1, cond1, ct, cell_type)]
            interaction_b = res[(datas2, cond2, ct, cell_type)]
        
        gene_as = []
        gene_bs = []
        for _,row in interaction_a.iterrows(): 
            gene_as.append([row['gene_name_a'],row['gene_name_b']])
        for _,row in interaction_b.iterrows(): 
            gene_bs.append([row['gene_name_a'],row['gene_name_b']])

        for _,row in interaction_a.iterrows(): 
            gene_a1 = row['gene_name_a']
            gene_a2 = row['gene_name_b']           
            if [gene_a1,gene_a2] not in gene_bs:
                only_a_inters[cell_type].append(row.tolist()) 

        for _,row in interaction_b.iterrows(): 
            gene_b1 = row['gene_name_a']
            gene_b2 = row['gene_name_b']
            if [gene_b1,gene_b2] not in gene_as:
                    only_b_inters[cell_type].append(row.tolist())
        
        for _,row1 in interaction_a.iterrows():
            gene_a1 = row1['gene_name_a']
            gene_a2 = row1['gene_name_b']
            for _,row2 in interaction_b.iterrows():
                gene_b1 = row2['gene_name_a']
                gene_b2 = row2['gene_name_b']
                if gene_a1 == gene_b1 and gene_a2 == gene_b2:
                    same_inters[cell_type].append(row1.tolist() + row2.tolist())

    for cell_type in cell_types:
        same_inters[cell_type]=pd.DataFrame(
            same_inters[cell_type],
            columns=['dataset', 'Condition',
                     'cell_type1','cell_type2', 
                     'gene_name_a', 'gene_name_b',
                     'frac1', 'frac2',
                     'avg1','avg2', 
                     'fra_sum', 'av_sum',
                     'dataset', 'Condition',
                     'cell_type1','cell_type2', 
                     'gene_name_a', 'gene_name_b',
                     'frac1', 'frac2',
                     'avg1','avg2', 
                     'fra_sum', 'av_sum'])
        only_a_inters[cell_type]=pd.DataFrame(
            only_a_inters[cell_type],
            columns=['dataset', 'Condition',
                     'cell_type1','cell_type2', 
                     'gene_name_a', 'gene_name_b',
                     'frac1', 'frac2',
                     'avg1','avg2', 
                     'fra_sum', 'av_sum'])
        only_b_inters[cell_type]=pd.DataFrame(
            only_b_inters[cell_type],
            columns=['dataset', 'Condition',
                     'cell_type1','cell_type2', 
                     'gene_name_a', 'gene_name_b',
                     'frac1', 'frac2',
                     'avg1','avg2', 
                     'fra_sum', 'av_sum'])
    
    return [same_inters, only_a_inters, only_b_inters]

def special_interactions(datas1, datas2, datas3, cond1, cond2, cond3, ct, res, save=None):
    '''analysis of the interactions between ct and other cell_types,
    Args:
        datas1, cond1 --> requiremets for interaction_1
        datas2, cond2 --> requiremets for interaction_2
        datas3, cond3 --> requiremets for interaction_3
        ct --> analysis of ct and other cell_types between interaction_a and interaction_b
        res --> the dict of dataframe containing columns we want to analyze, eg. 
                res = average_fraction(adatag, ['dataset', 'Condition', 'cell_type'], 
                                       axis='obs', greater_than=0, log=False)
    Returns:
        dict containing DataFrame with interactions expressed both in 
        [only in interaction_1 compared with interaction_2] & [only in interaction_1 compared with interaction_3]
    '''
    from collections import defaultdict
    
    cell_types = ['B_cells','NK_cells' , 'T_cells', 'Monocytes', 'Plasmablasts', 'pDCs', 'cDCs']
    same_inters = []    

    for cell_type in cell_types:
        interaction_a = interaction_analysis_only_a(datas1, datas2, cond1, cond2, ct, res)[cell_type]
        interaction_b = interaction_analysis_only_a(datas1,datas3, cond1, cond3, ct, res)[cell_type]
        
        for _,row1 in interaction_a.iterrows():
            gene_a1 = row1['gene_name_a']
            gene_a2 = row1['gene_name_b']
            for _,row2 in interaction_b.iterrows():
                gene_b1 = row2['gene_name_a']
                gene_b2 = row2['gene_name_b']
                if gene_a1 == gene_b1 and gene_a2 == gene_b2:
                    same_inters.append(row1.tolist())

    same_inters = pd.DataFrame(
            same_inters,
            columns = ['dataset', 'Condition',
                      'cell_type1','cell_type2', 
                      'gene_name_a', 'gene_name_b',
                      'frac1', 'frac2',
                      'avg1','avg2', 
                      'fra_sum', 'av_sum'])

    if save is not None:
        same_inters.to_excel(save)
        
    #fdn = '/home/yike/phd/dengue/data/excels/interaction_analysis/'
    #same_inters.to_excel(fdn + ct + '_special_inters.xls')

    return same_inters

def interaction_by_ct_kdeplot(adata, interactions, save=None):
    '''kdeplot of gene expression in interactions, the plots shown as:
    subplots by different cell_types,
    in each subplot: 'density of cells' vs 'gene expression [log10(cpm+0.1)]' in different samples (IDs),
    Args:
        adata: --> the data that want to analyze, eg. adata_sick, adata_healthy
        interactions --> a list containing interaction pairs, eg. [['TNFRSF14', 'BTLA'],['TNFSF13B', 'TNFRSF13B']]
    Returns:
        Return nothing
        will save plots into '/home/yike/phd/dengue/figures/express_by_sick' with the name as gene1+'-'+gene2+'_by_ct_kde.png'
    '''
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    cell_types = list(adata.obs['cell_type'].astype('category').cat.categories)
    IDs = list(adata.obs['ID'].astype('category').cat.categories)
    n_column = len(cell_types) +1 
    
    colors = sns.color_palette('hls', len(IDs))
    sns.set_palette(colors)

    from collections import defaultdict
    # adata_dic = defaultdic(list)
    adata_dic = {}
    for cell_type in cell_types:
        adata_ct = adata[adata.obs['cell_type'] == cell_type]
        for ID in IDs:
            adata_dic[(cell_type, ID)] = adata_ct[adata_ct.obs['ID'] == ID]
            adata_dic[('all', ID)] = adata[adata.obs['ID'] == ID]

    for interaction in interactions:
        gene1 = interaction[0]
        gene2 = interaction[1]
        fig, axs = plt.subplots(n_column, 2, figsize = (8, 22), dpi=80, facecolor='white', sharex=True, sharey=True)
        for i in range(n_column):
            
            for ID in IDs:
                if i == 0: # epression in all samples
                    cell_type = 'all'
                else:      # expression in different samples (IDs)
                    cell_type = cell_types[i-1]
                
                gene1_ct = adata_dic[(cell_type, ID)][:, gene1].X.toarray()[:, 0]
                np_gene1_ct = np.log10(0.1 + gene1_ct)
                sns.kdeplot(np_gene1_ct, bw_method=0.5, ax=axs[i, 0], bw=0.1, label='_nolegend_')
                
                gene2_ct = adata_dic[(cell_type, ID)][:, gene2].X.toarray()[:, 0]
                np_gene2_ct = np.log10(0.1 + gene2_ct)
                if i == 0:
                    sns.kdeplot(np_gene2_ct, bw_method=0.5, ax=axs[i, 1], bw=0.1, label=ID)
                else:
                    sns.kdeplot(np_gene2_ct, bw_method=0.5, ax=axs[i, 1], bw=0.1, label='_nolegend_')
                
                gene1_ct_avg = np.log10(0.1 + np.mean(gene1_ct))
                gene2_ct_avg = np.log10(0.1 + np.mean(gene2_ct))
                sns.scatterplot(x=[gene1_ct_avg], y=[4], ax=axs[i, 0])
                sns.scatterplot(x=[gene2_ct_avg], y=[4], ax=axs[i, 1])
            
            axs[i, 0].set_ylim(0, 5)
            axs[i, 1].set_ylim(0, 5)
            axs[i, 0].tick_params(labelsize=15)
            axs[i, 1].tick_params(labelsize=15)
            axs[i, 0].set_ylabel(None)
            axs[i, 1].set_ylabel(None)
            axs[0, 1].legend(loc='upper right', ncol=2, markerscale=0.1, frameon=False, framealpha=0, fontsize=8, columnspacing=12.5)
            #axs[0, 1].legend(loc='upper left', ncol=2, bbox_to_anchor=(1, 1, 0, 0), bbox_transform=axs[0, 1].transAxes, fontsize=6)
            axs[i, 0].set_title(gene1 + ' in ' + cell_type, fontsize=15)
            axs[i, 1].set_title(gene2 + ' in ' + cell_type, fontsize=15)
        
        fig.text(0, 0.5, 'density of cells', va='center', rotation='vertical', fontsize=20)
        fig.text(0.5, 0, 'gene expression [log10(cpm+0.1)]', ha='center', fontsize=20)
        
        fig.tight_layout(rect=[0.05,0.02,1,1], pad=0, h_pad=0, w_pad=0)
        
        if save is not None:
            fig.savefig(save + gene1 + '-' + gene2 + '_by_ct_kde.png')
             

def interaction_by_ID_kdeplot(adata, interactions):
    '''kdeplot of gene expression in interactions, the plots shown as:
    subplots by different samples (IDs),
    in each subplot: 'density of cells' vs 'gene expression [log10(cpm+0.1)]' in different cell_types,
    Args:
        adata: --> the data that want to analyze, eg. adata_sick, adata_healthy
        interactions --> a list containing interaction pairs, eg. [['TNFRSF14', 'BTLA'],['TNFSF13B', 'TNFRSF13B']]
    Returns:
        Return nothing
        will save plots into '/home/yike/phd/dengue/figures/express_by_sick' with the name as gene1+'-'+gene2+'_by_ID_kde.png'
    '''
    import matplotlib.pyplot as plt
    import seaborn as sns

    cell_types = list(adata.obs['cell_type'].astype('category').cat.categories)
    IDs = list(adata.obs['ID'].astype('category').cat.categories)
    n_column = len(IDs) +1 
    
    colors = sns.color_palette('hls', len(cell_types))
    sns.set_palette(colors)

    from collections import defaultdict
    adata_dic = {}
    for ID in IDs:
        adata_ID = adata[adata.obs['ID'] == ID]
        for cell_type in cell_types:
            adata_dic[(ID, cell_type)] = adata_ID[adata_ID.obs['cell_type'] == cell_type]
            adata_dic[('all', cell_type)] = adata[adata.obs['cell_type'] == cell_type]

    for interaction in interactions:
        gene1 = interaction[0]
        gene2 = interaction[1]
        
        fig, axs = plt.subplots(n_column, 2, figsize=(8, 20), dpi=80, facecolor='white')
        for i in range(n_column):
            for cell_type in cell_types:
                if i == 0: # epression in all samples
                    ID = 'all'
                else:      # expression in different samples (IDs)
                    ID = IDs[i-1] 
                
                gene1_ct = adata_dic[(ID, cell_type)][:, gene1].X.toarray()[:, 0]
                np_gene1_ct = np.log10(0.1 + gene1_ct)
                sns.kdeplot(np_gene1_ct, bw_method=0.5, ax=axs[i, 0], bw=0.1, label='_nolegend_')
                
                gene2_ct = adata_dic[(ID, cell_type)][:, gene2].X.toarray()[:, 0]
                np_gene2_ct = np.log10(0.1 + gene2_ct)
                
                if i == 0:
                    sns.kdeplot(np_gene2_ct, bw_method=0.5, ax=axs[i, 1], bw=0.1, label=cell_type)
                else:
                    sns.kdeplot(np_gene2_ct, bw_method=0.5, ax=axs[i, 1], bw=0.1, label='_nolegend_')
                    
                gene1_ct_avg = np.log10(0.1 + np.mean(gene1_ct))
                gene2_ct_avg = np.log10(0.1 + np.mean(gene2_ct))
                sns.scatterplot(x=[gene1_ct_avg], y=[4], ax=axs[i, 0])
                sns.scatterplot(x=[gene2_ct_avg], y=[4], ax=axs[i, 1])
                    
            axs[i, 0].set_ylim(0, 5)
            axs[i, 1].set_ylim(0, 5)
            axs[i, 0].tick_params(labelsize=15)
            axs[i, 1].tick_params(labelsize=15)
            axs[i, 0].set_ylabel(None)
            axs[i, 1].set_ylabel(None)
            axs[0, 1].legend(loc='upper right', markerscale=0.1, frameon=False, framealpha=0, fontsize=8)
            #axs[0, 1].legend(loc='upper left', ncol=2, bbox_to_anchor=(1, 1, 0, 0), bbox_transform=axs[0, 1].transAxes, fontsize=6)
            axs[i, 0].set_title(gene1 + ' in ' + ID, fontsize=15)
            axs[i, 1].set_title(gene2 + ' in ' + ID, fontsize=15)
        
        fig.text(0, 0.5, 'density of cells', va='center', rotation='vertical', fontsize=20)
        fig.text(0.5, 0, 'gene expression [log10(cpm+0.1)]', ha='center', fontsize=20)
        
        fig.tight_layout(rect=[0.05,0.02,1,1], pad=0, h_pad=0, w_pad=0)
        
        if save is not None:
            fig.savefig(save + gene1 + '-' + gene2 + '_by_ID_kde.png')



