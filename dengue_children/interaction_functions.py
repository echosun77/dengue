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
        
        gene_as = []
        gene_bs = []
        for _,row in interaction_a.iterrows(): 
            gene_as.append([row['gene_name_a'],row['gene_name_b']]) # 'gene_name_a', 'gene_name_b'
        for _,row in interaction_b.iterrows(): 
            gene_bs.append([row['gene_name_a'],row['gene_name_b']]) # 'gene_name_a', 'gene_name_b'

        for _,row in interaction_a.iterrows(): 
            gene_a1 = row[4]
            gene_a2 = row[5]           
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
        gene_bs = []
        for _,row in interaction_a.iterrows(): 
            gene_as.append([row['gene_name_a'],row['gene_name_b']]) 
        for _,row in interaction_b.iterrows(): 
            gene_bs.append([row['gene_name_a'],row['gene_name_b']]) 

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

