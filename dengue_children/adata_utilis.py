import pandas as pd 
import numpy as np 

def split(adata, columns, axis='obs'):
    '''Split AnnData by metadata column
    
    Args:
        adata (AnnData): The object to split
        columns (str or list): column(s) of metadata table to split by
        axis (str): 'obs' (default) or 'var'
        
    Returns:
        dict of AnnData with keys equal to the unique values found in that
        columns of the appropriate metadata (adata.obs or adata.var)
    '''
    if axis not in ('obs','var'):
        raise ValueError('axis must be "obs" or "var"')
    
    # getattr()用于返回一个对象属性值，eg.getattr(a, 'bar')获取a的‘bar’属性值
    meta = getattr(adata, axis)
    
    # isinstance() 函数来判断一个对象是否是一个已知的类型，类似type()
    if isinstance(columns, str):
        column = columns
        if column not in meta.columns:
            raise ValueError('column not found')
        metac = meta[column]
    else:
        if len(columns) < 1:
            return adata
        elif len(columns) == 1:
            return split(adata, columns[0], axis=axis)
        
        # Merge with @ (a hack, but alright in most cases)
        metac = meta[columns[0]].astype(str)
        for i, col in enumerate(columns[1:]):
            metac = metac + '@' + meta[col].astype(str)
            
    # Unique values
    metau = metac.unique()
    
    # Construct dict
    d = {}
    for key in metau:
        idx = metac.index[metac == key]
        if axis == 'obs':
            asub = adata[idx]
        else:
            asub = adata[:, idx]
    
        if not isinstance(columns, str) and (len(columns) > 1):
            key = key.split('@')
            for i, col in enumerate(columns):
                key[i] = type(meta[col].iloc[0])(key[i])
            key = tuple(key)

        d[key] = asub
    
    return d    

def average(adata, columns, axis='obs', log=False):
    '''Average expression by metadata column
    
    Args:
        adata (AnnData): The project to split
        columns (str or list): column(s) of metadata table to split by
        axis (str): 'obs' (default) or 'var'
        log: if False, return arithmetic average, if True, returns geometric average
    Returns:
        pandas DataFrame with rows equal to the non-integrated axis names and columns
        equal to the keys of the AnnData split.
    '''
    
    adatad = split(adata, columns, axis=axis)
    
    if axis == 'var':
        iax = 1
        index = adata.obs_names
    else:
        iax = 0
        index = adata.var_names
    
    ave_expre = {}
    for key, adatai in adatad.items():
        matrix = adatai.X
        # if log=False, get arithmetric average, elif log=True, get geometric average
        if log:
            matrix = matrix.copy()
            matrix.data = np.log(matrix.data + 0.1) # loge
        av_ex = np.asarray(adatai.X .mean(axis=iax))[0]
        # av_ex = np.asarray(adatai.X.sum(axis=iax) / (adatai.X > 0).sum(axis=iax))[0] # number of only expressed genes
        if log:
            av_ex = np.exp(av_ex)
        
        ave_expre[key] = av_ex
    
    ave_expre = pd.DataFrame(ave_expre, index=index)
    ave_expre.name = 'Average expression'
    
    if isinstance(columns, str):
        ave_expre.columns.name = columns
    else:
        ave_expre.columns.names = columns
    
    return ave_expre

def expressing_fractions(adata, columns, axis='obs', greater_than=0):
    '''Fraction of expressors by metadata column
    Args:
        adata (AnnData): The object to split
        columns (str or list): column(s) of metadata table to split by
        axis (str): 'obs' (default) or 'var'
        greater_than (float): only expressors stricly above this threshold are
          counted
    Returns:
        pandas DataFrame with rows equal to the non-integrated axis names
        and columns equal to the keys of the AnnData split.
    '''
    adatad = split(adata, columns, axis=axis)

    if axis == 'var':
        iax = 1
        index = adata.obs_names
    else:
        iax = 0
        index = adata.var_names

    fracd = {}
    for key, adatai in adatad.items():
        fr = np.asarray((adatai.X > greater_than).mean(axis=iax))[0]
        fracd[key] = fr
    fracd = pd.DataFrame(fracd, index=index)
    fracd.name = 'Fraction of expressing cells'

    if isinstance(columns, str):
        fracd.columns.name = columns
    else:
        fracd.columns.names = columns

    return fracd

def average_fraction(adata, columns, axis='obs', greater_than=0, log=False):
    '''Fraction of expressors and average expression by metadata column
    Args:
        adata (AnnData): The object to split
        columns (str or list): column(s) of metadata table to split by
                               from left to right [datas, cond, cell_type]
        axis (str): 'obs' (default) or 'var'
        greater_than (float): only expressors stricly above this threshold are
          counted
        log: if False, return arithmetic average, if True, returns geometric average
    Returns:
        pandas DataFrame with rows equal to the non-integrated axis names
        and columns equal to the keys of the AnnData split.'''
    
    print('Load interaction')
    fn_int = '/home/yike/phd/dengue/data/interaction_unpacked_mouse.tsv'
    interactions = pd.read_csv(fn_int, sep='\t')[['gene_name_a', 'gene_name_b']]
    ga, gb = interactions['gene_name_a'], interactions['gene_name_b']
    
    fracd = expressing_fractions(adata, columns, axis='obs', greater_than=0)
    avgd = average(adata, columns, axis='obs', log=False)
    stats = {
        'frac_exp': fracd,
        'avg_exp': avgd,
    }

    # Flexible criterion
    criterion = {'key': 'frac_exp', 'threshold': 0.1}
    criterion = {'key': 'avg_exp', 'threshold': {'child': 60, 'adult': 35}}

    from collections import defaultdict
    th = criterion['threshold']
    cell_types = ['B_cells','NK_cells' , 'T_cells', 'Monocytes', 'Plasmablasts', 'pDCs', 'cDCs']
    res = {}
    for col in fracd.columns:
        datas, cond, cell_type1 = col
        for cell_type2 in cell_types:
            col2 = (datas, cond, cell_type2)
            res[(datas, cond, cell_type1, cell_type2)] = []
            fra = fracd.loc[ga, col].values
            frb = fracd.loc[gb, col2].values
            avga = avgd.loc[ga, col].values
            avgb = avgd.loc[gb, col2].values
            key = criterion['key']
            if isinstance(th, dict):
                th1 = th[col[0]]
            else:
                th1 = th

            ind = (stats[key].loc[ga, col].values > th1) & (stats[key].loc[gb, col2].values > th1)
            ind = ind.nonzero()[0]
            for i in ind:
                resi = {
                    'dataset': datas,
                    'Condition': cond,
                    'cell_type1': cell_type1,
                    'cell_type2': cell_type2,
                    'gene_name_a': interactions.iloc[i]['gene_name_a'],
                    'gene_name_b': interactions.iloc[i]['gene_name_b'],
                    'frac1': fra[i],
                    'frac2': frb[i],
                    'avg1': avga[i],
                    'avg2': avgb[i],
                }
                res[(datas, cond, cell_type1, cell_type2)].append(resi)
            res[(datas, cond, cell_type1, cell_type2)] = pd.DataFrame(res[(datas, cond, cell_type1, cell_type2)])     
    
    for value in res.values():
        value['fra_sum'] = value['frac1'] + value['frac2']
        value['av_sum'] = value['avg1'] + value['avg2']
        
    # combine dataframe of T_cell, B_cell and dataframe of B_cell, T_cell
    merge_res_list = []
    datasets = list(adata.obs[columns[0]].astype('category').cat.categories)
    conds = list(adata.obs[columns[1]].astype('category').cat.categories)
    cell_types= list(adata.obs[columns[2]].astype('category').cat.categories)
    for datas in datasets:
        for cond in conds:
            for i, cell_type1 in enumerate(cell_types):
                for j, cell_type2 in enumerate(cell_types[: i+1]):
                    merge_res_list.append((datas, cond, cell_type1, cell_type2))
 
    res2={}
    for key, value in res.items():
        (datas, cond, cell_type1, cell_type2) = key
        key2 = (datas, cond, cell_type2, cell_type1)
        if key in merge_res_list:
            res[key2] = res[key2] [['dataset', 'Condition',
                                    'cell_type2','cell_type1', 
                                    'gene_name_b', 'gene_name_a',
                                    'frac2','frac1', 
                                    'avg2', 'avg1', 
                                    'fra_sum', 'av_sum']]
            res[key2].columns=['dataset', 'Condition',
                               'cell_type1','cell_type2', 
                               'gene_name_a', 'gene_name_b',
                               'frac1', 'frac2',
                               'avg1','avg2', 
                               'fra_sum', 'av_sum']
            res2[key]=pd.merge(res[key], res[key2],on=['dataset', 'Condition',
                                                       'cell_type1','cell_type2', 
                                                       'gene_name_a', 'gene_name_b',
                                                       'frac1', 'frac2',
                                                       'avg1','avg2', 
                                                       'fra_sum', 'av_sum'],
                               how='outer')
            res2[key] = res2[key].sort_values('av_sum',ascending=False)
    
    return res2