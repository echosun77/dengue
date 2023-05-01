import pandas as pd 
import numpy as np

import anndata
import scanpy as sc

import matplotlib.pyplot as plt
import seaborn as sns
import anndataks

def pair_comparison(adata, ct_obs, cd_SD, cd_D, cell_types, gene_cut_off, log1p=False):
    ress = pd.DataFrame([])
    log2FC = pd.DataFrame([])
    
    for cell_type in cell_types:
        adata_ct = adata[adata.obs[ct_obs].isin(cell_type)]
        
        ####### filter out genes expressed less than gene_cut_off in all patients
        IDs_o = list(adata_ct.obs['ID'].unique())
        fra = [np.asarray((adata_ct[adata_ct.obs['ID'] == ID].X > 0).mean(axis=0))[0] for ID in IDs_o]
        
        fra=pd.DataFrame(fra, index=IDs_o, columns=adata_ct.var_names).T
        gene_list = fra.index[(fra.values > gene_cut_off).sum(axis=1) > 0]
        adata_ct = adata_ct[:, gene_list]
        
        ####### filter out patient with less than 5 cells for the cell type
        IDs = []
        for ID in IDs_o:
            n = adata_ct[adata_ct.obs['ID'] == ID].obs.shape[0]
            n_thresh = [3, 5][ct_obs == 'cell_type_new']
            if n >= n_thresh:
                IDs.append(ID)
                
        ####### 
        adata_ct = adata_ct[adata_ct.obs['ID'].isin(IDs)]
        
        adata_S_ct = adata_ct[adata_ct.obs['Condition'].isin(cd_SD)]
        adata_NS_ct = adata_ct[adata_ct.obs['Condition'].isin(cd_D)]
        
        IDs_S = list(adata_S_ct.obs['ID'].unique())
        IDs_NS = list(adata_NS_ct.obs['ID'].unique())
                
        ####### pair comparison
        log2_fc = []
        for ID_S in IDs_S:
            adata_S_ID = adata_S_ct[adata_S_ct.obs['ID'] == ID_S]
            
            for ID_NS in IDs_NS:
                adata_NS_ID = adata_NS_ct[adata_NS_ct.obs['ID'] == ID_NS]

                X_S = adata_S_ID.X
                X_NS = adata_NS_ID.X
                avg_S = np.asarray(X_S.mean(axis=0))[0]
                avg_NS = np.asarray(X_NS.mean(axis=0))[0]

                if log1p is False:
                    log2_fc.append(np.log2(avg_S + 1) - np.log2(avg_NS + 1))
                elif log1p not in (True, 2):
                    log2_fc.append((avg_S - avg_NS) / np.log2(log1p))
                else:
                    log2_fc.append(avg_S - avg_NS)

        log2_fc = np.asarray(log2_fc)
        
        if log2_fc.size == 0:
            continue
        med_pair = np.median(log2_fc, axis=0)
        pos_fra_pair = [len(log2_fc[:, i][log2_fc[:, i] > 0])/log2_fc.shape[0] for i in range(log2_fc.shape[1])]
        neg_fra_pair = [len(log2_fc[:, i][log2_fc[:, i] < 0])/log2_fc.shape[0] for i in range(log2_fc.shape[1])]

        index_name = ''.join(str(item) + '_' for item in cell_type)
        res = pd.DataFrame([], index=adata_ct.var_names)
        res['med_pair'] = med_pair
        res['fra_pair'] = pos_fra_pair
        res['neg_fra_pair'] = neg_fra_pair
        res['cell_subtype'] = index_name
        ress = pd.concat([ress, res], join='outer')
        ress['gene'] = ress.index.tolist()

        FCs = pd.DataFrame(log2_fc, columns = adata_ct.var_names, index=[index_name] * log2_fc.shape[0])
        log2FC = pd.concat([log2FC, FCs], join='outer')
        
    return {'pair_res': ress, 'log_FCs': log2FC}


def fra_avg(adata, ct_obs, cell_types, cds, log1p=False):    

    fra = pd.DataFrame([])
    avg = pd.DataFrame([])

    # cds = {'SD': ['S_dengue'], 'N': ['dengue', 'DWS']}

    # if 'S_dengue' in adata.obs['Condition'].unique():
    #     cds = ['dengue', 'S_dengue', 'Healthy', 'DWS']
    # elif 'Severe' in adata.obs['Condition'].unique():
    #     cds = ['Moderate', 'Severe']

    # if 'donor' in adata.obs.columns:
    #     p_ID = 'donor'
    # elif 'ID' in adata.obs.columns:
    #     p_ID = 'ID'

    # for cell_type in cell_types:
    #     ct_obs = 'cell_subtype_new'
    #     if cell_type in adata.obs['cell_type_new'].unique():
    #         ct_obs = 'cell_type_new'

    adata_cd = {cd: (adata[adata.obs['Condition'].isin(cd_ls)]) for cd, cd_ls in cds.items()}

    for cell_type in cell_types:
        for cd in cds:
            if cell_type not in adata_cd[cd].obs[ct_obs].unique():
                continue
            adata_ct = adata_cd[cd][adata_cd[cd].obs[ct_obs] == cell_type]
            
            if adata_ct.obs.shape[0] == 0:
                continue
            fra_ct = np.asarray((adata_ct.X > 0).mean(axis=0))[0]
            fra_ct = pd.DataFrame(fra_ct, columns=['fra'], index=adata_ct.var_names)
            fra_ct[ct_obs] = cell_type
            fra_ct['condition'] = cd
            fra = pd.concat([fra, fra_ct])
            fra['gene'] = fra.index.tolist()

            avg_ct = np.asarray(adata_ct.X.mean(axis=0))[0]
                        
            if log1p == 2:
                avg_ct = np.exp2(avg_ct) -1
            elif log1p not in (True, 2):
                avg_ct = np.exp2(np.log2(log1p) * avg_ct) - 1

            avg_ct = pd.DataFrame(avg_ct, columns=['avg'], index=adata_ct.var_names)
            avg_ct[ct_obs] = cell_type
            avg_ct['condition'] = cd
            avg = pd.concat([avg, avg_ct])
            avg['gene'] = avg.index.tolist()
            
    return {'fra': fra, 'avg': avg}

def combination(pair, fra, avg):
    '''
    fra.index = ['cell_subtype_2', 'gene']
    '''
    pair = pair.set_index(['cell_subtype', 'gene'])
    
    ct = 'cell_type_new'
    if 'cell_subtype_new' in fra.columns:
        ct = 'cell_subtype_new'
    fra = fra.set_index([ct, 'gene', 'condition'])
    avg = avg.set_index([ct, 'gene', 'condition'])
        
    NS_idx = [(i[0], i[1], 'N') for i in pair.index]
    S_idx = [(i[0], i[1], 'SD') for i in pair.index]
    
    pair['S_fra'] = (fra.loc[S_idx])['fra'].tolist()
    pair['NS_fra'] = (fra.loc[NS_idx])['fra'].tolist()
    
    pair['S_avg'] = (avg.loc[S_idx])['avg'].tolist()
    pair['NS_avg'] = (avg.loc[NS_idx])['avg'].tolist()
    
    return pair

def ks(adata, cell_types, ct_obs, cd_SD, cd_D, log1p=2):

    results = {}
    for cell_type in cell_types:
        adata_ct = adata[adata.obs[ct_obs] == cell_type]
        if cell_type == 'cDCs':
            adata_ct = adata_ct[~adata_ct.obs['ID'].isin(['1_140_01', '5_193_01'])]

        adata_SD = adata_ct[adata_ct.obs['Condition'].isin(cd_SD)]
        adata_D = adata_ct[adata_ct.obs['Condition'].isin(cd_D)]
        results[cell_type] = anndataks.compare(adata_D, adata_SD, log1p=log1p)

    res = pd.DataFrame([])
    for cell_type in cell_types:
        results[cell_type]['cell_type'] = [cell_type] * results[cell_type].shape[0]
        res = pd.concat([res, results[cell_type]] )
    res.index.name = 'gene'
    
    return res

def random_pair_comparison(adata, cell_types, gene_cut_off, cell_cut_off, compare_number, log1p=False):    
    import random
    ress = pd.DataFrame([])
    log2FC = pd.DataFrame([])
    
    cds = ['dengue', 'S_dengue']
    p_ID = 'ID'
    
#     adata_cd = {cd: (adata[adata.obs['Condition'] == cd]) for cd in cds}
    for cell_type in cell_types:
        ct_obs = 'cell_subtype_2'
        if cell_type in adata.obs['cell_type'].unique():
            ct_obs = 'cell_type'
            
        ####### filter out genes expressed less than gene_cut_off in all patients
        adata_ct = adata[adata.obs[ct_obs] == cell_type]
        IDs = adata_ct.obs['ID'].unique()
        fra = [np.asarray((adata_ct[adata_ct.obs['ID'] == ID].X > 0).mean(axis=0))[0] for ID in IDs]
        fra=pd.DataFrame(fra, index=IDs, columns=adata.var_names).T
        gene_list = fra.index.tolist()
        
        if gene_cut_off is not False:
            for idx, row in fra.iterrows():
                n = 0
                for i in row.tolist():
                    if i <= gene_cut_off:
                        n += 1
                if n == len(IDs):
                    gene_list.remove(idx)
        
        adata_ct = adata_ct[:, gene_list]
        
        ####### filter out patient with less than 5 cells for the cell type
        adata_S_ct = adata_ct[adata_ct.obs['Condition'] == cds[1]]
        adata_NS_ct = adata_ct[adata_ct.obs['Condition'] == cds[0]]
        
        IDs_S = list(adata_S_ct.obs[p_ID].unique())
        for ID_S in IDs_S:
            if adata_S_ct[adata_S_ct.obs[p_ID] == ID_S].obs.shape[0] < cell_cut_off:
                IDs_S.remove(ID_S)
        
        IDs_NS = list(adata_NS_ct.obs[p_ID].unique())
        for ID_NS in IDs_NS:
            if adata_NS_ct[adata_NS_ct.obs[p_ID] == ID_NS].obs.shape[0] < cell_cut_off:
                IDs_NS.remove(ID_NS)
                
        ####### random pair comparison
        log2_fc = []
        for i in range(compare_number):
            ID_S = random.sample(IDs_S, 1)[0]
            adata_S_ID = adata_S_ct[adata_S_ct.obs[p_ID] == ID_S]

            ID_NS = random.sample(IDs_NS, 1)[0]
            adata_NS_ID = adata_NS_ct[adata_NS_ct.obs[p_ID] == ID_NS]

            X_S = adata_S_ID.X
            X_NS = adata_NS_ID.X
            avg_S = np.asarray(X_S.mean(axis=0))[0]
            avg_NS = np.asarray(X_NS.mean(axis=0))[0]

            if log1p is False:
                log2_fc.append(np.log2(avg_S + 1) - np.log2(avg_NS + 1))
            elif log1p not in (True, 2):
                log2_fc.append((avg_S - avg_NS) / np.log2(log1p))
            else:
                log2_fc.append(avg_S - avg_NS)

        log2_fc = np.asarray(log2_fc)
        
        if log2_fc.size == 0:
            continue
        med_pair = np.median(log2_fc, axis=0)
        fra_pair = [len(log2_fc[:, i][log2_fc[:, i] > 0])/log2_fc.shape[0] for i in range(log2_fc.shape[1])]

        res = pd.DataFrame([], index=adata_ct.var_names)
        res['med_pair'] = med_pair
        res['fra_pair'] = fra_pair
        res['cell_subtype'] = cell_type
        ress = pd.concat([ress, res], join='outer')
        ress['gene'] = ress.index.tolist()

        FCs = pd.DataFrame(log2_fc, columns = adata_ct.var_names, index=[cell_type] * log2_fc.shape[0])
        log2FC = pd.concat([log2FC, FCs], join='outer')
        
    return {'pair_res': ress, 'log_FCs': log2FC}