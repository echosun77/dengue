import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict

def get_sig_inters(gene_dic1, gene_dic2):
    # {'T_cells': ['TYROBP', 'XCL1', 'KIR3DL1', 'GPR25'], 'Monocytes': ['CD8A','HAVCR2'], 'cDCs': ['CSF2RB','HLA-G'],
    #'B_cells': ['CCL4'], 'pDCs': ['NAMPT','MS4A4A'], 'NK_cells': ['KLRB1'], 'Plasmablasts': ['TIMP1']}
     
    cta = []
    ctb = []
    gene_a = []
    gene_b = []

    print('Load interaction')
    fn_int = '/home/yike/phd/dengue/data/interaction_source_file/interactions_DB.tsv'
    interactions = pd.read_csv(fn_int, sep=',')[['gene_name_a', 'gene_name_b']]

    print('Get interactions')
    for _, row in interactions.iterrows():
        ga = row['gene_name_a']
        gb = row['gene_name_b']
        for ct1 in gene_dic1.keys():
            for ct2 in gene_dic2.keys():
                if (ga in gene_dic1[ct1]) & (gb in gene_dic2[ct2]):
                    cta.append(ct1)
                    ctb.append(ct2)
                    gene_a.append(ga)
                    gene_b.append(gb)
                if (gb in gene_dic1[ct1]) & (ga in gene_dic2[ct2]):
                    cta.append(ct1)
                    ctb.append(ct2)
                    gene_a.append(gb)
                    gene_b.append(ga)
                
    inters = pd.DataFrame([])
    inters['cta'] = cta
    inters['ctb'] = ctb
    inters['ga'] = gene_a
    inters['gb'] = gene_b

    idx = inters[inters['cta'] == inters['ctb']].index[::2]
    inters.drop(idx, inplace=True)
    return inters

def inter_number(inters_df, vmax, trend):
    it_n = inters_df.groupby(['cta', 'ctb']).size().unstack(fill_value=0)
    cell_types = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
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
    sns.heatmap(it_n.T, ax=ax, cmap='plasma', linecolor='w', linewidths=1, vmin=0, vmax=vmax)

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

    ax.axvline(7, c='black')
    ax.axhline(7, c='black')
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    xlabels = [i.get_text() for i in ax.get_xticklabels()]
    for i, ct in enumerate(xlabels):
        if ct in ['B_cells', 'T_cells', 'NK_cells']:
            xlabels[i] = ct.replace('_', ' ')
    ax.set_xticklabels(xlabels)
    #xlabels.reverse()
    ax.set_yticklabels(xlabels)
    ax.text(9.05, 3.8, 'Number of interactions', verticalalignment='center', rotation=90)
    ax.set_title('\n %sregulated'%trend)
    return {'figure': fig, 'ax': ax}

def inter_mix_number(inters_df, vmax):

    cell_types = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
    it_n = inters_df.groupby(['cta', 'ctb']).size().unstack(fill_value=0)

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

    if it_n.loc[best_idx[-1]][best_idx[-1]] > it_n.loc[best_idx[0]][best_idx[0]]:
        best_idx = best_idx
    else:
        best_idx.reverse()

    it_n = it_n.loc[best_idx]
    it_n = it_n[best_idx]
    #it_n = it_n.iloc[best_idx]
    #lkg_col = linkage(pdist(it_n.T.values), optimal_ordering=True)
    #best_col = leaves_list(lkg_col)
    #best_col = it_n.columns[best_col]
    #it_n = it_n[best_col]

    fig, ax = plt.subplots(figsize=[3, 2], dpi=300)
    sns.heatmap(it_n.T, ax=ax, cmap='plasma', linecolor='w', linewidths=1, vmin=0, vmax=vmax)
    ax.set_xlabel(None)
    ax.set_ylabel(None)

    ax.axvline(0, c='black')
    ax.axhline(0, c='black')

    ax.axvline(7, c='black')
    ax.axhline(7, c='black')

    xlabels = [i.get_text() for i in ax.get_xticklabels()]
    for i, ct in enumerate(xlabels):
        if ct in ['B_cells', 'T_cells', 'NK_cells']:
            xlabels[i] = ct.replace('_', ' ')
    ax.text(9.05, 3.8, 'Number of interactions', verticalalignment='center', rotation=90)
    ax.set_title('Mixregulated')

    ax.set_xticklabels(xlabels, c='r')
    ax.set_yticklabels(xlabels, c='b')

    return {'figure': fig, 'ax': ax}

# define a function to collect interactions meeting both log2_FC and exp fraction
def get_inters(fdn, genes, frac_n):
    ####################################### interactions both up --> head(50), both down --> tail(50)
    log2_fc_ave = pd.read_csv(fdn, index_col=['cell_type', 'gene'])
    
    cell_types = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
    both_up = defaultdict(list)
    both_down = defaultdict(list)

    fn_int = '/home/yike/phd/dengue/data/interaction_source_file/interactions_DB.tsv'
    interactions = pd.read_csv(fn_int, sep=',')[['gene_name_a', 'gene_name_b']]
    for _, row in interactions.iterrows():
        ga = row['gene_name_a']
        gb = row['gene_name_b']
        if (ga not in genes) | (gb not in genes):
            continue
        for ct1 in cell_types:
            for ct2 in cell_types:
                ga_ct1_log2FC = log2_fc_ave.loc[ct1, ga]['fold_2_change']
                ga_ct2_log2FC = log2_fc_ave.loc[ct2, ga]['fold_2_change']
                gb_ct1_log2FC = log2_fc_ave.loc[ct1, gb]['fold_2_change']
                gb_ct2_log2FC = log2_fc_ave.loc[ct2, gb]['fold_2_change']

                ga_ct1_comfrac = log2_fc_ave.loc[ct1, ga]['comp_frac']
                ga_ct2_comfrac = log2_fc_ave.loc[ct2, ga]['comp_frac']
                gb_ct1_comfrac = log2_fc_ave.loc[ct1, gb]['comp_frac']
                gb_ct2_comfrac = log2_fc_ave.loc[ct2, gb]['comp_frac']
                
                ct1_head = log2_fc_ave.loc[ct1][log2_fc_ave.loc[ct1]['fold_2_change'] > 0].index.tolist()
                ct2_head = log2_fc_ave.loc[ct2][log2_fc_ave.loc[ct2]['fold_2_change'] > 0].index.tolist()
                ct1_tail = log2_fc_ave.loc[ct1][log2_fc_ave.loc[ct1]['fold_2_change'] < 0].index.tolist()
                ct2_tail = log2_fc_ave.loc[ct2][log2_fc_ave.loc[ct2]['fold_2_change'] < 0].index.tolist()

                if (ga in ct1_head) & (gb in ct2_head):
                    both_up[(ct1, ct2)].append([ga, ga_ct1_log2FC, ga_ct1_comfrac, gb, gb_ct2_log2FC, gb_ct2_comfrac])
                if (gb in ct1_head) & (ga in ct2_head):
                    both_up[(ct1, ct2)].append([gb, gb_ct1_log2FC, gb_ct1_comfrac, ga, ga_ct2_log2FC, ga_ct2_comfrac])
                if (ga in ct1_tail) & (gb in ct2_tail):
                    both_down[(ct1, ct2)].append([ga, ga_ct1_log2FC, ga_ct1_comfrac, gb, gb_ct2_log2FC, gb_ct2_comfrac])
                if (gb in ct1_tail) & (ga in ct2_tail):
                    both_down[(ct1, ct2)].append([gb, gb_ct1_log2FC, gb_ct1_comfrac, ga, ga_ct2_log2FC, ga_ct2_comfrac])
    
    ####################################### both log2_FC and exp fraction
    # log2_FC in head(50), exp fraction >= 5 % in S_dengue
    # log2_FC in tail(50), exp fraction >= 5 % in dengue
    inter_final = defaultdict(list)

    for dic, cd, k in zip([both_up, both_down], ['SD', 'D'], ['up', 'down']):
        for key in dic.keys(): 
            ct1 = key[0]
            ct2 = key[1]
            for i in range(len(dic[key])): # eg., [ga, ga_ct1_log2FC, ga_ct1_comfrac, gb, gb_ct2_log2FC, gb_ct2_comfrac]
                ga = dic[key][i][0] # in ct1
                gb = dic[key][i][3] # in ct2
                
                ga_SD_exp_frac = log2_fc_ave.loc[ct1, ga]['SD_exp_frac']
                ga_D_exp_frac = log2_fc_ave.loc[ct1, ga]['D_exp_frac']
                gb_SD_exp_frac = log2_fc_ave.loc[ct2, gb]['SD_exp_frac']
                gb_D_exp_frac = log2_fc_ave.loc[ct2, gb]['D_exp_frac']

                #############################################
                ga_in_ct1 = log2_fc_ave.loc[ct1, ga][cd + '_exp_frac']
                gb_in_ct2 = log2_fc_ave.loc[ct2, gb][cd + '_exp_frac']   
                #############################################
                inf = dic[key][i]
                if (ga_in_ct1 >= frac_n) & (gb_in_ct2 >= frac_n):
                    inter_final[k].append([inf[0], inf[3], key[0], key[1], inf[1], inf[2], ga_D_exp_frac, ga_SD_exp_frac, inf[4], inf[5], gb_D_exp_frac, gb_SD_exp_frac])
                else:
                    pass

    inter_f = pd.DataFrame([])
    for key in inter_final.keys():
        inters = pd.DataFrame(inter_final[key], columns=['ga', 'gb', 'cta', 'ctb', 'ga_log2FC', 'ga_comfrac', 'ga_D_exp_frac', 'ga_SD_exp_frac', 'gb_log2FC', 'gb_comfrac', 'gb_D_exp_frac', 'gb_SD_exp_frac'], 
                             index=pd.Index([key]*len(inter_final[key])))
        inter_f =  pd.concat([inter_f, inters])
    
    inter_f['up/down'] = inter_f.index
    inter_f.index = pd.Index(list(range(inter_f.shape[0])))
    
    pairs = []
    for i, cell_type in enumerate(cell_types):
        cts = cell_types[(i+1):]
        for ct in cts:
            pairs.append([cell_type, ct])

    for i in inter_f.index:
        row = inter_f.loc[i]
        if [row['cta'], row['ctb']] not in pairs:
            inter_f.loc[i] = [row[1], row[0], row[3], row[2], row[8], row[9], row[10], row[11], row[4], row[5], row[6], row[7], row[12]]
            
    idx = inter_f[inter_f['cta'] == inter_f['ctb']].index[::2]
    inter_f.drop(idx, inplace=True)
    
    inter_f = inter_f[~inter_f.duplicated()]
        
    return inter_f

def get_mix_inters(fdn, genes, frac_n):
    log2_fc_ave = pd.read_csv(fdn, index_col=['cell_type', 'gene'])

    cell_types = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
    mix = []
    
    fn_int = '/home/yike/phd/dengue/data/interaction_source_file/interactions_DB.tsv'
    interactions = pd.read_csv(fn_int, sep=',')[['gene_name_a', 'gene_name_b']]
    for _, row in interactions.iterrows():
        ga = row['gene_name_a']
        gb = row['gene_name_b']
        if (ga not in genes) | (gb not in genes):
            continue
        for ct1 in cell_types:
            for ct2 in cell_types:

                ct1_head = log2_fc_ave.loc[ct1][log2_fc_ave.loc[ct1]['fold_2_change'] > 0].index.tolist()
                ct2_head = log2_fc_ave.loc[ct2][log2_fc_ave.loc[ct2]['fold_2_change'] > 0].index.tolist()
                ct1_tail = log2_fc_ave.loc[ct1][log2_fc_ave.loc[ct1]['fold_2_change'] < 0].index.tolist()
                ct2_tail = log2_fc_ave.loc[ct2][log2_fc_ave.loc[ct2]['fold_2_change'] < 0].index.tolist()
                #########
                log2_FC = {(gene, ct): log2_fc_ave.loc[ct, gene]['fold_2_change'] for gene in [ga, gb] for ct in [ct1, ct2]}
                com_fra = {(gene, ct): log2_fc_ave.loc[ct, gene]['comp_frac'] for gene in [ga, gb] for ct in [ct1, ct2]}
                exp = {(gene, ct, cd): log2_fc_ave.loc[ct, gene][cd + '_exp_frac']  for gene in [ga, gb] for ct in [ct1, ct2] for cd in ['SD', 'D']}
                ################### interaction == (up, down), up_exp_fra_SD >= 5%, down_exp_fra_D >= 5%
                if (ga in ct1_head) & (exp[(ga, ct1, 'SD')] >= frac_n) & (gb in ct2_tail) & (exp[(gb, ct2, 'D')] >= frac_n):
                    mix.append([ga, gb, ct1, ct2, 
                                log2_FC[(ga, ct1)], com_fra[(ga, ct1)], exp[(ga, ct1, 'D')], exp[(ga, ct1, 'SD')],
                                log2_FC[(gb, ct2)], com_fra[(gb, ct2)], exp[(gb, ct2, 'D')], exp[(gb, ct2, 'SD')],
                                ])
                if (ga in ct2_head) & (exp[(ga, ct2, 'SD')] >= frac_n) & (gb in ct1_tail) & (exp[(gb, ct1, 'D')] >= frac_n):
                    mix.append([ga, gb, ct2, ct1, 
                                log2_FC[(ga, ct2)], com_fra[(ga, ct2)], exp[(ga, ct2, 'D')], exp[(ga, ct2, 'SD')],
                                log2_FC[(gb, ct1)], com_fra[(gb, ct1)], exp[(gb, ct1, 'D')], exp[(gb, ct1, 'SD')],
                                ])
                if (gb in ct1_head) & (exp[(gb, ct1, 'SD')] >= frac_n) & (ga in ct2_tail) & (exp[(ga, ct2, 'D')] >= frac_n):
                    mix.append([gb, ga, ct1, ct2, 
                                log2_FC[(gb, ct1)], com_fra[(gb, ct1)], exp[(gb, ct1, 'D')], exp[(gb, ct1, 'SD')],
                                log2_FC[(ga, ct2)], com_fra[(ga, ct2)], exp[(ga, ct2, 'D')], exp[(ga, ct2, 'SD')],
                                ])
                if (gb in ct2_head) & (exp[(gb, ct2, 'SD')] >= frac_n) & (ga in ct1_tail) & (exp[(ga, ct1, 'D')] >= frac_n):
                    mix.append([gb, ga, ct2, ct1, 
                                log2_FC[(gb, ct2)], com_fra[(gb, ct2)], exp[(gb, ct2, 'D')], exp[(gb, ct2, 'SD')],
                                log2_FC[(ga, ct1)], com_fra[(ga, ct1)], exp[(ga, ct1, 'D')], exp[(ga, ct1, 'SD')],
                                ])
    inter_f = pd.DataFrame(mix, columns=['ga', 'gb', 'cta', 'ctb', 'ga_log2FC', 'ga_comfrac', 'ga_D_exp_frac', 'ga_SD_exp_frac', 'gb_log2FC', 'gb_comfrac', 'gb_D_exp_frac', 'gb_SD_exp_frac'])
    inter_f = inter_f[~inter_f.duplicated()]
    return inter_f
