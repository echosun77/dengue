import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict

def get_inters(gene_dic1, gene_dic2):
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
    best_idx = it_n.index[best_idx]
    #best_idx = ['T_cells', 'B_cells', 'NK_cells', 'Plasmablasts', 'Monocytes', 'pDCs', 'cDCs', ]

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
    ax.text(9.05, 3.8, '          Number of \n %sregulated interactions'%trend, verticalalignment='center', rotation=90)
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
    #best_idx = ['T_cells', 'B_cells', 'NK_cells', 'Plasmablasts', 'Monocytes', 'pDCs', 'cDCs', ]
    it_n = it_n.iloc[best_idx]

    lkg_col = linkage(pdist(it_n.T.values), optimal_ordering=True)
    best_col = leaves_list(lkg_col)
    best_col = it_n.columns[best_col]
    it_n = it_n[best_col]

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

    ax.set_xticklabels(xlabels, c='r')
    ax.set_yticklabels(xlabels, c='b')

    return {'figure': fig, 'ax': ax}