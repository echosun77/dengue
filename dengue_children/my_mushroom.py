def dotplot_inters(inter_type, partner1, partner2, height):
    '''
    the plot showing gene expression and fraction of interactions, i.e., up-up, down-down, up-down
    Args:
        inter_type: the type of interactions, including three types, i.e., 'upregulated', 'downregulated', 'mixregulated'
        partner1/partner2: a dict with genes as keys and cell types as values
            e.g., partner1 = {'CXCL10': 'Monocytes', 'DPP4': ['B_cells', 'NK_cells']}
        height: the figure height, inceasing with the increase of interaction number
    returns: {'figure': fig, 'ax': axs}
    '''
    from matplotlib.patches import Wedge
    import matplotlib as mpl
    import math

    cell_types = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
    
    l_genes = list(set(partner1.keys()))
    r_genes = list(set(partner2.keys()))
    all_genes = list(set(l_genes + r_genes))

    exp_fra = {(ct, cd, gene): gene_exp_ave.loc[(ct, cd)][gene_exp_ave.loc[(ct, cd)]['gene'] == gene]['gene_expre'][0] for ct in cell_types for cd in ['S_dengue', 'dengue'] for gene in all_genes}
    ave_exp = {(ct, cd, gene): np.log2(2 + ave_exp_ave.loc[(ct, cd)][ave_exp_ave.loc[(ct, cd)]['gene'] == gene]['ave_exp'][0]) for ct in cell_types for cd in ['S_dengue', 'dengue'] for gene in all_genes}

    rs = {key: pow(exp_fra[key]/5, 0.5) for key in exp_fra.keys()}
    colors = {key: plt.cm.get_cmap('viridis')((ave_exp[key] - min(ave_exp.values()))/max(ave_exp.values())) for key in ave_exp.keys()}

    fig,axs = plt.subplots(1, 2, figsize=[8, height], dpi=300)
    plt.subplots_adjust(wspace=1)
    
    if inter_type == 'upregulated':
        edge_c = [['r', 'none'], ['r', 'none']]
    elif inter_type == 'downregulated':
        edge_c = [['none' , 'r'], ['none' , 'r']]
    elif inter_type == 'mixregulated':
        edge_c = [['r', 'none'], ['none', 'r']]
            
    for x, ct in enumerate(cell_types):
        for genes, dic, i, tt, p_e in zip([l_genes, r_genes], [partner1, partner2], [0, 1], ['Partner 1', 'Partner 2'], edge_c):
            for y, gene in enumerate(genes):
                if ct in dic[gene]: 
                    e1 = Wedge((x,y+0.05), rs[(ct, 'S_dengue', gene)], 0, 180, facecolor=colors[(ct, 'S_dengue', gene)], edgecolor=p_e[0])
                    e2 = Wedge((x,y-0.05), rs[(ct, 'dengue', gene)], 180, 360, facecolor=colors[(ct, 'dengue', gene)], edgecolor=p_e[1])
                else:
                    e1 = Wedge((x,y+0.05), rs[(ct, 'S_dengue', gene)], 0, 180, facecolor=colors[(ct, 'S_dengue', gene)])
                    e2 = Wedge((x,y-0.05), rs[(ct, 'dengue', gene)], 180, 360, facecolor=colors[(ct, 'dengue', gene)])
                axs[i].add_patch(e1)
                axs[i].add_patch(e2)

                axs[i].set_xlim([-0.7, len(cell_types)-0.3])
                axs[i].set_ylim([-0.7, len(genes)-0.3])

                axs[i].set_xticks([i for i in range(len(cell_types))])
                axs[i].set_yticks([i for i in range(len(genes))])

                axs[i].axhline(y, c='w', lw=0.1, zorder=-1)

            axs[i].set_xticklabels(['B cells', 'Monocytes', 'NK cells', 'Plasmablasts', 'T cells', 'cDCs', 'pDCs'], rotation=90)
            axs[i].set_yticklabels(genes)
            axs[i].set_title(tt)
            axs[i].set_aspect(1)
    axs[0].yaxis.tick_right()
    axs[0].yaxis.set_label_position('right')


    norm = mpl.colors.Normalize(vmin=0, vmax=10) # (max(ave_exp.values())+1)
    cmap = plt.cm.get_cmap('viridis')
    position = fig.add_axes([0.91, 0.1, 0.02, 0.7]) # [left, bottom, width, height]
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=position, ax=axs[1], fraction=0.01, pad=0.01, orientation='vertical')
    axs[1].text(9.3, 3, 'Gene exp (log2[2+cpm])', verticalalignment='center', rotation=90)
    
    return {'figure': fig, 'ax': axs}

#################################################### Fabio's mushroom
def mushrooms(genes):

    from matplotlib.patches import Wedge
    import matplotlib as mpl
    import math
    cell_types = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
    conditions = ['S_dengue', 'dengue']
    cmap = plt.cm.get_cmap('viridis')
    vmin, vmax = -1, 3
    threshold = 0.1
    frac = pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/exp_fra_ave.tsv', index_col=['cell_type', 'condition', 'gene'], squeeze=True)
    avg = pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/ave_exp_ave.tsv', index_col=['cell_type', 'condition', 'gene'], squeeze=True)

    fig, axs = plt.subplots(
        2, 1, sharex=True, sharey=False,
        figsize=((1 + 0.8 * len(cell_types)) * 0.5, (1 + len(genes[0]) + len(genes[1])) * 0.5),
        gridspec_kw={'height_ratios': [len(genes[0]), len(genes[1])]},
        dpi=300
    )
    datap = []
    for iax, (genesi, ax) in enumerate(zip(genes, axs)):
        for i, cst in enumerate(cell_types):
            for j, gene in enumerate(genesi):
                avgs = []
                for k, cond in enumerate(conditions):
                    fr = frac.loc[(cst, cond, gene)]
                    av = np.log10(avg.loc[(cst, cond, gene)] + 0.1)
                    avgs.append(av)
                                    
                    r = 0.5 * fr**0.3
                    color = cmap((min(vmax, av) - vmin) / (vmax - vmin))
                    theta0, theta1 = 180 * (k > 0), 180 + 180 * (k > 0)
                    datap.append({
                        'r': r,
                        'facecolor': color,
                        'center': (i, j),
                        'theta': (theta0, theta1),
                        'ax': iax,
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
                
        ax.set_yticks(np.arange(len(genesi)))
        ax.set_yticklabels(genesi)
        ax.set_ylim(-0.6, len(genesi) - 0.4)        
    ax.set_xticks(np.arange(len(cell_types)))
    ax.set_xticklabels([x.replace('_', ' ') for x in cell_types], rotation=90)
    ax.set_xlim(-0.6, len(cell_types) - 0.4)
    
    for datum in datap:
        ax = axs[datum['ax']]
        r = datum['r']
        color = datum['facecolor']
        center = datum['center']
        theta0, theta1 = datum['theta']
        ec = datum['edgecolor']

        h = Wedge(
            center, r, theta0, theta1, edgecolor=ec, facecolor=color,
        )
        ax.add_artist(h)
        ax.set_aspect(1)
    
    e = Wedge(
            (7, 0), 0.5 * 0.05**0.3, 0, 180, facecolor='gray',
        )
    axs[0].add_artist(e)
    
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax) # (max(ave_exp.values())+1)
    cmap = plt.cm.get_cmap('viridis')
    position = fig.add_axes([1, 0.2, 0.02, 0.25]) # [left, bottom, width, height]
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=position, ax=axs[1], label='Gene exp (log10[cpm+0.1])')

    fig.tight_layout()

    return {'fig': fig, 'ax': ax}


def s_mushrooms(genes):
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
    vmin, vmax = -1, 3
    threshold = 0.1
    frac = pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/20210625_figure_4_code/fra.tsv', index_col=['cell_type', 'condition', 'gene'], squeeze=True)
    avg = pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/20210625_figure_4_code/exp.tsv', index_col=['cell_type', 'condition', 'gene'], squeeze=True)

    yl = sum([len(list(itertools.chain.from_iterable(genesi.values()))) for genesi in genes])
    fig = plt.figure(figsize=((1 + 0.8 * 2) * 0.6, (1 + yl)* 0.6), dpi=300)

    grid = plt.GridSpec(yl , 2, wspace=0.1, hspace=0.1)

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
            avgs = []
            for cst in csts:
                yticklabels.append(gene + ' in\n' + cst.replace('_', ' '))
                for k, cond in enumerate(conditions):
                    fr = frac.loc[(cst, cond, gene)]
                    av = np.log10(avg.loc[(cst, cond, gene)] + 0.1)
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
    position = fig.add_axes([0.7, 0.02*yl, 0.05, 2/yl])
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=position, ax=axs[-1], label='Gene exp \n(log10[cpm+0.1])')

    fig.tight_layout()
    return {'fig': fig, 'ax': axs}