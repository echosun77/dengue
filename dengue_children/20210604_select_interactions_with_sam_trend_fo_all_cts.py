log2_fc_ave = pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/log2_fc_ave_kids.tsv', index_col=['cell_type', 'gene'])
log2_fc =pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/log2_fc_ave_kids.tsv', index_col=['gene','cell_type', ])

genes = log2_fc_ave.loc['B_cells'].index.tolist()
cell_types = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
same_trend = {}
same_trend['up'] = defaultdict(list)
same_trend['down'] = defaultdict(list)
for gene in genes:
    if (log2_fc.loc[gene]['fold_2_change'] > 0).sum() >= 7:
        ct_ls = log2_fc.loc[gene][log2_fc.loc[gene]['fold_2_change'] > 0].index.tolist()
        for ct in ct_ls:
            same_trend['up'][ct].append(gene)
    elif (log2_fc.loc[gene]['fold_2_change'] < 0).sum() >= 7:
        ct_ls = log2_fc.loc[gene][log2_fc.loc[gene]['fold_2_change'] < 0].index.tolist()
        for ct in ct_ls:
            same_trend['down'][ct].append(gene)
####################################################
####################################################
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
                
    up_inters = pd.DataFrame([])
    up_inters['cta'] = cta
    up_inters['ctb'] = ctb
    up_inters['ga'] = gene_a
    up_inters['gb'] = gene_b

    idx = up_inters[up_inters['cta'] == up_inters['ctb']].index[::2]
    up_inters.drop(idx, inplace=True)
    return up_inters
####################################################    

##### same trend for 7 cell types
#{'up': ['PLD2','COPA','GAS6','SIRPB1','IL1R2','ADM','CD1D','NPW','IL17RA','SORT1','IFNGR1','ENTPD1','ARF1','SLC7A1'],
#'down': ['HLA-E','CD24','TNFRSF14','CD48','COL6A3','B2M','CD40','LTB','CD200','ICAM2','KLRG1','TYROBP','BST2','TNFSF10','XCL2']}
####################################################
same_trend_up = {ct:same_trend['up'] for ct in cell_types}
same_trend_down = {ct:same_trend['down'] for ct in cell_types}
#up & up
a = get_inters(cell_types, same_trend_up, cell_types, same_trend_up)
a[~a[['ga', 'gb']].duplicated()]
#######['SORT1', 'COPA'], ['ARF1', 'PLD2']

#down & down
get_inters(cell_types, same_trend_down, cell_types, same_trend_down)
####### None

# up & down
get_inters(cell_types, same_trend_up, cell_types, same_trend_down)
####### None

##### same trend for 6 cell types

up_up.to_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/inters/same_trend_inters/up_up_6cts.tsv')
down_down.to_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/inters/same_trend_inters/down_down_6cts.tsv')
up_down.to_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/inters/same_trend_inters/up_down_6cts.tsv')

# up_up
a = ['SLC7A1', 'SORT1', 'ARF1', 'CADM1', 'IL10', 'AMH', 'ITGB7', 'ITGA2']
b = ['CSF1', 'COPA', 'PLD2', 'CADM1', 'IL10RA', 'BMPR1A', 'CDH1', 'CDH1']

dotplot_inters(a, b, 5, ['Partner 1', 'Partner 2'])

#down_down
a = ['TNFRSF14',
 'CD244',
 'TNFRSF14',
 'TNFRSF14',
 'KLRK1',
 'KIR3DL1',
 'BST2',
 'DAG1',
 'CLEC2D',
 'HLA-E']

 gb = ['TNFSF14',
 'CD48',
 'BTLA',
 'CD160',
 'HLA-E',
 'HLA-F',
 'LILRA4',
 'LGALS9',
 'KLRB1',
 'KLRD1']

 