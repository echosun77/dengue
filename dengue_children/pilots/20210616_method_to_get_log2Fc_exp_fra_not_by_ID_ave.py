from log2_FC_functions import log2_FC_all

gbs = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']

adata_gb = {}
log2_fc_all = {}
log2_fc_ave = pd.DataFrame([])

for gb in gbs:
    adata_gb[gb] = adata[adata.obs['cell_type'] == gb]
    if gb == 'cDCs':
        adata_gb[gb] = adata_gb[gb][ ~ adata_gb[gb].obs['ID'].isin(['1_140_01', '5_193_01'])]
    
    log2_fc_all[gb] = log2_FC_all(genes, adata_gb[gb], 'S_dengue', 'dengue', 'child')[1]
    log2_fc_all[gb]['cell_type']=[gb]*log2_fc_all[gb].shape[0]
    log2_fc_ave = pd.concat([log2_fc_ave, log2_fc_all[gb]])

log2_fc_ave['gene'] = log2_fc_ave.index.values
log2_fc_ave = log2_fc_ave.set_index('cell_type')

##################################################################
adatag_children = adatag[adatag.obs['dataset'] == 'child']
cell_types = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
conditions = ['S_dengue', 'dengue', 'DWS', 'Healthy']

adatag_ch_ct_cd = {}
for cell_type in cell_types:
    for condition in conditions:
        adatag_ch_ct_cd[(cell_type, condition)] = adatag_children[adatag_children.obs['cell_type'] == cell_type][adatag_children[adatag_children.obs['cell_type'] == cell_type].obs['Condition'] == condition] 

exp_fra = {}
for key in adatag_ch_ct_cd.keys():
    fra = (adatag_ch_ct_cd[key].X > 0).toarray().mean(axis=0)
    exp_fra[key] = pd.DataFrame(fra, index=adatag_ch_ct_cd[key].var.index, columns=[['exp_fra']])

ave_exp = {}
for key in adatag_ch_ct_cd.keys():
    exp = adatag_ch_ct_cd[key].X.toarray().mean(axis=0)
    ave_exp[key] = pd.DataFrame(exp, index=adatag_ch_ct_cd[key].var.index, columns=[['ave_exp']])

fra = pd.DataFrame([])
exp = pd.DataFrame([])
for key in exp_fra.keys():
    exp_fra[key]['cell_type'] = [key[0]]*(exp_fra[key].shape[0])
    exp_fra[key]['condition'] = [key[1]]*(exp_fra[key].shape[0])
    fra = pd.concat([fra, exp_fra[key]])
    fra['gene'] = fra.index
    fra.to_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/not_ID_ave/fra.tsv')
    
    ave_exp[key]['cell_type'] = [key[0]]*(ave_exp[key].shape[0])
    ave_exp[key]['condition'] = [key[1]]*(ave_exp[key].shape[0])
    exp = pd.concat([exp, ave_exp[key]])
    exp['gene'] = exp.index
    exp.to_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/not_ID_ave/exp.tsv')

########################################################################
log2_fc = log2_fc_ave.set_index([log2_fc_ave.index, 'gene'])
pvalue = pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/inters/inter_genes_pvalue.tsv', index_col=['cell_type', 'Unnamed: 0'])

pvalue = pvalue.loc[log2_fc.index]

log2_fc['statistic'] = pvalue['statistic']
log2_fc['pvalue'] = pvalue['pvalue']
######################################################################
fra = pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/not_ID_ave/fra.tsv', index_col=['cell_type', 'gene', 'condition'])
exp = pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/not_ID_ave/exp.tsv', index_col=['cell_type', 'gene', 'condition'])

fra_SD = fra.loc[[(idx[0], idx[1], 'S_dengue') for idx in log2_fc.index]]
fra_D = fra.loc[[(idx[0], idx[1], 'dengue') for idx in log2_fc.index]]

exp_SD = exp.loc[[(idx[0], idx[1], 'S_dengue') for idx in log2_fc.index]]
exp_D = exp.loc[[(idx[0], idx[1], 'dengue') for idx in log2_fc.index]]

log2_fc['SD_exp_frac'] = fra_SD['exp_fra'].tolist()
log2_fc['D_exp_frac'] = fra_D['exp_fra'].tolist()
log2_fc['SD_ave_exp'] = exp_SD['ave_exp'].tolist()
log2_fc['D_ave_exp'] = exp_D['ave_exp'].tolist()

log2_fc.to_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/not_ID_ave/log2_fc.tsv')
