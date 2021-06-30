from numpy import * # 调用numpy所有函数

adatag_children = adatag[adatag.obs['dataset'] == 'child']

D_IDs = list(adata_children[adata_children.obs['Condition'].isin(['dengue'])].obs['ID'].astype('category').cat.categories)
SD_IDs = list(adata_children[adata_children.obs['Condition'].isin(['S_dengue'])].obs['ID'].astype('category').cat.categories)
IDs = list(adata_children.obs['ID'].astype('category').cat.categories)

cell_types = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
conditions = ['S_dengue', 'dengue', 'DWS', 'Healthy']
from collections import defaultdict
adatag_ch_ct_cd_ID = {}
adatag_ch_ct_cd = {}
for cell_type in cell_types:
    for condition in conditions:
        adatag_ch_ct_cd[(cell_type, condition)] = adatag_children[adatag_children.obs['cell_type'] == cell_type][adatag_children[adatag_children.obs['cell_type'] == cell_type].obs['Condition'] == condition]
        for ID in IDs:
            adatag_ch_ct_cd_ID[(cell_type, condition, ID)] = adatag_ch_ct_cd[(cell_type, condition)][adatag_ch_ct_cd[(cell_type, condition)].obs['ID'] == ID]

gene_exp_ave = {}
for key in adatag_ch_ct_cd.keys():
    condition = key[1]
    exp = []
    for ID in list(adata_children[adata_children.obs['Condition'].isin([condition])].obs['ID'].astype('category').cat.categories):
        exp.append((adatag_ch_ct_cd_ID[(key[0], key[1], ID)].X > 0).toarray().mean(axis=0))
    exp_ave = np.array(exp).mean(axis=0)
    gene_exp_ave[key] = pd.DataFrame(exp_ave, index=adatag_ch_ct_cd[key].var.index, columns=['gene_expre'])

    
exp_fra_ave = pd.DataFrame([])
for key in gene_exp_ave.keys():
    gene_exp_ave[key]['cell_type'] = [key[0]]*(gene_exp_ave[key].shape[0])
    gene_exp_ave[key]['condition'] = [key[1]]*(gene_exp_ave[key].shape[0])
    #exp_fra = pd.DataFrame([], columns=['gene', 'exp_fra_ave'], index=key)
    exp_fra_ave = pd.concat([exp_fra_ave, gene_exp_ave[key]])
    exp_fra_ave['gene'] = exp_fra_ave.index
    exp_fra_ave.to_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/exp_fra_ave.tsv')  