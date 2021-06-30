import anndataks

cell_types = ['B_cells', 'Monocytes', 'NK_cells', 'Plasmablasts', 'T_cells', 'cDCs', 'pDCs']
conditions = ['S_dengue', 'dengue']

#sc.pp.log1p(adatag)

adata_kids = adatag[adatag.obs['dataset'] == 'child']
results = {}
for cell_type in cell_types:
    adata_ct = adata_kids[adata_kids.obs['cell_type'] == cell_type]
    if cell_type == 'cDCs':
        adata_ct = adata_ct[~adata_ct.obs['ID'].isin(['1_140_01', '5_193_01'])]
    
    adata_SD = adata_ct[adata_ct.obs['Condition'] == 'S_dengue']
    adata_D = adata_ct[adata_ct.obs['Condition'] == 'dengue']
    results[cell_type] = anndataks.compare(adata_D, adata_SD) # log1p=False # log2_fold_change: adata_Sd vs adata_D

res = pd.DataFrame([])
for cell_type in cell_types:
    results[cell_type]['cell_type'] = [cell_type] * results[cell_type].shape[0]
    res = pd.concat([res, results[cell_type]] )
    
res.to_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/inters/inter_genes_pvalue.tsv')