fn_int = '/home/yike/phd/dengue/data/interaction_source_file/inters_YK_20220324.tsv'
interactions = pd.read_csv(fn_int, sep='\t')[['genesymbol_intercell_source', 'genesymbol_intercell_target']]
genes = np.unique(interactions)
genes = [gene for gene in genes if gene in adata_kid.var_names]

up_DEGs = {}
down_DEGs = {}

for ct in cell_types:
    data_ct = ct_pair.loc[ct].loc[genes] ## interacting genes
    up_ct = ((data_ct['med_pair'] >= 1) & (data_ct['fra_pair'] >= 39/56) & (data_ct['S_fra'] >= 0.02))
    up_DEGs[ct] = up_ct
    down_ct = ((data_ct['med_pair'] <= -1) & (data_ct['neg_fra_pair'] >= 39/56) & (data_ct['NS_fra'] >= 0.02))
    down_DEGs[ct] = down_ct

res = []

data_ct = ct_pair
for _, row in interactions.iterrows():
    ga = row['genesymbol_intercell_source']
    gb = row['genesymbol_intercell_target']
    if (ga not in genes) | (gb not in genes):
        continue

    for cst in cts:
        cst_med = {gene: data_ct.loc[cst, gene]['med_pair'] for gene in [ga, gb]} 
        cst_fra_pair = {gene: data_ct.loc[cst, gene]['fra_pair'] for gene in [ga, gb]} 
        cst_neg_fra_pair = {gene: data_ct.loc[cst, gene]['neg_fra_pair'] for gene in [ga, gb]}
        cst_fra = {(gene, cd): data_ct.loc[cst, gene][cd+'_fra'] for gene in [ga, gb] for cd in ['S', 'NS']}
        for ct in cts:
            ct_med = {gene: data_ct.loc[ct, gene]['med_pair'] for gene in [ga, gb]}
            ct_fra_pair = {gene: data_ct.loc[ct, gene]['fra_pair'] for gene in [ga, gb]}
            ct_neg_fra_pair = {gene: data_ct.loc[cst, gene]['neg_fra_pair'] for gene in [ga, gb]}
            ct_fra = {(gene, cd): data_ct.loc[ct, gene][cd+'_fra'] for gene in [ga, gb] for cd in ['S', 'NS']}

            if (ga in up_DEGs[cst].index) & (gb in up_DEGs[ct].index):

                res.append([
                        ga, cst, cst_med[ga], cst_fra_pair[ga], cst_neg_fra_pair[ga], cst_fra[ga, 'S'], cst_fra[ga, 'NS'],
                        gb, ct, ct_med[gb], ct_fra_pair[gb], ct_neg_fra_pair[gb], ct_fra[gb, 'S'], ct_fra[gb, 'NS'], 'up'
                    ])

            if (ga in down_DEGs[cst].index) & (gb in down_DEGs[ct].index):
                res.append([
                        ga, cst, cst_med[ga], cst_fra_pair[ga], cst_neg_fra_pair[ga], cst_fra[ga, 'S'], cst_fra[ga, 'NS'],
                        gb, ct, ct_med[gb], ct_fra_pair[gb], ct_neg_fra_pair[gb], ct_fra[gb, 'S'], ct_fra[gb, 'NS'], 'down'
                    ])

            if (ga in up_DEGs[cst].index) & (gb in down_DEGs[ct].index):
                res.append([
                        ga, cst, cst_med[ga], cst_fra_pair[ga], cst_neg_fra_pair[ga], cst_fra[ga, 'S'], cst_fra[ga, 'NS'],
                        gb, ct, ct_med[gb], ct_fra_pair[gb], ct_neg_fra_pair[gb], ct_fra[gb, 'S'], ct_fra[gb, 'NS'], 'mix'
                    ])

            if (ga in down_DEGs[cst].index) & (gb in up_DEGs[ct].index):
                res.append([
                        gb, ct, ct_med[gb], ct_fra_pair[gb], ct_neg_fra_pair[gb], ct_fra[gb, 'S'], ct_fra[gb, 'NS'], 
                        ga, cst, cst_med[ga], cst_fra_pair[ga], cst_neg_fra_pair[ga], cst_fra[ga, 'S'], cst_fra[ga, 'NS'], 'mix'
                    ])

res = pd.DataFrame(res, columns=[
'ga', 'csta', 'ga_med_pair', 'ga_fra_pair','ga_neg_fra_pair','ga_SD_fra', 'ga_D_fra', 
'gb', 'cstb', 'gb_med_pair', 'gb_fra_pair','gb_neg_fra_pair','gb_SD_fra', 'gb_D_fra',
'inter_type'])
            
res_inters.to_csv(save_tabels + 'inters_med_pair1_pair39_54_exp002.tsv', sep='\t', index=False)




