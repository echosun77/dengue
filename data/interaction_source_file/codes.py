fn_int = '/home/yike/phd/dengue/data/interaction_source_file/interaction_unpacked_mouse_original.tsv'
interactions = pd.read_csv(fn_int, sep='\t')
#[['gene_name_a', 'gene_name_b']]
ga, gb = interactions['gene_name_a'], interactions['gene_name_b']

cellphoneDB = pd.read_csv('/home/yike/phd/dengue/data/interaction_source_file/cellphoneDB.csv', sep=',')
interactions_DB = cellphoneDB.iloc[466:]

human_gene_name = pd.read_excel('/home/yike/phd/dengue/data/interaction_source_file/Human_gene_name.xlsx')
mouse_gene_name = pd.read_excel('/home/yike/phd/dengue/data/interaction_source_file/Mouse_gene_name.xlsx')

gene_name_a = []
gene_name_b = []

for i in interactions_DB['protein_name_a'].tolist():
    if i in human_gene_name['Entry_name'].tolist():
        gene = human_gene_name[human_gene_name['Entry_name'] == i]['Gene_names_primary'].tolist()[0]
        gene_name_a.append(gene)
    else:
        gene_name_a.append('-')

for i in interactions_DB['protein_name_b'].tolist():
    if i in human_gene_name['Entry_name'].tolist():
        gene = human_gene_name[human_gene_name['Entry_name'] == i]['Gene_names_primary'].tolist()[0]
        gene_name_b.append(gene)
    else:
        gene_name_b.append('-')

gene_name_mouse_a = [str.capitalize(i) for i in gene_name_a]
gene_name_mouse_b = [str.capitalize(i) for i in gene_name_b]

interactions_DB['gene_name_a'] = gene_name_a
interactions_DB['gene_name_b'] = gene_name_b

interactions_DB['gene_name_mouse_a'] = gene_name_mouse_a
interactions_DB['gene_name_mouse_b'] = gene_name_mouse_b

interactions_DB['index'] = interactions_DB.index.tolist()
interactions_DB = interactions_DB[(interactions_DB['gene_name_a'] != '-') & (interactions_DB['gene_name_b'] != '-')]

inters = interactions.iloc[:732]
inters.index = pd.Index(inters['index'].tolist())

for n in interactions_DB['index'].tolist():
    if n not in inters['index'].tolist():
        inters.loc[n] = interactions_DB.loc[n]

sup_1 = pd.read_excel('/home/yike/phd/dengue/data/interaction_source_file/sup_1.xlsx').fillna(0)
sup_1 = sup_1[(sup_1['gene_name_a'] != 0) & (sup_1['gene_name_b'] != 0)]

ID_ls = []
for ID in sup_1['id_cp_interaction'].tolist():
    if ID not in inters['id_cp_interaction'].tolist():
        ID_ls.append(ID)

sup1_ap = sup_1[sup_1['id_cp_interaction'].isin(ID_ls)]
inters = pd.concat([inters, sup1_ap], join='outer')


sup_2 = pd.read_excel('/home/yike/phd/dengue/data/interaction_source_file/sup_2.xlsx').fillna(0)
#sup_2 = sup_2[(sup_1['gene_name_a'] != 0) & (sup_1['gene_name_b'] != 0)]
sup_2

partner_a = [i.split(':')[-1] for i in sup_2['partner_a'].tolist()]
partner_b = [i.split(':')[-1] for i in sup_2['partner_b'].tolist()]

sup_2['partner_a'] = partner_a
sup_2['partner_b'] = partner_b

gene_name_a = []
gene_name_b = []

for i in sup_2['partner_a'].tolist():
    if i in human_gene_name['Entry'].tolist():
        gene = human_gene_name[human_gene_name['Entry'] == i]['Gene_names_primary'].tolist()[0]
        gene_name_a.append(gene)
    else:
        gene_name_a.append('-')

for i in sup_2['partner_b'].tolist():
    if i in human_gene_name['Entry'].tolist():
        gene = human_gene_name[human_gene_name['Entry'] == i]['Gene_names_primary'].tolist()[0]
        gene_name_b.append(gene)
    else:
        gene_name_b.append('-')

gene_name_mouse_a = [str.capitalize(i) for i in gene_name_a]
gene_name_mouse_b = [str.capitalize(i) for i in gene_name_b]

sup_2['gene_name_a'] = gene_name_a
sup_2['gene_name_b'] = gene_name_b
sup_2['gene_name_mouse_a'] = gene_name_mouse_a
sup_2['gene_name_mouse_b'] = gene_name_mouse_b

sup_2.drop(['interacting_pair', 'gene_a', 'gene_b'], axis=1, inplace=True)
sup_2 = sup_2[(sup_2['gene_name_a'] != '-') & (sup_2['gene_name_b'] != '-')]

ID_ls = []
for ID in sup_2['id_cp_interaction'].tolist():
    if ID not in inters['id_cp_interaction'].tolist():
        ID_ls.append(ID)

sup2_ap = sup_2[sup_2['id_cp_interaction'].isin(ID_ls)]

inters = pd.concat([inters, sup2_ap], join='outer')

inter_ls = []
for _,row in inters[['gene_name_a', 'gene_name_b']].iterrows():
    ga = row['gene_name_a']
    gb = row['gene_name_b']
    inter_ls.append([ga, gb])

interactions = interactions.iloc[732:]
ga, gb, idx = interactions['gene_name_a'].tolist(), interactions['gene_name_b'].tolist(), interactions.index.tolist()
idx_ap = []
for gene_a, gene_b, i in zip(ga, gb, idx):
    if ([gene_a, gene_b] not in inter_ls) & ([gene_b, gene_a] not in inter_ls):
        idx_ap.append(i)
inter_ap = interactions.loc[idx_ap]
inters = pd.concat([inters, inter_ap], join='outer')
inters.to_csv('/home/yike/phd/dengue/data/interaction_source_file/interactions_DB.tsv')

