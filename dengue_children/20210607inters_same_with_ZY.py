inters = pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/inters/inters_im_75_002ave.tsv', index_col=0)

genes_zy = pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/inters/all_deg_febe_list_Zhiyuan.tsv', sep='\t', index_col=0)

inters_zy = pd.DataFrame([])
for _, row1 in genes_zy.iterrows():
    ct = row1['celltype']
    gene = row1['gene_name']
    for i,row2 in inters.iterrows():
        ga = row2['ga']
        gb = row2['gb']
        cta = row2['cta']
        ctb = row2['ctb']
        if ((ga == gene) & (cta == ct)) | ((gb == gene) & (ctb == ct)):
            row2 = row2.append(row1)
            #inters_zy.iloc[i] = row2.tolist()
            #inters_zy.columns=row2.index.tolist() 
            it = pd.DataFrame(row2.tolist(), index=row2.index.tolist())
            inters_zy = pd.concat([inters_zy, it], axis=1)

inters_zy.T.to_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/inters/inters_zy_75_002ave.tsv')
inters_zy.T


inter_mix = pd.read_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/inters/inters_mix_75_002ave.tsv', index_col=0)
inter_mix_zy = pd.DataFrame([])
for _, row1 in genes_zy.iterrows():
    ct = row1['celltype']
    gene = row1['gene_name']
    for i,row2 in inter_mix.iterrows():
        ga = row2['ga']
        gb = row2['gb']
        cta = row2['cta']
        ctb = row2['ctb']
        if ((ga == gene) & (cta == ct)) | ((gb == gene) & (ctb == ct)):
            row2 = row2.append(row1)
            #inters_zy.iloc[i] = row2.tolist()
            #inters_zy.columns=row2.index.tolist() 
            it = pd.DataFrame(row2.tolist(), index=row2.index.tolist())
            inter_mix_zy = pd.concat([inter_mix_zy, it], axis=1)
inter_mix_zy.T.to_csv('/home/yike/phd/dengue/data/excels/log2_fc/S_dengue_vs_dengue/inters/inter_mix_zy_75_002ave.tsv')