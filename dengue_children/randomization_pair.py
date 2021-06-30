def randomization_pair(genes):
    
    log2fc = defaultdict(list)
    r = defaultdict(list)
    pvalue = {}

    for inter in genes:
        ga = list(inter.keys())[0]
        cta = inter[ga]
        gb = list(inter.keys())[1]
        ctb = inter[gb]

        adata_g_ct = {}
        adata_g = {}
        adata_kids = adata_children[adata_children.obs['Condition'].isin(['S_dengue', 'dengue'])]
        for gene, ct in zip([ga, gb], [cta, ctb]):
            if ct == 'cDCs':
                adata_kids = adata_kids[~adata_kids.obs['ID'].isin(['1_140_01', '5_193_01'])]
            adata_g_ct[gene]= adata_kids[adata_kids.obs['cell_type'] == ct][:, gene]
            for cd in ['S_dengue', 'dengue']:
                adata_g[(gene, cd)] = adata_g_ct[gene][adata_g_ct[gene].obs['Condition'] == cd]
        IDs = {key: adata_g[key].obs['ID'].unique().tolist() for key in adata_g.keys()}
        avg = {(key[0], key[1], ID): adata_g[key][adata_g[key].obs['ID'] == ID].X.toarray().mean() for key in IDs.keys() for ID in IDs[key]}

        pair = defaultdict(list)
        for gene in [ga, gb]:
            for ID_sd in IDs[gene, 'S_dengue']:
                for ID_d in IDs[gene, 'dengue']:
                    pair[gene].append(np.log2(avg[gene, 'S_dengue', ID_sd] + 0.1) - np.log2(avg[gene, 'dengue', ID_d] + 0.1))

        lfc = {gene: np.median(pair[gene]) for gene in [ga, gb]}
        log2fc[(ga, cta, gb, ctb)].append([lfc[ga], lfc[gb]])
        r0 = (float(lfc[ga])**2 + float(lfc[gb])**2)**0.5
        r[(ga, cta, gb, ctb)].append(r0)

        adata_i = {}
        for i in range(1000):
            adata_g_ct_i = adata_g_ct
            raw = {gene: adata_g_ct_i[gene].obs['Condition'].tolist() for gene in [ga, gb]}
            for gene in [ga, gb]:
                random.shuffle(raw[gene])
                adata_g_ct_i[gene].obs['Condition'] = raw[gene]

                for cd in ['S_dengue', 'dengue']:
                    adata_i[(gene, cd)] = adata_g_ct_i[gene][adata_g_ct_i[gene].obs['Condition'] == cd]
            avg_i = {(key[0], key[1], ID): adata_i[key][adata_i[key].obs['ID'] == ID].X.toarray().mean() for key in IDs.keys() for ID in IDs[key]}

            pair_i = defaultdict(list)
            for gene in [ga, gb]:
                for ID_sd in IDs[gene, 'S_dengue']:
                    for ID_d in IDs[gene, 'dengue']:
                        pair_i[gene].append(np.log2(avg_i[gene, 'S_dengue', ID_sd] + 0.1) - np.log2(avg_i[gene, 'dengue', ID_d] + 0.1))

            lfc_i = {gene:  np.median(pair_i[gene]) for gene in [ga, gb]}
            log2fc[(ga, cta, gb, ctb)].append([lfc_i[ga], lfc_i[gb]])
            ri = (float(lfc_i[ga])**2 + float(lfc_i[gb])**2)**0.5
            r[(ga, cta, gb, ctb)].append(ri)
            pvalue[(ga, cta, gb, ctb)] = 0
            if ri >= r0:
                pvalue[(ga, cta, gb, ctb)] += 1
            pvalue[(ga, cta, gb, ctb)] = pvalue[(ga, cta, gb, ctb)]/1000

            res = pd.DataFrame([])
    for key in log2fc.keys():
        log2fc[key] = pd.DataFrame(log2fc[key], columns = ['log2fc_ga', 'log2fc_gb'])
        log2fc[key]['r'] = r[key]
        log2fc[key]['pvalue'] = pvalue[key]
        for i, s in enumerate(['ga', 'cta', 'gb', 'ctb']):
            log2fc[key][s] = key[i]
        res = pd.concat([res, log2fc[key]])
    res = res.set_index(['ga', 'cta', 'gb', 'ctb'])
    return res