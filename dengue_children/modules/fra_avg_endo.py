def fra_avg_endo(adata):    

    fra = pd.DataFrame([])
    avg = pd.DataFrame([])
    
    cts = adata.obs['Annotation'].astype('category').cat.categories    
    for ct in cts:
        adata_ct = adata[adata.obs['Annotation'] == ct]

        fra_ct = np.asarray((adata_ct.X > 0).mean(axis=0))[0]
        fra_ct = pd.DataFrame(fra_ct, columns=['fra'], index=adata_ct.var_names)
        fra_ct['cell_subtype'] = ct
        fra = pd.concat([fra, fra_ct])
        fra['gene'] = fra.index.tolist()

        avg_ct = np.asarray(adata_ct.X.mean(axis=0))[0]
        avg_ct = pd.DataFrame(avg_ct, columns=['avg'], index=adata_ct.var_names)
        avg_ct['cell_subtype'] = ct
        avg = pd.concat([avg, avg_ct])
        avg['gene'] = avg.index.tolist()
            
    return {'fra': fra, 'avg': avg}
