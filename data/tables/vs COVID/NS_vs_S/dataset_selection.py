### Italy data

print('load COVID Italy dataset')
COVID_I = sc.read_h5ad('/home/yike/phd/dengue/data/dataset_from_google/COVID PBMC Italy/50__subtype__ID/adata_filter.h5ad')
print("select 'timepoint' == 'T0'")
COVID_I = COVID_I[COVID_I.obs['timepoint'] == 'T0']
print("select 'days_since_hospitalized' <= 10")
COVID_I = COVID_I[COVID_I.obs['days_since_hospitalized'].isin(['1.0', '3.0', '0.0', '6.0', '9.0', '2.0', '7.0'])]
sc.pp.normalize_total(COVID_I, target_sum=1e6) 


### Germany data
print('load COVID Germany dataset')
COVID_G = sc.read_h5ad('/home/yike/phd/dengue/data/dataset_from_google/COVID PBMC Germany/filtered_cellXgene.h5ad')
print("select 'disease_stage' != 'late'")
COVID_G = COVID_G[COVID_G.obs['disease_stage'] != 'late'] # select control, early