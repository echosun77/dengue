%config Completer.use_jedi = False

import numpy as np
import pandas as pd

import anndata
import scanpy as sc

data = pd.read_csv('/home/yike/phd/dengue/data/dataset_from_google/human_lung/mouse/mouse_droplet_TMS_UMIs.csv', index_col=0)
dataT = data.T
matrix = dataT.values

obsname = pd.read_csv('/home/yike/phd/dengue/data/dataset_from_google//human_lung/mouse/mouse_droplet_TMS_metadata.csv', index_col=0)
obsname.index.name = None

varname = pd.DataFrame([], index=dataT.columns)
varname.index.name = None

adata_SS2 = anndata.AnnData(X=matrix, obs=obsname, var=varname)

save_file = '/home/yike/phd/dengue/data/dataset_from_google/human_lung/mouse/mouse_droplet.h5ad'
adata_SS2.write(save_file)

adata_droplet = sc.read_h5ad('/home/yike/phd/dengue/data/dataset_from_google/human_lung/mouse/mouse_droplet.h5ad')
adata_facs = sc.read_h5ad('/home/yike/phd/dengue/data/dataset_from_google/human_lung/mouse/mouse_facs.h5ad')
adata = anndata.concat([adata_droplet, adata_facs], join='outer')
adata.write('/home/yike/phd/dengue/data/dataset_from_google/human_lung/mouse/mouse_lung.h5ad')