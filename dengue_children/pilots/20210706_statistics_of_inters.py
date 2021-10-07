
import h5py
import numpy as np
import pandas as pd
import anndata
import scipy.sparse as sp
import scanpy as sc
import random
import scipy.sparse
from scipy.stats import ks_2samp

adata_T = adata_children[adata_children.obs['cell_type'] == 'T_cells']
adata_S = adata_T[(adata_T.obs['Condition'] == 'S_dengue')] 

frac = (adata_S[:, ['CTLA4', 'CD80']].X > 0).toarray()
from collections import Counter
c = Counter((frac * [1, 10]).sum(axis=1))

a = [[c[0], c[1]],
     [c[11], c[10]]]

from scipy.stats import fisher_exact
fisher_exact(a)
# there is a correction between the two genes, CTLA4, CD86

np.median(adata_S.obs['n_genes'].tolist())

b = [
    [1197, 7],
    [4384, 148]
]
fisher_exact(b)
# fraction of double positive T cells in S_dengue is more significant than dnegue

adata_T = adata_children[adata_children.obs['cell_type'] == 'Plasmablasts']
adata_S = adata_T[(adata_T.obs['Condition'] != 'Healthy') & (adata_T.obs['n_genes'] > 0)] 

c = {}
from scipy.stats import spearmanr, pearsonr

ar = adata_S.X.toarray()
ast = adata_S[:, 'IL6ST'].X.toarray()[:, 0]

for i, gene in enumerate(adata_S.var_names):
    if i%100 == 0:
        print(i)
    c[gene] = spearmanr(ar[:,i], ast)[0]

c = pd.Series(c).sort_values( ascending=False)

c = c.fillna(0).sort_values(ascending=False)

c = c.to_frame(name='r')

c['rank'] = np.arange(c.shape[0])