%config Completer.use_jedi = False
import os

import numpy as np
import pandas as pd

import anndata
import scanpy as sc

import matplotlib.pyplot as plt
import seaborn as sns

fn_loom = '/home/yike/phd/dengue/data/mergedata_20200930.loom'
adata = anndata.read_loom(fn_loom)

fn_obs = '/home/yike/phd/dengue/data/obs_2021_02_26.tsv'
obs_ZY = pd.read_csv(fn_obs, sep='\t', index_col=0)

adata_c = adata.copy()

obs_ZY = obs_ZY.loc[adata_c.obs.index]

adata_c.obs = obs_ZY

results_file = '/home/yike/phd/dengue/data/mergedata_new_ct_3_3_2021.h5ad'
adata_c.write(results_file)