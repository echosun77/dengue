import numpy as np
import pandas as pd

save_tables = '/home/yike/projects/dengue/data/'
DEIs_fn = '/home/yike/projects/dengue/data/cst_inters_med_pair1_pair39_56_exp002.tsv'
pair_fn = '/home/yike/projects/dengue/data/cst_pair.tsv'

cst_pair = pd.read_csv(pair_fn, sep='\t', index_col=['cell_subtype', 'gene'])

df = pd.read_csv(DEIs_fn, sep='\t')

ga_info = cst_pair.loc[np.array(df[['csta', 'ga']]).tolist()][['med_pair', 'fra_pair', 'neg_fra_pair','S_fra', 'NS_fra']]
ga_info.columns = ['ga_' + col for col in ga_info.columns]

gb_info = cst_pair.loc[np.array(df[['cstb', 'gb']]).tolist()][['med_pair', 'fra_pair', 'neg_fra_pair','S_fra', 'NS_fra']]
gb_info.columns = ['gb_' + col for col in gb_info.columns]

df2 = pd.merge(df, ga_info, left_on = ['csta', 'ga'], right_index=True)
df2 = pd.merge(df2, gb_info, left_on = ['cstb', 'gb'], right_index=True)

df2 = df2[~ df2[['ga', 'csta', 'gb', 'cstb']].duplicated()]
df2.to_csv('/home/yike/projects/dengue/data/cst_inters_med_pair1_pair39_56_exp002_full.tsv', sep='\t', index=False)
           


