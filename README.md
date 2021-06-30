# cell-cell communications in dengue children

To map the complex interactions between distinct immune cell types, we extracted the list of known potential interaction partners from CellPhoneDB1and complemented it with series additional interactions from recent immunology literature. 

The interactions for kids with SD and D were counted if both ligand and receptor are expressed in more than 2 % of cells within the respective cell type.  Differential interactions between D and SD were classified as upregulated, downregulated, or mixregulated according to which category either of the genes belongs to.

A nonparametric label randomization test was utilized to check the significance of differential interactions. The P-value was calculated as the fraction of randomizations for which the distance from the origin was larger (i.e., differential expressions were more extreme) than in the real data. 

## Folder structure

- `dengue_children`: this folder contains the code used for the analysis

- `data`: this folder contains some Excel/TSV files after analyzing the data and associated figures
  - `paper_figure`: this folder contains the cell-cell communication panels for the dengue children paper
  - `tables`: this folder contains useful cell-cell communication tables for the dengue children paper, e.g., differential gene expression, differential interactions
