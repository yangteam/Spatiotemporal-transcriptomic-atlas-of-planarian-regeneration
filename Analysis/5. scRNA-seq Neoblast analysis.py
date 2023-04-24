import scanpy as sc
import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.utils import shuffle

adata = sc.read('adata_scRNA_Annotated.h5ad')

subadata = adata[adata.obs['BigCellType']=='Neoblast',]

sc.pp.normalize_total(subadata, target_sum=1e4)
sc.pp.log1p(subadata)

sc.pp.highly_variable_genes(subadata, n_top_genes=3000)
subadata = subadata[:, subadata.var.highly_variable]

sc.pp.scale(subadata, max_value=10)
sc.tl.pca(subadata, svd_solver='arpack')

sc.pp.neighbors(subadata, n_neighbors=30, metric='cosine')
sc.tl.umap(subadata)
sc.tl.leiden(subadata,resolution=1.2)

subadata.write('adata_Neoblast.h5ad')