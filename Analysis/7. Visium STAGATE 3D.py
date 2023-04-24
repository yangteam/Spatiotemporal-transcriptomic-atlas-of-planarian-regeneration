import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys

import warnings
warnings.filterwarnings("ignore")

import STAGATE

temp_timepoint = 'cut0h'

adata = sc.read(os.path.join('ST_scanpy', temp_timepoint+'_3D.h5ad'))

#Normalization
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


# Running PCA (cell type-aware module of STAGATE)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata)
sc.tl.umap(adata)

sc.tl.louvain(adata, resolution=0.1)

adata.obsm['spatial'] = adata.obs.loc[:, ['X_3D', 'Y_3D']].values

# Running STAGATE
STAGATE.Cal_Spatial_Net_3D(adata, rad_cutoff_2D=150, rad_cutoff_Zaxis=100,
                           key_section='Section_id', section_order = [str(it+1) for it in range(17)],
                           verbose=True)

STAGATE.Stats_Spatial_Net(adata)

adata = STAGATE.train_STAGATE(adata, alpha=0.5, pre_labels='louvain',
                              n_epochs=1000)


# mclust clustering
adata = STAGATE.mclust_R(adata, used_obsm='STAGATE', num_cluster=10)