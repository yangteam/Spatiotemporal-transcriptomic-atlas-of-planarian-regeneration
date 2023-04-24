import sys

IN_COLAB = "google.colab" in sys.modules
#if IN_COLAB and branch == "stable":
#    !pip install --quiet git+https://github.com/BayraktarLab/cell2location#egg=cell2location[tutorials]

import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location
import scvi

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
import seaborn as sns
import os

from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel

import torch
torch.cuda.set_device('cuda:0')

TimePoint_list = ['0h', '6h', '12h', '24h', '3d', '7d']

# Load Visium Data
slides = []
for temp_timepoint in TimePoint_list:
    input_file = os.path.join('ST_scanpy',  'cut'+temp_timepoint+'_3D.h5ad')
    adata_vis = sc.read(input_file)
    
    #adata_ref = sc.read('../../Data/adata_Consensus_scRNA.h5ad')
    slides.append(adata_vis.copy())

# Combine anndata objects together
adata_vis = slides[0].concatenate(
    slides[1:],
    batch_key="sample",
    #uns_merge="unique",
    batch_categories=['cut'+x for x in TimePoint_list],
    index_unique=None
)

# Load scRNA-seq Data
adata_ref = sc.read('adata_scRNA_Annotated.h5ad')


# Running Cell2location
    
Overlap_Gene = np.intersect1d(adata_vis.var_names, adata_ref.var_names)

print(adata_vis.shape, adata_ref.shape, Overlap_Gene.shape)
    
    
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
# filter the object
adata_ref = adata_ref[:, selected].copy()
    
# prepare anndata for the regression model
scvi.data.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        #batch_key='L3', 
                        labels_key='BigCellType',
                        # cell type, covariate used for constructing signatures
                        #labels_key='Merge_Annotation',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        #categorical_covariate_keys=['Condition']
                        )
scvi.data.view_anndata_setup(adata_ref)
    
# create and train the regression model
mod = RegressionModel(adata_ref)
# Use all data for training (validation not implemented yet, train_size=1)
mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, use_gpu=True)
    
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
#inf_aver.iloc[0:5, 0:5]
    
    
#### Spatial Mapping
# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()
    
# prepare anndata for cell2location model
scvi.data.setup_anndata(adata=adata_vis, batch_key="sample")
scvi.data.view_anndata_setup(adata_vis)
    
# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection (using default here):
    detection_alpha=200
)

mod.train(max_epochs=30000,
            # train using full data (batch_size=None)
            batch_size=10000,
            # use all data points in training because
            # we need to estimate cell abundance at all locations
            train_size=1,
            use_gpu=True)

# plot ELBO loss history during training, removing first 100 epochs from the plot
#mod.plot_history(1000)
#plt.legend(labels=['full data training']);
    
    
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': False}
)
# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
    
adata_vis.obs.to_csv('cell2loc_L3.txt', sep='\t')
    
adata_vis.write('adata_cell2location.h5ad')