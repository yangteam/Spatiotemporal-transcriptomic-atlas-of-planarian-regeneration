import pandas as pd
import numpy as np
import scanpy as sc
import os
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.utils import shuffle
import scrublet as scr

sc.set_figure_params(dpi=200, dpi_save=300)

TimePoint = ['cut0d', 'cut6h', 'cut12h', 'cut24h', 'cut2d', 'cut3d', 'cut5d', 'cut7d']


# Load scRNA-seq data
adata_list = {}
for temp_time in TimePoint:
    print('--------------', temp_time)
    adata_list[temp_time] = sc.read_10x_mtx('../Data/scRNA/'+temp_time, cache=False)
    adata_list[temp_time].obs['TimePoint'] = temp_time
    sc.pp.calculate_qc_metrics(adata_list[temp_time], inplace=True)


num_df = pd.DataFrame()
for temp_time in TimePoint:
    num_df.loc[temp_time, 'Cell Number'] = adata_list[temp_time].shape[0]
fig, ax = plt.subplots(figsize=[3.8,2.8])
g=sns.barplot(x=num_df.index, y=num_df['Cell Number'], palette="deep", ax=ax)
it = 0
for index, row in num_df.iterrows():
    g.text(it,row['Cell Number']+100, round(row['Cell Number']), color='black', ha="center", size=7)
    it+=1
plt.title('#Cells before filtering')
plt.ylim([0,10000])
plt.yticks([0,2500,5000,7500,10000])
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('figures/Num_Cells_Before_filtering.png',dpi=300)


# Data preprocessing
fig, ax = plt.subplots(3,3, figsize=[15, 10])
for it, temp_time in enumerate(TimePoint):
    sc.pl.scatter(adata_list[temp_time], x='total_counts', y='n_genes_by_counts', title='scRNA-seq ('+temp_time+')', ax=ax[it//3,it%3], show=False)
plt.tight_layout()
plt.savefig('figures/scatter_before_filter.png',dpi=300)
plt.show()

for temp_time in TimePoint:
    print('-----', temp_time)
    print(adata_list[temp_time].shape)
    adata_list[temp_time] = adata_list[temp_time][adata_list[temp_time].obs.n_genes_by_counts < 6000, :]
    print(adata_list[temp_time].shape)
    adata_list[temp_time] = adata_list[temp_time][adata_list[temp_time].obs.n_genes_by_counts > 250, :]
    print(adata_list[temp_time].shape)

fig, ax = plt.subplots(3,3, figsize=[15, 10])
for it, temp_time in enumerate(TimePoint):
    sc.pl.scatter(adata_list[temp_time], x='total_counts', y='n_genes_by_counts', title='scRNA-seq ('+temp_time+')', ax=ax[it//3,it%3], show=False)
plt.tight_layout()
plt.savefig('figures/scatter_After_filter.png',dpi=300)
plt.show()

num_df = pd.DataFrame()
for temp_time in TimePoint:
    num_df.loc[temp_time, 'Cell Number'] = adata_list[temp_time].shape[0]

fig, ax = plt.subplots(figsize=[3.8,2.8])
g=sns.barplot(x=num_df.index, y=num_df['Cell Number'], palette="deep", ax=ax)
it = 0
for index, row in num_df.iterrows():
    g.text(it,row['Cell Number']+100, round(row['Cell Number']), color='black', ha="center", size=7)
    it+=1
plt.title('#Cells after filtering')
plt.ylim([0,10000])
plt.yticks([0,2500,5000,7500,10000])
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('figures/Num_Cells_After_filtering.png',dpi=300)


# Get merged scRNA-seq SCANPY object
adata_merge = adata_list['cut0d'].concatenate(adata_list['cut6h'], adata_list['cut12h'],  adata_list['cut24h'], 
                                              adata_list['cut2d'], adata_list['cut3d'], adata_list['cut5d'], adata_list['cut7d'])


# Filteing doublets using scrublet
adata_hvg = {}
for temp_time in TimePoint:
    print('------------', temp_time)
    sc.pp.normalize_total(adata_list[temp_time], target_sum=1e4)
    sc.pp.log1p(adata_list[temp_time])
    sc.pp.highly_variable_genes(adata_list[temp_time])
    adata_hvg[temp_time] = adata_list[temp_time][:, adata_list[temp_time].var.highly_variable]
    sc.pp.regress_out(adata_hvg[temp_time], ['total_counts'])

for temp_time in TimePoint:
    print('------------', temp_time)
    counts_matrix = pd.DataFrame(adata_list[temp_time].X.toarray(), index=adata_list[temp_time].obs_names, columns=adata_list[temp_time].var_names)
    scrub = scr.Scrublet(counts_matrix)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    scr_res = pd.DataFrame({'scrublet_scores': doublet_scores, 'scr_res': predicted_doublets}, index=counts_matrix.index)
    
    adata_hvg[temp_time].obs['scrublet_scores'] = scr_res.loc[adata_hvg[temp_time].obs_names, 'scrublet_scores']
    adata_hvg[temp_time].obs['scr_res'] = scr_res.loc[adata_hvg[temp_time].obs_names, 'scr_res']

adata_merge_hvg = adata_list['cut0d'].concatenate(adata_hvg['cut6h'], adata_hvg['cut12h'],  adata_hvg['cut24h'], 
                                              adata_hvg['cut2d'], adata_hvg['cut3d'], adata_hvg['cut5d'], adata_hvg['cut7d'])

adata_merge.obs['scrublet_score'] = adata_merge_hvg.obs.loc[adata_merge.obs_names, 'scrublet_score']
adata_merge = adata_merge[adata_merge.obs['scrublet_score']<0.5,]

# Saving scanpy object
adata_merge.write('adata_scRNA_merged.h5ad')