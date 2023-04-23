
import numpy as np
import pandas as pd
import hotspot
import matplotlib.pyplot as plt
import matplotlib.colors
import seaborn as sns
import os
import pickle


counts_file = 'cut12h_counts.txt'
pos_file = 'top8_slices_position.txt'

pos = pd.read_csv(pos_file, index_col=0, sep=" ",header=None)
counts = pd.read_csv(counts_file, index_col=0,sep=" ") # Takes a while, ~10min

# Align the indices
counts = counts.loc[:, pos.index]
barcodes = pos.index.values
num_umi = counts.sum(axis=0)

# Filter genes
# gene_counts = (counts > 0).sum(axis=1)
# valid_genes = gene_counts >= 50
# counts = counts.loc[valid_genes]

hs = hotspot.Hotspot(counts, model='danb', latent=pos, umi_counts=num_umi)

hs.create_knn_graph(
    weighted_graph=True, n_neighbors=30,
)

hs_results = hs.compute_autocorrelations(jobs=10)

hs.results.to_csv(path_or_buf="gene_results.txt",sep=" ")

hs_genes = hs_results.head(6000).index

# hs_genes = hs_results.index[hs_results.FDR < 0.01]

lcz = hs.compute_local_correlations(hs_genes, jobs=10)


modules = hs.create_modules(
    min_gene_threshold=50, core_only=False, fdr_threshold=0.05
)

modules.to_csv(path_or_buf="cluster.txt",sep=" ")

fig=plt.figure(figsize=(6,6))
fig=hs.plot_local_correlations()
plt.savefig('cor.png',dpi=200)
