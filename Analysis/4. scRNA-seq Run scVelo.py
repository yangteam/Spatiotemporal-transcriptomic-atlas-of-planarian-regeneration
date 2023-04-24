import scvelo as scv
import scanpy as sc
scv.logging.print_version()
import numpy as np
import pandas as pd
import os

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization


# Load Data
adata = sc.read('adata_scRNA_clutered.h5ad')

file_list = [x for x in os.listdir('../Data/loom_adata') if x.endswith('h5ad')]
adata_loom = sc.read(os.path.join('loom_adata', file_list[0]))
for it in file_list[1:]:
    print(it)
    temp_adata = sc.read(os.path.join('loom_adata', it))
    adata_loom = adata_loom.concatenate(temp_adata)

adata_loom.obs_names = [x[:20] for x in adata_loom.obs_names]
adata_Velo = adata_loom[np.intersect1d(adata.obs_names, adata_loom.obs_names),]
UMAP_df = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names)
adata_Velo.obsm['X_umap'] = UMAP_df.loc[adata_Velo.obs_names,].values
adata_Velo.obs['Annotation'] = adata.obs.loc[adata_Velo.obs_names, 'Annotation']


# Run scVelo
scv.pp.filter_genes(adata_Velo, min_shared_counts=20)
scv.pp.normalize_per_cell(adata_Velo)
scv.pp.filter_genes_dispersion(adata_Velo, n_top_genes=2000)
scv.pp.log1p(adata_Velo)

scv.pp.filter_and_normalize(adata_Velo, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata_Velo, n_pcs=30, n_neighbors=30)

scv.tl.velocity(adata_Velo)
scv.tl.velocity_graph(adata_Velo)

sc.set_figure_params(figsize=[8, 8])
scv.pl.velocity_embedding_stream(adata_Velo, basis='umap', color='Annotation', density=5, arrow_style='->',
                                 legend_loc=False, max_length=10, save='Velo.png')