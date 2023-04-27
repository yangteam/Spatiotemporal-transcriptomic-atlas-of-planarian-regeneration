from doctest import FAIL_FAST
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sc.settings.verbosity = 1 
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')

adata=sc.read_10x_mtx('./PLK1/outs/filtered_feature_bc_matrix')
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata, save="variable_genes.pdf")
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, 'total_counts')
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='SMED30007406',save="SMED30007406.pdf")

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata, init_pos='paga')
sc.tl.louvain(adata)

adata_ref = sc.read_h5ad("adata_scRNA_Annotated.h5ad")

var_names = adata_ref.var_names.intersection(adata.var_names)
adata_ref = adata_ref[:, var_names]
adata = adata[:, var_names]

sc.pp.log1p(adata_ref)
sc.tl.pca(adata_ref, svd_solver='arpack')
sc.pp.highly_variable_genes(adata_ref, min_mean=0.0125, max_mean=3, min_disp=0.5)

# adata_ref.uns['pca']['params']['use_highly_variable']
# adata_ref.varm['PCs'][adata.var['highly_variable']]

anno = pd.read_csv("adata_CellType.txt",sep="\t",header=0, index_col=0)

sc.tl.ingest(adata, adata_ref, obs='Annotation',embedding_method="umap")

adata.obs["TimePoint"]="plk1"
adata_concat = adata_ref.concatenate(adata, batch_categories=['ref', 'plk1'])
adata_concat.obs.Annotation = adata_concat.obs.Annotation.astype('category')

adata_concat_plk1 = adata_concat[adata_concat.obs['TimePoint'].isin(['cut5d', 'plk1'])]
sc.pl.umap(adata_concat_plk1, color='TimePoint',save="TimePoint.pdf", palette=["#8993D5","#B34C65"])

