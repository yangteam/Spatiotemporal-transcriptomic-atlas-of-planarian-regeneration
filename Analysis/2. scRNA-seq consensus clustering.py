import pandas as pd
import numpy as np
import scanpy as sc
import os
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.utils import shuffle

sc.set_figure_params(dpi=200, dpi_save=300)


# Load scRNA-seq data
adata_merge = sc.read('adata_scRNA_merged.h5ad')

sc.pp.normalize_total(adata_merge, target_sum=1e4)
sc.pp.log1p(adata_merge)
sc.pp.highly_variable_genes(adata_merge)

adata = adata_merge[:, adata_merge.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=30, metric='cosine')
sc.tl.umap(adata)
sc.tl.louvain(adata, resolution=2)

sc.set_figure_params(figsize=[3,3])
sc.pl.pca(adata[shuffle(adata.obs.index),], color='TimePoint')

sc.set_figure_params(figsize=[6,5.5])
sc.pl.umap(adata, color='louvain', legend_loc='on data', title='Louvain (resolution=2)')


# Run Louvain (resolution=2) 100 times with different random seeds
for it in range(100):
    sc.tl.louvain(adata, resolution=2, random_state=it)
    sc.set_figure_params(figsize=[6,5.5])
    sc.pl.umap(adata, color='louvain', legend_loc='on data', title='Louvain_Iter_'+str(it), 
               save='_Iter_'+str(it)+'.png', show=False)
    adata.obs.loc[:, ['TimePoint', 'louvain']].to_csv(os.path.join('Louvain', 'Iter_'+str(it)+'.txt'), sep='\t')
    print(it, '----', np.unique(adata.obs['louvain']).shape[0])
    plt.close()


# Run consensus clustering (DBSCAN algorithm)
file_list = ['Iter_'+str(x)+'.txt' for x in range(100)]
Total_df = pd.DataFrame()
for temp_file in file_list:
    temp_df = pd.read_csv('Louvain/'+temp_file, sep='\t', index_col=0)
    Total_df[temp_file.split('.txt')[0]] = temp_df['louvain']

from sklearn.cluster import DBSCAN

temp_epi = 0.1
min_s = 100

y_pred = DBSCAN(min_samples=min_s, metric='hamming', eps=temp_epi).fit_predict(Total_df)
res_df = pd.DataFrame(y_pred)
res_df.index = Total_df.index

adata.obs['DBSCAN'] = res_df.loc[adata.obs_names, 0]
adata_merge.obs['DBSCAN'] = res_df.loc[adata.obs_names, 0]


# Dealing with outliers using random forest
from imblearn.ensemble import BalancedRandomForestClassifier
from sklearn.metrics import confusion_matrix, pairwise_distances, balanced_accuracy_score, adjusted_rand_score

n_estimators = 1000
clf = BalancedRandomForestClassifier(n_estimators=n_estimators,
                                     criterion='gini',
                                     max_depth=None,
                                     min_samples_split=2,
                                     min_samples_leaf=2,
                                     min_weight_fraction_leaf=0.0,
                                     max_features='auto',
                                     max_leaf_nodes=None,
                                     min_impurity_decrease=0.0,
                                     bootstrap=True,
                                     oob_score=False,
                                     sampling_strategy='auto',
                                     replacement=False,
                                     random_state=0,
                                     verbose=0,
                                     warm_start=False,
                                     class_weight=None)

DBSCAN_df = adata.obs.copy()
PCA_df = pd.DataFrame(adata.obsm['X_pca'], index=adata.obs_names)

##### Outliers
Train_df = DBSCAN_df.loc[DBSCAN_df['DBSCAN']!=-1,]
Train_df = PCA_df.loc[Train_df.index,]
Test_df = DBSCAN_df.loc[DBSCAN_df['DBSCAN']==-1,]
Test_df = PCA_df.loc[Test_df.index,]

X_train = Train_df.values
Y_train = DBSCAN_df.loc[Train_df.index, 'DBSCAN'].values
X_test = Test_df.values

clf.fit(X_train, Y_train)

y_train_pred = clf.predict(X_train)
accuracy_train = balanced_accuracy_score(y_true=Y_train,
                                         y_pred=y_train_pred)

print('Random Forest accuracy:', accuracy_train)

y_test_pred = clf.predict(X_test)
y_test_pred_proba = clf.predict_proba(X_test)
Test_Max_proba = pd.DataFrame(np.max(y_test_pred_proba, axis=1))
Test_Max_proba.index = Test_df.index
Test_res = pd.DataFrame(y_test_pred)
Test_res.index = Test_df.index

Outlier_list = np.array(Test_proba.loc[Test_proba[1]<0.5,].index)
Outlier_list.shape

Test_res.loc[Outlier_list, 1] = -1

adata.obs.loc[Test_res.index, 'DBSCAN'] = Test_res[1]
adata_merge.obs.loc[Test_res.index, 'DBSCAN'] = Test_res[1]


# Filtering outliers
adata = adata[adata.obs['DBSCAN']!=-1,]
adata_merge = adata_merge[adata_merge.obs['DBSCAN']!=-1,]


# Saving scanpy object
adata.write('adata_scRNA_clutered.h5ad')
adata_merge.write('adata_scRNA_clutered-AllGene.h5ad')