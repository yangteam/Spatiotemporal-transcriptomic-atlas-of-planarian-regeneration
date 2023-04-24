import pandas as pd
import numpy as np
import scanpy as sc
import os
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.utils import shuffle

sc.set_figure_params(dpi=200, dpi_save=300)

# Load Data
adata = sc.read('adata_scRNA_clutered.h5ad')


# Annotation
Annotation_dict = {}
Annotation_dict["0"]="neural progenitors"
Annotation_dict["1"]="neural progenitors"
Annotation_dict["2"]="early epidermal progenitors 2"
Annotation_dict["3"]="GABA neurons"
Annotation_dict["4"]="parenchymal progenitors"
Annotation_dict["5"]="muscle DV"
Annotation_dict["25"]="secretory 3"
Annotation_dict["6"]="mecom+ neurons"
Annotation_dict["7"]="muscle body 1"
Annotation_dict["8"]="psap+ parenchymal cells"
Annotation_dict["9"]="phagocyte progenitors"
Annotation_dict["10"]="brn+ parenchymal cells"
Annotation_dict["11"]="gut neoblast"
Annotation_dict["12"]="gnmt+ goblet cells"
Annotation_dict["13"]="otf+ cells 1"
Annotation_dict["14"]="aqp+ parenchymal cells"
Annotation_dict["15"]="secretory 6"
Annotation_dict["16"]="cpa-2+ goblet cells"
Annotation_dict["17"]="secretory neoblast"
Annotation_dict["18"]="tspan-1+ neoblast"
Annotation_dict["19"]="pgrn+ parenchymal cells"
Annotation_dict["32"]="early epidermal progenitors 3"
Annotation_dict["20"]="novel neoblast"
Annotation_dict["21"]="neuron neoblast"
Annotation_dict["22"]="secretory 2"
Annotation_dict["23"]="muscle/parenchymal neoblast"
Annotation_dict["24"]="npp-4+ neurons"
Annotation_dict["26"]="zic-3+ neurons"
Annotation_dict["27"]="early epidermal progenitors 1"
Annotation_dict["28"]="pharynx cell type"
Annotation_dict["29"]="muscle pharynx progenitors"
Annotation_dict["30"]="pigment"
Annotation_dict["31"]="epidermal neoblast"
Annotation_dict["33"]="secretory 5"
Annotation_dict["34"]="psd+ cells"
Annotation_dict["35"]="phagocytes"
Annotation_dict["36"]="otf+ cells progenitors"
Annotation_dict["37"]="protonephridia"
Annotation_dict["38"]="Idlrr-1+ parenchymal cells"
Annotation_dict["39"]="protonephridia neoblast"
Annotation_dict["40"]="muscle DV progenitors"
Annotation_dict["41"]="goblet cell progenitors"
Annotation_dict["42"]="muscle body progenitors"
Annotation_dict["43"]="late epidermal progenitors"
Annotation_dict["44"]="secretory 4"
Annotation_dict["45"]="muscle body 2"
Annotation_dict["46"]="secretory 8"
Annotation_dict["47"]="ChAT neurons 1"
Annotation_dict["48"]="secretory 9"
Annotation_dict["49"]="muscle pharynx"
Annotation_dict["50"]="muscle genital"
Annotation_dict["51"]="epidermis DVb"
Annotation_dict["52"]="secretory 1"
Annotation_dict["53"]="parapharyngeal neoblast"
Annotation_dict["54"]="neural progenitors"
Annotation_dict["55"]="secretory 7"
Annotation_dict["56"]="cav-1+ neurons"
Annotation_dict["57"]="epidermis DVb progenitors"
Annotation_dict["58"]="slc6a-3+ neurons"
Annotation_dict["59"]="scna-1+ neurons"
Annotation_dict["60"]="ChAT neurons 2"
Annotation_dict["61"]="epidermis"
Annotation_dict["62"]="otf+ cells 2"


# Plot Annotation Colors
Color_dict = {}
Color_dict["ChAT neurons 1"]="#f16913"
Color_dict["ChAT neurons 2"]="#7a0177"
Color_dict["GABA neurons"]="#fd8d3c"
Color_dict["Idlrr-1+ parenchymal cells"]="#ae017e"
Color_dict["aqp+ parenchymal cells"]="#C71585"
Color_dict["slc6a-3+ neurons"]="#7f2704"
Color_dict["early epidermal progenitors 1"]="#9e9ac8"
Color_dict["early epidermal progenitors 2"]="#807dba"
Color_dict["early epidermal progenitors 3"]="#6a51a3"
Color_dict["epidermal neoblast"]="#dadaeb"
Color_dict["epidermis"]="#54278f"
Color_dict["epidermis DVb"]="#A289AA"
Color_dict["epidermis DVb progenitors"]="#bcbddc"
Color_dict["gnmt+ goblet cells"]="#006d2c"
Color_dict["cpa-2+ goblet cells"]="#238b45"
Color_dict["gut neoblast"]="#e5f5e0"
Color_dict["goblet cell progenitors"]="#a1d99b"
Color_dict["phagocyte progenitors"]="#41ab5d"
Color_dict["late epidermal progenitors"]="#885578"
Color_dict["muscle body 1"]="#67000d"
Color_dict["muscle body 2"]="#cb181d"
Color_dict["muscle DV"]="#ef3b2c"
Color_dict["muscle/parenchymal neoblast"]="#fff5f0"
Color_dict["muscle pharynx progenitors"]="#fee0d2"
Color_dict["muscle genital"]="#fb6a4a"
Color_dict["muscle pharynx"]="#a50f15"
Color_dict["muscle body progenitors"]="#fb6a4a"
Color_dict["muscle DV progenitors"]="#fcbba1"
Color_dict["neural progenitors"]="#fdd0a2"
Color_dict["neuron neoblast"]="#F4DFC8"
Color_dict["npp-4+ neurons"]="#FFA500"
Color_dict["zic-3+ neurons"]="#f16913"
Color_dict["cav-1+ neurons"]="#67001f"
Color_dict["scna-1+ neurons"]="#7a0177"
Color_dict["mecom+ neurons"]="#BDB76B"
Color_dict["tspan-1+ neoblast"]="#FFF0F5"
Color_dict["novel neoblast"]="#737373"
Color_dict["otf+ cells 1"]="#b91d73"
Color_dict["otf+ cells 2"]="#b91d73"
Color_dict["otf+ cells progenitors"]="#DDA0DD"
Color_dict["parapharyngeal neoblast"]="#F4DFBA"
Color_dict["parenchymal progenitors"]="#fde0dd"
Color_dict["pgrn+ parenchymal cells"]="#f768a1"
Color_dict["phagocytes"]="#00441b"
Color_dict["pharynx cell type"]="#3f007d"
Color_dict["pigment"]="#7a0177"
Color_dict["protonephridia"]="#A0522D"
Color_dict["protonephridia neoblast"]="#DEB887"
Color_dict["brn+ parenchymal cells"]="#fa9fb5"
Color_dict["psap+ parenchymal cells"]="#dd3497"
Color_dict["psd+ cells"]="#74c476"
Color_dict["secretory 1"]="#6baed6"
Color_dict["secretory 2"]="#08519c"
Color_dict["secretory 3"]="#2171b5"
Color_dict["secretory 4"]="#9ecae1"
Color_dict["secretory neoblast"]="#F0F8FF"
Color_dict["secretory 5"]="#dadaeb"
Color_dict["secretory 7"]="#08306b"
Color_dict["secretory 6"]="#c6dbef"
Color_dict["secretory 8"]="#deebf7"
Color_dict["secretory 9"]="#4292c6"

# BigCellType 
BigCellType_dict = {}
BigCellType_dict["neural progenitors"]="Neuronal"
BigCellType_dict["early epidermal progenitors 2"]="Epidermal"
BigCellType_dict["GABA neurons"]="Neuronal"
BigCellType_dict["parenchymal progenitors"]="Parenchymal"
BigCellType_dict["muscle DV"]="Muscle"
BigCellType_dict["secretory 3"]="Secretory"
BigCellType_dict["mecom+ neurons"]="Neuronal"
BigCellType_dict["muscle body 1"]="Muscle"
BigCellType_dict["psap+ parenchymal cells"]="Parenchymal"
BigCellType_dict["phagocyte progenitors"]="Gut"
BigCellType_dict["brn+ parenchymal cells"]="Parenchymal"
BigCellType_dict["gut neoblast"]="Neoblast"
BigCellType_dict["gnmt+ goblet cells"]="Gut"
BigCellType_dict["otf+ cells 1"]="Neuronal"
BigCellType_dict["aqp+ parenchymal cells"]="Parenchymal"
BigCellType_dict["secretory 6"]="Secretory"
BigCellType_dict["cpa-2+ goblet cells"]="Gut"
BigCellType_dict["secretory neoblast"]="Neoblast"
BigCellType_dict["tspan-1+ neoblast"]="Neoblast"
BigCellType_dict["pgrn+ parenchymal cells"]="Parenchymal"
BigCellType_dict["early epidermal progenitors 3"]="Epidermal"
BigCellType_dict["novel neoblast"]="Neoblast"
BigCellType_dict["neuron neoblast"]="Neoblast"
BigCellType_dict["secretory 2"]="Secretory"
BigCellType_dict["muscle/parenchymal neoblast"]="Neoblast"
BigCellType_dict["npp-4+ neurons"]="Neuronal"
BigCellType_dict["zic-3+ neurons"]="Neuronal"
BigCellType_dict["early epidermal progenitors 1"]="Epidermal"
BigCellType_dict["pharynx cell type"]="Epidermal"
BigCellType_dict["muscle pharynx progenitors"]="Muscle"
BigCellType_dict["pigment"]="Parenchymal"
BigCellType_dict["epidermal neoblast"]="Neoblast"
BigCellType_dict["secretory 5"]="Secretory"
BigCellType_dict["psd+ cells"]="Gut"
BigCellType_dict["phagocytes"]="Gut"
BigCellType_dict["otf+ cells progenitors"]="Neuronal"
BigCellType_dict["protonephridia"]="Protonephridia"
BigCellType_dict["Idlrr-1+ parenchymal cells"]="Parenchymal"
BigCellType_dict["protonephridia neoblast"]="Neoblast"
BigCellType_dict["muscle DV progenitors"]="Muscle"
BigCellType_dict["goblet cell progenitors"]="Gut"
BigCellType_dict["muscle body progenitors"]="Muscle"
BigCellType_dict["late epidermal progenitors"]="Epidermal"
BigCellType_dict["secretory 4"]="Secretory"
BigCellType_dict["muscle body 2"]="Muscle"
BigCellType_dict["secretory 8"]="Secretory"
BigCellType_dict["ChAT neurons 1"]="Neuronal"
BigCellType_dict["secretory 9"]="Secretory"
BigCellType_dict["muscle pharynx"]="Muscle"
BigCellType_dict["muscle genital"]="Muscle"
BigCellType_dict["epidermis DVb"]="Epidermal"
BigCellType_dict["secretory 1"]="Secretory"
BigCellType_dict["parapharyngeal neoblast"]="Neoblast"
BigCellType_dict["secretory 7"]="Secretory"
BigCellType_dict["cav-1+ neurons"]="Neuronal"
BigCellType_dict["epidermis DVb progenitors"]="Epidermal"
BigCellType_dict["slc6a-3+ neurons"]="Neuronal"
BigCellType_dict["scna-1+ neurons"]="Neuronal"
BigCellType_dict["ChAT neurons 2"]="Neuronal"
BigCellType_dict["epidermis"]="Epidermal"
BigCellType_dict["otf+ cells 2"]="Neuronal"


adata.obs['Annotation'] = list(adata.obs['DBSCAN'].map(Annotation_dict))
adata.obs['Annotation'] = adata.obs['Annotation'].astype('category')

adata.uns['Annotation_colors'] = [Color_dict[x] for x in adata.obs['Annotation'].cat.categories]

adata.obs['BigCellType'] = list(adata.obs['Annotation'].map(BigCellType_dict))
adata.obs['BigCellType'] = adata.obs['BigCellType'].astype('category')


sc.set_figure_params(figsize=[10, 10])
sc.pl.umap(adata, color='BigCellType', legend_loc='on data', s=12,
           add_outline=True, legend_fontsize=12, 
           legend_fontoutline=2)

sc.set_figure_params(figsize=[10, 10])
sc.pl.umap(adata, color='Annotation', legend_loc='on data', s=12,
           add_outline=True, legend_fontsize=12, 
           legend_fontoutline=2)


# Saving scanpy object
adata.write('adata_scRNA_Annotated.h5ad')