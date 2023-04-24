import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

def sc_read_visium_mtx(path, library_id, count_file='filtered_feature_bc_matrix'):
    from pathlib import Path, PurePath
    from matplotlib.image import imread
    import json
    
    adata = sc.read_10x_mtx(os.path.join(path, count_file))
    
    path = Path(path)
    adata.uns["spatial"] = dict()
    adata.uns["spatial"][library_id] = dict()
    
    files = dict(tissue_positions_file=path / 'spatial/tissue_positions_list.csv',
                 scalefactors_json_file=path / 'spatial/scalefactors_json.json',
                 hires_image=path / 'spatial/tissue_hires_image.png',
                 lowres_image=path / 'spatial/tissue_lowres_image.png')
    
    
    for f in files.values():
        if not f.exists():
            if any(x in str(f) for x in ["hires_image", "lowres_image"]):
                logg.warning(
                    f"You seem to be missing an image file.\n"
                    f"Could not find '{f}'.")
            else:
                raise OSError(f"Could not find '{f}'")
            
    adata.uns["spatial"][library_id]['images'] = dict()
    for res in ['hires', 'lowres']:
        try:
            adata.uns["spatial"][library_id]['images'][res] = imread(str(files[f'{res}_image']))
        except Exception:
            raise OSError(f"Could not find '{res}_image'")

            adata.uns["spatial"][library_id]['scalefactors'] = json.loads(
                files['scalefactors_json_file'].read_bytes()
            )
            
            
    adata.uns["spatial"][library_id]['images'] = dict()
    for res in ['hires', 'lowres']:
        try:
            adata.uns["spatial"][library_id]['images'][res] = imread(str(files[f'{res}_image']))
        except Exception:
            raise OSError(f"Could not find '{res}_image'")

        # read json scalefactors
    adata.uns["spatial"][library_id]['scalefactors'] = json.loads(
        files['scalefactors_json_file'].read_bytes()
    )

    #adata.uns["spatial"][library_id]["metadata"] = {
        #k: (str(attrs[k], "utf-8") if isinstance(attrs[k], bytes) else attrs[k])
        #for k in ("chemistry_description", "software_version")
        #if k in attrs}
        
    positions = pd.read_csv(files['tissue_positions_file'], header=None)
    positions.columns = [
        'barcode',
        'in_tissue',
        'array_row',
        'array_col',
        'pxl_col_in_fullres',
        'pxl_row_in_fullres',
    ]
    positions.index = positions['barcode']

    adata.obs = adata.obs.join(positions, how="left")

    adata.obsm['spatial'] = adata.obs[
        ['pxl_row_in_fullres', 'pxl_col_in_fullres']
    ].to_numpy()
    adata.obs.drop(
        columns=['barcode', 'pxl_row_in_fullres', 'pxl_col_in_fullres'],
        inplace=True,
    )

    return adata


TimePoint_list = ['cut0h', 'cut6h', 'cut12h', 'cut24h', 'cut3d', 'cut7d']
section_list = ['A1', 'B1', 'C1', 'D1']

# Load cut0h data
section_list = ['A1', 'B1', 'C1', 'D1']
adata_list = dict()
for section_it in section_list:
    input_dir = os.path.join('../Data/ST', 'cut0h', section_it)
    
    adata_list[section_it] = sc.read_visium(path=input_dir, library_id='cut0h_'+section_it)
    
    adata_list[section_it].var_names_make_unique()
    sc.pp.calculate_qc_metrics(adata_list[section_it], inplace=True)

adata_spatial = adata_list['A1'].concatenate([adata_list['B1'], adata_list['C1'], adata_list['D1']],
                                  batch_key="library_id",
                                  uns_merge="unique",
                                  batch_categories=[k 
                                                    for d in [adata_list['A1'].uns["spatial"], adata_list['B1'].uns["spatial"], adata_list['C1'].uns["spatial"], adata_list['D1'].uns["spatial"]]
                                                    for k, v in d.items()
                                                   ],)

ST_aligned_dir = '../Data/ST_aligned/cut0d'

Total_spot = 0
trans_dict = {}
trans_dict['_1'] = '-cut0h_A1'
trans_dict['_2'] = '-cut0h_B1'
trans_dict['_3'] = '-cut0h_C1'
trans_dict['_4'] = '-cut0h_D1'
file_list = [x for x in os.listdir(ST_aligned_dir) if x.endswith('.txt')]
for temp_file in file_list:
    temp = pd.read_csv(os.path.join(ST_aligned_dir, temp_file), sep='\t', header=None)
    temp.columns = ['X_3D', 'Y_3D', 'spot']
    temp['spot'] = temp['spot'].map(lambda x: x.split('_')[0] + trans_dict['_'+x.split('_')[1]])
    temp.index = temp['spot']
    overlap_spot = np.intersect1d(temp['spot'], adata_spatial.obs_names)
    print(temp.shape, overlap_spot.shape)
    
    Total_spot+=temp.shape[0]
    
    adata_spatial.obs.loc[overlap_spot, 'X_3D'] = temp.loc[overlap_spot, 'X_3D']
    adata_spatial.obs.loc[overlap_spot, 'Y_3D'] = temp.loc[overlap_spot, 'Y_3D']
    adata_spatial.obs.loc[overlap_spot, 'Z_3D'] = int(temp_file.split('_')[0][4:])
    
print('-----', Total_spot, '---', adata_spatial.shape[0])

# Saving cut0h scanpy object
adata_spatial.write('ST_scanpy/cut0h_merge.h5ad')


# cut6h
section_list = ['A1', 'B1', 'C1', 'D1']
adata_list = dict()
for section_it in section_list:
    input_dir = os.path.join('../Data/ST', 'cut6h', section_it)
    adata_list[section_it] = sc_read_visium_mtx(path=input_dir, library_id='cut6h_'+section_it)
    adata_list[section_it].obs['library_id'] = 'cut6h_'+section_it
    
    adata_list[section_it].var_names_make_unique()
    sc.pp.calculate_qc_metrics(adata_list[section_it], inplace=True)

adata_spatial = adata_list['A1'].concatenate([adata_list['B1'], adata_list['C1'], adata_list['D1']],
                                  batch_key="library_id",
                                  uns_merge="unique",
                                  batch_categories=[k 
                                                    for d in [adata_list['A1'].uns["spatial"], adata_list['B1'].uns["spatial"], adata_list['C1'].uns["spatial"], adata_list['D1'].uns["spatial"]]
                                                    for k, v in d.items()
                                                   ],)

ST_aligned_dir = '../Data/ST_aligned/cut6h'

Total_spot = 0
trans_dict = {}
trans_dict['_1'] = '-cut6h_A1'
trans_dict['_2'] = '-cut6h_B1'
trans_dict['_3'] = '-cut6h_C1'
trans_dict['_4'] = '-cut6h_D1'
file_list = [x for x in os.listdir(ST_aligned_dir) if x.endswith('.txt')]
for temp_file in file_list:
    temp = pd.read_csv(os.path.join(ST_aligned_dir, temp_file), sep='\t', header=None)
    temp.columns = ['X_3D', 'Y_3D', 'spot']
    temp['spot'] = temp['spot'].map(lambda x: x.split('_')[0] + trans_dict['_'+x.split('_')[1]])
    temp.index = temp['spot']
    overlap_spot = np.intersect1d(temp['spot'], adata_spatial.obs_names)
    print(temp.shape, overlap_spot.shape)
    
    Total_spot+=temp.shape[0]
    
    adata_spatial.obs.loc[overlap_spot, 'X_3D'] = temp.loc[overlap_spot, 'X_3D']
    adata_spatial.obs.loc[overlap_spot, 'Y_3D'] = temp.loc[overlap_spot, 'Y_3D']
    adata_spatial.obs.loc[overlap_spot, 'Z_3D'] = int(temp_file.split('_')[0][4:])
    
print('-----', Total_spot, '---', adata_spatial.shape[0])

adata_spatial.write('ST_scanpy/cut6h_merge.h5ad')


# cut12h
section_list = ['A1', 'B1', 'C1', 'D1']
adata_list = dict()
for section_it in section_list:
    input_dir = os.path.join('../Data/ST', 'cut12h', section_it)
    
    adata_list[section_it] = sc.read_visium(path=input_dir, library_id='cut12h_'+section_it)
    
    adata_list[section_it].var_names_make_unique()
    sc.pp.calculate_qc_metrics(adata_list[section_it], inplace=True)

adata_spatial = adata_list['A1'].concatenate([adata_list['B1'], adata_list['C1'], adata_list['D1']],
                                  batch_key="library_id",
                                  uns_merge="unique",
                                  batch_categories=[k 
                                                    for d in [adata_list['A1'].uns["spatial"], adata_list['B1'].uns["spatial"], adata_list['C1'].uns["spatial"], adata_list['D1'].uns["spatial"]]
                                                    for k, v in d.items()
                                                   ],)

ST_aligned_dir = '../Data/ST_aligned/cut12h'

Total_spot = 0
trans_dict = {}
trans_dict['_1'] = '-cut12h_A1'
trans_dict['_2'] = '-cut12h_B1'
trans_dict['_3'] = '-cut12h_C1'
trans_dict['_4'] = '-cut12h_D1'
file_list = [x for x in os.listdir(ST_aligned_dir) if x.endswith('.txt')]
for temp_file in file_list:
    temp = pd.read_csv(os.path.join(ST_aligned_dir, temp_file), sep='\t', header=None)
    temp.columns = ['spot', 'X_3D', 'Y_3D']
    temp['spot'] = temp['spot'].map(lambda x: x.split('_')[0] + trans_dict['_'+x.split('_')[1]])
    temp.index = temp['spot']
    overlap_spot = np.intersect1d(temp['spot'], adata_spatial.obs_names)
    print(temp.shape, overlap_spot.shape)
    
    Total_spot+=temp.shape[0]
    
    adata_spatial.obs.loc[overlap_spot, 'X_3D'] = temp.loc[overlap_spot, 'X_3D']
    adata_spatial.obs.loc[overlap_spot, 'Y_3D'] = temp.loc[overlap_spot, 'Y_3D']
    adata_spatial.obs.loc[overlap_spot, 'Z_3D'] = int(temp_file.split('.txt')[0][4:])
    
print('-----', Total_spot, '---', adata_spatial.shape[0])

adata_spatial.write('ST_scanpy/cut12h_merge.h5ad')


# cut24h
section_list = ['A1', 'B1', 'C1', 'D1']
adata_list = dict()
for section_it in section_list:
    input_dir = os.path.join('../Data/ST', 'cut24h', section_it)
    adata_list[section_it] = sc_read_visium_mtx(path=input_dir, library_id='cut24h_'+section_it)
    adata_list[section_it].obs['library_id'] = 'cut24h_'+section_it
    
    adata_list[section_it].var_names_make_unique()
    sc.pp.calculate_qc_metrics(adata_list[section_it], inplace=True)

adata_spatial = adata_list['A1'].concatenate([adata_list['B1'], adata_list['C1'], adata_list['D1']],
                                  batch_key="library_id",
                                  uns_merge="unique",
                                  batch_categories=[k 
                                                    for d in [adata_list['A1'].uns["spatial"], adata_list['B1'].uns["spatial"], adata_list['C1'].uns["spatial"], adata_list['D1'].uns["spatial"]]
                                                    for k, v in d.items()
                                                   ],)

ST_aligned_dir = '../Data/ST_aligned/cut24h'

Total_spot = 0
trans_dict = {}
trans_dict['_1'] = '-cut24h_A1'
trans_dict['_2'] = '-cut24h_B1'
trans_dict['_3'] = '-cut24h_C1'
trans_dict['_4'] = '-cut24h_D1'
file_list = [x for x in os.listdir(ST_aligned_dir) if x.endswith('.txt')]
for temp_file in file_list:
    temp = pd.read_csv(os.path.join(ST_aligned_dir, temp_file), sep='\t', header=None)
    temp.columns = ['X_3D', 'Y_3D', 'spot']
    temp['spot'] = temp['spot'].map(lambda x: x.split('_')[0] + trans_dict['_'+x.split('_')[1]])
    temp.index = temp['spot']
    overlap_spot = np.intersect1d(temp['spot'], adata_spatial.obs_names)
    print(temp.shape, overlap_spot.shape)
    
    Total_spot+=temp.shape[0]
    
    adata_spatial.obs.loc[overlap_spot, 'X_3D'] = temp.loc[overlap_spot, 'X_3D']
    adata_spatial.obs.loc[overlap_spot, 'Y_3D'] = temp.loc[overlap_spot, 'Y_3D']
    adata_spatial.obs.loc[overlap_spot, 'Z_3D'] = int(temp_file.split('_')[0][4:])
    
print('-----', Total_spot, '---', adata_spatial.shape[0])

adata_spatial.write('ST_scanpy/cut24h_merge.h5ad')


# cut3d
section_list = ['A1', 'B1', 'C1', 'D1']
adata_list = dict()
for section_it in section_list:
    input_dir = os.path.join('../Data/ST', 'cut3d', section_it)
    
    adata_list[section_it] = sc.read_visium(path=input_dir, library_id='cut3d_'+section_it)
    
    adata_list[section_it].var_names_make_unique()
    sc.pp.calculate_qc_metrics(adata_list[section_it], inplace=True)

adata_spatial = adata_list['A1'].concatenate([adata_list['B1'], adata_list['C1'], adata_list['D1']],
                                  batch_key="library_id",
                                  uns_merge="unique",
                                  batch_categories=[k 
                                                    for d in [adata_list['A1'].uns["spatial"], adata_list['B1'].uns["spatial"], adata_list['C1'].uns["spatial"], adata_list['D1'].uns["spatial"]]
                                                    for k, v in d.items()
                                                   ],)

ST_aligned_dir = '../Data/ST_aligned/cut3d'

Total_spot = 0
trans_dict = {}
trans_dict['_1'] = '-cut3d_A1'
trans_dict['_2'] = '-cut3d_B1'
trans_dict['_3'] = '-cut3d_C1'
trans_dict['_4'] = '-cut3d_D1'
file_list = [x for x in os.listdir(ST_aligned_dir) if x.endswith('.txt')]
for temp_file in file_list:
    temp = pd.read_csv(os.path.join(ST_aligned_dir, temp_file), sep='\t', header=None)
    temp.columns = ['X_3D', 'Y_3D', 'spot']
    temp['spot'] = temp['spot'].map(lambda x: x.split('_')[0] + trans_dict['_'+x.split('_')[1]])
    temp.index = temp['spot']
    overlap_spot = np.intersect1d(temp['spot'], adata_spatial.obs_names)
    print(temp.shape, overlap_spot.shape)
    
    Total_spot+=temp.shape[0]
    
    adata_spatial.obs.loc[overlap_spot, 'X_3D'] = temp.loc[overlap_spot, 'X_3D']
    adata_spatial.obs.loc[overlap_spot, 'Y_3D'] = temp.loc[overlap_spot, 'Y_3D']
    adata_spatial.obs.loc[overlap_spot, 'Z_3D'] = int(temp_file.split('_')[0][4:])
    
print('-----', Total_spot, '---', adata_spatial.shape[0])

adata_spatial.write('ST_scanpy/cut3d_merge.h5ad')


# cut7d
section_list = ['A1', 'B1', 'C1', 'D1']
adata_list = dict()
for section_it in section_list:
    input_dir = os.path.join('../Data/ST', 'cut7d', section_it)
    
    adata_list[section_it] = sc.read_visium(path=input_dir, library_id='cut7d_'+section_it)
    
    adata_list[section_it].var_names_make_unique()
    sc.pp.calculate_qc_metrics(adata_list[section_it], inplace=True)

adata_spatial = adata_list['A1'].concatenate([adata_list['B1'], adata_list['C1'], adata_list['D1']],
                                  batch_key="library_id",
                                  uns_merge="unique",
                                  batch_categories=[k 
                                                    for d in [adata_list['A1'].uns["spatial"], adata_list['B1'].uns["spatial"], adata_list['C1'].uns["spatial"], adata_list['D1'].uns["spatial"]]
                                                    for k, v in d.items()
                                                   ],)

ST_aligned_dir = '../Data/ST_aligned/cut7d'

Total_spot = 0
trans_dict = {}
trans_dict['_1'] = '-cut7d_A1'
trans_dict['_2'] = '-cut7d_B1'
trans_dict['_3'] = '-cut7d_C1'
trans_dict['_4'] = '-cut7d_D1'
file_list = [x for x in os.listdir(ST_aligned_dir) if x.endswith('.txt')]
for temp_file in file_list:
    temp = pd.read_csv(os.path.join(ST_aligned_dir, temp_file), sep='\t', header=None)
    temp.columns = ['X_3D', 'Y_3D', 'spot']
    temp['spot'] = temp['spot'].map(lambda x: x.split('_')[0] + trans_dict['_'+x.split('_')[1]])
    temp.index = temp['spot']
    overlap_spot = np.intersect1d(temp['spot'], adata_spatial.obs_names)
    print(temp.shape, overlap_spot.shape)
    
    Total_spot+=temp.shape[0]
    
    adata_spatial.obs.loc[overlap_spot, 'X_3D'] = temp.loc[overlap_spot, 'X_3D']
    adata_spatial.obs.loc[overlap_spot, 'Y_3D'] = temp.loc[overlap_spot, 'Y_3D']
    adata_spatial.obs.loc[overlap_spot, 'Z_3D'] = int(temp_file.split('_')[0][4:])
    
print('-----', Total_spot, '---', adata_spatial.shape[0])

adata_spatial.write('ST_scanpy/cut7d_merge.h5ad')


# Selected used sections
TimePoint_list = ['cut0h', 'cut6h', 'cut12h', 'cut24h', 'cut3d', 'cut7d']

z_axis = {}
z_axis['cut0h'] = [60, 80, 110, 130, 160, 180, 210, 230, 250, 280, 320, 360]
z_axis['cut6h'] = [6, 12, 18, 26, 34, 38, 42, 46, 54, 60, 64, 72]
z_axis['cut12h'] = [40, 60, 80, 100, 140, 180, 220, 280, 310, 340, 370, 390]
z_axis['cut24h'] = [50, 70, 110, 130, 190, 210, 240, 300, 320, 340, 360, 380]
z_axis['cut3d'] = [131, 196, 327, 393, 436, 480, 589, 633, 676, 720, 764, 807]
z_axis['cut7d'] = [60, 120, 180, 240, 280, 340, 380, 440, 520, 580, 660, 780]

for temp_timepoint in TimePoint_list:
    adata = sc.read('ST_scanpy/'+temp_timepoint+'_merge.h5ad')
    adata = adata[adata.obs['Z_3D'].isin(range(1, 13)),]
    
    adata.obs['Aligned_X'] = adata.obs['X_3D']
    adata.obs['Aligned_Y'] = adata.obs['Y_3D']
    adata.obs['Aligned_Z'] = adata.obs['Z_3D'].map(dict(zip(range(1,13), z_axis['cut0h'])))
    
    adata.write('ST_scanpy/'+temp_timepoint+'_3D.h5ad')