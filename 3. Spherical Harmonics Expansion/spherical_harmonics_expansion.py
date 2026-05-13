#Author: Lucas Philipp
#Code heavily relies on: https://github.com/AllenCell/aics-shparam

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from aicsshparam import shtools, shparam
from tifffile import imread
import vtk
from vtk.util import numpy_support as vtknp
import tifffile
from skimage import feature
from skimage import morphology
import matplotlib.colors as mcolors
from matplotlib.font_manager import FontProperties
import napari
from PIL import Image
from matplotlib.colors import ListedColormap, to_rgb
import random
from math import radians, cos, sin
from scipy import stats
from scipy.ndimage import binary_erosion

def PCA_coords_to_shape(pca: PCA, PCA_coords: np.ndarray, lmax: int, save_directory: str, save_as_tiff: bool = False, save_cross_sections: bool = False):
    '''
    Visualize PCA coordinates as shape. See Fig. S5 b.
    
    Parameters:
    - PCA_coords (np.ndarray) (N, 2) First column is PC1 coordinate. Second column is PC2 coordinate. N rows are different shapes to be created.
    - save_directory
    - Set save_as_tiff = True to save a 3D binary z-stack (TIFF) of shape. 
    - Set save_cross_sections = True to save xy, yz, and xz cross-sections of shape.
    - 2*(lmax+1)^2 is the number of spherical harmonics terms used for shape reconstruction
    '''
    
    def normalize_to_uint8(img):
        img = img.astype(np.float32)
        img -= img.min()
        if img.max() != 0:
            img /= img.max()
        return (img * 255).astype(np.uint8)
    
    num_sample_shapes = PCA_coords.shape[0]
    
    shape_params = np.zeros((num_sample_shapes, 2*(lmax+1)**2))
    
    for i in range(num_sample_shapes):
        shape_params[i, :] = pca.inverse_transform(tuple(PCA_coords[i]))
    
    shape_params_reshape = np.zeros((2, lmax, lmax, num_sample_shapes))
    
    for d in range(num_sample_shapes):
        count = 0
        for i in range(2):
            for j in range(lmax):
                for k in range(lmax):
                    shape_params_reshape[i, j, k, d] = shape_params[d, count]
                    count = count + 1
                
    for i in range(num_sample_shapes):
        mesh, _ = shtools.get_reconstruction_from_coeffs(shape_params_reshape[:, :, :, i])
        coords = vtknp.vtk_to_numpy(mesh.GetPoints().GetData())
    
        rmin = (coords.min(axis=0) - 0.5).astype(int)
        rmax = (coords.max(axis=0) + 0.5).astype(int)

        pad = 6
    
        w = int(2 * pad + (rmax[0] - rmin[0]))
        h = int(2 * pad + (rmax[1] - rmin[1]))
        d = int(2 * pad + (rmax[2] - rmin[2]))
    
        imagedata = vtk.vtkImageData()
        imagedata.SetDimensions([w, h, d])
        imagedata.SetExtent(0, w - 1, 0, h - 1, 0, d - 1)
        imagedata.SetOrigin(rmin - pad)
        imagedata.AllocateScalars(vtk.VTK_UNSIGNED_CHAR, 1)
    
        imagedata.GetPointData().GetScalars().FillComponent(0, 1)
        
        voxelized = shtools.voxelize_mesh(imagedata=imagedata, shape=(d, h, w), mesh=mesh, origin=rmin - pad)
    
        voxelized = voxelized.astype('int8')
        if save_as_tiff == True:
            tifffile.imwrite(save_directory + '/shape {}.tiff'.format(i), voxelized, bigtiff=True)

        z, y, x = voxelized.shape
        xy_section = voxelized[z // 2, :, :]
        yz_section = voxelized[:, :, x // 2]
        xz_section = voxelized[:, y // 2, :]
        
        xy_section = feature.canny(normalize_to_uint8(xy_section), sigma=1.0)
        yz_section = feature.canny(normalize_to_uint8(yz_section), sigma=1.0)
        xz_section = feature.canny(normalize_to_uint8(xz_section), sigma=1.0)
        
        xy_section = morphology.dilation(xy_section, morphology.square(pad))
        yz_section = morphology.dilation(yz_section, morphology.square(pad))
        xz_section = morphology.dilation(xz_section, morphology.square(pad))
        
        xy_img = Image.fromarray(normalize_to_uint8(xy_section))
        yz_img = Image.fromarray(normalize_to_uint8(yz_section))
        xz_img = Image.fromarray(normalize_to_uint8(xz_section))
    
        if save_cross_sections == True:
            xy_img.save(save_directory + '/xy/{}.png'.format(i))
            yz_img.save(save_directory + '/yz/{}.png'.format(i))
            xz_img.save(save_directory + '/xz/{}.png'.format(i))
       
        print(i)
        
def estimate_surface_area_from_binary(volume: np.ndarray, voxel_size_nm: float) -> float:
    """
    Estimate surface area of a binary 3D volume by counting exposed voxel faces.

    Parameters:
    - volume: 3D numpy array of bool or 0/1 (True=object, False=background)
    - voxel_size_nm: size of one voxel edge in nanometers (assumes isotropic voxels)

    Returns:
    - surface area in square nanometers (nm²)
    """
    if volume.dtype != bool:
        volume = volume.astype(bool)

    # Calculate face area (each voxel face is square)
    face_area = voxel_size_nm ** 2

    # Erode to find interior voxels
    eroded = binary_erosion(volume)

    # Surface voxels: object voxels that get removed by erosion
    surface_voxels = volume & (~eroded)

    # Pad volume to check neighbors at edges safely
    padded = np.pad(volume, pad_width=1, mode='constant', constant_values=0)

    # Directions to check neighbors (6-connected)
    neighbors_offsets = [(1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)]

    exposed_faces = 0
    for dx, dy, dz in neighbors_offsets:
        neighbors = padded[
            1 + dx : 1 + dx + volume.shape[0],
            1 + dy : 1 + dy + volume.shape[1],
            1 + dz : 1 + dz + volume.shape[2]
        ]
        # Count faces where surface voxel has background neighbor
        exposed_faces += np.sum(surface_voxels & (~neighbors))

    surface_area = exposed_faces * face_area
    return surface_area


#get images and angles and perform spherical harmonic expansion
def spherical_harmonic_expansion(ROI_start,ROI_end,image_path,angle_csv_path,cell_ID,species_name,lmax):
    """
    Represent shape in terms of spherical harmonic coefficients. Store other annotations alongside coefficients in dataframe.
    
    Parameters:
    - ROI_start: int, starting index/name of image file. e.g. ROI 1.tiff
    - ROI_end: int, ending index/name of image file. e.g. ROI 100.tiff
    - image_path: folder with segmented chromosome TIFF files
    - angle_csv_path: location of surface ridge angle measurements
    - cell_ID
    - species_name
    - 2*(lmax+1)^2 is the number of spherical harmonics terms used for shape reconstruction
    
    assumes 4x4x4 nm voxel size. 
    rejects chromosome if center of mass is outside chromosome.
    rejects chromosome if mean square reconstruction error > 10.

    Returns:
    - df_coeffs: pandas dataframe of spherical harmonic coefficients for all shapes

    - surface_ridges: boolean TRUE or FALSE
    - cell_ID
    - species
    - ROI: int
    - volume: [nm^3]
    - surface_area: [nm^2]
    - length: longest axis[nm]
    - width: longest orthogonal axis to length [nm}
    - delta_theta: difference of surface ridge angles: \theta_{back} - \theta_{front}
    - LH_or_RH: helical handedness if surface ridges present. NaN if not. either 'left handed', 'right handed', or '|Δθ| < 5°' (flat discs)
    - aspect_ratio: []
    """
    
    images = []
    ROI = []
    flag = []
    for i in range(ROI_start,ROI_end+1):
        images.append(imread(image_path+'ROI {}.tiff'.format(i)))
        ROI.append(i)
    #load whether surface ridge angle was extracted
    angles = pd.read_csv(angle_csv_path)
    angles.drop(angles.columns[[0,2,3,4,5]],axis=1, inplace=True)
    angles['Label'] = values(angles.loc[:,"Label"]).astype(str).astype(int)
    surface_ridges = []
    delta_theta = []
    LH_or_RH = []
    for i in range(ROI_start,ROI_end+1):
        surface_ridges.append(angles.Label.eq(i).any()) #was an angle extracted for this chromosome?        
        temp = angles.index[angles['Label'] == i]
        if len(temp)>0:
            delta_theta.append(angles.Angle[temp[0]]-angles.Angle[temp[1]])
            if (angles.Angle[temp[0]]-angles.Angle[temp[1]])<-5:
                LH_or_RH.append('left handed')
            elif (angles.Angle[temp[0]]-angles.Angle[temp[1]])>5:
                LH_or_RH.append('right handed')
            else:
                LH_or_RH.append('|Δθ| < 5°')
        else:
            delta_theta.append(np.nan)
            LH_or_RH.append(np.nan)
            
    #spherical harmonic expansion
    df_coeffs_list = []
    volume = [] #[nm^3]
    aspect_ratio = []
    surface_area = [] #[nm^2]
    length = [] #[nm]
    width = [] #[nm]
    
    violation = [] #list of ROIs to reject
    for i in range(len(images)):
        (coeffs, grid_rec), (_, _, grid_input, _), flag = shparam.get_shcoeffs(image=images[i], lmax=lmax)
        if flag == True:
            violation.append(i)
        error = shtools.get_reconstruction_error(grid_rec,grid_input)
        if error > 25:
            violation.append(i)
        volume.append((4**3)*np.sum(images[i]/np.max(images[i]))) #[nm^3] assuming 4nmx4nm4xnm voxels 
        surface_area.append(estimate_surface_area_from_binary(images[i]/np.max(images[i]), 4)) #[nm^2] assuming 4nmx4nm4xnm voxels
        aspect_ratio.append(np.size(images[i],axis=0)/np.size(images[i],axis=1))     
        length.append(images[i].shape[0]*4) #[nm] assuming 4nmx4nm4xnm voxels
        width.append(images[i].shape[1]*4) #[nm] assuming 4nmx4nm4xnm voxels
        df_coeffs_list.append(coeffs)
        print(i)
    df_coeffs = pd.DataFrame(df_coeffs_list)
    #add surface_ridges label
    df_coeffs['surface_ridges'] = surface_ridges
    #add cell label
    df_coeffs['cell'] = np.ones(len(images))*cell_ID
    #add species label
    species = list()
    for i in range(len(images)):
        species.append(species_name)
    df_coeffs['species'] = species
    #add ROI label
    df_coeffs['ROI']=ROI
    #add volume label
    df_coeffs['volume']=volume
    #add surface area
    df_coeffs['surface_area']=surface_area
    #add length
    df_coeffs['length']=length
    #add width
    df_coeffs['width']=width
    #add delta_theta label
    df_coeffs['delta_theta']=delta_theta
    #add LH_or_RH label
    df_coeffs['LH_or_RH']=LH_or_RH
    #add aspect ratio label
    df_coeffs['aspect_ratio']=aspect_ratio
    #return df_coeffs
    
    #remove duplicate ROI indices
    violation = list(dict.fromkeys(violation))
    df_coeffs = df_coeffs.drop(index=violation)
    print(violation)
    return df_coeffs

def values(x):
    return x.str.extract(r'([0-9]+)')

def generate_shades(base_color, n_shades):
    """Generate n_shades from light to dark of a base RGB color."""
    shades = []
    for i in range(n_shades):
        factor = 0.5 + 1.5 * (i / (n_shades - 1))  # avoid pure white
        shaded = tuple(np.clip(factor * np.array(base_color), 0, 1))
        shades.append(shaded)
    return shades

# Compute spherical harmonics coefficients of shape and store them in a pandas dataframe.
lmax = 40 #2*(lmax+1)^2 is the number of expansion terms

microadriaticum_cell1 = spherical_harmonic_expansion(2,114,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/Cell 1/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/SM Cell 1.csv',1,'microadriaticum',lmax)
microadriaticum_cell2 = spherical_harmonic_expansion(1,100,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/Cell 2/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/SM Cell 2.csv',2,'microadriaticum',lmax)
microadriaticum_cell3 = spherical_harmonic_expansion(2,105,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/Cell 3/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/SM Cell 3.csv',3,'microadriaticum',lmax)

pilosum_cell1 = spherical_harmonic_expansion(4,102,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 1/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/SP Cell 1.csv',1,'pilosum',lmax)
pilosum_cell2 = spherical_harmonic_expansion(22,124,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 2/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/SP Cell 2.csv',2,'pilosum',lmax)
pilosum_cell3 = spherical_harmonic_expansion(27,124,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 3/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/SP Cell 3.csv',3,'pilosum',lmax)

minutum_cell1 = spherical_harmonic_expansion(1,26,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Breviolum minutum/Cell 1/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Breviolum minutum/SMin Cell 1.csv',1,'minutum',lmax)
minutum_cell2 = spherical_harmonic_expansion(2,35,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Breviolum minutum/Cell 2/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Breviolum minutum/SMin Cell 2.csv',2,'minutum',lmax)
minutum_cell3 = spherical_harmonic_expansion(1,33,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Breviolum minutum/Cell 3/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Breviolum minutum/SMin Cell 3.csv',3,'minutum',lmax)

nutricula_cell1 = spherical_harmonic_expansion(2,162,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/brandtodinium cell 1 chromosomes 4nm sampling/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/BN Cell 1.csv',1,'nutricula',lmax)
nutricula_cell2 = spherical_harmonic_expansion(1,100,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/brandtodinium cell 2 chromosomes 4nm sampling/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/BN Cell 2.csv',2,'nutricula',lmax)
nutricula_cell3 = spherical_harmonic_expansion(4,327,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/brandtodinium cell 3 chromosomes 4nm sampling/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/BN Cell 3.csv',3,'nutricula',lmax)

cohnii_cell1 = spherical_harmonic_expansion(2,118,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/Cell 1/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/CC cell 1.csv',1,'cohnii',lmax)
cohnii_cell2 = spherical_harmonic_expansion(2,187,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/Cell 2/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/CC cell 2.csv',2,'cohnii',lmax)
cohnii_cell3 = spherical_harmonic_expansion(2,135,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/Cell 3/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/CC cell 3.csv',3,'cohnii',lmax)

tyrrhenica_cell1 = spherical_harmonic_expansion(2,113,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Ensiculifera tyrrhenica/ensiculifera tyrrhenica chromosomes 4nm sampling/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Ensiculifera tyrrhenica/ET Cell 1.csv',1,'tyrrhenica',lmax)

ross_sea_dinoflagellate_cell2 = spherical_harmonic_expansion(1,165,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Ross Sea Dinoflagellate/Ross Sea Dinoflagellate Cell 2 Chromosomes 4nm voxels/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Ross Sea Dinoflagellate/RSD Cell 2.csv',2,'ross sea dinoflagellate',lmax)
ross_sea_dinoflagellate_cell3 = spherical_harmonic_expansion(1,130,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Ross Sea Dinoflagellate/Ross Sea Dinoflagellate Cell 3 Chromosomes 4nm voxels/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Ross Sea Dinoflagellate/RSD Cell 3.csv',3,'ross sea dinoflagellate',lmax)

SHE_all = pd.concat([microadriaticum_cell1, microadriaticum_cell2, microadriaticum_cell3, pilosum_cell1, pilosum_cell2, pilosum_cell3, minutum_cell1, minutum_cell2, minutum_cell3, nutricula_cell1, nutricula_cell2, nutricula_cell3, cohnii_cell1, cohnii_cell2, cohnii_cell3, tyrrhenica_cell1, ross_sea_dinoflagellate_cell2, ross_sea_dinoflagellate_cell3], ignore_index=True, sort=False)

#uncomment below for species specific PCA (e.g. B. nutricula)
#SHE_all = pd.concat([nutricula_cell1, nutricula_cell2, nutricula_cell3], ignore_index=True, sort=False)

print("Total Number of Chromosomes Reconstructed")
print(SHE_all.shape[0])

###Uncomment code block below to ensure PCA is fed an equal # of chromosomes per species
### start code block
# num_chroms = []
# num_chroms.append((SHE_all['species'] == 'microadriaticum').sum())
# num_chroms.append((SHE_all['species'] == 'pilosum').sum())
# num_chroms.append((SHE_all['species'] == 'minutum').sum())
# num_chroms.append((SHE_all['species'] == 'nutricula').sum())
# num_chroms.append((SHE_all['species'] == 'cohnii').sum())
# num_chroms.append((SHE_all['species'] == 'tyrrhenica').sum())
# num_chroms.append((SHE_all['species'] == 'ross sea dinoflagellate').sum())

# match_rows = SHE_all[SHE_all['species'] == 'microadriaticum']
# rows_to_drop = match_rows.sample(n=num_chroms[0]-min(num_chroms))
# SHE_all = SHE_all.drop(rows_to_drop.index)

# match_rows = SHE_all[SHE_all['species'] == 'pilosum']
# rows_to_drop = match_rows.sample(n=num_chroms[1]-min(num_chroms))
# SHE_all = SHE_all.drop(rows_to_drop.index)

# match_rows = SHE_all[SHE_all['species'] == 'minutum']
# rows_to_drop = match_rows.sample(n=num_chroms[2]-min(num_chroms))
# SHE_all = SHE_all.drop(rows_to_drop.index)

# match_rows = SHE_all[SHE_all['species'] == 'nutricula']
# rows_to_drop = match_rows.sample(n=num_chroms[3]-min(num_chroms))
# SHE_all = SHE_all.drop(rows_to_drop.index)

# match_rows = SHE_all[SHE_all['species'] == 'cohnii']
# rows_to_drop = match_rows.sample(n=num_chroms[4]-min(num_chroms))
# SHE_all = SHE_all.drop(rows_to_drop.index)

# match_rows = SHE_all[SHE_all['species'] == 'tyrrhenica']
# rows_to_drop = match_rows.sample(n=num_chroms[5]-min(num_chroms))
# SHE_all = SHE_all.drop(rows_to_drop.index)

# match_rows = SHE_all[SHE_all['species'] == 'ross sea dinoflagellate']
# rows_to_drop = match_rows.sample(n=num_chroms[6]-min(num_chroms))
# SHE_all = SHE_all.drop(rows_to_drop.index)

# SHE_all = SHE_all.reset_index(drop=True)
### end code block

#write chromosome lengths/widths to csv file
selected_data = SHE_all[['length', 'width']]
selected_data.to_csv('chromosome_length_width.csv', index=False)

# Vizualize spherical harmonic coefficients dataframe
with pd.option_context('display.max_rows', 5, 'display.max_columns', 5):
    print(SHE_all.shape[0])
    
# initialize PCA object
pca_all = PCA(n_components=2)
#perform PCA on SHE coefficients only. annotations are dropped
trans = pca_all.fit_transform(SHE_all.drop(columns=['surface_ridges','cell','species','ROI','volume','length','width','surface_area','delta_theta','LH_or_RH','aspect_ratio']))

#percentage of total variance explained by top PCs
print("Variance Explained by PC1")
print(f"{pca_all.explained_variance_ratio_[0] * 100:.1f}%")

print("Variance Explained by PC2")
print(f"{pca_all.explained_variance_ratio_[1] * 100:.1f}%")

df_trans = pd.DataFrame(trans)
df_trans.columns = ['PC1', 'PC2']

#add annotations to PCA result. chromosome order is preserved.
df_trans['surface_ridges'] = SHE_all.surface_ridges
df_trans['cell'] = SHE_all.cell
df_trans['species'] = SHE_all.species
df_trans['ROI'] = SHE_all.ROI
df_trans['volume'] = SHE_all.volume
df_trans['surface_area'] = SHE_all.surface_area
df_trans['delta_theta'] = SHE_all.delta_theta
df_trans['LH_or_RH'] = SHE_all.LH_or_RH
df_trans['aspect_ratio'] = SHE_all.aspect_ratio
df_trans['width'] = SHE_all.width
df_trans['length'] = SHE_all.length

#Main PCA plot
#no color (light grey dots, black outline)
with pd.option_context('display.max_rows', 5, 'display.max_columns', 5):
    print(df_trans)
    # Scatter plot to show how similar shapes are grouped together.
fig, ax = plt.subplots(1,1, figsize=(8,8))
ax.scatter(df_trans.PC1, df_trans.PC2, facecolors='lightgrey', edgecolors='black', s=200)
ax.tick_params(axis='both', which='major', labelsize=22.5)
plt.xlabel(f'PC1 ({pca_all.explained_variance_ratio_[0]:.1%} Variance Explained)', fontsize=25)
plt.ylabel(f'PC2 ({pca_all.explained_variance_ratio_[1]:.1%} Variance Explained)', fontsize=25)
left, right = plt.xlim()
up, down = plt.ylim()
#plt.xlim((left,100)) 
left, right = plt.xlim()
plt.show()

#plot PC2 vs PC3
# fig, ax = plt.subplots(1,1, figsize=(8,8))
# ax.scatter(df_trans.PC2, df_trans.PC3, facecolors='lightgrey', edgecolors='black', s=200)
# ax.tick_params(axis='both', which='major', labelsize=22.5)
# plt.xlabel('PC2', fontsize=25)
# plt.ylabel('PC3', fontsize=25)
# left, right = plt.xlim()
# up, down = plt.ylim()
# #plt.xlim((left,100)) 
# left, right = plt.xlim()
# plt.show()

#Surface ridges? Yes/No
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
# Define manual colors for the two groups
colors = {True: '#333333', False: 'white'}
for label, df_label in df_trans.groupby('surface_ridges'):
    ax.scatter(
        df_label.PC1, 
        df_label.PC2, 
        label=label, 
        s=200,
        color=colors[label],
        edgecolor='black',   # Ensures white dots are visible
        linewidth=0.5
    )
ax.tick_params(axis='both', which='major', labelsize=22.5)

plt.xlabel('PC1', fontsize=30)
plt.ylabel('PC2', fontsize=30)
plt.legend(['No Surface Ridges', 'Yes Surface Ridges'], loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=30, frameon=False)
plt.show()

#color by volume
fig, ax = plt.subplots(1,1, figsize=(8,8))
sc = ax.scatter(df_trans.PC1, df_trans.PC2, c=df_trans['volume'], cmap='hsv', edgecolors='black', norm='log', s=200)
ax.tick_params(axis='both', which='major', labelsize=22.5)
#cbar = fig.colorbar(sc, ax=ax)
#cbar.set_label('Volume [nm³]', fontsize=30)
#cbar.ax.tick_params(labelsize=16)
plt.xlabel('PC1', fontsize=30)
plt.ylabel('PC2', fontsize=30)
plt.xlim((left, right))
plt.ylim((up, down))
plt.show()

slope, intercept, r_value, p_value, std_err = stats.linregress(df_trans.PC1, df_trans['volume'])
print(f'Correlation between PC1 and volume: R² = {r_value**2:.4f}')

#setup aspect-ratio plot colormap
vmin = 1
vmax = 3
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

#define base colors
base_hues = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
]
base_colors = [to_rgb(color) for color in base_hues]

#generate 5 shades per base color
custom_colors = []
for base in base_colors:
    custom_colors.extend(generate_shades(base, 5))

# Create the new colormap
custom_cmap = ListedColormap(custom_colors)

#color by aspect ratio
fig, ax = plt.subplots(1,1, figsize=(8,8))
sc = ax.scatter(df_trans.PC1, df_trans.PC2, c=df_trans['aspect_ratio'], norm=norm, cmap=custom_cmap, edgecolors='black', s=200)
ax.tick_params(axis='both', which='major', labelsize=22.5)
#cbar = fig.colorbar(sc, ax=ax)
#cbar.set_label('Aspect ratio', fontsize = 22.5)
plt.xlabel('PC1', fontsize=30)
plt.ylabel('PC2', fontsize=30)
plt.xlim((left, right))
plt.ylim((up, down))
#cbar.set_ticks([1, 1.5, 2, 2.5, 3])
#cbar.set_ticklabels(['1', '1.5', '2', '2.5', '>3'])
plt.show()

slope, intercept, r_value, p_value, std_err = stats.linregress(df_trans.PC2, df_trans['aspect_ratio'])
print(f'Correlation between PC2 and aspect-ratio: R² = {r_value**2:.4f}')

#color by delta_theta
vmin = -90
vmax = 90
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

fig, ax = plt.subplots(1,1, figsize=(8,8))
sc = ax.scatter(df_trans.PC1, df_trans.PC2, c=df_trans['delta_theta'], norm=norm, cmap='coolwarm', edgecolors='black', s=200)
ax.tick_params(axis='both', which='major', labelsize=22.5)
cbar = fig.colorbar(sc, ax=ax, ticks=[-90, -60, -30, 0, 30, 60, 90])
cbar.set_label('Δθ [°]', fontsize=22.5)
cbar.ax.tick_params(labelsize=22.5)
plt.xlabel('PC1', fontsize=30)
plt.ylabel('PC2', fontsize=30)
plt.xlim((left, right))
plt.ylim((up, down))
plt.show()

#color by LH_or_RH
with pd.option_context('display.max_rows', 5, 'display.max_columns', 5):
    print(df_trans)
    # Scatter plot to show how similar shapes are grouped together.
fig, ax = plt.subplots(1,1, figsize=(8,8))
for label, df_label in df_trans.groupby('LH_or_RH'):
    ax.scatter(df_label.PC1, df_label.PC2, label=label, edgecolors='black', s=200)
ax.tick_params(axis='both', which='major', labelsize=22.5)
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=22.5, frameon=False)
plt.xlabel('PC1', fontsize=30)
plt.ylabel('PC2', fontsize=30)
plt.xlim((left, right))
plt.ylim((up, down))
plt.show()

#color by species
colors = ['magenta', 'blue', 'orange', 'red', 'purple','grey', 'black']
count = 0
fig, ax = plt.subplots(1,1, figsize=(8,8))
for label, df_label in df_trans.groupby('species'):
    ax.scatter(df_label.PC1, df_label.PC2, label=label, s=100, color = colors[count], alpha = 0.35)
    count = count + 1
ax.tick_params(axis='both', which='major', labelsize=22.5)
ax.set_xticks([-50, -25, 0, 25, 50, 75])
italic_font = FontProperties(style='italic')
legend = plt.legend(['C. cohnii','S. microadriaticum','S. minutum','B. nutricula','S. pilosum','E. tyrrhenica', 'K. sp.'],
                    loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=22.5, prop=italic_font, frameon=False)
for handle in legend.legend_handles:
    handle.set_alpha(1)
plt.xlabel('PC1', fontsize=30)
plt.ylabel('PC2', fontsize=30)
x_lim = ax.get_xlim()
y_lim = ax.get_ylim()
left, right = plt.xlim()
ax.grid(False)
plt.savefig('specis_PCA.svg', format='svg', bbox_inches='tight')
plt.show()

#species centroids
from matplotlib.font_manager import FontProperties

colors = ['magenta', 'blue', 'orange', 'red', 'purple', 'grey', 'black']
italic_font = FontProperties(style='italic')

fig, ax = plt.subplots(1, 1, figsize=(8, 8))

# Loop through each species and plot its centroid
for count, (label, df_label) in enumerate(df_trans.groupby('species')):
    centroid_x = df_label.PC1.mean()
    centroid_y = df_label.PC2.mean()

    ax.scatter(
        centroid_x, centroid_y,
        s=250,
        color=colors[count],
        edgecolors='black',
        linewidth=1.5,
        alpha=1.0,
        label=label
    )

ax.tick_params(axis='both', which='major', labelsize=22.5)
ax.set_xticks([-50, -25, 0, 25, 50, 75])
ax.set_xlim(x_lim)
ax.set_ylim(y_lim)

# plt.legend(['C. cohnii', 'S. microadriaticum', 'S. minutum', 'B. nutricula', 'S. pilosum', 'E. tyrrhenica', 'K. sp.'],
#     loc='upper left',
#     bbox_to_anchor=(1.05, 1),
#     fontsize=22.5,
#     prop=italic_font,
#     frameon=False
# )

plt.xlabel('PC1', fontsize=30)
plt.ylabel('PC2', fontsize=30)
plt.show()

#color by surface_area/volume
fig, ax = plt.subplots(1, 1, figsize=(8, 8))
for label, df_label in df_trans.groupby('LH_or_RH'):
    ratio = df_label['surface_area'] / df_label['volume']
    
    scatter = ax.scatter(
        df_label.PC1,
        df_label.PC2,
        label=label,
        c=ratio,            # color by surface_area / volume
        cmap='viridis',     # choose colormap you like
        edgecolors='black',
        s=200
    )

# Add colorbar for the ratio
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('Surface Area / Volume [nm]', fontsize=22.5)

ax.tick_params(axis='both', which='major', labelsize=22.5)
plt.xlabel('PC1', fontsize=30)
plt.ylabel('PC2', fontsize=30)
plt.show()

##############################################################################
#repeat PCA but only using data from a single species. compare species specific PC1 to global PC1

#get global eigenvector
pc1_all = pca_all.components_[0]
pc2_all = pca_all.components_[1]

#microadriaticum
SHE_microadriaticum = pd.concat([microadriaticum_cell1, microadriaticum_cell2, microadriaticum_cell3], ignore_index=True, sort=False)
pca_microadriaticum = PCA(n_components=2)
pca_microadriaticum.fit_transform(SHE_microadriaticum.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_microadriaticum = pca_microadriaticum.components_[0]
pc2_microadriaticum = pca_microadriaticum.components_[1]

#pilosum
SHE_pilosum = pd.concat([pilosum_cell1, pilosum_cell2, pilosum_cell3], ignore_index=True, sort=False)
pca_pilosum = PCA(n_components=2)
pca_pilosum.fit_transform(SHE_pilosum.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_pilosum = pca_pilosum.components_[0]
pc2_pilosum = pca_pilosum.components_[1]

#minutum
SHE_minutum = pd.concat([minutum_cell1, minutum_cell2, minutum_cell3], ignore_index=True, sort=False)
pca_minutum = PCA(n_components=2)
pca_minutum.fit_transform(SHE_minutum.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_minutum = pca_minutum.components_[0]
pc2_minutum = pca_minutum.components_[1]

#cohnii
SHE_cohnii = pd.concat([cohnii_cell1, cohnii_cell2, cohnii_cell3], ignore_index=True, sort=False)
pca_cohnii = PCA(n_components=2)
pca_cohnii.fit_transform(SHE_cohnii.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_cohnii = pca_cohnii.components_[0]
pc2_cohnii = pca_cohnii.components_[1]

#nutricula
SHE_nutricula = pd.concat([nutricula_cell1, nutricula_cell2, nutricula_cell3], ignore_index=True, sort=False)
pca_nutricula = PCA(n_components=2)
pca_nutricula.fit_transform(SHE_nutricula.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_nutricula = pca_nutricula.components_[0]
pc2_nutricula = pca_nutricula.components_[1]

#tyrrhenica
SHE_tyrrhenica = pd.concat([tyrrhenica_cell1], ignore_index=True, sort=False)
pca_tyrrhenica = PCA(n_components=2)
pca_tyrrhenica.fit_transform(SHE_tyrrhenica.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_tyrrhenica = pca_tyrrhenica.components_[0]
pc2_tyrrhenica = pca_tyrrhenica.components_[1]

#ross sea dinoflagellate
SHE_ross_sea_dinoflagellate = pd.concat([ross_sea_dinoflagellate_cell2, ross_sea_dinoflagellate_cell3], ignore_index=True, sort=False)
pca_ross_sea_dinoflagellate = PCA(n_components=2)
pca_ross_sea_dinoflagellate.fit_transform(SHE_ross_sea_dinoflagellate.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_ross_sea_dinoflagellate = pca_ross_sea_dinoflagellate.components_[0]
pc2_ross_sea_dinoflagellate = pca_ross_sea_dinoflagellate.components_[1]

#compute angle between species-specific and global eigenvectors
all_angles_PC1 = []
all_angles_PC2 = []

#microadriaticum
cos_angle = np.dot(pc1_all, pc1_microadriaticum) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_microadriaticum))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 microadriaticum: {angle_deg:.2f} degrees")
all_angles_PC1.append(angle_deg)

cos_angle = np.dot(pc2_all, pc2_microadriaticum) / (np.linalg.norm(pc2_all) * np.linalg.norm(pc2_microadriaticum))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC2 Global and PC2 microadriaticum: {angle_deg:.2f} degrees")
all_angles_PC2.append(angle_deg)

#pilosum
cos_angle = np.dot(pc1_all, pc1_pilosum) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_pilosum))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 pilosum: {angle_deg:.2f} degrees")
all_angles_PC1.append(angle_deg)

cos_angle = np.dot(pc2_all, pc2_pilosum) / (np.linalg.norm(pc2_all) * np.linalg.norm(pc2_pilosum))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC2 Global and PC2 pilosum: {angle_deg:.2f} degrees")
all_angles_PC2.append(angle_deg)

#minutum
cos_angle = np.dot(pc1_all, pc1_minutum) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_minutum))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 minutum: {angle_deg:.2f} degrees")
all_angles_PC1.append(angle_deg)

cos_angle = np.dot(pc2_all, pc2_minutum) / (np.linalg.norm(pc2_all) * np.linalg.norm(pc2_minutum))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC2 Global and PC2 minutum: {angle_deg:.2f} degrees")
all_angles_PC2.append(angle_deg)

#cohnii
cos_angle = np.dot(pc1_all, pc1_cohnii) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_cohnii))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 cohnii: {angle_deg:.2f} degrees")
all_angles_PC1.append(angle_deg)

cos_angle = np.dot(pc2_all, pc2_cohnii) / (np.linalg.norm(pc2_all) * np.linalg.norm(pc2_cohnii))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC2 Global and PC2 cohnii: {angle_deg:.2f} degrees")
all_angles_PC2.append(angle_deg)

#nutricula
cos_angle = np.dot(pc1_all, pc1_nutricula) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_nutricula))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 nutricula: {angle_deg:.2f} degrees")
all_angles_PC1.append(angle_deg)

cos_angle = np.dot(pc2_all, pc2_nutricula) / (np.linalg.norm(pc2_all) * np.linalg.norm(pc2_nutricula))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC2 Global and PC2 nutricula: {angle_deg:.2f} degrees")
all_angles_PC2.append(angle_deg)

#tyrrhenica
cos_angle = np.dot(pc1_all, pc1_tyrrhenica) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_tyrrhenica))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 tyrrhenica: {angle_deg:.2f} degrees")
all_angles_PC1.append(angle_deg)

cos_angle = np.dot(pc2_all, pc2_tyrrhenica) / (np.linalg.norm(pc2_all) * np.linalg.norm(pc2_tyrrhenica))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC2 Global and PC2 tyrrhenica: {angle_deg:.2f} degrees")
all_angles_PC2.append(angle_deg)

#ross sea dinoflagellate
cos_angle = np.dot(pc1_all, pc1_ross_sea_dinoflagellate) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_ross_sea_dinoflagellate))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 K. sp.: {angle_deg:.2f} degrees")
all_angles_PC1.append(angle_deg)

cos_angle = np.dot(pc2_all, pc2_ross_sea_dinoflagellate) / (np.linalg.norm(pc2_all) * np.linalg.norm(pc2_ross_sea_dinoflagellate))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC2 Global and PC2 K. sp.: {angle_deg:.2f} degrees")
all_angles_PC2.append(angle_deg)

###########################################
#plot \theta_{species}
angle_rad = [radians(a) for a in all_angles_PC1]

# Calculate unit vector components
x = np.cos(angle_rad)
y = np.sin(angle_rad)

# Plot setup
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_aspect('equal')

# Limit to top-right quadrant
ax.set_xlim(0, 1.1)
ax.set_ylim(0, 1.2)

# Color list (extend as needed)
colors = ['blue', 'purple', 'orange', 'pink', 'red', 'grey', 'black']

# Draw unit vectors as arrows
for i in range(len(angle_rad)):
    ax.arrow(0, 0, x[i], y[i],
             head_width=0.03, head_length=0.05,
             fc=colors[i], ec=colors[i])

# Radial guide lines for 0°, 90°, 30°, and 60°
def draw_angle_line(deg, ax, length=1.1, label=None):
    rad = radians(deg)
    ax.plot([0, cos(rad)], [0, sin(rad)], color='black', linestyle='dotted')
    if label:
        ax.text(cos(rad) * length, sin(rad) * length, label,
                fontsize=16, ha='left', va='bottom')

draw_angle_line(0, ax, label='0°')
draw_angle_line(30, ax, label='30°')
draw_angle_line(60, ax, label='60°')
draw_angle_line(90, ax, label='90°')

for spine in ax.spines.values():
    spine.set_visible(False)

ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
plt.title('Angle between species-specific PC1 and global PC1', fontsize=14)
plt.ylabel('Orthogonal shape variation', fontsize=16)
plt.xlabel('Identical shape variation', fontsize=16)
plt.show()

######################################################################################
#repeat PCA but only using data from a single cells. compare cell specific PC1 to global PC1

#get global eigenvector
pc1_all = pca_all.components_[0]
pc2_all = pca_all.components_[1]

#microadriaticum
SHE_microadriaticum_cell1 = pd.concat([microadriaticum_cell1], ignore_index=True, sort=False)
pca_microadriaticum_cell1 = PCA(n_components=2)
pca_microadriaticum_cell1.fit_transform(SHE_microadriaticum_cell1.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_microadriaticum_cell1 = pca_microadriaticum_cell1.components_[0]
pc2_microadriaticum_cell1 = pca_microadriaticum_cell1.components_[1]

SHE_microadriaticum_cell2 = pd.concat([microadriaticum_cell2], ignore_index=True, sort=False)
pca_microadriaticum_cell2 = PCA(n_components=2)
pca_microadriaticum_cell2.fit_transform(SHE_microadriaticum_cell2.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_microadriaticum_cell2 = pca_microadriaticum_cell2.components_[0]
pc2_microadriaticum_cell2 = pca_microadriaticum_cell2.components_[1]

SHE_microadriaticum_cell3 = pd.concat([microadriaticum_cell3], ignore_index=True, sort=False)
pca_microadriaticum_cell3 = PCA(n_components=2)
pca_microadriaticum_cell3.fit_transform(SHE_microadriaticum_cell3.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_microadriaticum_cell3 = pca_microadriaticum_cell3.components_[0]
pc2_microadriaticum_cell3 = pca_microadriaticum_cell3.components_[1]

#pilosum
SHE_pilosum_cell1 = pd.concat([pilosum_cell1], ignore_index=True, sort=False)
pca_pilosum_cell1 = PCA(n_components=2)
pca_pilosum_cell1.fit_transform(SHE_pilosum_cell1.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_pilosum_cell1 = pca_pilosum_cell1.components_[0]
pc2_pilosum_cell1 = pca_pilosum_cell1.components_[1]

SHE_pilosum_cell2 = pd.concat([pilosum_cell2], ignore_index=True, sort=False)
pca_pilosum_cell2 = PCA(n_components=2)
pca_pilosum_cell2.fit_transform(SHE_pilosum_cell2.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_pilosum_cell2 = pca_pilosum_cell2.components_[0]
pc2_pilosum_cell2 = pca_pilosum_cell2.components_[1]

SHE_pilosum_cell3 = pd.concat([pilosum_cell3], ignore_index=True, sort=False)
pca_pilosum_cell3 = PCA(n_components=2)
pca_pilosum_cell3.fit_transform(SHE_pilosum_cell3.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_pilosum_cell3 = pca_pilosum_cell3.components_[0]
pc2_pilosum_cell3 = pca_pilosum_cell3.components_[1]

#minutum
SHE_minutum_cell1 = pd.concat([minutum_cell1], ignore_index=True, sort=False)
pca_minutum_cell1 = PCA(n_components=2)
pca_minutum_cell1.fit_transform(SHE_minutum_cell1.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_minutum_cell1 = pca_minutum_cell1.components_[0]
pc2_minutum_cell1 = pca_minutum_cell1.components_[1]

SHE_minutum_cell2 = pd.concat([minutum_cell2], ignore_index=True, sort=False)
pca_minutum_cell2 = PCA(n_components=2)
pca_minutum_cell2.fit_transform(SHE_minutum_cell2.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_minutum_cell2 = pca_minutum_cell2.components_[0]
pc2_minutum_cell2 = pca_minutum_cell2.components_[1]

SHE_minutum_cell3 = pd.concat([minutum_cell3], ignore_index=True, sort=False)
pca_minutum_cell3 = PCA(n_components=2)
pca_minutum_cell3.fit_transform(SHE_minutum_cell3.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_minutum_cell3 = pca_minutum_cell3.components_[0]
pc2_minutum_cell3 = pca_minutum_cell3.components_[1]

#cohnii
SHE_cohnii_cell1 = pd.concat([cohnii_cell1], ignore_index=True, sort=False)
pca_cohnii_cell1 = PCA(n_components=2)
pca_cohnii_cell1.fit_transform(SHE_cohnii_cell1.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_cohnii_cell1 = pca_cohnii_cell1.components_[0]
pc2_cohnii_cell1 = pca_cohnii_cell1.components_[1]

SHE_cohnii_cell2 = pd.concat([cohnii_cell2], ignore_index=True, sort=False)
pca_cohnii_cell2 = PCA(n_components=2)
pca_cohnii_cell2.fit_transform(SHE_cohnii_cell2.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_cohnii_cell2 = pca_cohnii_cell2.components_[0]
pc2_cohnii_cell2 = pca_cohnii_cell2.components_[1]

SHE_cohnii_cell3 = pd.concat([cohnii_cell3], ignore_index=True, sort=False)
pca_cohnii_cell3 = PCA(n_components=2)
pca_cohnii_cell3.fit_transform(SHE_cohnii_cell3.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_cohnii_cell3 = pca_cohnii_cell3.components_[0]
pc2_cohnii_cell3 = pca_cohnii_cell3.components_[1]

#nutricula
SHE_nutricula_cell1 = pd.concat([nutricula_cell1], ignore_index=True, sort=False)
pca_nutricula_cell1 = PCA(n_components=2)
pca_nutricula_cell1.fit_transform(SHE_nutricula_cell1.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_nutricula_cell1 = pca_nutricula_cell1.components_[0]
pc2_nutricula_cell1 = pca_nutricula_cell1.components_[1]

SHE_nutricula_cell2 = pd.concat([nutricula_cell2], ignore_index=True, sort=False)
pca_nutricula_cell2 = PCA(n_components=2)
pca_nutricula_cell2.fit_transform(SHE_nutricula_cell2.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_nutricula_cell2 = pca_nutricula_cell2.components_[0]
pc2_nutricula_cell2 = pca_nutricula_cell2.components_[1]

SHE_nutricula_cell3 = pd.concat([nutricula_cell3], ignore_index=True, sort=False)
pca_nutricula_cell3 = PCA(n_components=2)
pca_nutricula_cell3.fit_transform(SHE_nutricula_cell3.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_nutricula_cell3 = pca_nutricula_cell3.components_[0]
pc2_nutricula_cell3 = pca_nutricula_cell3.components_[1]

#tyrrhenica
SHE_tyrrhenica_cell1 = pd.concat([tyrrhenica_cell1], ignore_index=True, sort=False)
pca_tyrrhenica_cell1 = PCA(n_components=2)
pca_tyrrhenica_cell1.fit_transform(SHE_tyrrhenica_cell1.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_tyrrhenica_cell1 = pca_tyrrhenica_cell1.components_[0]
pc2_tyrrhenica_cell1 = pca_tyrrhenica_cell1.components_[1]

#ross sea dinoflagellate
SHE_ross_sea_dinoflagellate_cell2 = pd.concat([ross_sea_dinoflagellate_cell2], ignore_index=True, sort=False)
pca_ross_sea_dinoflagellate_cell2 = PCA(n_components=2)
pca_ross_sea_dinoflagellate_cell2.fit_transform(SHE_ross_sea_dinoflagellate_cell2.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_ross_sea_dinoflagellate_cell2 = pca_ross_sea_dinoflagellate_cell2.components_[0]
pc2_ross_sea_dinoflagellate_cell2 = pca_ross_sea_dinoflagellate_cell2.components_[1]

SHE_ross_sea_dinoflagellate_cell3 = pd.concat([ross_sea_dinoflagellate_cell3], ignore_index=True, sort=False)
pca_ross_sea_dinoflagellate_cell3 = PCA(n_components=2)
pca_ross_sea_dinoflagellate_cell3.fit_transform(SHE_ross_sea_dinoflagellate_cell3.drop(columns=['surface_ridges','cell','length','width','surface_area','species','ROI','volume','delta_theta','LH_or_RH','aspect_ratio']))
pc1_ross_sea_dinoflagellate_cell3 = pca_ross_sea_dinoflagellate_cell3.components_[0]
pc2_ross_sea_dinoflagellate_cell3 = pca_ross_sea_dinoflagellate_cell3.components_[1]

#compute angle between cell-specific and global eigenvectors
cells_angles_PC1 = []
cells_angles_PC2 = []

#microadriaticum
cos_angle = np.dot(pc1_all, pc1_microadriaticum_cell1) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_microadriaticum_cell1))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 microadriaticum Cell 1: {angle_deg:.2f} degrees")
cells_angles_PC1.append(angle_deg)

cos_angle = np.dot(pc1_all, pc1_microadriaticum_cell2) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_microadriaticum_cell2))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 microadriaticum Cell 2: {angle_deg:.2f} degrees")
cells_angles_PC1.append(angle_deg)

cos_angle = np.dot(pc1_all, pc1_microadriaticum_cell3) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_microadriaticum_cell3))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 microadriaticum Cell 3: {angle_deg:.2f} degrees")
cells_angles_PC1.append(angle_deg)

#pilosum
cos_angle = np.dot(pc1_all, pc1_pilosum_cell1) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_pilosum_cell1))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 pilosum Cell 1: {angle_deg:.2f} degrees")
cells_angles_PC1.append(angle_deg)

cos_angle = np.dot(pc1_all, pc1_pilosum_cell2) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_pilosum_cell2))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 pilosum Cell 2: {angle_deg:.2f} degrees")
cells_angles_PC1.append(angle_deg)

cos_angle = np.dot(pc1_all, pc1_pilosum_cell3) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_pilosum_cell3))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 pilosum Cell 3: {angle_deg:.2f} degrees")
cells_angles_PC1.append(angle_deg)

#minutum
cos_angle = np.dot(pc1_all, pc1_minutum_cell1) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_minutum_cell1))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 minutum Cell 1: {angle_deg:.2f} degrees")
cells_angles_PC1.append(angle_deg)

cos_angle = np.dot(pc1_all, pc1_minutum_cell2) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_minutum_cell2))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 minutum Cell 2: {angle_deg:.2f} degrees")
cells_angles_PC1.append(angle_deg)

cos_angle = np.dot(pc1_all, pc1_minutum_cell3) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_minutum_cell3))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 minutum Cell 3: {angle_deg:.2f} degrees")
cells_angles_PC1.append(angle_deg)

#cohnii
cos_angle = np.dot(pc1_all, pc1_cohnii_cell1) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_cohnii_cell1))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 cohnii Cell 1: {angle_deg:.2f} degrees")
cells_angles_PC1.append(angle_deg)

cos_angle = np.dot(pc1_all, pc1_cohnii_cell2) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_cohnii_cell2))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 cohnii Cell 2: {angle_deg:.2f} degrees")
cells_angles_PC1.append(angle_deg)

cos_angle = np.dot(pc1_all, pc1_cohnii_cell3) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_cohnii_cell3))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 cohnii Cell 3: {angle_deg:.2f} degrees")
cells_angles_PC1.append(angle_deg)

#nutricula
cos_angle = np.dot(pc1_all, pc1_nutricula_cell1) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_nutricula_cell1))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 nutricula Cell 1: {angle_deg:.2f} degrees")
cells_angles_PC1.append(angle_deg)

cos_angle = np.dot(pc1_all, pc1_nutricula_cell2) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_nutricula_cell2))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 nutricula Cell 2: {angle_deg:.2f} degrees")
cells_angles_PC1.append(angle_deg)

cos_angle = np.dot(pc1_all, pc1_nutricula_cell3) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_nutricula_cell3))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 nutricula Cell 3: {angle_deg:.2f} degrees")
cells_angles_PC1.append(angle_deg)

#tyrrhenica
cos_angle = np.dot(pc1_all, pc1_tyrrhenica_cell1) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_tyrrhenica_cell1))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 tyrrhenica Cell 1: {angle_deg:.2f} degrees")
cells_angles_PC1.append(angle_deg)

#ross sea dinoflagellate
cos_angle = np.dot(pc1_all, pc1_ross_sea_dinoflagellate_cell2) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_ross_sea_dinoflagellate_cell2))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 K. sp. Cell 2: {angle_deg:.2f} degrees")
cells_angles_PC1.append(angle_deg)

cos_angle = np.dot(pc1_all, pc1_ross_sea_dinoflagellate_cell3) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_ross_sea_dinoflagellate_cell3))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 K. sp. Cell 3: {angle_deg:.2f} degrees")
cells_angles_PC1.append(angle_deg)

######################################################################################
#plot \theta_{cell}
angle_rad = [radians(a) for a in cells_angles_PC1]

# Calculate unit vector components
x = np.cos(angle_rad)
y = np.sin(angle_rad)

# Plot setup
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_aspect('equal')

# Limit to top-right quadrant
ax.set_xlim(0, 1.1)
ax.set_ylim(0, 1.2)

# Color list (extend as needed)
colors = ['blue', 'blue', 'blue', 'purple', 'purple', 'purple', 'orange', 'orange', 'orange', 'pink', 'pink', 'pink', 'red', 'red', 'red', 'grey', 'black', 'black']

# Draw unit vectors as arrows
for i in range(len(angle_rad)):
    arrow_length = np.linspace(1,0.2,len(angle_rad))
    Nchrs = [len(microadriaticum_cell1), len(microadriaticum_cell2), len(microadriaticum_cell3), len(pilosum_cell1), len(pilosum_cell2), len(pilosum_cell3),
             len(minutum_cell1), len(minutum_cell2), len(minutum_cell3), len(cohnii_cell1), len(cohnii_cell2), len(cohnii_cell3),
             len(nutricula_cell1), len(nutricula_cell2), len(nutricula_cell3), len(tyrrhenica_cell1), len(ross_sea_dinoflagellate_cell2), len(ross_sea_dinoflagellate_cell3),
             ]
    Nchrs = [x / max(Nchrs) for x in Nchrs]
    ax.arrow(0, 0, x[i]*Nchrs[i], y[i]*Nchrs[i],
             head_width=0.03, head_length=0.05,
             fc=colors[i], ec=colors[i])

# Radial guide lines for 0°, 90°, 30°, and 60°
def draw_angle_line(deg, ax, length=1.1, label=None):
    rad = radians(deg)
    ax.plot([0, cos(rad)], [0, sin(rad)], color='black', linestyle='dotted')
    if label:
        ax.text(cos(rad) * length, sin(rad) * length, label,
                fontsize=16, ha='left', va='bottom')

# Draw labeled guide lines
draw_angle_line(0, ax, label='0°')
draw_angle_line(30, ax, label='30°')
draw_angle_line(60, ax, label='60°')
draw_angle_line(90, ax, label='90°')

# Remove axis ticks and frame
# Turn off spines
for spine in ax.spines.values():
    spine.set_visible(False)

# Turn off ticks and tick labels
ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
plt.title('Angle between cell-specific PC1 and global PC1', fontsize=14)
plt.ylabel('Orthogonal shape variation', fontsize=16)
plt.xlabel('Identical shape variation', fontsize=16)
plt.show()

######################################################################################

#inperpolate across PC1 and save 3D shapes
coords_sample_along_PC1 = np.column_stack((np.linspace(min(trans[:,0]), max(trans[:,0]), 20), np.zeros(20))) #create 20 shapes
PCA_coords_to_shape(pca=pca_all, PCA_coords=coords_sample_along_PC1, lmax=lmax, save_directory='/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Spherical Harmonics Expansion/all/PC1/', save_as_tiff=True, save_cross_sections=True)

#inperpolate across PC2 and save 3D shapes
coords_sample_along_PC2 = np.column_stack((np.zeros(20), np.linspace(min(trans[:,0]), max(trans[:,0]), 20))) #create 20 shapes
PCA_coords_to_shape(pca=pca_all, PCA_coords=coords_sample_along_PC2, lmax=lmax, save_directory='/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Spherical Harmonics Expansion/all/PC2/', save_as_tiff=True, save_cross_sections=True)

#optional: view these shapes
#napari.view_image(voxelized)

# # Get shapes corresponding to corner points in PC1/PC2 space
# corner_points = np.array([
#     [-60,   0],   #bottom left triangle
#     [-60, -15],   #bottom left triangle
#     [  0, -15],   #bottom left triangle
#     [-60,  10],   #top left triangle
#     [-60,  30],   #top left triangle
#     [-20,  30]    #top left triangle
# ])

# PCA_coords_to_shape(pca=pca_all, PCA_coords=corner_points, lmax=lmax, save_directory='/Users/lucasphilipp/Downloads/corners/', save_as_tiff=False, save_cross_sections=False)

random.seed(42)  # Set seed for reproducibility.
#compute reconstruction error for 50 random chromosomes across all datasets
i=0
num_chroms = 3
max_lmax = 40
error = np.zeros((num_chroms, max_lmax))
colors = [None] * num_chroms
while i < num_chroms:
    skip_outer = False
    dataset_num = random.randint(1, 19)
    if dataset_num == 1:
        chr_num = random.randint(2,114)
        image_file = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/Cell 1/ROI '+str(chr_num)+'.tiff'
        colors[i] = 'blue'
    elif dataset_num == 2:
        chr_num = random.randint(1,100)
        image_file = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/Cell 2/ROI '+str(chr_num)+'.tiff'
        colors[i] = 'blue'
    elif dataset_num == 3:
        chr_num = random.randint(2,105)
        image_file = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/Cell 3/ROI '+str(chr_num)+'.tiff'
        colors[i] = 'blue'
    elif dataset_num == 4:
        chr_num = random.randint(4,102)
        image_file = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 1/ROI '+str(chr_num)+'.tiff'
        colors[i] = 'purple'
    elif dataset_num == 5:
        chr_num = random.randint(22,124)
        image_file = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 2/ROI '+str(chr_num)+'.tiff'
        colors[i] = 'purple'
    elif dataset_num == 6:
        chr_num = random.randint(27,124)
        image_file = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 3/ROI '+str(chr_num)+'.tiff'
        colors[i] = 'purple'
    elif dataset_num == 7:
        chr_num = random.randint(1,26)
        image_file = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Breviolum minutum/Cell 1/ROI '+str(chr_num)+'.tiff'
        colors[i] = 'orange'
    elif dataset_num == 8:
        chr_num = random.randint(2,35)
        image_file = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Breviolum minutum/Cell 2/ROI '+str(chr_num)+'.tiff'
        colors[i] = 'orange'
    elif dataset_num == 9:
        chr_num = random.randint(1,33)
        image_file = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Breviolum minutum/Cell 3/ROI '+str(chr_num)+'.tiff'
        colors[i] = 'orange'
    elif dataset_num == 10:
        chr_num = random.randint(2,162)
        image_file = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/brandtodinium cell 1 chromosomes 4nm sampling/ROI '+str(chr_num)+'.tiff'
        colors[i] = 'red'
    elif dataset_num == 11:
        chr_num = random.randint(1,100)
        image_file = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/brandtodinium cell 2 chromosomes 4nm sampling/ROI '+str(chr_num)+'.tiff'
        colors[i] = 'red'
    elif dataset_num == 12:
        chr_num = random.randint(4,327)
        image_file = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/brandtodinium cell 3 chromosomes 4nm sampling/ROI '+str(chr_num)+'.tiff'
        colors[i] = 'red'
    elif dataset_num == 14:
        chr_num = random.randint(2,118)
        image_file = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/Cell 1/ROI '+str(chr_num)+'.tiff'
        colors[i] = 'pink'
    elif dataset_num == 15:
        chr_num = random.randint(2,187)
        image_file = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/Cell 2/ROI '+str(chr_num)+'.tiff'
        colors[i] = 'pink'
    elif dataset_num == 16:
        chr_num = random.randint(2,135)
        image_file = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/Cell 3/ROI '+str(chr_num)+'.tiff'
        colors[i] = 'pink'
    elif dataset_num == 17:
        chr_num = random.randint(2,113)
        image_file = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Ensiculifera tyrrhenica/ensiculifera tyrrhenica chromosomes 4nm sampling/ROI '+str(chr_num)+'.tiff'
        colors[i] = 'grey'
    elif dataset_num == 18:
        chr_num = random.randint(1,165)
        image_file = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Ross Sea Dinoflagellate/Ross Sea Dinoflagellate Cell 2 Chromosomes 4nm voxels/ROI '+str(chr_num)+'.tiff'
        colors[i] = 'black'
    elif dataset_num == 19:
        chr_num = random.randint(1,130)
        image_file = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Ross Sea Dinoflagellate/Ross Sea Dinoflagellate Cell 3 Chromosomes 4nm voxels/ROI '+str(chr_num)+'.tiff'
        colors[i] = 'black'
    for j in range(1,max_lmax+1): #sample lmax 1-50
        (coeffs, grid_rec), (_, _, grid_input, _), flag = shparam.get_shcoeffs(image=imread(image_file), lmax=j)
        if flag == True:
            skip_outer = True
            print(flag)
            break
        error[i,j-1] = shtools.get_reconstruction_error(grid_rec,grid_input)
        
        mat=shtools.convert_coeffs_dict_to_matrix(coeffs,j)
        mesh_rec, _ = shtools.get_reconstruction_from_coeffs(mat)
        coords = vtknp.vtk_to_numpy(mesh_rec.GetPoints().GetData())

        # Find bounds of the mesh
        rmin = (coords.min(axis=0) - 0.5).astype(int)
        rmax = (coords.max(axis=0) + 0.5).astype(int)

        # Width, height and depth
        w = int(2 + (rmax[0] - rmin[0]))
        h = int(2 + (rmax[1] - rmin[1]))
        d = int(2 + (rmax[2] - rmin[2]))

        # Create image data
        imagedata = vtk.vtkImageData()
        imagedata.SetDimensions([w, h, d])
        imagedata.SetExtent(0, w - 1, 0, h - 1, 0, d - 1)
        imagedata.SetOrigin(rmin)
        imagedata.AllocateScalars(vtk.VTK_UNSIGNED_CHAR, 1)

        # Set all values to 1
        imagedata.GetPointData().GetScalars().FillComponent(0, 1)

        # Create an empty 3D numpy array to sum up
        # voxelization of all meshes
        img = np.zeros((d, h, w), dtype=np.uint8)

        voxelized = shtools.voxelize_mesh(imagedata=imagedata, shape=(d, h, w), mesh=mesh_rec, origin=rmin)

        voxelized = voxelized.astype('int8')
        #tifffile.imwrite('lmax {}.tiff'.format(i), voxelized, bigtiff=True) #uncomment if you want to write the reconstructed chromosome volume as a tiff file
    if skip_outer:
        continue
    print(i)
    i = i+1

for i in range(num_chroms):
    plt.plot(np.linspace(1,max_lmax, max_lmax), error[i,:], color=colors[i])
plt.xlabel('lmax', fontsize=16)
plt.ylabel('Reconstruction Error', fontsize=16)
plt.yscale('log')
plt.tick_params(axis='both', which='major', labelsize=14)
plt.show()
        
#lmax = 40 looks good

