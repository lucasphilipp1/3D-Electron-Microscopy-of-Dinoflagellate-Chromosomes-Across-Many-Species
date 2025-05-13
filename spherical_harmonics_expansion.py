"""
@author: lucasphilipp
code adapted from: https://github.com/AllenCell/aics-shparam
"""

# Import required packages
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
import matplotlib.colors as mcolors
from matplotlib.font_manager import FontProperties

import napari
#viewer = napari.view_image(numpy_data)

from PIL import Image
import warnings

from matplotlib.colors import ListedColormap, to_rgb

def values(x):
    return x.str.extract(r'([0-9]+)')

def normalize_to_uint8(img):
    img = img.astype(np.float32)
    img -= img.min()
    if img.max() != 0:
        img /= img.max()
    return (img * 255).astype(np.uint8)

def generate_shades(base_color, n_shades):
    """Generate n_shades from light to dark of a base RGB color."""
    shades = []
    for i in range(n_shades):
        factor = 0.5 + 1.5 * (i / (n_shades - 1))  # avoid pure white
        shaded = tuple(np.clip(factor * np.array(base_color), 0, 1))
        shades.append(shaded)
    return shades

#get images and angles and perform spherical harmonic expansion
def spherical_harmonic_expansion(ROI_start,ROI_end,image_path,angle_csv_path,cell_num,species_name,lmax):
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
    is_angle_ext = []
    delta_angle = []
    LH_or_RH = []
    for i in range(ROI_start,ROI_end+1):
        is_angle_ext.append(angles.Label.eq(i).any()) #was an angle extracted for this chromosome?        
        temp = angles.index[angles['Label'] == i]
        if len(temp)>0:
            delta_angle.append(angles.Angle[temp[0]]-angles.Angle[temp[1]])
            if (angles.Angle[temp[0]]-angles.Angle[temp[1]])<-5:
                LH_or_RH.append('left handed')
            elif (angles.Angle[temp[0]]-angles.Angle[temp[1]])>5:
                LH_or_RH.append('right handed')
            else:
                LH_or_RH.append('|Δθ|<5 deg')
        else:
            delta_angle.append(np.NAN)
            LH_or_RH.append(np.NAN)
    #spherical harmonic expansion
    df_coeffs_list = []
    volume = []
    aspect_ratio = []
    warnings_iter = []
    for i in range(len(images)):
        # with warnings.catch_warnings(record=True) as caught_warnings:
        #     warnings.simplefilter("always")
        (coeffs, _), _ , flag = shparam.get_shcoeffs(image=images[i], lmax=lmax)
        if flag == True:
            warnings_iter.append(i)
        volume.append((4**3)*np.sum(images[i]/np.max(images[i]))) #volume in units [nm^3] assuming 4nmx4nm4xnm voxels 
        aspect_ratio.append(np.size(images[i],axis=0)/np.size(images[i],axis=1))      
        df_coeffs_list.append(coeffs)
        print(i)
    df_coeffs = pd.DataFrame(df_coeffs_list)
    #add is_angle_ext label
    df_coeffs['is_angle_ext'] = is_angle_ext
    #add cell label
    df_coeffs['cell'] = np.ones(len(images))*cell_num
    #add species label
    species = list()
    for i in range(len(images)):
        species.append(species_name)
    df_coeffs['species'] = species
    #add ROI label
    df_coeffs['ROI']=ROI
    #add volume label
    df_coeffs['volume']=volume
    #add delta_angle label
    df_coeffs['delta_angle']=delta_angle
    #add LH_or_RH label
    df_coeffs['LH_or_RH']=LH_or_RH
    #add aspect ratio label
    df_coeffs['aspect_ratio']=aspect_ratio
    #add toroid/rod label
    #...
    #return df_coeffs
    df_coeffs = df_coeffs.drop(index=warnings_iter)
    print(warnings_iter)
    return df_coeffs

# Compute spherical harmonics coefficients of shape and store them in a pandas dataframe.
lmax = 40 #number of expansion terms (sort of)

microadriaticum_cell1 = spherical_harmonic_expansion(2,114,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/Cell 1/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/SM Cell 1.csv',1,'microadriaticum',lmax)
microadriaticum_cell2 = spherical_harmonic_expansion(1,100,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/Cell 2/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/SM Cell 2.csv',2,'microadriaticum',lmax)
microadriaticum_cell3 = spherical_harmonic_expansion(2,105,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/Cell 3/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/SM Cell 3.csv',3,'microadriaticum',lmax)

pilosum_cell1 = spherical_harmonic_expansion(4,102,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 1/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/SP Cell 1.csv',1,'pilosum',lmax)
pilosum_cell2 = spherical_harmonic_expansion(22,124,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 2/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/SP Cell 2.csv',2,'pilosum',lmax)
pilosum_cell3 = spherical_harmonic_expansion(27,124,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 3/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/SP Cell 3.csv',3,'pilosum',lmax)

minutum_cell1 = spherical_harmonic_expansion(1,26,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/Cell 1/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/SMin Cell 1.csv',1,'minutum',lmax)
minutum_cell2 = spherical_harmonic_expansion(2,35,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/Cell 2/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/SMin Cell 2.csv',2,'minutum',lmax)
minutum_cell3 = spherical_harmonic_expansion(1,33,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/Cell 3/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/SMin Cell 3.csv',3,'minutum',lmax)

nutricula_cell1 = spherical_harmonic_expansion(2,162,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/brandtodinium cell 1 chromosomes 4nm sampling/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/BN Cell 1.csv',1,'nutricula',lmax)
nutricula_cell2 = spherical_harmonic_expansion(1,100,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/brandtodinium cell 2 chromosomes 4nm sampling/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/BN Cell 2.csv',2,'nutricula',lmax)
nutricula_cell3 = spherical_harmonic_expansion(4,327,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/brandtodinium cell 3 chromosomes 4nm sampling/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/BN Cell 3.csv',3,'nutricula',lmax)

cohnii_cell1 = spherical_harmonic_expansion(2,118,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/Cell 1/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/CC cell 1.csv',1,'cohnii',lmax)
cohnii_cell2 = spherical_harmonic_expansion(2,187,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/Cell 2/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/CC cell 2.csv',2,'cohnii',lmax)
cohnii_cell3 = spherical_harmonic_expansion(2,135,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/Cell 3/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/CC cell 3.csv',3,'cohnii',lmax)

tyrrhenica_cell1 = spherical_harmonic_expansion(2,113,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Ensiculifera tyrrhenica/ensiculifera tyrrhenica chromosomes 4nm sampling/','/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Ensiculifera tyrrhenica/ET Cell 1.csv',1,'tyrrhenica',lmax)

SHE_all = pd.concat([microadriaticum_cell1, microadriaticum_cell2, microadriaticum_cell3, pilosum_cell1, pilosum_cell2, pilosum_cell3, minutum_cell1, minutum_cell2, minutum_cell3, nutricula_cell1, nutricula_cell2, nutricula_cell3, cohnii_cell1, cohnii_cell2, cohnii_cell3, tyrrhenica_cell1], ignore_index=True, sort=False)

#ensure PCA is fed an equal # of chromosomes per species
num_chroms = []
num_chroms.append((SHE_all['species'] == 'microadriaticum').sum())
num_chroms.append((SHE_all['species'] == 'pilosum').sum())
num_chroms.append((SHE_all['species'] == 'minutum').sum())
num_chroms.append((SHE_all['species'] == 'nutricula').sum())
num_chroms.append((SHE_all['species'] == 'cohnii').sum())
num_chroms.append((SHE_all['species'] == 'tyrrhenica').sum())

match_rows = SHE_all[SHE_all['species'] == 'microadriaticum']
rows_to_drop = match_rows.sample(n=num_chroms[0]-min(num_chroms))
SHE_all = SHE_all.drop(rows_to_drop.index)

match_rows = SHE_all[SHE_all['species'] == 'pilosum']
rows_to_drop = match_rows.sample(n=num_chroms[1]-min(num_chroms))
SHE_all = SHE_all.drop(rows_to_drop.index)

match_rows = SHE_all[SHE_all['species'] == 'minutum']
rows_to_drop = match_rows.sample(n=num_chroms[2]-min(num_chroms))
SHE_all = SHE_all.drop(rows_to_drop.index)

match_rows = SHE_all[SHE_all['species'] == 'nutricula']
rows_to_drop = match_rows.sample(n=num_chroms[3]-min(num_chroms))
SHE_all = SHE_all.drop(rows_to_drop.index)


match_rows = SHE_all[SHE_all['species'] == 'cohnii']
rows_to_drop = match_rows.sample(n=num_chroms[4]-min(num_chroms))
SHE_all = SHE_all.drop(rows_to_drop.index)

match_rows = SHE_all[SHE_all['species'] == 'tyrrhenica']
rows_to_drop = match_rows.sample(n=num_chroms[5]-min(num_chroms))
SHE_all = SHE_all.drop(rows_to_drop.index)

SHE_all = SHE_all.reset_index(drop=True)

# Vizualize the resulting dataframe
with pd.option_context('display.max_rows', 5, 'display.max_columns', 5):
    print(SHE_all)
    # Let's use PCA to reduce the dimensionality of the coefficients
# dataframe from 51 down to 2.
pca_all = PCA(n_components=2)

trans = pca_all.fit_transform(SHE_all.drop(columns=['is_angle_ext','cell','species','ROI','volume','delta_angle','LH_or_RH','aspect_ratio']))

#% variance explained by top PCs?
print(pca_all.explained_variance_ratio_)

num_dim_pre_PCA = SHE_all.drop(columns=['is_angle_ext','cell','species','ROI','volume','delta_angle','LH_or_RH','aspect_ratio']).shape[1]
df_trans = pd.DataFrame(trans)
df_trans.columns = ['PC1', 'PC2']

#CARRY OVER FROM PRE-PCA
df_trans['is_angle_ext'] = SHE_all.is_angle_ext
df_trans['cell'] = SHE_all.cell
df_trans['species'] = SHE_all.species
df_trans['ROI'] = SHE_all.ROI
df_trans['volume'] = SHE_all.volume
df_trans['delta_angle'] = SHE_all.delta_angle
df_trans['LH_or_RH'] = SHE_all.LH_or_RH
df_trans['aspect_ratio'] = SHE_all.aspect_ratio

#color by is_angle_ext
with pd.option_context('display.max_rows', 5, 'display.max_columns', 5):
    print(df_trans)
    # Scatter plot to show how similar shapes are grouped together.
fig, ax = plt.subplots(1,1, figsize=(8,8))
for label, df_label in df_trans.groupby('is_angle_ext'):
    ax.scatter(df_label.PC1, df_label.PC2, label=label, s=30)
ax.tick_params(axis='both', which='major', labelsize=16)
plt.legend(['No Surface Ridges', 'Surface Ridges'], loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=20)
plt.xlabel('PC1', fontsize=20)
plt.ylabel('PC2', fontsize=20)
left, right = plt.xlim()
up, down = plt.ylim()
#plt.xlim((left,100)) 
left, right = plt.xlim()
plt.show()

#color by volume
fig, ax = plt.subplots(1,1, figsize=(8,8))
sc = ax.scatter(df_trans.PC1, df_trans.PC2, c=df_trans['volume'], cmap='hsv', norm='log', s=20)
ax.tick_params(axis='both', which='major', labelsize=16)
cbar = fig.colorbar(sc, ax=ax)
cbar.set_label('Volume [nm³]', fontsize=20)
cbar.ax.tick_params(labelsize=16)
plt.xlabel('PC1', fontsize=20)
plt.ylabel('PC2', fontsize=20)
plt.xlim((left, right))
plt.ylim((up, down))
plt.show()

vmin = 1
vmax = 3
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

# Define 8 base hues (manually picked or derived from a colormap like 'tab20')
base_hues = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
]
base_colors = [to_rgb(color) for color in base_hues]

# Generate 5 shades per base color
custom_colors = []
for base in base_colors:
    custom_colors.extend(generate_shades(base, 5))

# Create the new colormap
custom_cmap = ListedColormap(custom_colors)

#color by aspect ratio
fig, ax = plt.subplots(1,1, figsize=(8,8))
sc = ax.scatter(df_trans.PC1, df_trans.PC2, c=df_trans['aspect_ratio'], norm=norm, cmap=custom_cmap, s=20)
ax.tick_params(axis='both', which='major', labelsize=16)
cbar = fig.colorbar(sc, ax=ax)
cbar.set_label('Aspect ratio')
plt.xlabel('PC1', fontsize=20)
plt.ylabel('PC2', fontsize=20)
plt.xlim((left, right))
plt.ylim((up, down))
cbar.set_ticks([1, 1.5, 2, 2.5, 3])
cbar.set_ticklabels(['1', '1.5', '2', '2.5', '>3'])
plt.show()

#color by delta_angle
vmin = -90
vmax = 90
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

fig, ax = plt.subplots(1,1, figsize=(8,8))
sc = ax.scatter(df_trans.PC1, df_trans.PC2, c=df_trans['delta_angle'], norm=norm, cmap='coolwarm', s=20)
ax.tick_params(axis='both', which='major', labelsize=16)
cbar = fig.colorbar(sc, ax=ax)
cbar.set_label('Δθ [degrees]', fontsize=20)
cbar.ax.tick_params(labelsize=16)
plt.xlabel('PC1', fontsize=20)
plt.ylabel('PC2', fontsize=20)
plt.xlim((left, right))
plt.ylim((up, down))
plt.show()

#color by LH_or_RH
with pd.option_context('display.max_rows', 5, 'display.max_columns', 5):
    print(df_trans)
    # Scatter plot to show how similar shapes are grouped together.
fig, ax = plt.subplots(1,1, figsize=(8,8))
for label, df_label in df_trans.groupby('LH_or_RH'):
    ax.scatter(df_label.PC1, df_label.PC2, label=label, s=30)
ax.tick_params(axis='both', which='major', labelsize=16)
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=20)
plt.xlabel('PC1', fontsize=20)
plt.ylabel('PC2', fontsize=20)
plt.xlim((left, right))
plt.ylim((up, down))
plt.show()

#color by species
colors = ['pink', 'blue', 'orange', 'red', 'purple','grey']
count = 0
fig, ax = plt.subplots(1,1, figsize=(8,8))
for label, df_label in df_trans.groupby('species'):
    ax.scatter(df_label.PC1, df_label.PC2, label=label, s=20, color = colors[count])
    count = count + 1
ax.tick_params(axis='both', which='major', labelsize=16)
italic_font = FontProperties(style='italic')
plt.legend(['C. cohnii','S. microadriaticum','S. minutum','B. nutricula','S. pilosum','E. ttyrrhenica'],loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=20, prop=italic_font)
plt.xlabel('PC1', fontsize=20)
plt.ylabel('PC2', fontsize=20)
left, right = plt.xlim()
plt.xlim((left,100)) 
plt.show()

#aspect ratio vs volume
colors = ['pink', 'blue', 'orange', 'red', 'purple','grey']
count = 0
fig, ax = plt.subplots(1,1, figsize=(8,8))
for label, df_label in df_trans.groupby('species'):
    ax.scatter(df_label.volume, df_label.aspect_ratio, label=label, s=20, color = colors[count])
    count = count + 1
ax.tick_params(axis='both', which='major', labelsize=16)
italic_font = FontProperties(style='italic')
plt.legend(['C. cohnii','S. microadriaticum','S. minutum','B. nutricula','S. pilosum','E. ttyrrhenica'],loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=20, prop=italic_font)
plt.xlabel('Volume [nm³]', fontsize=20)
plt.ylabel('Aspect ratio', fontsize=20)
plt.xlim((left,1e9)) 
plt.show()

##############################################################################
#do species-specific shape analayis using PCA

#get global eigenvector
pc1_all = pca_all.components_[0]
pc2_all = pca_all.components_[1]

#microadriaticum
SHE_microadriaticum = pd.concat([microadriaticum_cell1, microadriaticum_cell2, microadriaticum_cell3], ignore_index=True, sort=False)
pca_microadriaticum = PCA(n_components=2)
pca_microadriaticum.fit_transform(SHE_microadriaticum.drop(columns=['is_angle_ext','cell','species','ROI','volume','delta_angle','LH_or_RH','aspect_ratio']))
pc1_microadriaticum = pca_microadriaticum.components_[0]
pc2_microadriaticum = pca_microadriaticum.components_[1]

#pilosum
SHE_pilosum = pd.concat([pilosum_cell1, pilosum_cell2, pilosum_cell3], ignore_index=True, sort=False)
pca_pilosum = PCA(n_components=2)
pca_pilosum.fit_transform(SHE_pilosum.drop(columns=['is_angle_ext','cell','species','ROI','volume','delta_angle','LH_or_RH','aspect_ratio']))
pc1_pilosum = pca_pilosum.components_[0]
pc2_pilosum = pca_pilosum.components_[1]

#minutum
SHE_minutum = pd.concat([minutum_cell1, minutum_cell2, minutum_cell3], ignore_index=True, sort=False)
pca_minutum = PCA(n_components=2)
pca_minutum.fit_transform(SHE_minutum.drop(columns=['is_angle_ext','cell','species','ROI','volume','delta_angle','LH_or_RH','aspect_ratio']))
pc1_minutum = pca_minutum.components_[0]
pc2_minutum = pca_minutum.components_[1]

#cohnii
SHE_cohnii = pd.concat([cohnii_cell1, cohnii_cell2, cohnii_cell3], ignore_index=True, sort=False)
pca_cohnii = PCA(n_components=2)
pca_cohnii.fit_transform(SHE_cohnii.drop(columns=['is_angle_ext','cell','species','ROI','volume','delta_angle','LH_or_RH','aspect_ratio']))
pc1_cohnii = pca_cohnii.components_[0]
pc2_cohnii = pca_cohnii.components_[1]

#nutricula
SHE_nutricula = pd.concat([nutricula_cell1, nutricula_cell2, nutricula_cell3], ignore_index=True, sort=False)
pca_nutricula = PCA(n_components=2)
pca_nutricula.fit_transform(SHE_nutricula.drop(columns=['is_angle_ext','cell','species','ROI','volume','delta_angle','LH_or_RH','aspect_ratio']))
pc1_nutricula = pca_nutricula.components_[0]
pc2_nutricula = pca_nutricula.components_[1]

#tyrrhenica
SHE_tyrrhenica = pd.concat([tyrrhenica_cell1], ignore_index=True, sort=False)
pca_tyrrhenica = PCA(n_components=2)
pca_tyrrhenica.fit_transform(SHE_tyrrhenica.drop(columns=['is_angle_ext','cell','species','ROI','volume','delta_angle','LH_or_RH','aspect_ratio']))
pc1_tyrrhenica = pca_tyrrhenica.components_[0]
pc2_tyrrhenica = pca_tyrrhenica.components_[1]

#compute angle between species-specific and global eigenvectors
#microadriaticum
cos_angle = np.dot(pc1_all, pc1_microadriaticum) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_microadriaticum))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 microadriaticum: {angle_deg:.2f} degrees")

cos_angle = np.dot(pc2_all, pc2_microadriaticum) / (np.linalg.norm(pc2_all) * np.linalg.norm(pc2_microadriaticum))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC2 Global and PC2 microadriaticum: {angle_deg:.2f} degrees")

#pilosum
cos_angle = np.dot(pc1_all, pc1_pilosum) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_pilosum))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 pilosum: {angle_deg:.2f} degrees")

cos_angle = np.dot(pc2_all, pc2_pilosum) / (np.linalg.norm(pc2_all) * np.linalg.norm(pc2_pilosum))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC2 Global and PC2 pilosum: {angle_deg:.2f} degrees")

#minutum
cos_angle = np.dot(pc1_all, pc1_minutum) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_minutum))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 minutum: {angle_deg:.2f} degrees")

cos_angle = np.dot(pc2_all, pc2_minutum) / (np.linalg.norm(pc2_all) * np.linalg.norm(pc2_minutum))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC2 Global and PC2 minutum: {angle_deg:.2f} degrees")

#cohnii
cos_angle = np.dot(pc1_all, pc1_cohnii) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_cohnii))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 cohnii: {angle_deg:.2f} degrees")

cos_angle = np.dot(pc2_all, pc2_cohnii) / (np.linalg.norm(pc2_all) * np.linalg.norm(pc2_cohnii))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC2 Global and PC2 cohnii: {angle_deg:.2f} degrees")

#nutricula
cos_angle = np.dot(pc1_all, pc1_nutricula) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_nutricula))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 nutricula: {angle_deg:.2f} degrees")

cos_angle = np.dot(pc2_all, pc2_nutricula) / (np.linalg.norm(pc2_all) * np.linalg.norm(pc2_nutricula))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC2 Global and PC2 nutricula: {angle_deg:.2f} degrees")

#tyrrhenica
cos_angle = np.dot(pc1_all, pc1_tyrrhenica) / (np.linalg.norm(pc1_all) * np.linalg.norm(pc1_tyrrhenica))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC1 Global and PC1 tyrrhenica: {angle_deg:.2f} degrees")

cos_angle = np.dot(pc2_all, pc2_tyrrhenica) / (np.linalg.norm(pc2_all) * np.linalg.norm(pc2_tyrrhenica))
angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
if angle_deg > 90:
    angle_deg=180-angle_deg #get the acute angle
print(f"Angle between PC2 Global and PC2 tyrrhenica: {angle_deg:.2f} degrees")

#inperpolate across PC1 and show 3D shape
pca_all.components_

#project all data onto PC1 and extract variance
fig, ax = plt.subplots(1,1, figsize=(3,3))
ax.hist(trans[:,0])
plt.show()

num_sample_shapes = 20;

shape_params_along_PC1 = np.zeros((num_sample_shapes,num_dim_pre_PCA))
shape_params_along_PC2 = np.zeros((num_sample_shapes,num_dim_pre_PCA))

for i in range(num_sample_shapes):
    PC1_coord = np.linspace(min(trans[:,0]),max(trans[:,0]), num_sample_shapes)
    PC2_coord = np.linspace(min(trans[:,1]),max(trans[:,1]), num_sample_shapes)
    shape_params_along_PC1[i,:] = pca_all.inverse_transform((PC1_coord[i],0))
    shape_params_along_PC2[i,:] = pca_all.inverse_transform((0,PC2_coord[i]))

shape_params_along_PC1_reshape=np.zeros((2,lmax,lmax,num_sample_shapes))
shape_params_along_PC2_reshape=np.zeros((2,lmax,lmax,num_sample_shapes))

for d in range(num_sample_shapes):
    count = 0
    for i in range(2):
        for j in range(lmax):
            for k in range(lmax):
                #CHANGE shape_params_along_PC1[0
                shape_params_along_PC1_reshape[i,j,k,d] = shape_params_along_PC1[d,count]
                shape_params_along_PC2_reshape[i,j,k,d] = shape_params_along_PC2[d,count]
                #print(i,j,k,d)
                count = count + 1
            
for i in range(num_sample_shapes):
    mesh_PC1, _ = shtools.get_reconstruction_from_coeffs(shape_params_along_PC1_reshape[:,:,:,i])
    coords = vtknp.vtk_to_numpy(mesh_PC1.GetPoints().GetData())

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
    
    voxelized = shtools.voxelize_mesh(imagedata=imagedata, shape=(d, h, w), mesh=mesh_PC1, origin=rmin)

    voxelized = voxelized.astype('int8')
    #tifffile.imsave('/Users/lucasphilipp/Downloads/Spherical Harmonics Expansion/all/PC1/PC1 {}.tiff'.format(i), voxelized, bigtiff=True)
    
    z, y, x = voxelized.shape
    xy_section = voxelized[z // 2, :, :]
    yz_section = voxelized[:, :, x // 2]
    xz_section = voxelized[:, y // 2, :]
    
    xy_section = feature.canny(normalize_to_uint8(xy_section), sigma=1.0)
    yz_section = feature.canny(normalize_to_uint8(yz_section), sigma=1.0)
    xz_section = feature.canny(normalize_to_uint8(xz_section), sigma=1.0)
    
    xy_img = Image.fromarray(normalize_to_uint8(xy_section))
    yz_img = Image.fromarray(normalize_to_uint8(yz_section))
    xz_img = Image.fromarray(normalize_to_uint8(xz_section))

    # Save as PNGs (or change extension to .tiff if preferred)
    #xy_img.save('/Users/lucasphilipp/Downloads/Spherical Harmonics Expansion/all/PC1/xy/{}.png'.format(i))
    #yz_img.save('/Users/lucasphilipp/Downloads/Spherical Harmonics Expansion/all/PC1/yz/{}.png'.format(i))
    #xz_img.save('/Users/lucasphilipp/Downloads/Spherical Harmonics Expansion/all/PC1/xz/{}.png'.format(i))

    mesh_PC2, _ = shtools.get_reconstruction_from_coeffs(shape_params_along_PC2_reshape[:,:,:,i])
    coords = vtknp.vtk_to_numpy(mesh_PC2.GetPoints().GetData())

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
    
    voxelized = shtools.voxelize_mesh(imagedata=imagedata, shape=(d, h, w), mesh=mesh_PC2, origin=rmin)

    voxelized = voxelized.astype('int8')
    #tifffile.imsave('/Users/lucasphilipp/Downloads/Spherical Harmonics Expansion/all/PC2/PC2 {}.tiff'.format(i), voxelized, bigtiff=True)
    
    z, y, x = voxelized.shape
    xy_section = voxelized[z // 2, :, :]
    yz_section = voxelized[:, :, x // 2]
    xz_section = voxelized[:, y // 2, :]
    
    xy_section = feature.canny(normalize_to_uint8(xy_section), sigma=1.0)
    yz_section = feature.canny(normalize_to_uint8(yz_section), sigma=1.0)
    xz_section = feature.canny(normalize_to_uint8(xz_section), sigma=1.0)
    
    xy_img = Image.fromarray(normalize_to_uint8(xy_section))
    yz_img = Image.fromarray(normalize_to_uint8(yz_section))
    xz_img = Image.fromarray(normalize_to_uint8(xz_section))

    # Save as PNGs (or change extension to .tiff if preferred)
    #xy_img.save('/Users/lucasphilipp/Downloads/Spherical Harmonics Expansion/all/PC2/xy/{}.png'.format(i))
    #yz_img.save('/Users/lucasphilipp/Downloads/Spherical Harmonics Expansion/all/PC2/yz/{}.png'.format(i))
    #xz_img.save('/Users/lucasphilipp/Downloads/Spherical Harmonics Expansion/all/PC2/xz/{}.png'.format(i))
    
    #print(i)
    
#napari.view_image(voxelized)

#check accuracy of spherical harmonics expansion
# error = []
# for i in range(1,51):
#     (coeffs, grid_rec), (_, _, grid_input, _) = shparam.get_shcoeffs(image=imread('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/Cell 2/ROI 95.tiff'), lmax=i)
#     error.append(shtools.get_reconstruction_error(grid_rec,grid_input))
    
#     mat=shtools.convert_coeffs_dict_to_matrix(coeffs,i)
#     mesh_rec, _ =shtools.get_reconstruction_from_coeffs(mat)
#     coords = vtknp.vtk_to_numpy(mesh_rec.GetPoints().GetData())

#     # Find bounds of the mesh
#     rmin = (coords.min(axis=0) - 0.5).astype(int)
#     rmax = (coords.max(axis=0) + 0.5).astype(int)

#     # Width, height and depth
#     w = int(2 + (rmax[0] - rmin[0]))
#     h = int(2 + (rmax[1] - rmin[1]))
#     d = int(2 + (rmax[2] - rmin[2]))

#     # Create image data
#     imagedata = vtk.vtkImageData()
#     imagedata.SetDimensions([w, h, d])
#     imagedata.SetExtent(0, w - 1, 0, h - 1, 0, d - 1)
#     imagedata.SetOrigin(rmin)
#     imagedata.AllocateScalars(vtk.VTK_UNSIGNED_CHAR, 1)

#     # Set all values to 1
#     imagedata.GetPointData().GetScalars().FillComponent(0, 1)

#     # Create an empty 3D numpy array to sum up
#     # voxelization of all meshes
#     img = np.zeros((d, h, w), dtype=np.uint8)

#     voxelized = shtools.voxelize_mesh(imagedata=imagedata, shape=(d, h, w), mesh=mesh_rec, origin=rmin)

#     voxelized = voxelized.astype('int8')
#     tifffile.imsave('lmax {}.tiff'.format(i), voxelized, bigtiff=True)


# fig, ax = plt.subplots(1,1, figsize=(3,3))
# ax.plot(np.linspace(1,50), error)
# plt.xlabel('lmax')
# plt.ylabel('Reconstruction Error')
# plt.show()

#lmax = 40 looks good

