import numpy as np
import matplotlib.pyplot as plt
from tifffile import imread
from scipy.ndimage import binary_erosion

#get images and angles and perform spherical harmonic expansion
def spherical_harmonic_expansion(ROI_start,ROI_end,image_path):
    contours = []
    for i in range(ROI_start,ROI_end+1):
        image=imread(image_path+'ROI {}.tiff'.format(i))
        middle_slice = image[:, image.shape[1] // 2, :]
        binary = (middle_slice > 0).astype(np.uint8)
        boundary = binary ^ binary_erosion(binary)
        upper_half_middle_slice = boundary[:,boundary.shape[1] // 2:]
        z, x = np.where(upper_half_middle_slice == 1)
        unique_z = np.unique(z)
        result = []
        for j in unique_z:
            indices = np.where(z == j)[0]
            max_x = np.max(x[indices])
            result.append([j,max_x])
        result = np.array(result)
        contours.append(result)
        print(i)
    return contours

microadriaticum_cell1 = spherical_harmonic_expansion(2,114,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/Cell 1/')
microadriaticum_cell2 = spherical_harmonic_expansion(1,100,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/Cell 2/')
microadriaticum_cell3 = spherical_harmonic_expansion(2,105,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/Cell 3/')

pilosum_cell1 = spherical_harmonic_expansion(4,102,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 1/')
pilosum_cell2 = spherical_harmonic_expansion(22,124,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 2/')
pilosum_cell3 = spherical_harmonic_expansion(27,124,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 3/')

minutum_cell1 = spherical_harmonic_expansion(1,26,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/Cell 1/')
minutum_cell2 = spherical_harmonic_expansion(2,35,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/Cell 2/')
minutum_cell3 = spherical_harmonic_expansion(1,33,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/Cell 3/')

nutricula_cell1 = spherical_harmonic_expansion(2,162,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/brandtodinium cell 1 chromosomes 4nm sampling/')
nutricula_cell2 = spherical_harmonic_expansion(1,100,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/brandtodinium cell 2 chromosomes 4nm sampling/')
nutricula_cell3 = spherical_harmonic_expansion(4,327,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/brandtodinium cell 3 chromosomes 4nm sampling/')

cohnii_cell1 = spherical_harmonic_expansion(2,118,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/Cell 1/')
cohnii_cell2 = spherical_harmonic_expansion(2,187,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/Cell 2/')
cohnii_cell3 = spherical_harmonic_expansion(2,135,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/Cell 3/')

tyrrhenica_cell1 = spherical_harmonic_expansion(2,113,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Ensiculifera tyrrhenica/ensiculifera tyrrhenica chromosomes 4nm sampling/')

microadriaticum = microadriaticum_cell1+microadriaticum_cell2+microadriaticum_cell3
pilosum = pilosum_cell1+pilosum_cell2+pilosum_cell3
minutum = minutum_cell1+minutum_cell2+minutum_cell3
cohnii = cohnii_cell1+cohnii_cell2+cohnii_cell3
nutricula = nutricula_cell1+nutricula_cell2+nutricula_cell3
tyrrhenica = tyrrhenica_cell1

# Plot
plt.figure(figsize=(8, 6))
for i, pair in enumerate(microadriaticum):
    z = pair[:, 0]  # Horizontal axis
    x = pair[:, 1]  # Vertical axis
    z = z-np.mean(z)
    plt.plot(z*4, x*4, 'blue', alpha=0.1)
plt.xlabel('Long axis [nm]')
plt.ylabel('Perpendicular axis [nm]')
plt.ylim((0,300))
plt.xlim((-500,500)) 
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

plt.figure(figsize=(8, 6))
for i, pair in enumerate(pilosum):
    z = pair[:, 0]  # Horizontal axis
    x = pair[:, 1]  # Vertical axis
    z = z-np.mean(z)
    plt.plot(z*4, x*4, 'purple', alpha=0.1)
plt.xlabel('Long axis [nm]')
plt.ylabel('Perpendicular axis [nm]')
plt.ylim((0,300))
plt.xlim((-500,500)) 
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

plt.figure(figsize=(8, 6))
for i, pair in enumerate(minutum):
    z = pair[:, 0]  # Horizontal axis
    x = pair[:, 1]  # Vertical axis
    z = z-np.mean(z)
    plt.plot(z*4, x*4, 'orange', alpha=0.2)
plt.xlabel('Long axis [nm]')
plt.ylabel('Perpendicular axis [nm]')
plt.ylim((0,550)) 
plt.xlim((-750,750)) 
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

plt.figure(figsize=(8, 6))
for i, pair in enumerate(cohnii):
    z = pair[:, 0]  # Horizontal axis
    x = pair[:, 1]  # Vertical axis
    z = z-np.mean(z)
    plt.plot(z*4, x*4, 'pink', alpha=0.2)
plt.xlabel('Long axis [nm]')
plt.ylabel('Perpendicular axis [nm]')
plt.ylim((0,600)) 
plt.xlim((-1000,1000)) 
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

plt.figure(figsize=(8, 6))
for i, pair in enumerate(nutricula):
    z = pair[:, 0]  # Horizontal axis
    x = pair[:, 1]  # Vertical axis
    z = z-np.mean(z)
    plt.plot(z*4, x*4, 'red', alpha=0.075)
plt.xlabel('Long axis [nm]')
plt.ylabel('Perpendicular axis [nm]')
plt.ylim((0,600)) 
plt.xlim((-1000,1000)) 
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

plt.figure(figsize=(8, 6))
for i, pair in enumerate(tyrrhenica):
    z = pair[:, 0]  # Horizontal axis
    x = pair[:, 1]  # Vertical axis
    z = z-np.mean(z)
    plt.plot(z*4, x*4, 'grey', alpha=0.075)
plt.xlabel('Long axis [nm]')
plt.ylabel('Perpendicular axis [nm]')
plt.ylim((0,600))
plt.xlim((-1000,1000)) 
plt.gca().set_aspect('equal', adjustable='box')
plt.show()


