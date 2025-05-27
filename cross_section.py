import numpy as np
import matplotlib.pyplot as plt
from tifffile import imread
from scipy.ndimage import binary_erosion
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize
import csv, os

def get_cross_section(ROI_start,ROI_end,image_path):
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

def plot_rectellipse(z_all,x_all,zlim_start,z_lim_end,xlim_end,a,b,p):
    # Function to compute rectellipse coordinates from parameters
    def generate_rectellipse(theta, a, b, p):
        denom = (np.abs(np.cos(theta))**p + np.abs(np.sin(theta))**p)**(1/p)
        x = a * np.cos(theta) / denom
        y = b * np.sin(theta) / denom
        return np.vstack((x, y)).T

    z_all = np.array(z_all)
    x_all = np.array(x_all)
    xbins = np.linspace(zlim_start, z_lim_end, 200)
    ybins = np.linspace(0, xlim_end, 60)
    heatmap, zedges, xedges = np.histogram2d(z_all, x_all, bins=[xbins, ybins])
    
    # Generate the full rectellipse using the fitted parameters
    theta_full = np.linspace(0, np.pi, 200)  # Full range of theta

    rectellipse = generate_rectellipse(theta_full, *[a, b, p])
    
    return rectellipse

microadriaticum_cell1 = get_cross_section(2,114,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/Cell 1/')
microadriaticum_cell2 = get_cross_section(1,100,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/Cell 2/')
microadriaticum_cell3 = get_cross_section(2,105,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium microadriaticum/Cell 3/')

pilosum_cell1 = get_cross_section(4,102,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 1/')
pilosum_cell2 = get_cross_section(22,124,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 2/')
pilosum_cell3 = get_cross_section(27,124,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 3/')

minutum_cell1 = get_cross_section(1,26,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/Cell 1/')
minutum_cell2 = get_cross_section(2,35,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/Cell 2/')
minutum_cell3 = get_cross_section(1,33,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/Cell 3/')

nutricula_cell1 = get_cross_section(2,162,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/brandtodinium cell 1 chromosomes 4nm sampling/')
nutricula_cell2 = get_cross_section(1,100,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/brandtodinium cell 2 chromosomes 4nm sampling/')
nutricula_cell3 = get_cross_section(4,327,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/brandtodinium cell 3 chromosomes 4nm sampling/')

cohnii_cell1 = get_cross_section(2,118,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/Cell 1/')
cohnii_cell2 = get_cross_section(2,187,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/Cell 2/')
cohnii_cell3 = get_cross_section(2,135,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/Cell 3/')

tyrrhenica_cell1 = get_cross_section(2,113,'/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Ensiculifera tyrrhenica/ensiculifera tyrrhenica chromosomes 4nm sampling/')

microadriaticum = microadriaticum_cell1+microadriaticum_cell2+microadriaticum_cell3
pilosum = pilosum_cell1+pilosum_cell2+pilosum_cell3
minutum = minutum_cell1+minutum_cell2+minutum_cell3
cohnii = cohnii_cell1+cohnii_cell2+cohnii_cell3
nutricula = nutricula_cell1+nutricula_cell2+nutricula_cell3
tyrrhenica = tyrrhenica_cell1

all_fit_params=[]

z_all = []
x_all = []
plt.figure(figsize=(8, 6))
for i, pair in enumerate(microadriaticum):
    z = pair[:, 0]  # Horizontal axis
    x = pair[:, 1]  # Vertical axis
    z = z-np.mean(z)
    z_all.extend(z * 4)
    x_all.extend(x * 4)
    plt.plot(z*4, x*4, 'blue', alpha=0.1)
rectellipse = plot_rectellipse(z_all,x_all,-500,500,300,220,200,1.8)
plt.plot(rectellipse[:, 0], rectellipse[:, 1], 'k--', linewidth=2)
plt.xlabel('Long axis [nm]')
plt.ylabel('Perpendicular axis [nm]')
plt.ylim((0,300))
plt.xlim((-500,500))
plt.ylim((0,600))
plt.xlim((-1000,1000)) 
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

all_fit_params.append([220,200,1.8,'S. microadriaticum'])

z_all = []
x_all = []
plt.figure(figsize=(8, 6))
for i, pair in enumerate(pilosum):
    z = pair[:, 0]  # Horizontal axis
    x = pair[:, 1]  # Vertical axis
    z = z-np.mean(z)
    z_all.extend(z * 4)
    x_all.extend(x * 4)
    plt.plot(z*4, x*4, 'purple', alpha=0.1)
rectellipse=plot_rectellipse(z_all,x_all,-500,500,300,240,180,2)
plt.plot(rectellipse[:, 0], rectellipse[:, 1], 'k--', linewidth=2)
plt.xlabel('Long axis [nm]')
plt.ylabel('Perpendicular axis [nm]')
plt.ylim((0,300))
plt.xlim((-500,500))
plt.ylim((0,600))
plt.xlim((-1000,1000)) 
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

all_fit_params.append([240,180,2,'S. pilosum'])

z_all = []
x_all = []
plt.figure(figsize=(8, 6))
for i, pair in enumerate(minutum):
    z = pair[:, 0]  # Horizontal axis
    x = pair[:, 1]  # Vertical axis
    z = z-np.mean(z)
    z_all.extend(z * 4)
    x_all.extend(x * 4)
    plt.plot(z*4, x*4, 'orange', alpha=0.2)
rectellipse=plot_rectellipse(z_all,x_all,-750,750,550,320,230,2.2)
plt.plot(rectellipse[:, 0], rectellipse[:, 1], 'k--', linewidth=2)
plt.xlabel('Long axis [nm]')
plt.ylabel('Perpendicular axis [nm]')
plt.ylim((0,550)) 
plt.xlim((-750,750)) 
plt.ylim((0,600))
plt.xlim((-1000,1000)) 
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

all_fit_params.append([320,230,2.2,'S. minutum'])

z_all = []
x_all = []
plt.figure(figsize=(8, 6))
for i, pair in enumerate(cohnii):
    z = pair[:, 0]  # Horizontal axis
    x = pair[:, 1]  # Vertical axis
    z = z-np.mean(z)
    z_all.extend(z * 4)
    x_all.extend(x * 4)
    plt.plot(z*4, x*4, 'pink', alpha=0.2)
rectellipse=plot_rectellipse(z_all,x_all,-1000,1000,600,450,290,2.1)
plt.plot(rectellipse[:, 0], rectellipse[:, 1], 'k--', linewidth=2)
plt.xlabel('Long axis [nm]')
plt.ylabel('Perpendicular axis [nm]')
plt.ylim((0,600)) 
plt.xlim((-1000,1000)) 
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

all_fit_params.append([450,290,2.1,'C. cohnii'])

z_all_short = []
x_all_short = []
z_all_tall = []
x_all_tall = []
plt.figure(figsize=(8, 6))
for i, pair in enumerate(nutricula):
    z = pair[:, 0]  # Horizontal axis
    x = pair[:, 1]  # Vertical axis
    z = z-np.mean(z)
    if np.max(x*4)<200:
        z_all_short.extend(z * 4)
        x_all_short.extend(x * 4)
    else:
        z_all_tall.extend(z * 4)
        x_all_tall.extend(x * 4)
    plt.plot(z*4, x*4, 'red', alpha=0.075)
rectellipse=plot_rectellipse(z_all_short,x_all_short,-1000,1000,600,155,100,1.7)
plt.plot(rectellipse[:, 0], rectellipse[:, 1], 'k--', linewidth=2)
rectellipse=plot_rectellipse(z_all_tall,x_all_tall,-1000,1000,600,440,290,2.2)
plt.plot(rectellipse[:, 0], rectellipse[:, 1], 'k--', linewidth=2)
plt.xlabel('Long axis [nm]')
plt.ylabel('Perpendicular axis [nm]')
plt.ylim((0,600)) 
plt.xlim((-1000,1000)) 
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

all_fit_params.append([155,100,1.7,'B. nutricula inner ring'])
all_fit_params.append([440,290,2.2,'B. nutricula outer ring'])

z_all = []
x_all = []
plt.figure(figsize=(8, 6))
for i, pair in enumerate(tyrrhenica):
    z = pair[:, 0]  # Horizontal axis
    x = pair[:, 1]  # Vertical axis
    z = z-np.mean(z)
    z_all.extend(z * 4)
    x_all.extend(x * 4)
    plt.plot(z*4, x*4, 'grey', alpha=0.075)
rectellipse=plot_rectellipse(z_all,x_all,-1000,1000,600,720,410,2.7)
plt.plot(rectellipse[:, 0], rectellipse[:, 1], 'k--', linewidth=2)
plt.xlabel('Long axis [nm]')
plt.ylabel('Perpendicular axis [nm]')
plt.ylim((0,600))
plt.xlim((-1000,1000)) 
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

all_fit_params.append([720,410,2.7,'E. tyrrhenica'])

file_path = r'/Users/lucasphilipp/Downloads/rectellipse_fit_params.csv'
os.makedirs(os.path.dirname(file_path), exist_ok=True)
with open(file_path, 'w', newline='') as f:
    csv.writer(f).writerows([['semi major axis: a [nm]', 'semi minor axis: b [nm]', 'curvature: p', 'species name']] + all_fit_params)

#parameter sweep to illustrate what a,b,p do
p_sweep = np.logspace(0, 1, 10)
colors = plt.cm.rainbow(np.linspace(0, 1, len(p_sweep)))  # Choose colormap here (e.g., viridis)

fig, ax = plt.subplots()
for i, color in zip(p_sweep, colors):
    rectellipse = plot_rectellipse(z_all, x_all, -500, 500, 300, 400, 250, i)
    ax.plot(rectellipse[:, 0], rectellipse[:, 1], color=color, linewidth=2)

ax.set_xlabel('Long axis [nm]')
ax.set_ylabel('Perpendicular axis [nm]')
ax.set_ylim((0, 300))
ax.set_xlim((-500, 500))
ax.set_aspect('equal', adjustable='box')

#add colorbar
sm = plt.cm.ScalarMappable(cmap="rainbow", norm=plt.Normalize(vmin=np.min(p_sweep), vmax=np.max(p_sweep)))
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, label='Curvature p', shrink=0.5, ticks=[1,2,3,4,5,6,7,8,9,10])
plt.rcParams.update({'font.size': 12})

plt.show()

#relate curvature and aspect ratio
a = [row[0] for row in all_fit_params]
b = [row[1] for row in all_fit_params]
p = [row[2] for row in all_fit_params]
labels = [row[3] for row in all_fit_params]
colors = ['blue', 'purple', 'yellow', 'pink', 'red', 'red', 'grey']

aspect_ratio = np.array(a) / np.array(b)

plt.rcParams.update({'font.size': 20})
plt.figure(figsize=(8,6))
plt.figure(figsize=(8,6))
for i in range(len(aspect_ratio)):
    plt.scatter(aspect_ratio[i], p[i], color=colors[i], s=100)
    plt.text(aspect_ratio[i], p[i], labels[i], fontsize=12, ha='left', va='bottom')

plt.xlabel('Aspect ratio a/b')
plt.ylabel('Curvature p')
plt.show()



