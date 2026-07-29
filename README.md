# This GitHub repository contains content related to:
<img width="1792" height="1006" alt="Screenshot 2026-02-26 at 1 18 15 PM" src="https://github.com/user-attachments/assets/1f78d268-63a3-4d3c-9e4c-33ab8ed80cda" />

bioRxiv link to pre-print: __________ (coming soon) <br />

# Data Availability:
In this study we used a roughly 50:50 mix of original image data, and image data available from the literature. <br />

### Original Image Data:
*Symbiodinium microadriaticum* (3 cells) 4x4x4 nm voxels <br />
*Breviolum minutum* (3 cells) 4x4x4 nm voxels <br />
*Fugacium kawagutii* (3 cells) 4x4x4 nm voxels <br />
*Crypthecodinium cohnii* (3 cells) 4x4x4 nm voxels <br />
Raw image files have been uploaded to EMPIAR under accession code: EMPIAR-13824

### Image Data From The Literature:
*Brandtodinium nutricula* (3 cells) 8x8x8 nm voxels <br />
raw data: https://doi.org/10.1111/1462-2920.15766, EMPIAR-47483651

*Symbiodinium pilosum* (3 cells) 8x8x8 nm voxels <br />
raw data: https://doi.org/10.1038/s41467-021-21314-0, https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BSST575

*Ensiculifera tyrrhenica* (1 cell) 8x8x8 nm voxels <br />
raw data: https://doi.org/10.1242/jcs.261355, EMPIAR-11399

*Kareniaceae sp.* (2 cells) 8x8x8 nm voxels <br />
raw data: https://doi.org/10.1016/j.cub.2025.03.076, EMPIAR-12627

<img width="2520" height="4252" alt="figure_S3" src="https://github.com/user-attachments/assets/601b2069-0a2b-40c8-8bd1-fe03bf1363f4" />

See also our Zenodo data repository: https://zenodo.org/records/21478140, which includes: a spreadsheet of volume measurements; 3D binary TIFF files of vertically oriented chromosomes; 3D binary TIFF files of rods, crescents, and toroids; image galleries of extracted surface ridge angles (SI Data 1-21); 3D videos of *F. kawagutii* chromosome network segmentations (SI Videos 1-3); and *F. kawagutii* skeletonization PLY files.
# ATLAS Browser Based Viewer:
This online infrastructure is developed and maintained by Fibics Incorporated. For more information visit: https://www.petapixelproject.com/about.html <br />

Segmented images and 3D reconstructions of dinoflagellate nuclei are available in interactive format online at: <br /> 
👉👉👉 https://petapixelproject.com/urls/philipp2026 👈👈👈 <br />

No download or installation required. Use the K and L keys to scroll through the z stacks. <br />
Click the play button to view a slide-show highlighting interesting features in the data.

<img width="1792" height="1120" alt="BBV screenshot" src="https://github.com/user-attachments/assets/26b8ea76-92f7-403c-86cd-70cd78aa6f8b" />

# Descriptions of code:
Scripts are organized into separate folders according to their functionality.

## 1. Dragonfly Workflows:
Dragonfly is a free software for academics: https://dragonfly.comet.tech/

<img width="1790" height="737" alt="Dragonfly Workflow" src="https://github.com/user-attachments/assets/38f4d324-1184-4d5f-8711-ae89edd2b50b" />

### CropAndRotateAndExportChromosomes.py
Description: A custom menu option to re-orient each ROI (chromosome) in a multi-ROI so the chromosome long axis corresponds to the z-axis. Saves everychromosome as binary 3D .tiff files in a folder. There are two version of this script depending on which Dragonfly version you are using (2022.2 or 2024.1).

1. To load the menu option: Close Dragonfly. Navigate to either of:
C:\Users\Username\AppData\Local/ORS/Dragonfly(version)/pythonUserExtensions/GenericMenuItems
C:\Users\Username\AppData\Local/ORS/Dragonfly(version)/pythonAllUsersExtensions/GenericMenuItems
C:\Users\Username\AppData\Local/ORS/Dragonfly(version)/python/OrsPythonPlugins/OrsGenericMenuItems/menuItems <br />
(it varies from system to system)

To help locate the folder, in Dragonfly try: Utilities -> Open All Users Folder in File Browser.

2. Place CropAndRotateAndExportChromosomes.py in this folder.
3. Start Dragonfly. The menu option should be accessible after right-clicking on a multi-ROI.

### 3D mask to extract raw EM intensities within volume:
Input the following commands in Dragonfly's python console:
```python
Image = #drag and drop image object here from object list (top right of screen)
inverted_ROI = #drag and drop inverted ROI here from object list (top right of screen)   
Image.overwriteValueWithROI(inverted_ROI,0)
Image.setDataDirty()
```
<!-- ADD FIGURE COMPARING TOROID AND ROD EM INTENSITIES TO CHROMOSOME INTENSITIES AND TO NUCLEOLUS INTENSITIES -->

### Voronoi simulation using Dragonfly's python console:
Quantifies whether chromosomes are spatially clustered or located homogenously throughout the nucleus. <br />

1. Downsample Nucleus ROI, Nucleolus ROI, chromosome multiROI to 40nmx40nmx40nm (from 4nmx4nmx4nm).
2. Initialize an empty ROI, call it ROI_COM (center of mass) with the geometry of Nucleus ROI (downsampled). The voronoi simulation will continue until all voxels in this geometry are painted so check that the extent of this ROI is a tight box crop around the nucleus and not the whole cell or a larger volume.
3. Extract ROIs from chromosome multiROI (downsampled).
4. Input the following commands in Dragonfly's python console:
```
ROI_list = #drag and drop chromosome ROIs here from drop down menu (top right of screen)
ROI_COM = #drag and drop empty ROI here from drop down menu (top right of screen)
for i in range(len(ROI_list)):
        ROI = ROI_list[i]
	center_of_mass = ROI.getCenterOfMass(0)
	voxel_coords = ROI.getWorldToVoxelCoordinates(center_of_mass)
	ROI_COM.paintSubset(voxel_coords[0],voxel_coords[1],voxel_coords[2],voxel_coords[0],voxel_coords[1],voxel_coords[2],1,0)
ROI_COM.setDataDirty()
```
<img src="https://github.com/user-attachments/assets/e8c62466-8134-470e-94d2-3348fea41876" width="700"/>

5. Use connected components -> new multiROI 6-connected on ROI_COM to make a multiROI. Name it COM_multiROI.
6. Execute the following commands (credit: Dr. Joseph Deering. See https://pubs.acs.org/doi/10.1021/acsnano.5c11922) in Dragonfly's python console:   
<details>
  <summary><strong>Click to expand code block</strong></summary>

  ```python
# Your MultiROI should contain a cloud of single-voxel points (each as their own unique ROI in the MultiROI) to serve as the seeds for Voronoi tessellation
COM_multiROI = #drag and drop COM_multiROI here from drop down menu

from ORSModel import ConvolutionKernel

# Initialize the kernel for normal isotropic Voronoi tessellation
uniConvolutionKernel = ConvolutionKernel()
uniConvolutionKernel.initializeAs3DKernel(5,5,5)  # Uniform dilation 5x5x5 kernel

# Initialize kernel values to form a 5x5x5 spherical kernel for dilation
uniConvolutionKernel.setValueAt(0,0,0,0)
uniConvolutionKernel.setValueAt(1,0,0,0)
uniConvolutionKernel.setValueAt(2,0,0,0)
uniConvolutionKernel.setValueAt(3,0,0,0)
uniConvolutionKernel.setValueAt(4,0,0,0)
uniConvolutionKernel.setValueAt(0,1,0,0)
uniConvolutionKernel.setValueAt(1,1,0,0)
uniConvolutionKernel.setValueAt(2,1,0,0)
uniConvolutionKernel.setValueAt(3,1,0,0)
uniConvolutionKernel.setValueAt(4,1,0,0)
uniConvolutionKernel.setValueAt(0,2,0,0)
uniConvolutionKernel.setValueAt(1,2,0,0)
uniConvolutionKernel.setValueAt(2,2,0,1)
uniConvolutionKernel.setValueAt(3,2,0,0)
uniConvolutionKernel.setValueAt(4,2,0,0)
uniConvolutionKernel.setValueAt(0,3,0,0)
uniConvolutionKernel.setValueAt(1,3,0,0)
uniConvolutionKernel.setValueAt(2,3,0,0)
uniConvolutionKernel.setValueAt(3,3,0,0)
uniConvolutionKernel.setValueAt(4,3,0,0)
uniConvolutionKernel.setValueAt(0,4,0,0)
uniConvolutionKernel.setValueAt(1,4,0,0)
uniConvolutionKernel.setValueAt(2,4,0,0)
uniConvolutionKernel.setValueAt(3,4,0,0)
uniConvolutionKernel.setValueAt(4,4,0,0)

uniConvolutionKernel.setValueAt(0,0,1,0)
uniConvolutionKernel.setValueAt(1,0,1,0)
uniConvolutionKernel.setValueAt(2,0,1,0)
uniConvolutionKernel.setValueAt(3,0,1,0)
uniConvolutionKernel.setValueAt(4,0,1,0)
uniConvolutionKernel.setValueAt(0,1,1,0)
uniConvolutionKernel.setValueAt(1,1,1,1)
uniConvolutionKernel.setValueAt(2,1,1,1)
uniConvolutionKernel.setValueAt(3,1,1,1)
uniConvolutionKernel.setValueAt(4,1,1,0)
uniConvolutionKernel.setValueAt(0,2,1,0)
uniConvolutionKernel.setValueAt(1,2,1,1)
uniConvolutionKernel.setValueAt(2,2,1,1)
uniConvolutionKernel.setValueAt(3,2,1,1)
uniConvolutionKernel.setValueAt(4,2,1,0)
uniConvolutionKernel.setValueAt(0,3,1,0)
uniConvolutionKernel.setValueAt(1,3,1,1)
uniConvolutionKernel.setValueAt(2,3,1,1)
uniConvolutionKernel.setValueAt(3,3,1,1)
uniConvolutionKernel.setValueAt(4,3,1,0)
uniConvolutionKernel.setValueAt(0,4,1,0)
uniConvolutionKernel.setValueAt(1,4,1,0)
uniConvolutionKernel.setValueAt(2,4,1,0)
uniConvolutionKernel.setValueAt(3,4,1,0)
uniConvolutionKernel.setValueAt(4,4,1,0)

uniConvolutionKernel.setValueAt(0,0,2,0)
uniConvolutionKernel.setValueAt(1,0,2,0)
uniConvolutionKernel.setValueAt(2,0,2,1)
uniConvolutionKernel.setValueAt(3,0,2,0)
uniConvolutionKernel.setValueAt(4,0,2,0)
uniConvolutionKernel.setValueAt(0,1,2,0)
uniConvolutionKernel.setValueAt(1,1,2,1)
uniConvolutionKernel.setValueAt(2,1,2,1)
uniConvolutionKernel.setValueAt(3,1,2,1)
uniConvolutionKernel.setValueAt(4,1,2,0)
uniConvolutionKernel.setValueAt(0,2,2,1)
uniConvolutionKernel.setValueAt(1,2,2,1)
uniConvolutionKernel.setValueAt(2,2,2,1)
uniConvolutionKernel.setValueAt(3,2,2,1)
uniConvolutionKernel.setValueAt(4,2,2,1)
uniConvolutionKernel.setValueAt(0,3,2,0)
uniConvolutionKernel.setValueAt(1,3,2,1)
uniConvolutionKernel.setValueAt(2,3,2,1)
uniConvolutionKernel.setValueAt(3,3,2,1)
uniConvolutionKernel.setValueAt(4,3,2,0)
uniConvolutionKernel.setValueAt(0,4,2,0)
uniConvolutionKernel.setValueAt(1,4,2,0)
uniConvolutionKernel.setValueAt(2,4,2,1)
uniConvolutionKernel.setValueAt(3,4,2,0)
uniConvolutionKernel.setValueAt(4,4,2,0)

uniConvolutionKernel.setValueAt(0,0,3,0)
uniConvolutionKernel.setValueAt(1,0,3,0)
uniConvolutionKernel.setValueAt(2,0,3,0)
uniConvolutionKernel.setValueAt(3,0,3,0)
uniConvolutionKernel.setValueAt(4,0,3,0)
uniConvolutionKernel.setValueAt(0,1,3,0)
uniConvolutionKernel.setValueAt(1,1,3,1)
uniConvolutionKernel.setValueAt(2,1,3,1)
uniConvolutionKernel.setValueAt(3,1,3,1)
uniConvolutionKernel.setValueAt(4,1,3,0)
uniConvolutionKernel.setValueAt(0,2,3,0)
uniConvolutionKernel.setValueAt(1,2,3,1)
uniConvolutionKernel.setValueAt(2,2,3,1)
uniConvolutionKernel.setValueAt(3,2,3,1)
uniConvolutionKernel.setValueAt(4,2,3,0)
uniConvolutionKernel.setValueAt(0,3,3,0)
uniConvolutionKernel.setValueAt(1,3,3,1)
uniConvolutionKernel.setValueAt(2,3,3,1)
uniConvolutionKernel.setValueAt(3,3,3,1)
uniConvolutionKernel.setValueAt(4,3,3,0)
uniConvolutionKernel.setValueAt(0,4,3,0)
uniConvolutionKernel.setValueAt(1,4,3,0)
uniConvolutionKernel.setValueAt(2,4,3,0)
uniConvolutionKernel.setValueAt(3,4,3,0)
uniConvolutionKernel.setValueAt(4,4,3,0)

uniConvolutionKernel.setValueAt(0,0,4,0)
uniConvolutionKernel.setValueAt(1,0,4,0)
uniConvolutionKernel.setValueAt(2,0,4,0)
uniConvolutionKernel.setValueAt(3,0,4,0)
uniConvolutionKernel.setValueAt(4,0,4,0)
uniConvolutionKernel.setValueAt(0,1,4,0)
uniConvolutionKernel.setValueAt(1,1,4,0)
uniConvolutionKernel.setValueAt(2,1,4,0)
uniConvolutionKernel.setValueAt(3,1,4,0)
uniConvolutionKernel.setValueAt(4,1,4,0)
uniConvolutionKernel.setValueAt(0,2,4,0)
uniConvolutionKernel.setValueAt(1,2,4,0)
uniConvolutionKernel.setValueAt(2,2,4,1)
uniConvolutionKernel.setValueAt(3,2,4,0)
uniConvolutionKernel.setValueAt(4,2,4,0)
uniConvolutionKernel.setValueAt(0,3,4,0)
uniConvolutionKernel.setValueAt(1,3,4,0)
uniConvolutionKernel.setValueAt(2,3,4,0)
uniConvolutionKernel.setValueAt(3,3,4,0)
uniConvolutionKernel.setValueAt(4,3,4,0)
uniConvolutionKernel.setValueAt(0,4,4,0)
uniConvolutionKernel.setValueAt(1,4,4,0)
uniConvolutionKernel.setValueAt(2,4,4,0)
uniConvolutionKernel.setValueAt(3,4,4,0)
uniConvolutionKernel.setValueAt(4,4,4,0)

aProgress = None  # No progress bar
theOutputChannel = None  #overwrites the existing MultiROI
labelRange = range(1,COM_multiROI.getLabelCount()) 
aSetOfLabels = None
aSetOfLabels = COM_multiROI.getNonEmptyLabels(aSetOfLabels)  # Get the label names that Dragonfly uses for each ROI in the MultiROI

for i in range(0,25):  # Repeat a sufficient number of times to fill the volume
	COM_multiROI.getDilatedWithKernel(uniConvolutionKernel, aSetOfLabels, 0, aProgress, COM_multiROI)

COM_multiROI.setDataDirty()
```
</details> 

<img src="https://github.com/user-attachments/assets/2a4d371a-949f-486c-94cc-13de31539afc" width="700"/>

Voronoi cells should not extend outside nucleus or overlap with nucleolus.

8. A-B of COM_multiROI (A) = Nucleus ROI (downsampled) (B) = outside_nucleus -> (save to new).
9. A-B of COM_multiROI (A) = outside_nucleus (B) = COM_multiROI ->(overwrite).
10. A-B of COM_multiROI (A) = Nucleolus ROI (downsampled) (B) = COM_multiROI -> (overwrite).
<img src="https://github.com/user-attachments/assets/1a3fb76d-593f-46e7-94ce-f31c537e4dca" width="700"/>

11. Compute volume of voronoi cells.
Voronoi cell volume distribution is sharply peaked -> chromosomes are evenly spaced. Voronoi cell volume distribution is flat/broad -> chromosomes are clustered.

### Write coloured segmentation into image stack:
1. Make 3 copies of the image stacks, name them R, G, B respectively.
2. Export all ROIs from multiROI.
3. Select all ROIs from drop down menu, drag and drop into Dragonfly's python console and initialize as ROI_list.
4. Execute the following commands in python console:
```python
ROI_list = *drag and drop ROIs from drop down menu*
R = #drag image copy 1 here from object list (top right of screen)
B = #drag image copy 2 here from object list (top right of screen)
G = #drag image copy 3 here from object list (top right of screen)
for i in range(len(ROI_list)):
        ROI = ROI_list[i]
	color = ROI.getInitialColor()
	R.overwriteValueWithROI(ROI,round(color.getRed()*255))
	R.setDataDirty()
	G.overwriteValueWithROI(ROI,round(color.getGreen()*255))
	G.setDataDirty()
	B.overwriteValueWithROI(ROI,round(color.getBlue()*255))
	B.setDataDirty()
```
6. Select R,G,B images in dropdown. Right click. Export->as RGB->.png
<img width="1792" alt="RGB" src="https://github.com/user-attachments/assets/1820a837-e22a-439e-b5bb-b9a921aa8b98" />
<img width="800" alt="writing_segmentation_into_image" src="https://github.com/user-attachments/assets/825b3bd0-ea48-4792-951e-b6d1d374b500" />

## 2. Analysis of Chromosome Surface Ridges:
<img width="1769" height="582" alt="Screenshot 2026-05-17 at 1 09 06 PM" src="https://github.com/user-attachments/assets/aa98072c-06bf-4be3-8ed5-b0566cb62b81" />

### project_ROI_front_and_save.m
Description: Used to extract surface ridge angles from dinoflagellate chromosomes. Input here is output from CropAndRotateAndExportChromosomes.py. Projects 3D binary ROI onto a middle bisecting plane parallel to the z-axis. Front hemi-volume is projected inwards onto the middle plane and saved as an image.

### project_ROI_back_and_save.m
Description: Used to extract surface ridge angles from dinoflagellate chromosomes. Input is output from CropAndRotateAndExportChromosomes.py. Projects 3D binary ROI onto a middle bisecting plane parallel to the z-axis. Back hemi-volume is projected inwards onto the middle plane and saved as an image.

### save_clock_images.m
Description: Used to corroborate manually extracted surface ridge angles. Superimposes surface ridge angles onto projection images.

### make_clock_pdf.m
Description: Synthesizes output from save_clock_images.m into a single .pdf document (one .pdf per cell).

<img width="1737" height="626" alt="Screenshot 2026-04-13 at 3 16 53 PM" src="https://github.com/user-attachments/assets/1c612a78-d9aa-44fc-ac04-7971cef73093" />

## 3. Spherical Harmonics Expansion:

After installing the aicsshparam python package. A small modification was made to shparam.py located at /opt/anaconda3/lib/python3.13/site-packages/aicsshparam/. <br>
A flag is now returned when the spherical harmonics reconstruction fails.

Replace the "if statement" after line 132 in shparam.py with:
```
flag = False
   if not image_[tuple([int(u) for u in centroid[::-1]])]:
       warnings.warn(
           "Mesh centroid seems to fall outside the object. This indicates\
       the mesh may not be a manifold suitable for spherical harmonics\
       parameterization."
       )
       flag = True

   return (coeffs_dict, grid_rec), (image_, mesh, grid_down, transform), flag
```

### spherical_harmonics_expansion.py
Description: An adaptation of the pipeline originally developed by Viana et al. Paper:  https://doi.org/10.1038/s41586-022-05563-7. Github: https://github.com/AllenCell/aics-shparam. <br>
<br>
<img width="1612" height="382" alt="Screenshot 2026-07-08 at 5 44 51 AM" src="https://github.com/user-attachments/assets/cc9ec604-ec17-4e51-acd2-c42eca67c2ab" />
<br>
<br>

1. PCA clustering based on shape similarity.
2. Quantitative analysis of dinoflagellate chromosome shape characteristics.
<img width="1725" height="523" alt="Screenshot 2026-07-08 at 1 37 07 PM" src="https://github.com/user-attachments/assets/370069b9-9e18-477d-a39b-dac39b0bf265" />
<br>
<br>

3. Comparison of chromosome shape variation across species.
<img width="1698" height="480" alt="Screenshot 2026-07-08 at 1 37 37 PM" src="https://github.com/user-attachments/assets/874a943c-172b-4150-82a2-16affef4269e" />

## 4. Analysis of Chromosome Cross-sections:
### cross_section.py
Description: Take 3D binary .tiff of segmented dinoflagellate chromosome and view their cross-sectional profiles.
<img width="929" height="500" alt="figure_S6_b" src="https://github.com/user-attachments/assets/447a4959-0963-433c-8806-57e628581948" />

### make_gif_from_cross_sections.py
Description: Take 3D binary .tiff volumes for shapes along PC1 or PC2 and produce a gif showing how its XY, XZ, YZ cross-sections vary with principal component coordinate. <br>
<img src="https://github.com/user-attachments/assets/099b7f07-7e86-42c6-a08b-7fedd0cc2157" width="270"/> <img src="https://github.com/user-attachments/assets/555a01f3-9053-48d9-9b38-f88ec9914007" width="448"/>

## 5. Skeletonization Post-Processing:
### plot_skeleton.py
Description: Plot output skeleton from Voxel Cores (https://github.com/danielyan86129/voxel_ma) and Erosion Thickness (https://github.com/danielyan86129/ET). Colors skeleton according to maximally inscribed ball radii, r, centred at points along the skeleton. See Methods section in paper for more details.

Schematic (left) modified with permission from: Tagliasacchi, A., et al. (2016). 3D Skeletons: A State‐of‐the‐Art Report. In Computer Graphics Forum (Vol. 35, No. 2, pp. 573-597).
<img width="1769" height="739" alt="Screenshot 2026-05-17 at 1 24 05 PM" src="https://github.com/user-attachments/assets/40cce40b-8383-4d23-b4c7-fa8803115b3f" />
<img width="1768" height="739" alt="Screenshot 2026-05-17 at 1 24 28 PM" src="https://github.com/user-attachments/assets/aa69e184-0e63-44ac-ae0e-7f6cc5d42e6e" />

### plot_skeleton_connected_components.py
Description: Subdivide skeleton into thick tubes and thin bridges by thresholding the max inscribed ball radii. Uses connected components to further subdivide skeleton at junction nodes with connectivity degree >= 3.
<img width="1774" height="566" alt="Screenshot 2026-05-17 at 1 25 12 PM" src="https://github.com/user-attachments/assets/1250f1a2-61f0-4fa6-b395-f0570acc2649" />

# Questions:
If you have questions about this repository please contact Lucas Philipp (lucas.philipp@mail.mcgill.ca).
