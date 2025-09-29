# This GitHub repository contains content related to:
<img width="1792" height="1005" alt="Title Page" src="https://github.com/user-attachments/assets/c271f475-1a78-4c24-8690-013589a63f57" />

bioRxiv link to pre-print: __________

<img width="1779" alt="Screenshot 2025-06-03 at 3 35 26 PM" src="https://github.com/user-attachments/assets/5db5a207-3f23-4f1b-ace6-e3443c9b967c" />

# Data Availability:

Raw data is available from EMPIAR: https://www.ebi.ac.uk/bioimage-archive/submit/ 

Datasets are available in an interactive format online at: https://petapixelproject.com/mosaics/biology/dinoflagellates/philipp2025/
use the K and L keys to scroll through the z stacks

# Descriptions of code:
Scripts are organized into separate folders according to their functionality.

## 1. Dragonfly Workflows:
Dragonfly 2024.1 See documentation: https://dev.theobjects.com/dragonfly_2024_1_release/contents.html

<img width="1790" height="737" alt="Dragonfly Workflow" src="https://github.com/user-attachments/assets/38f4d324-1184-4d5f-8711-ae89edd2b50b" />

### CropAndRotateAndExportChromosomes.py

Description: A custom menu option to re-orient each ROI (chromosome) in a multi-ROI so the chromosome long axis corresponds to the z-axis. Saves everychromosome as binary 3D .tiff files in a folder.

To load the menu option: Close Dragonfly. Navigate to either of:

C:\Users\Username\AppData\Local/ORS/Dragonfly(version)/pythonUserExtensions/GenericMenuItems
C:\Users\Username\AppData\Local/ORS/Dragonfly(version)/pythonAllUsersExtensions/GenericMenuItems
C:\Users\Username\AppData\Local/ORS/Dragonfly(version)/python/OrsPythonPlugins/OrsGenericMenuItems/menuItems

(it varies from system to system)

To help locate the folder, in Dragonfly try: Utilities -> Open All Users Folder in File Browser.

Place CropAndRotateAndExportChromosomes.py in this folder. Start Dragonfly. The menu option should be accessible after right-clicking on a multi-ROI.

### 3D mask to extract raw EM intensities within volume:
Input the following commands in Dragonfly's python console:

```
Image = *drag and drop image object from object list (top right of screen)*
inverted_ROI = *drag and drop inverted ROI from object list (top right of screen)*   
Image.overwriteValueWithROI(inverted_ROI,0)
Image.setDataDirty()
```
ADD FIGURE OF KAWAGUTII GENOME
ADD FIGURE COMPARING TOROID AND ROD EM INTENSITIES TO CHROMOSOME INTENSITIES AND TO NUCLEOLUS INTENSITIES

### Voronoi simulation using Dragonfly's python console:
Are chromosomes spatially clustered or homogenous? <br />
CITE Dr. Deering's Paper!!!  <br />
Downsample Nucleus ROI, Nucleolus ROI, chromosome multiROI to 40nmx40nmx40nm (from 4nmx4nmx4nm).
Initialize an empty ROI, call it ROI_COM (center of mass) with the geometry of Nucleus ROI (downsampled). The voronoi simulation will continue until all voxels in this geometry are painted so check that the extent of this ROI is a tight box crop around the nucleus and not the whole cell or a larger volume. Extract ROIs from chromosome multiROI (downsampled).
Input the following commands in Dragonfly's python console:
```
ROI_list = *drag and drop chromosome ROIs from drop down menu (top right of screen)*
ROI_COM = *drag and drop empty ROI from drop down menu (top right of screen)*
for i in range(len(ROI_list)):
        ROI = ROI_list[i]
	center_of_mass = ROI.getCenterOfMass(0)
	voxel_coords = ROI.getWorldToVoxelCoordinates(center_of_mass)
	ROI_COM.paintSubset(voxel_coords[0],voxel_coords[1],voxel_coords[2],voxel_coords[0],voxel_coords[1],voxel_coords[2],1,0)
ROI_COM.setDataDirty()
```
Use connected components -> new multiROI 6-connected on ROI_COM to make a multiROI. Name it COM_multiROI.

![seed points](https://github.com/user-attachments/assets/e8c62466-8134-470e-94d2-3348fea41876)

```
# Your MultiROI should contain a cloud of single-voxel points (each as their own unique ROI in the MultiROI) to serve as the seeds for Voronoi tessellation
COM_multiROI = *drag and drop COM_multiROI from drop down menu to set object identifier*

from ORSModel import ConvolutionKernel #I don't think this import statement is necessary, Dragonfly's python console seems to have everything already imported

# Initialize the kernel for normal isotropic Voronoi tessellation
uniConvolutionKernel = ConvolutionKernel()
uniConvolutionKernel.initializeAs3DKernel(5,5,5)  # Uniform dilation 5x5x5 kernel

# Initialize kernel values to form a 5x5x5 spherical kernel for dilation/erosion
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
theOutputChannel = None  # My script overwrites the existing MultiROI, no need to create a new one but it is required as an argument for the method
labelRange = range(1,COM_multiROI.getLabelCount()) 
aSetOfLabels = None  # Initialize the label names within the MultiROI
aSetOfLabels = COM_multiROI.getNonEmptyLabels(aSetOfLabels)  # Get the label names that Dragonfly uses for each ROI in the MultiROI

for i in range(0,25):  # Repeat a sufficient number of times to fill the volume
	COM_multiROI.getDilatedWithKernel(uniConvolutionKernel, aSetOfLabels, 0, aProgress, COM_multiROI)

COM_multiROI.setDataDirty()
```
![cube](https://github.com/user-attachments/assets/2a4d371a-949f-486c-94cc-13de31539afc)

Voronoi cells should not extend outside nucleus or overlap with nucleolus. <br /> 1. A-B of COM_multiROI (A) - Nucleus ROI (downsampled) (B) -> outside_nucleus (save to new). <br /> 2. A-B of COM_multiROI (A) - outside_nucleus (B) = COM_multiROI (overwrite). <br /> 3. A-B of COM_multiROI (A) - Nucleolus ROI (downsampled) (B) = COM_multiROI (overwrite). <br /> 4. Compute volume of voronoi cells. <br /> Voronoi cell volume distribution is sharply peaked -> chromosomes are evenly spaced. Voronoi cell volume distribution is flat/broad -> chromosomes are clustered.

![nucleus crop](https://github.com/user-attachments/assets/74c102f9-c997-45ce-a69a-054084d4836f)
![nucleolus crop with chromosomes and nucleolus](https://github.com/user-attachments/assets/4505bd72-b6a3-4106-924d-5bbcf75de952)
![2D view](https://github.com/user-attachments/assets/1a3fb76d-593f-46e7-94ce-f31c537e4dca)

### Write coloured segmentation into image stack:
Make 3 copies of the image stacks, name them R, G, B respectively
Export all ROIs from multiROI
Select all ROIs from drop down menu
In dragonfly python console
```
ROI_list = *drag and drop ROIs from drop down menu*
R = *drag image copy 1*
B = *drag image copy 2*
G = *drag image copy 3*
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
Select R,G,B images in dropdown. Export->as RGB->.png
<img width="1792" alt="RGB" src="https://github.com/user-attachments/assets/1820a837-e22a-439e-b5bb-b9a921aa8b98" />
<img width="1055" alt="writing_segmentation_into_image" src="https://github.com/user-attachments/assets/825b3bd0-ea48-4792-951e-b6d1d374b500" />
### Protocol for measuring distance of DNA toroids/rods or chromosomes to nuclear membrane/nucleolus/nearest distance to both:
You can identify toroids in multiROIs by computing the Euler characteristic for all objects: https://dev.theobjects.com/dragonfly_2024_1_release/ORSModel/sphinxIndexORSModelClasses/sphinxIndexORSModelMesh.html#ORSModel.ors.Mesh.getEulerCharacteristicNumber getEulerCharacteristicNumber(self, iTIndex: int) → int. 0 is toroid. 2 is topolically equilvalent to a sphere. <br />
segment nucleus <br />
create contour mesh <br />
create ROI from contour mesh (boundary). 1. create empty ROI 2. export mesh to empty ROI <br />
create distance map of boundary ROI (make sure geometry of ROI is of the tight crop around the nucleus, or the same as the geometry of the chromosome multi-ROI and DNA rods and toroids multi-ROI) <br />
on multi-ROI: compute measurements, basic measurements with dataset, min-intensity values, select distance map of nuclear membrane <br />
export as csv

<img width="1784" alt="distance maps" src="https://github.com/user-attachments/assets/b4b29acd-a39e-4c95-ac22-e8631b34467e" />
<img width="1783" alt="Toroid Positioning" src="https://github.com/user-attachments/assets/4a56d90f-f7ad-443a-9879-667c5ae6d42b" />


## 2. Analysis of Chromosome Surface Ridges:
### project_ROI_front_and_save.m

Description: Used to extract surface ridge angles from dinoflagellate chromosomes. Input here is output from CropAndRotateAndExportChromosomes.py. Projects 3D binary ROI onto a middle bisecting plane parallel to the z-axis. Front hemi-volume is projected inwards onto the middle plane and saved as an image (see Fig. 3 a).

### project_ROI_back_and_save.m

Description: Used to extract surface ridge angles from dinoflagellate chromosomes. Input is output from CropAndRotateAndExportChromosomes.py. Projects 3D binary ROI onto a middle bisecting plane parallel to the z-axis. Back hemi-volume is projected inwards onto the middle plane and saved as an image (see Fig. 3 a).

### save_clock_images.m

Description: Used to corroborate manually extracted surface ridge angles. Superimposes surface ridge angles onto projection images.

### make_clock_pdf.m

Description: Synthesizes output from save_clock_images.m into a single .pdf document (one .pdf per cell).

<img width="1769" height="838" alt="angle verification" src="https://github.com/user-attachments/assets/c12db561-cc1e-4bfc-b05b-e9485ee22b55" />
<img width="1721" height="556" alt="angle analysis" src="https://github.com/user-attachments/assets/03d45e98-ae1d-4e22-bc5c-edb6243efca8" />
<img width="1461" height="1101" alt="tilted discs diagram" src="https://github.com/user-attachments/assets/0b28c932-0695-4a09-809f-ad93342dc0c1" />

## 3. Spherical Harmonics Expansion:

### spherical_harmonics_expansion.py
Description: An adaptation of the pipeline originally developed by Viana et al. Paper:  https://doi.org/10.1038/s41586-022-05563-7 \& Github: https://github.com/AllenCell/aics-shparam. <br>
Quantitative analysis of dinoflagellate chromosome shape characteristics. PCA clustering based on shape similarity. Comparison of chromosome shape variation across species.

<img width="1207" height="647" alt="SHE Schematic" src="https://github.com/user-attachments/assets/e727b33f-3e09-4d1b-a7ae-1bb0be5a4846" />


<img width="1790" height="499" alt="SHE PCA" src="https://github.com/user-attachments/assets/9fc170e1-f14d-4c5b-bf15-a69c45452af0" />



<img width="1776" height="849" alt="SHE species comparison" src="https://github.com/user-attachments/assets/7ec8b4da-ce85-4429-869d-df7797df0fb1" />

<img width="1776" height="849" alt="SHE species comparison" src="https://github.com/user-attachments/assets/9d7340ea-97bc-409b-b4dd-f931bb63d7da" />

## 4. Analysis of Chromosome Cross-sections:
### cross_section_um.py
Description: Take 3D binary .tiff dinoflagellate chromosome segmentations and view their cross-sectional profiles.
<img width="1636" height="823" alt="chromosome cross sections" src="https://github.com/user-attachments/assets/e3b3052b-f0f5-40e9-b0e7-5a294061bf4d" />
### make_gif_from_cross_sections.py
Description: Take 3D binary .tiff volumes for shapes along PC1 or PC2 and produce a gif showing how its cross-sections vary with principal component coordinate. <br>
<img src="https://github.com/user-attachments/assets/b2b271be-c362-439d-8402-0ed5299510aa" width="355"/> <img src="https://github.com/user-attachments/assets/738cfdef-5a85-49d4-9f3e-467094be817f" width="600"/>

ADD XY/YZ/XZ LABELS

# Questions:
If you have questions about this repository please contact Lucas Philipp (lucas.philipp@mail.mcgill.ca).
