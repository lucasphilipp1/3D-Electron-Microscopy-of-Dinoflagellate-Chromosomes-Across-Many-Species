# This GitHub repository contains code related to:
Philipp, L.*, Ittah, E.*, Schumann D., de Fourestier J., Reznikov N., Weber S. C., 3D Electron Microscopy of Dinoflagellate Chromosomes Across Many Species. In review. (* means equal contribution).

bioRxiv link to pre-print: __________

# Data Availability:

Raw data is available from EMPIAR: https://www.ebi.ac.uk/bioimage-archive/submit/ 

Datasets are available in an interactive format online at: https://www.petapixelproject.com/

use the K and L keys to scroll through the z stacks

_Symbiodinium microadriaticum_ <br />
https://petapixelproject.com/mosaics/biology/dinoflagellates/microadriaticum/cell-1/ROI
https://petapixelproject.com/mosaics/biology/dinoflagellates/microadriaticum/cell-1/Overview

https://petapixelproject.com/mosaics/biology/dinoflagellates/microadriaticum/cell-2/ROI
https://petapixelproject.com/mosaics/biology/dinoflagellates/microadriaticum/cell-2/Overview

https://petapixelproject.com/mosaics/biology/dinoflagellates/microadriaticum/cell-3/ROI
https://petapixelproject.com/mosaics/biology/dinoflagellates/microadriaticum/cell-3/Overview

_Symbiodinium (Breviolum) minutum_ <br />
https://petapixelproject.com/mosaics/biology/dinoflagellates/minutum/cell-1/ROI
https://petapixelproject.com/mosaics/biology/dinoflagellates/minutum/cell-1/Overview

https://petapixelproject.com/mosaics/biology/dinoflagellates/minutum/cell-2/ROI
https://petapixelproject.com/mosaics/biology/dinoflagellates/minutum/cell-2/Overview

https://petapixelproject.com/mosaics/biology/dinoflagellates/minutum/cell-3/ROI
https://petapixelproject.com/mosaics/biology/dinoflagellates/minutum/cell-3/Overview

_Symbiodinium (Fugacium) kawagutii_ <br />
https://petapixelproject.com/mosaics/biology/dinoflagellates/kawagutii/cell-1/ROI
https://petapixelproject.com/mosaics/biology/dinoflagellates/kawagutii/cell-1/Overview

https://petapixelproject.com/mosaics/biology/dinoflagellates/kawagutii/cell-2/ROI
https://petapixelproject.com/mosaics/biology/dinoflagellates/kawagutii/cell-2/Overview

https://petapixelproject.com/mosaics/biology/dinoflagellates/kawagutii/cell-3/ROI
https://petapixelproject.com/mosaics/biology/dinoflagellates/kawagutii/cell-3/Overview

_Crypthecodinium cohnii_ <br />
https://petapixelproject.com/mosaics/biology/dinoflagellates/cohnii/cell-1/ROI
https://petapixelproject.com/mosaics/biology/dinoflagellates/cohnii/cell-1/Overview

https://petapixelproject.com/mosaics/biology/dinoflagellates/cohnii/cell-2/ROI/
https://petapixelproject.com/mosaics/biology/dinoflagellates/cohnii/cell-2/Overview/

https://petapixelproject.com/mosaics/biology/dinoflagellates/cohnii/cell-3/ROI
https://petapixelproject.com/mosaics/biology/dinoflagellates/cohnii/cell-3/Overview

_Symbiodinium pilosum_ <br />
https://petapixelproject.com/mosaics/biology/dinoflagellates/pilosum/cell-1/ROI
https://petapixelproject.com/mosaics/biology/dinoflagellates/pilosum/cell-1/Overview

https://petapixelproject.com/mosaics/biology/dinoflagellates/pilosum/cell-2/ROI
https://petapixelproject.com/mosaics/biology/dinoflagellates/pilosum/cell-2/Overview

https://petapixelproject.com/mosaics/biology/dinoflagellates/pilosum/cell-3/ROI
https://petapixelproject.com/mosaics/biology/dinoflagellates/pilosum/cell-3/Overview

_Brandtodinium nutricula_ <br />
https://petapixelproject.com/mosaics/biology/dinoflagellates/brandtodinium/cell-1/ROI/
https://petapixelproject.com/mosaics/biology/dinoflagellates/brandtodinium/cell-1/Overview/

https://petapixelproject.com/mosaics/biology/dinoflagellates/brandtodinium/cell-2/ROI/
https://petapixelproject.com/mosaics/biology/dinoflagellates/brandtodinium/cell-2/Overview/

https://petapixelproject.com/mosaics/biology/dinoflagellates/brandtodinium/cell-3/ROI/
https://petapixelproject.com/mosaics/biology/dinoflagellates/brandtodinium/cell-3/Overview/

# Descriptions of code:

### CropAndRotateAndExportChromosomes.py

Description: A custom menu option to re-orient each ROI (chromosome) in a multi-ROI so the chromosome long axis corresponds to the z-axis. Saves everychromosome as binary 3D .tiff files in a folder.

To load the menu option: Close Dragonfly. Navigate to either of:

C:\Users\Username\AppData\Local/ORS/Dragonfly(version)/pythonUserExtensions/GenericMenuItems
C:\Users\Username\AppData\Local/ORS/Dragonfly(version)/pythonAllUsersExtensions/GenericMenuItems
C:\Users\Username\AppData\Local/ORS/Dragonfly(version)/python/OrsPythonPlugins/OrsGenericMenuItems/menuItems

(it varies from system to system)

To help locate the folder, in Dragonfly try: Utilities -> Open All Users Folder in File Browser.

Place CropAndRotateAndExportChromosomes.py in this folder. Start Dragonfly. The menu option should be accessible after right-clicking on a multi-ROI.

NEED TO UPDATE CropAndRotateAndExportChromosomes.py, PyQt5 -> PyQt6 import issue 

### helix.m

Description: Create 3D binary images of left- and right- handed helices. Used as a positive control to ensure CropAndRotateAndExportChromosomes.py does not alter the helical handedness of processed ROIs.

### project_ROI_front_and_save.m

Description: Used to extract surface ridge angles from dinoflagellate chromosomes. Input here is output from CropAndRotateAndExportChromosomes.py. Projects 3D binary ROI onto a middle bisecting plane parallel to the z-axis. Front hemi-volume is projected inwards onto the middle plane and saved as an image (see Fig. 3 a).

### project_ROI_back_and_save.m

Description: Used to extract surface ridge angles from dinoflagellate chromosomes. Input is output from CropAndRotateAndExportChromosomes.py. Projects 3D binary ROI onto a middle bisecting plane parallel to the z-axis. Back hemi-volume is projected inwards onto the middle plane and saved as an image (see Fig. 3 a).

### save_clock_images.m

Description: Used to corroborate manually extracted surface ridge angles. Superimposes surface ridge angles onto projection images.

### make_clock_pdf.m

Description: Synthesizes output from save_clock_images.m into a single .pdf document (one .pdf per cell).

Dragonfly 2024.1 See documentation: https://dev.theobjects.com/dragonfly_2024_1_release/contents.html
### 3D mask to extract raw EM intensities within volume:
Input the following commands in Dragonfly's python console:

```
Image = *drag and drop image object from object list (top right menu)*
inverted_ROI = *drag and drop inverted ROI from object list (top right menu)*   
Image.overwriteValueWithROI(inverted_ROI,0)
Image.setDataDirty()
```
### Voronoi simulation (Are chromosomes spatially clustered or homogenous?):
CITE Dr. Deering's Paper!!!  <br />
Downsample Nucleus ROI, Nucleolus ROI, chromosome multiROI to 40nmx40nmx40nm (from 4nmx4nmx4nm).
Initialize an empty ROI, call it ROI_COM (center of mass) with the geometry of Nucleus ROI (downsampled). The voronoi simulation will continue until all voxels in this geometry are painted so check that the extent of this ROI is a tight box crop around the nucleus and not the whole cell or a larger volume.
Input the following commands in Dragonfly's python console:
```
ROI_list = *drag and drop ROIs from drop down menu*
for i in range(len(ROI_list)):
        ROI = ROI_list[i]
	center_of_mass = ROI.getCenterOfMass(0)
	voxel_coords = ROI.getWorldToVoxelCoordinates(center_of_mass)
	ROI_COM.paintSubset(voxel_coords[0],voxel_coords[1],voxel_coords[2],voxel_coords[0],voxel_coords[1],voxel_coords[2],1,0)
```
Use connected components -> new multiROI 6-connected on ROI_COM to make a multiROI. Name it COM_multiROI.
```
from ORSModel import ConvolutionKernel

# Your MultiROI should contain a cloud of single-voxel points (each as their own unique ROI in the MultiROI) to serve as the seeds for Voronoi tessellation
from ORSModel import ConvolutionKernel

# Your MultiROI should contain a cloud of single-voxel points (each as their own unique ROI in the MultiROI) to serve as the seeds for Voronoi tessellation
COM_multiROI = *drag and drop COM_multiROI from drop down menu to set object identifier*

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

aProgress = None  # No progress bar for the process, have not figured out this aspect of Dragonfly yet
theOutputChannel = None  # My script overwrites the existing MultiROI, no need to create a new one but it is required as an argument for the method
labelRange = range(1,COM_multiROI.getLabelCount()) 
aSetOfLabels = None  # Initialize the label names within the MultiROI
aSetOfLabels = COM_multiROI.getNonEmptyLabels(aSetOfLabels)  # Get the label names that Dragonfly uses for each ROI in the MultiROI

for i in range(0,20):  # Repeat a sufficient number of times to fill the volume
	COM_multiROI.getDilatedWithKernel(uniConvolutionKernel, aSetOfLabels, 0, aProgress, COM_multiROI)

COM_multiROI.setDataDirty()
```
Then take the union of Nucleus ROI (downsampled) and the ROI_COM, remove the intersection of ROI_COM and Nucleolus ROI (downsampled).
Compute volume of voronoi cells. Voronoi cell volume distribution is sharply peaked -> chromosomes are evenly spaced. Voronoi cell volume distribution is flat/broad -> chromosomes are clustered.
### Write coloured segmentation into image:
Make 3 copies of the image, name them R, G, B respectively
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
# Questions:
If you have questions about this repository please contact Lucas Philipp (lucas.philipp@mail.mcgill.ca).
