# This GitHub repository contains the code to recreate the analyses of:
Philipp, L., Ittah E., Schumann D., Reznikov N., Weber S. C., 3D Electron Microscopy of Dinoflagellate Chromosomes Across Many Species. In review.

bioRxiv link to pre-print: __________

# Data Availability:

https://www.ebi.ac.uk/bioimage-archive/submit/ 

Datasets are available in an interactive format online at: https://www.petapixelproject.com/

# Descriptions of code:

### CropAndRotateAndExportChromosomes.py

Description: A custom menu option to re-orient each ROI (chromosome) in a multi-ROI so the chromosome long axis corresponds to the z-axis. Saves everychromosome as binary 3D .tiff files in a folder.

To load the menu option: Close Dragonfly. Navigate to: C:\Users\Username\AppData\Local\ORS\Dragonfly2022.2\pythonUserExtensions\GenericMenultems. Place CropAndRotateAndExportChromosomes.py in this folder. Start Dragonfly. The menu option should be accessible after right-clicking on a multi-ROI.

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

# Questions:
If you have questions about this repository please contact Lucas Philipp (lucas.philipp@mail.mcgill.ca).
