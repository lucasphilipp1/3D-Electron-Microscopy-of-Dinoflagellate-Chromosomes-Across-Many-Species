"""
Take a MultiROI containing chromosomes, crop them, rotate them to be aligned with the Z axis and export as 3D TIFF image.
 Works with any multi-ROI objects (not just chromosomes). Specific to Dragonfly version 2024.1

:author: Lucas Philipp, Benjamin Z. Rudski
:contact: 
:email: lucas.philipp@mail.mcgill.ca; benjamin.rudski@mail.mcgill.ca
:organization: McGill University
:address: 
:copyright: 
:date: Feb 22 2024 09:23
:dragonflyVersion: 2024.1
:UUID: ee6b6e04d18d11ee93ecc45ab1da15d2
"""

__version__ = '1.0.0'

import os.path
import numpy as np
import skimage
import skimage.io

from ORSServiceClass.menuItems.contextualMenuItem import ContextualMenuItem
from ORSServiceClass.actionAndMenu.menu import Menu
from ORSServiceClass.decorators.infrastructure import interfaceMethod
from ORSServiceClass.ORSWidget.orsfiledialog import OrsFileDialog
from OrsHelpers.ListHelper import ListHelper
from OrsHelpers.roihelper import ROIHelper
from ORSModel.ors import MultiROI, VisualBox
from OrsPythonPlugins.orsCreateABoxFromCurrentObjectBoxOrBoundingBox_df18ee80548711e8b181107b441a0862 import orsCreateABoxFromCurrentObjectBoxOrBoundingBox_df18ee80548711e8b181107b441a0862
from OrsPythonPlugins.OrsGenericMenuItems.menuItems.takeShapeOfMinimalBox_264ff9b8a66511e88120107b441a0862 import \
    TakeShapeOfMinimalBox_264ff9b8a66511e88120107b441a0862
from OrsPythonPlugins.OrsGenericMenuItems.menuItems.extractBoxFromChannel_bdfbc9a8110b11e8819130e37aed33ab import \
    extractBoxFromChannel_bdfbc9a8110b11e8819130e37aed33ab
from OrsPythonPlugins.OrsGenericMenuItems.menuItems.exportROIAsBinary_cfdc4c58867011e888c684a6c8f5618e import exportROIAsBinary_cfdc4c58867011e888c684a6c8f5618e
from OrsPythonPlugins.OrsDatasetInvertor.OrsDatasetInvertor import OrsDatasetInvertor
from OrsLibraries.workingcontext import WorkingContext
from OrsPlugins.orsimagesaver import OrsImageSaver


class CropAndRotateAndExportChromosomes_ee6b6e04d18d11ee93ecc45ab1da15d2(ContextualMenuItem):

    @classmethod
    def getIsSelectionValid(cls, aCollectionOfObjects, implementation):
        # We only want this menu item to run on a single MultiROI. 
	#First, check to see if the user has only selected one item.
        if len(aCollectionOfObjects) != 1:
            return False

        # Check to see if the selected item is a MultiROI
        selectedItem = aCollectionOfObjects[0]

        return isinstance(selectedItem, MultiROI)

    @classmethod
    def getMenuItemForSelection(cls, aCollectionOfObjects, implementation):
        collectionString = ListHelper.asPythonCollectionString(aCollectionOfObjects)
        myMenu = Menu(title='Crop, Rotate and Export Chromosomes',
                      id_='CropAndRotateAndExportChromosomes_ee6b6e04d18d11ee93ecc45ab1da15d2',
                      section='',
                      action=f'CropAndRotateAndExportChromosomes_ee6b6e04d18d11ee93ecc45ab1da15d2.menuItemEntryPoint({collectionString}, {implementation.getVarName()})')
        return myMenu

    @classmethod
    def menuItemEntryPoint(cls, collectionString, implementation):        
        # aCollectionOfObjects is a Python list of objects currently being selected
        aCollectionOfObjects = ListHelper.fromPythonCollection(collectionString, asPythonList=True)
       
        # Prompt the user for an output folder
        outputFolder = OrsFileDialog.getExistingDirectory(WorkingContext.getCurrentContextWindow(),
                                                         caption="Select output folder")
        
        if not outputFolder:
            return
        
        # Select the chromosome MultiROI
        chromosomeMultiROI = aCollectionOfObjects[0]

        # Extract ROIs from MultiROI
        chromosomeROIList = ROIHelper.extractROISFromMultiROI(source_multiROI=chromosomeMultiROI)

        # rotate and crop each chromosome in the list
        for i, chromosomeROI in enumerate(chromosomeROIList):
            chromosomeName = f"ROI {i + 1}"
            transformedChromosome = cls.rotateAndCropChromosome(
                chromosomeROI=chromosomeROI,
                chromosomeName=chromosomeName
            )
            # Convert the ROI to binary
            voxelValue = 255

            chromosomeArray = transformedChromosome.getNDArray() * voxelValue
            filename = f"{outputFolder}/{chromosomeName}.tiff"
            skimage.io.imsave(filename, chromosomeArray)

    @classmethod
    def rotateAndCropChromosome(cls, chromosomeROI, chromosomeName=None, deleteBox=False, publishDataset=True):
        """
        Rotate and crop a chromosome from an ROI.

        :param chromosomeROI: The ROI containing the chromosome.
        :returns: the cropped and rotated ROI.
        """

        # Get the bounding box of the chromosome
        box = chromosomeROI.getBoundingBox(0)
        visualBox = VisualBox()
        visualBox.setBox(box, 0)
        instance = orsCreateABoxFromCurrentObjectBoxOrBoundingBox_df18ee80548711e8b181107b441a0862()
        instance.createAndPublishBox(visualBox)

        # Get the minimal bounding box
        iTIndex = 0
        output = TakeShapeOfMinimalBox_264ff9b8a66511e88120107b441a0862.takeMinimalBoxShapeOfObject(
            visualBox=visualBox, object=chromosomeROI, iTIndex=iTIndex
        )

        # Crop the ROI to the minimal bounding box
        timestep = 0

        title = ''

        croppedChromosome = extractBoxFromChannel_bdfbc9a8110b11e8819130e37aed33ab.createROIFromVisualBox(
            sourceDataset=chromosomeROI, visualBox=visualBox, timestep=timestep, title=title
        )

        # Perform the data rotation / inversion
        invertX = False
        invertY = False
        invertZ = False
        invertData = False
        axisTransformation = 4 #0 -> XYZ (no transformation) 1 -> XZY 2 -> YXZ 3 -> YZX 4 -> ZXY 5 -> ZYX
        createNewDataset = True

        transformedChromosome = OrsDatasetInvertor.invert(
            structuredGrid=croppedChromosome,
            invertX=invertX,
            invertY=invertY,
            invertZ=invertZ,
            invertData=invertData,
            axisTransformation=axisTransformation,
            createNewDataset=createNewDataset
        )

        if publishDataset:
            transformedChromosome.publish()

        if chromosomeName is not None:
            transformedChromosome.setTitle(chromosomeName)

        if deleteBox:
            del visualBox

        return transformedChromosome