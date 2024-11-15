"""
Crop and rotate chromosomes and export the individual chromosomes as 3D Tiff files.

Take a MultiROI containing cleaned chromosomes, crop them, rotate them to be aligned with the Z axis and export as 3D TIFF image.

:author: Lucas Philipp, Benjamin Z. Rudski
:contact: 
:email: lucas.philipp@mail.mcgill.ca; benjamin.rudski@mail.mcgill.ca
:organization: McGill University
:address: 
:copyright: 
:date: Feb 22 2024 09:23
:dragonflyVersion: 2022.2.0.1409
:UUID: ee6b6e04d18d11ee93ecc45ab1da15d2
"""

__version__ = '1.0.0'

import os.path

from PyQt5.QtWidgets import QFileDialog

from ORSServiceClass.menuItems.contextualMenuItem import ContextualMenuItem
from ORSServiceClass.actionAndMenu.menu import Menu
from ORSServiceClass.decorators.infrastructure import interfaceMethod
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
        """
        :param aCollectionOfObjects: a list of objects currently being selected, i.e. on which the menu item could be applied.
        :param implementation: a subclass of AbstractPlugin in the current context.
        :return: if True, the menu item will be displayed.
        """
        
        """
        # Example: the selection should be made of exactly 2 Channels
        from ORSModel.ors import Channel
        if aCollectionOfObjects is None or len(aCollectionOfObjects) != 2:
            return False
        selectionIsOnlyChannels = all([isinstance(obj, Channel) for obj in aCollectionOfObjects])
        return selectionIsOnlyChannels
        """
        
        # Put your code here

        # We only want this menu item to run on a single MultiROI. First, check to see if the
        # user has only selected one item.
        if len(aCollectionOfObjects) != 1:
            return False

        # Check to see if the selected item is a MultiROI
        selectedItem = aCollectionOfObjects[0]

        return isinstance(selectedItem, MultiROI)

    @classmethod
    def getMenuItemForSelection(cls, aCollectionOfObjects, implementation):
        """
        Returns the menu item
        :param aCollectionOfObjects: a list of objects currently being selected, i.e. on which the menu item will be applied.
        :param implementation: a subclass of AbstractPlugin in the current context.
        :return: Menu
        """
        
        collectionString = ListHelper.asPythonCollectionString(aCollectionOfObjects)
        myMenu = Menu(title='Crop, Rotate and Export Chromosomes',
                      id_='CropAndRotateAndExportChromosomes_ee6b6e04d18d11ee93ecc45ab1da15d2',
                      section='',
                      action=f'CropAndRotateAndExportChromosomes_ee6b6e04d18d11ee93ecc45ab1da15d2.menuItemEntryPoint({collectionString}, {implementation.getVarName()})')
        return myMenu

    @classmethod
    def menuItemEntryPoint(cls, collectionString, implementation):
        """
        Will be executed when the menu item is selected.
        :param collectionString: a list of objects representation currently being selected, i.e. on which the menu item will be applied.
        :param implementation: a subclass of AbstractPlugin in the current context.
        """
        
        # aCollectionOfObjects is a Python list of objects currently being selected
        aCollectionOfObjects = ListHelper.fromPythonCollection(collectionString, asPythonList=True)
        # Put your code here
        
        # If logging is relevant (used for macro recording and playing), a call to an interface methods is mandatory.
        # A prototype of such an interface method is written below (method "execute").
        # A full example is given in the demonstration file:
        # %ORSPYTHON%/OrsPythonPlugins/OrsGenericMenuItems/menuItems/demo_c18408e83c9511e7a502448a5b5d70c0.py
        """
        # Example: calling an interface method to copy the first dataset into the second.
        firstDataset = aCollectionOfObjects[0]
        secondDataset = aCollectionOfObjects[1]
        cls.copyInto(firstDataset, secondDataset)
        """

        # Prompt the user for an output folder

        outputFolder = QFileDialog.getExistingDirectory(WorkingContext.getCurrentContextWindow(),
                                                         caption="Select output folder",
                                                         options=QFileDialog.ShowDirsOnly)

        # Select the chromosome MultiROI
        chromosomeMultiROI = aCollectionOfObjects[0]

        # Extract ROIs from MultiROI
        chromosomeROIList = ROIHelper.extractROISFromMultiROI(source_multiROI=chromosomeMultiROI)

        # Now we have a list of ROIs! We can run the macro code on each of them
        # to produce the transformed result.
        for i, chromosomeROI in enumerate(chromosomeROIList):
            chrosomeName = f"ROI {i + 1}"
            transformedChromosome = cls.rotateAndCropChromosome(
                chromosomeROI=chromosomeROI,
                chromosomeName=chrosomeName
            )

            # Convert the ROI to binary
            voxelValue = 255
            chromosomeChannel = transformedChromosome.convertToChannel(value=voxelValue)

            # Give the channel the same title as the ROI
            chromosomeChannel.setTitle(chrosomeName)

            # Save as 3D tiff in the specified output folder
            filename = os.path.join(outputFolder, f"{chrosomeName}.tiff")
            OrsImageSaver.exportDatasetToTiffFile(dataset=chromosomeChannel.getGUID(),
                                                  fileName=filename,
                                                  useLZWCompression=False, exportMultiFrameFile=True, showProgress=True,
                                                  wlw=0, wlc=0, gamma=0)

            # Delete the channel
            # chromosomeChannel.deleteObject()
            # del chromosomeChannel



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
        axisTransformation = 4
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

    # This is the prototype of an interface method:
    # @classmethod
    # @interfaceMethod
    # def execute(cls):
    #     """
    #     Doing this and that.
    #     """
    #     
    #     # Put your code here
    #     pass

    # This is an example of an interface method:
    # @classmethod
    # @interfaceMethod
    # def copyInto(cls, sourceDataset, destinationDataset):
    #     """
    #     Copies the first dataset into the second dataset.
    #     
    #     :param sourceDataset: source dataset
    #     :type sourceDataset: ORSModel.ors.Channel
    #     :param destinationDataset: destination dataset
    #     :type destinationDataset: ORSModel.ors.Channel
    #     """
    #     
    #     sourceDataset.copyInto(destinationDataset)
    #     destinationDataset.setDataDirty()
    #     destinationDataset.propagateDataDirty()

