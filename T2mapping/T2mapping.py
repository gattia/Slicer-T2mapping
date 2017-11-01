import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import numpy
import string 

#
# T2mapping
#

class T2mapping(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "T2 Mapping" 
    self.parent.categories = ["Quantification"]
    self.parent.dependencies = []
    self.parent.contributors = ["Anthony Gatti)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
    This program take a MultiVolume image as input, where, the 4th dimension is echo time. 
    The program fits a monoexponential decay to the echoes using a log transformed least squares regression method. 
    The resulting fit is used to determine the T2 relaxation time for every pixel in the image. 
    This program will produce 3 resulting images, if they are selected. 1) T2 relaxation map, 2) Proton Density (PD) map,
    3) R2 map of the resulting fit.  
    """
    self.parent.acknowledgementText = """
    We would like to acknowledge Jean-Christophe Fillion-Robin, and Steve Pieper for producing the original template python module used to creat this package. We would also like to Acknowledge Michael Noseworthy for aid in image processing and analyses that lead to production of this module. 
""" # replace with organization, grant and thanks.

#
# T2mappingWidget
#

class T2mappingWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    #
    # input volume selector
    #
    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ( ("vtkMRMLMultiVolumeNode"), "" )
    self.inputSelector.selectNodeUponCreation = True
    self.inputSelector.addEnabled = False
    self.inputSelector.removeEnabled = False
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Pick the MultiVolume or Multiple Echo image" )
    parametersFormLayout.addRow("Input Volume: ", self.inputSelector)

    #
    # output volume selector
    #
    self.outputSelector = slicer.qMRMLNodeComboBox()
    self.outputSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.outputSelector.selectNodeUponCreation = True
    self.outputSelector.addEnabled = True
    self.outputSelector.removeEnabled = True
    self.outputSelector.noneEnabled = True
    self.outputSelector.showHidden = False
    self.outputSelector.showChildNodeTypes = False
    self.outputSelector.setMRMLScene( slicer.mrmlScene )
    self.outputSelector.setToolTip( "Pick the volume to output the T2 map, or create a new one" )
    parametersFormLayout.addRow("Output Volume: ", self.outputSelector)

    #
    # R2 threshold
    #
    self.r2ThresholdSliderWidget = ctk.ctkSliderWidget()
    self.r2ThresholdSliderWidget.singleStep = 0.01
    self.r2ThresholdSliderWidget.minimum = 0
    self.r2ThresholdSliderWidget.maximum = 1
    self.r2ThresholdSliderWidget.value = 0.5
    self.r2ThresholdSliderWidget.setToolTip("Set threshold value for T2 fit (R^2). Voxels that have R^2 lower than this value will set to zero.")
    parametersFormLayout.addRow("R^2 threshold", self.r2ThresholdSliderWidget)

    #
    # Upper T2 threshold
    #
    self.t2ThresholdSliderWidget = ctk.ctkSliderWidget()
    self.t2ThresholdSliderWidget.singleStep = 0.1
    self.t2ThresholdSliderWidget.minimum = 0
    self.t2ThresholdSliderWidget.maximum = 500
    self.t2ThresholdSliderWidget.value = 100
    self.t2ThresholdSliderWidget.setToolTip("Set upper threshold value for maximum T2 value. Voxels that have T2 higher than this value will set to zero.")
    parametersFormLayout.addRow("T2 Upper Threshold", self.t2ThresholdSliderWidget)
    
    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = True
    parametersFormLayout.addRow(self.applyButton)

    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    # self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    # self.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  # def onSelect(self):
  #   self.applyButton.enabled = self.inputSelector.currentNode() and self.outputSelector.currentNode()

  def onApplyButton(self):
    logic = T2mappingLogic()
    r2Cutoff = self.r2ThresholdSliderWidget.value
    t2Cutoff = self.t2ThresholdSliderWidget.value
    logic.run(self.inputSelector.currentNode(), self.outputSelector.currentNode(), r2Cutoff, t2Cutoff)
    
    inputVolume = self.inputSelector.currentNode()
    outputVolume = self.outputSelector.currentNode()

    

#
# T2mappingLogic
#

class T2mappingLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def run(self, inputVolume, outputVolume, r2Cutoff, t2Cutoff):
    """
    Run the actual algorithm
    SI = Signal Intensity
    T2 = Transverse Relaxation Time (ms) 
    TE = Echo Time(ms)
    """
    if not (inputVolume and outputVolume):
      qt.QMessageBox.critical(slicer.util.mainWindow(), 'T2 Map', 'Input & output volumes needed!')
      return

    logging.info('Processing started')
    
    t2Node = slicer.util.getNode(inputVolume.GetID())
    t2Data = slicer.util.array(inputVolume.GetID())
    
    TE = numpy.array(map(float,string.split(t2Node.GetAttribute('MultiVolume.FrameLabels'),','))) #list of echoes
    fourDSize = t2Data.shape #Shape of the multidimensional array
    
    #log transform the T2 data to be able to do a linear least squares estimate. 
    t2Log = numpy.log(t2Data)
    noPixels = fourDSize[0] * fourDSize[1] * fourDSize[2] #number of pixels in the entire 3D MRI
    
    #get number of echoes and create the array needed to do the least squares estimate. 
    noEchoes = len(TE)
    regressTE = numpy.matrix([numpy.ones(noEchoes), TE]).T
    
    '''
    Calcualte beta coefficients from: beta =(X.T * X).I * X.T * Y
    As can be seen from the above equation, most of the math can be done from just the X(TE) data before multiplying the Y(SI) as the last step. 
    Because TEs are constant across pixels, this math can be done once and then be multiplied by each respective set of Y(SIs).
    '''
    
    regressTE = numpy.matrix(regressTE) # These are to prepare the Echoes... or the X variable from the above equation.
    regressPartOne = (regressTE.T*regressTE).I *regressTE.T #this is that part of the equation that can be done without SI (Y-variable).
    
    '''
    We now prepare the dependent variable. Turn the 3D volume for each TE into a single-dimension and concatenate. So we end up with a matrix that is 2 dimensions. 
    The first dimension is the number of echoes and the second dimension is the number of pixels in the entire MRI. 
    '''
    dependentVariable = numpy.empty([noEchoes, noPixels])
    for echo in range(noEchoes):
      dependentVariable[echo,:] = t2Log[:,:,:,echo].flatten()
    dependentVariable = numpy.matrix(dependentVariable)
    
    regressionResult = numpy.dot(regressPartOne, dependentVariable)# produces a 2xn matrix with row 1 = intercept, and row 2 = T2 fit or 1/T2. 
    t2Values = numpy.asarray(-1/regressionResult[1,:]) #these are the T2 values - still in a single dimension
    pdValues = numpy.asarray(regressionResult[0,:]) # these are the PD values - still in single dimension
    
    nanValues = numpy.where( ((numpy.isnan(t2Values)) | (numpy.isinf(t2Values)) )) #determines where we have NANs / Infinity
    t2Values[nanValues] = 0 #set nans = 0 
    outlierValues = numpy.where( ((t2Values>t2Cutoff) | (t2Values<0)) )  #Identified outliers based on the t2 cutoff. Take all T2 values that are between our upper cutoff determined in the GUI and 0. 0 is the lower-bound because this is the lowest value that makes logical sense, it shouldn't / cant be <0. 
    t2Values[outlierValues] = 0 #set outliers = 0
    
    
    '''
    Now we calculate the R2 for each pixel. We use the general notation that R2 = 1 - ssRes/ssTotal
    Where ssTotal = sum(yi-ymean)^2;
    ssRes = sum(ymeasured - ypredicted)^2
    '''
    determineR2For = numpy.where((t2Values[0,:] >0 )) #only calcualte R2 for non-zero numbers (outliers etc. removed)
    dvDoMath = numpy.squeeze(dependentVariable[:,determineR2For]) #only the Y (dependent) values that are not "outliers"
    #resultR2 = numpy.zeros(shape =(dvDoMath.shape[1]))
    r2WholeImage = numpy.zeros(shape=(t2Values.shape[1])) #pre-allocate Array that we will fill with R2 values later
    regressionCoefficientsCandidates = numpy.squeeze(regressionResult[:,determineR2For])
    
    ssTotalMean = numpy.mean(dvDoMath, axis=0) #calculate mean Log transfortmed SI for each pixel. 
    ssTotalDiff = dvDoMath-ssTotalMean #What is the difference between the SI of each pixel and the mean SI for that pixel
    ssTotalDiffSquared = numpy.square(ssTotalDiff) #Square the difference between the mean value and the actual SI
    ssTotal = numpy.sum(ssTotalDiffSquared, axis=0) #sum of squares total = sum of the squared differences. 
    
    intercept = regressionCoefficientsCandidates[0,:] #Get intercepts for only the pixels of interest. 
    slope = regressionCoefficientsCandidates[1,:] #Get the slopes for only the pixels of interest
    echoTime = numpy.zeros(shape=(dvDoMath.shape))
    for index in range(noEchoes):
      echoTime[index,:] = TE[index] #creates an array that essentially tells the TE (independent variable) for every collected dependent variable (Y or SI). 
    slopeTimesEchoTime = numpy.multiply(echoTime, slope) 
    ssResPredicted = numpy.add(intercept, slopeTimesEchoTime)# Add intercept to the above portion of the equation and determine the predicted SI (Y) for every pixel. 
    ssResDiff = numpy.subtract(dvDoMath, ssResPredicted) #Difference between the predicted values and the actual values
    ssResDiffSquared = numpy.square(ssResDiff) #Squared Difference
    ssRes = numpy.sum(ssResDiffSquared, axis=0) #Sum of the squared differeences, or the sum of the squared residuals. 
    
    unexplainedVariance = numpy.divide(ssRes, ssTotal) # Unexplained variance is determined by ssRes/ssTotal. 
    r2 = numpy.subtract(1, unexplainedVariance) #R2 is just 1- the unexplained variance. 
    
    #Now we place the calcualted r2 values back in those pixels where we wanted to determine it... i.e. those with an "acceptable" fit as chosen in the GUI
    r2WholeImage[determineR2For] = r2

    poorfit = numpy.where(r2WholeImage<r2Cutoff)
    t2Values = numpy.squeeze(t2Values)
    t2Values[poorfit] = 0 #setthe locations that had poor fit to be equal to 0.


    '''
    Create 3D images of T2/PD data
    '''
    t2Map = numpy.reshape(t2Values, fourDSize[0:3])
    pdMap = numpy.reshape(pdValues, fourDSize[0:3])
    
    imageSize = (fourDSize[2], fourDSize[1], fourDSize[0])
    imageSpacing = [1,1,1]
    voxelType = vtk.VTK_UNSIGNED_CHAR
    outputVolume.SetSpacing(imageSpacing)
#
    
    # ##V3 (http://massmail.spl.harvard.edu/public-archives/slicer-devel/2015/037087.html)
    # importer = vtk.vtkImageImport()
    # importer.CopyImportVoidPointer(t2Map, t2Map.nbytes)
    # setDataType = 'importer.SetDataScalarTypeTo' + 'UnsignedShort' + '()'
    # eval(setDataType)
    # importer.SetNumberOfScalarComponents(1)
    # importer.SetWholeExtent(0, t2Map.shape[2]-1, 0, t2Map.shape[1]-1, 0, t2Map.shape[0]-1)
    # importer.SetDataExtentToWholeExtent()
    # importer.Update()
    #
    # ijkToRAS = vtk.vtkMatrix4x4()
    # t2Node.GetIJKToRASMatrix(ijkToRAS)
    # outputVolume.SetIJKToRASMatrix(ijkToRAS)
    #
    # outputVolume.SetAndObserveImageData(importer.GetOutput())
    # slicer.mrmlScene.AddNode(outputVolume)
    # volumeDisplayNode = slicer.vtkMRMLScalarVolumeDisplayNode()
    # slicer.mrmlScene.AddNode(volumeDisplayNode)
    # volumeDisplayNode.SetAndObserveColorNodeID('vtkMRMLColorTableNodeFileHotToColdRainbow.txt')
    #
    # outputVolume.SetAndObserveDisplayNodeID(volumeDisplayNode.GetID())
    # selectionNode = slicer.app.applicationLogic().GetSelectionNode()
    # selectionNode.SetReferenceActiveVolumeID(outputVolume.GetID())
    # slicer.app.applicationLogic().PropagateVolumeSelection(0)
    #
     
  
  

    # # V2 (https://github.com/stevedaxiao/T1_Mapping/blob/master/T1_Mapping/T1_Mapping.py)
    # importer = vtk.vtkImageImport()
    # data_string = t2Map.tostring()
    # importer.CopyImportVoidPointer(data_string, len(data_string))
    # setDataType = 'importer.SetDataScalarTypeTo' + 'UnsignedShort' + '()'
    # eval(setDataType)
    # importer.SetNumberOfScalarComponents(1)
    # importer.SetWholeExtent(0, t2Map.shape[2]-1, 0, t2Map.shape[1]-1, 0, t2Map.shape[0]-1)
    # importer.SetDataExtentToWholeExtent()
    # print importer.GetDataExtent()
    # importer.Update()
    #
    # ijkToRAS = vtk.vtkMatrix4x4()
    # t2Node.GetIJKToRASMatrix(ijkToRAS)
    # outputVolume.SetIJKToRASMatrix(ijkToRAS)
    # outputVolume.SetAndObserveImageData(importer.GetOutput())
    # slicer.mrmlScene.AddNode(outputVolume)
    # volumeDisplayNode = slicer.vtkMRMLScalarVolumeDisplayNode()
    # slicer.mrmlScene.AddNode(volumeDisplayNode)
    # volumeDisplayNode.SetAndObserveColorNodeID('vtkMRMLColorTableNodeFileHotToColdRainbow.txt')
    #
    # selectionNode = slicer.app.applicationLogic().GetSelectionNode()
    # selectionNode.SetReferenceActiveVolumeID(outputVolume.GetID())
    # slicer.app.applicationLogic().PropagateVolumeSelection(0)
    #
    
    
    
    #V1
    #create empty image volume
    imageData = vtk.vtkImageData()
    imageData.SetDimensions(imageSize)
    imageData.AllocateScalars(voxelType,1)
    thresholder = vtk.vtkImageThreshold()
    thresholder.SetInputData(imageData)
    thresholder.SetInValue(0)
    thresholder.SetOutValue(0)

    #set volume information
    outputVolume.SetSpacing(imageSpacing)
    outputVolume.SetImageDataConnection(thresholder.GetOutputPort())
    #Add volume to scene
    slicer.mrmlScene.AddNode(outputVolume)
    displayNode = slicer.vtkMRMLScalarVolumeDisplayNode()
    slicer.mrmlScene.AddNode(displayNode)
    outputVolume.SetAndObserveDisplayNodeID(displayNode.GetID())
    outputVolume.CreateDefaultStorageNode()
    outputVolume.SetName('T2 Map')

    volumeArray = slicer.util.array(outputVolume.GetID())
    #volumeArray = slicer.util.array('T2 Map')
    volumeArray[:] = t2Map

    ijkToRAS = vtk.vtkMatrix4x4()
    t2Node.GetIJKToRASMatrix(ijkToRAS)
    outputVolume.SetIJKToRASMatrix(ijkToRAS)

    selectionNode = slicer.app.applicationLogic().GetSelectionNode()
    selectionNode.SetReferenceActiveVolumeID(outputVolume.GetID())
    slicer.app.applicationLogic().PropagateVolumeSelection(0)

    changeColorMapNode = slicer.util.getNode(outputVolume.GetID())
    displayNode = changeColorMapNode.GetDisplayNode()
    displayNode.AutoWindowLevelOff()
    displayNode.SetLevel(t2Cutoff/2)
    displayNode.SetWindow(t2Cutoff)
    # displayNode.SetAndObserveColorNodeID('vtkMRMLColorTableNodeFileHotToColdRainbow.txt')
    displayNode.SetAndObserveColorNodeID('vtkMRMLColorTableNodeFileViridis.txt')
    outputVolume.CreateDefaultStorageNode()
  #
  #
    logging.info('Processing completed')

    return True

