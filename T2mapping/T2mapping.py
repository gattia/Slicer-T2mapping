import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import numpy

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
    self.parent.contributors = ["Anthony Gatti"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
    This program take a MultiVolume image as input, where, the 4th dimension is echo time. 
    The program fits a monoexponential decay to the echoes using a log transformed least squares regression method. 
    The resulting fit is used to determine the T2 relaxation time for every pixel in the image. 
    This program will produce 3 resulting images, if they are selected. 1) T2 relaxation map, 2) Proton Density (PD) map,
    3) R2 map of the resulting fit.  
    """
    self.parent.acknowledgementText = """
    We would like to acknowledge Jean-Christophe Fillion-Robin, and Steve Pieper for producing the original template python module used to creat this package. We would also like to Acknowledge Michael Noseworthy for aid in image processing and analyses that lead to production of this module. 
""" 

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
    self.T2outputSelector = slicer.qMRMLNodeComboBox()
    self.T2outputSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.T2outputSelector.selectNodeUponCreation = True
    self.T2outputSelector.addEnabled = True
    self.T2outputSelector.removeEnabled = True
    self.T2outputSelector.noneEnabled = True
    self.T2outputSelector.showHidden = False
    self.T2outputSelector.showChildNodeTypes = False
    self.T2outputSelector.setMRMLScene( slicer.mrmlScene )
    self.T2outputSelector.setToolTip( "Pick the volume to output the T2 map, or create a new one" )
    parametersFormLayout.addRow("T2 Map Output Volume: ", self.T2outputSelector)
    #
    # output volume selector - IF THEY WANT A PD MAP
    #
    self.PDoutputSelector = slicer.qMRMLNodeComboBox()
    self.PDoutputSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.PDoutputSelector.selectNodeUponCreation = True
    self.PDoutputSelector.addEnabled = True
    self.PDoutputSelector.removeEnabled = True
    self.PDoutputSelector.noneEnabled = True
    self.PDoutputSelector.showHidden = False
    self.PDoutputSelector.showChildNodeTypes = False
    self.PDoutputSelector.setMRMLScene( slicer.mrmlScene )
    self.PDoutputSelector.setToolTip( "If you want a PD Map, select a volume or create a new one new one to write the Map to" )
    parametersFormLayout.addRow("Proton Density Output Volume: ", self.PDoutputSelector)
    #
    # output volume selector - IF THEY WANT AN R2 MAP
    #
    self.R2outputSelector = slicer.qMRMLNodeComboBox()
    self.R2outputSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.R2outputSelector.selectNodeUponCreation = True
    self.R2outputSelector.addEnabled = True
    self.R2outputSelector.removeEnabled = True
    self.R2outputSelector.noneEnabled = True
    self.R2outputSelector.showHidden = False
    self.R2outputSelector.showChildNodeTypes = False
    self.R2outputSelector.setMRMLScene( slicer.mrmlScene )
    self.R2outputSelector.setToolTip( "If you want an R^2 Map for the fitted T2 values select a volume or create a new one new one to write the Map to" )
    parametersFormLayout.addRow("R^2 Output Volume: ", self.R2outputSelector)

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
    self.t2ThresholdSliderWidget.maximum = 1000
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

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    # self.onSelect()

  def cleanup(self):
    pass

  def onApplyButton(self):
    logic = T2mappingLogic()
    r2Cutoff = self.r2ThresholdSliderWidget.value
    t2Cutoff = self.t2ThresholdSliderWidget.value
    logic.run(self.inputSelector.currentNode(), self.T2outputSelector.currentNode(), self.PDoutputSelector.currentNode(), self.R2outputSelector.currentNode(), r2Cutoff, t2Cutoff)
    
    inputVolume = self.inputSelector.currentNode()
    T2outputVolume = self.T2outputSelector.currentNode()
#
# T2mappingLogic
#

class T2mappingLogic(ScriptedLoadableModuleLogic):
  def createVolume(self, data, outputVolume, inputVolume, inputNode, imageSize, voxelType, name):
        dataMax = numpy.nanmax(data[numpy.isfinite(data)])
        #create empty image volume
        imageData = vtk.vtkImageData()
        imageData.SetDimensions(imageSize)
        imageData.AllocateScalars(voxelType,1)
        thresholder = vtk.vtkImageThreshold()
        thresholder.SetInputData(imageData)
        thresholder.SetInValue(0)
        thresholder.SetOutValue(0)

        #set volume information
        outputVolume.SetSpacing(inputVolume.GetSpacing())

        outputVolume.SetImageDataConnection(thresholder.GetOutputPort())
        #Add volume to scene
        slicer.mrmlScene.AddNode(outputVolume)
        displayNode = slicer.vtkMRMLScalarVolumeDisplayNode()
        slicer.mrmlScene.AddNode(displayNode)


        outputVolume.SetAndObserveDisplayNodeID(displayNode.GetID())
        outputVolume.CreateDefaultStorageNode()
        outputVolume.SetName(name)

        volumeArray = slicer.util.array(outputVolume.GetID())
        volumeArray[:] = data[:]

        ijkToRAS = vtk.vtkMatrix4x4()
        inputNode.GetIJKToRASMatrix(ijkToRAS)
        outputVolume.SetIJKToRASMatrix(ijkToRAS)

        displayNode.AutoWindowLevelOff()
        displayNode.SetLevel(dataMax/2)
        displayNode.SetWindow(dataMax)
        displayNode.SetAndObserveColorNodeID('vtkMRMLColorTableNodeFileViridis.txt')
        outputVolume.CreateDefaultStorageNode()
        return()

  def run(self, inputVolume, T2outputVolume, PDoutputVolume, R2outputVolume, r2Cutoff, t2Cutoff):
    """
    Run the actual algorithm
    SI = Signal Intensity
    T2 = Transverse Relaxation Time (ms) 
    TE = Echo Time(ms)
    """
    if not (inputVolume and T2outputVolume):
      qt.QMessageBox.critical(slicer.util.mainWindow(), 'T2 Map', 'Input & output volumes needed!')
      return

    logging.info('Processing started')
    
    t2Node = slicer.util.getNode(inputVolume.GetID())
    t2Data = slicer.util.array(inputVolume.GetID())
    
    TE = t2Node.GetAttribute('MultiVolume.FrameLabels').split(',') # list of echoes as strings
    TE = numpy.array(list(map(float, TE))) # convert strings to float and wrap in np.array
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
    pdValues = numpy.asarray(regressionResult[0,:]) # these are the PD values - still in single dimension - These are in log transformed Units
    pdValues = numpy.exp(pdValues)
    
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
    r2Map = numpy.reshape(r2WholeImage, fourDSize[0:3])
    
    imageSize = (fourDSize[2], fourDSize[1], fourDSize[0])
    imageSpacing = [1,1,1]
    voxelType = vtk.VTK_FLOAT
    # T2outputVolume.SetSpacing(imageSpacing)

    t2Name = 'T2 Map'
    pdName = 'PD Map'
    r2Name = 'R^2 Map'

    #These functions put the different images into the "containers" needed to show the data in Slicer. 
    self.createVolume(t2Map, T2outputVolume, inputVolume, t2Node, imageSize, voxelType, t2Name)
    if PDoutputVolume:
        self.createVolume(pdMap, PDoutputVolume, inputVolume, t2Node, imageSize, voxelType, pdName)
    if R2outputVolume:
        self.createVolume(r2Map, R2outputVolume, inputVolume, t2Node, imageSize, voxelType, r2Name)
        
    logging.info('Processing completed')

    # I called the result node ResultVolume in this code
    yellow_logic = slicer.app.layoutManager().sliceWidget("Yellow").sliceLogic()
    yellow_cn = yellow_logic.GetSliceCompositeNode()
    yellow_cn.SetBackgroundVolumeID(T2outputVolume.GetID())
    lm = slicer.app.layoutManager()
    lm.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUpYellowSliceView)

    # Fit slice to window
    sliceNodes = slicer.util.getNodes('vtkMRMLSliceNode*')
    layoutManager = slicer.app.layoutManager()
    for sliceNode in sliceNodes.values():
        sliceWidget = layoutManager.sliceWidget(sliceNode.GetLayoutName())
        if sliceWidget:
            sliceWidget.sliceLogic().FitSliceToAll()

    return True

