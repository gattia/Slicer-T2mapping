## T2 Mapping Slicer Extension 

### Inputs/Outputs
There are 4 selectors within the module. 
The first and second selectors are required. They are the Input 4D mri with multiple echos and selection of an output T2 map volume. Users are recommended to select the create new volume from the drop down. If users also want to save a coinciding Proton density (PD) map and/or a coinciding R^2 they are advised to create a new volume for each of these input parameters. 

### Parameter Selectors
Next, users are advised to determine a minimum R^2 value and an upper T2 threshold value. The R^2 threshold is used to remove pixels that have a poor fit to the T2 data. The T2 upper threshold is used to remove pixels that are above a specified value. This was created to allow users to remove pixels with T2 values that are not characteristic of their tissues of interest. 

On Apply, if the user has not selected an appropriate Input volume and a T2 Map output volume they will get an error. If no PD output or R2 outpur are selected, these maps will not be created. 

### Background
T2 is calculated using a linear least squares method. T2 is a mono-exponential decaying signal. Therefore, we log transform the pixel data so that it becomes linear. Using the log transformed data, we fit the data using linear least squares. From the fitted data, we get the slope and intercept. The slope (m) of the fitted equation is equivalent to 1/T2. T2 in ms is extracted as T2 = 1/m. The intercept of the fit is equivalent to the proton density. PD is calcualted every time, and PD maps are created if the user selects a PD output volume from the module. 

There are supporting documentation notes within the module repository that explain some of the linear algebra used. https://github.com/gattia/Slicer-T2mapping/tree/master/SupportingDocuments

### Examples/Screenshots
Included in the screenshots folder are examples of the program. Image 1_ shows an example of the module when first opened. Image 2_ shows creation of a new volume where the T2map will be saved. Image 3_ shows an example of the resulting segmentation created using provided example data (). Image 4_ shows an example of a resulting R^2 map. 

### Example Data. 
There is example data included at . This data is known to work with the module and successfully create T2/R^2/PD maps. This data was acquired using the GE Cartigram sequence, a multi-echo-spin-echo sequence with 8 echoes acquired at ~6ms, 12ms, 18ms, 24ms, 30ms, 36ms, 42ms, 48ms. 