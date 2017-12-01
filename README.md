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
Included in the screenshots folder are examples of the program (https://github.com/gattia/Slicer-T2mapping/tree/master/Screenshots). 

#### Image 1 shows an example of the module when first opened. 
![alt text][image_1] {:height="50%" width="50%"}
#### Image 2 shows creation of a new volume where the T2map will be saved. 
![alt text][image_2] {:height="50%" width="50%"} 
#### Image 3 shows an example of the resulting segmentation created using provided example data (https://github.com/gattia/Slicer-T2mapping/tree/master/Example_multi_echo_mri). 
![alt text][image_3] {:height="50%" width="50%"} 
#### Image 4 shows an example of a resulting R^2 map. 
![alt text][image_4] {:height="50%" width="50%"} 

### Example Data. 
There is example data included at https://github.com/gattia/Slicer-T2mapping/tree/master/Example_multi_echo_mri. This data is known to work with the module and successfully creates T2/R^2/PD maps. This data was acquired using the GE Cartigram sequence, a multi-echo-spin-echo sequence with 8 echoes acquired at ~6ms, 12ms, 18ms, 24ms, 30ms, 36ms, 42ms, 48ms. 



[image_1]: https://github.com/gattia/Slicer-T2mapping/blob/master/Screenshots/1_Module_and_multi_echo_t2_image.png "Image 1"
[image_2]: https://github.com/gattia/Slicer-T2mapping/blob/master/Screenshots/2_Module_create_new_volume_t2.png "Image 2"
[image_3]: https://github.com/gattia/Slicer-T2mapping/blob/master/Screenshots/3_Example_resulting_t2_map_r2_threshold_0.7_t2_upper_threshold_100.png "Image 3"
[image_4]: https://github.com/gattia/Slicer-T2mapping/blob/master/Screenshots/4_Example_resulting_R2_map_r2_threshold_0.7_t2_upper_threshold_100.png "Image 4"
