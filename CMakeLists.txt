cmake_minimum_required(VERSION 2.8.9)

project(T2mapping)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://github.com/gattia/Slicer-T2mapping")
set(EXTENSION_CATEGORY "Quantification")
set(EXTENSION_CONTRIBUTORS "Anthony Gatti")
set(EXTENSION_DESCRIPTION "Used to create T2 maps from multi echo images.")
set(EXTENSION_ICONURL "https://raw.githubusercontent.com/gattia/Slicer-T2mapping/master/T2mapping.png")
set(EXTENSION_SCREENSHOTURLS "https://raw.githubusercontent.com/gattia/Slicer-T2mapping/master/Screenshots/1_Module_and_multi_echo_t2_image.png https://raw.githubusercontent.com/gattia/Slicer-T2mapping/master/Screenshots/2_Module_create_new_volume_t2.png https://raw.githubusercontent.com/gattia/Slicer-T2mapping/master/Screenshots/3_Example_resulting_t2_map_r2_threshold_0.7_t2_upper_threshold_100.png https://raw.githubusercontent.com/gattia/Slicer-T2mapping/master/Screenshots/4_Example_resulting_R2_map_r2_threshold_0.7_t2_upper_threshold_100.png")

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(T2mapping)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
