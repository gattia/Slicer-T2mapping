cmake_minimum_required(VERSION 2.8.9)

project(T2mapping)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/T2mapping")
set(EXTENSION_CATEGORY "Quantification")
set(EXTENSION_CONTRIBUTORS "Anthony Gatti")
set(EXTENSION_DESCRIPTION "Used to create T2 maps from multi echo images. ")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/T2mapping.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/T2mapping/Screenshots/1.png")

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