cmake_minimum_required(VERSION 3.1)

project(GUDHI)

include(CMakeGUDHIVersion.txt)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules/")

set(GUDHI_MODULES "" CACHE INTERNAL "GUDHI_MODULES")
set(GUDHI_MISSING_MODULES "" CACHE INTERNAL "GUDHI_MISSING_MODULES")

# This variable is used by Cython CMakeLists.txt and by GUDHI_third_party_libraries to know its path
set(GUDHI_CYTHON_PATH "cython")

# For third parties libraries management - To be done last as CGAL updates CMAKE_MODULE_PATH
include(GUDHI_third_party_libraries NO_POLICY_SCOPE)

include(GUDHI_compilation_flags)

# Add your new module in the list, order is not important
include(GUDHI_modules)

add_gudhi_module(common)
add_gudhi_module(Alpha_complex)
add_gudhi_module(Bitmap_cubical_complex)
add_gudhi_module(Bottleneck_distance)
add_gudhi_module(Cech_complex)
add_gudhi_module(Contraction)
add_gudhi_module(Hasse_complex)
add_gudhi_module(Persistence_representations)
add_gudhi_module(Persistent_cohomology)
add_gudhi_module(Rips_complex)
add_gudhi_module(Simplex_tree)
add_gudhi_module(Skeleton_blocker)
add_gudhi_module(Spatial_searching)
add_gudhi_module(Strong_collapse)
add_gudhi_module(Subsampling)
add_gudhi_module(Tangential_complex)
add_gudhi_module(Witness_complex)
add_gudhi_module(Nerve_GIC)

# For "make doxygen" - Requires GUDHI_USER_VERSION_DIR to be set
set(GUDHI_USER_VERSION_DIR ${CMAKE_SOURCE_DIR})
include(GUDHI_doxygen_target)

#---------------------------------------------------------------------------------------
# Gudhi compilation part
include_directories(include)

# Include module CMake subdirectories
# GUDHI_SUB_DIRECTORIES is managed in CMAKE_MODULE_PATH/GUDHI_modules.cmake
foreach(GUDHI_MODULE ${GUDHI_MODULES})
  foreach(GUDHI_SUB_DIRECTORY ${GUDHI_SUB_DIRECTORIES})
    if(EXISTS ${CMAKE_SOURCE_DIR}/${GUDHI_SUB_DIRECTORY}/${GUDHI_MODULE}/CMakeLists.txt)
      add_subdirectory(${CMAKE_SOURCE_DIR}/${GUDHI_SUB_DIRECTORY}/${GUDHI_MODULE}/)
    endif()
  endforeach()
endforeach()

add_subdirectory(GudhUI)

if (WITH_GUDHI_PYTHON)
  # specific for cython module
  add_subdirectory(${GUDHI_CYTHON_PATH})
else()
  message("++ Python module will not be compiled because WITH_GUDHI_PYTHON is set to OFF")
  set(GUDHI_MISSING_MODULES ${GUDHI_MISSING_MODULES} "python")
endif()

message("++ GUDHI_MODULES list is:\"${GUDHI_MODULES}\"")
message("++ GUDHI_MISSING_MODULES list is:\"${GUDHI_MISSING_MODULES}\"")

#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# GUDHIConfig.cmake
# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE GUDHI)

message("++ make install will install ${PROJECT_NAME} in the following directory : ${CMAKE_INSTALL_PREFIX}")
# Create the GUDHIConfig.cmake and GUDHIConfigVersion files
set(CONF_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/include;${CMAKE_INSTALL_PREFIX}/include")
configure_file(GUDHIConfig.cmake.in "${PROJECT_BINARY_DIR}/GUDHIConfig.cmake" @ONLY)
configure_file(GUDHIConfigVersion.cmake.in "${PROJECT_BINARY_DIR}/GUDHIConfigVersion.cmake" @ONLY)

#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# Gudhi installation part

# Install the GUDHIConfig.cmake and GUDHIConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}/GUDHIConfig.cmake"
  "${PROJECT_BINARY_DIR}/GUDHIConfigVersion.cmake"
  DESTINATION share/gudhi)

# install the include file on "make install"
install(DIRECTORY include/gudhi DESTINATION include)
#---------------------------------------------------------------------------------------
