cmake_minimum_required(VERSION 2.6)
project(GUDHI)

include("CMakeGUDHIVersion.txt")

enable_testing()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules/")

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
add_gudhi_module(Subsampling)
add_gudhi_module(Tangential_complex)
add_gudhi_module(Witness_complex)
add_gudhi_module(Nerve_GIC)

message("++ GUDHI_MODULES list is:\"${GUDHI_MODULES}\"")

# For "make doxygen" - Requires GUDHI_USER_VERSION_DIR to be set
set(GUDHI_USER_VERSION_DIR ${CMAKE_SOURCE_DIR})
include(GUDHI_doxygen_target)

# This variable is used by Cython CMakeLists.txt to know its path
set(GUDHI_CYTHON_PATH "cython")
# For third parties libraries management - To be done last as CGAL updates CMAKE_MODULE_PATH
include(GUDHI_third_party_libraries)

if(MSVC)
  # Turn off some VC++ warnings
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4668 /wd4311 /wd4800 /wd4820 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -pedantic")
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug)
  message("++ Debug compilation flags are: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
else()
  message("++ Release compilation flags are: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

if (DEBUG_TRACES)
  message(STATUS "DEBUG_TRACES are activated")
  # For programs to be more verbose
  add_definitions(-DDEBUG_TRACES)
endif()

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
endif()
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# GUDHIConfig.cmake
# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE GUDHI)

message("++ make install will install ${PROJECT_NAME} in the following directory : ${CMAKE_INSTALL_PREFIX}")
# Create the GUDHIConfig.cmake and GUDHIConfigVersion files
set(CONF_INCLUDE_DIRS "${CMAKE_INSTALL_PREFIX}/include")
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
install(DIRECTORY include/Miniball DESTINATION include)
#---------------------------------------------------------------------------------------
