cmake_minimum_required(VERSION 2.6)
project(GUDHI)

include("CMakeGUDHIVersion.txt")
# Generate GUDHI official version file
configure_file(GUDHIVersion.cmake.in "${CMAKE_SOURCE_DIR}/GUDHIVersion.cmake" @ONLY)

enable_testing()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules/")

# For "make doxygen"
set(GUDHI_USER_VERSION_DIR ${CMAKE_SOURCE_DIR})
include(${CMAKE_MODULE_PATH}/GUDHI_doxygen_target.txt)

# This variable is used by Cython CMakeLists.txt to know its path
set(GUDHI_CYTHON_PATH "cython")
# For third parties libraries management - To be done last as CGAL updates CMAKE_MODULE_PATH
include("${CMAKE_MODULE_PATH}/GUDHI_third_party_libraries.txt")

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

add_subdirectory(example/common)
add_subdirectory(example/Simplex_tree)
add_subdirectory(example/Persistent_cohomology)
add_subdirectory(example/Skeleton_blocker)
add_subdirectory(example/Contraction)
add_subdirectory(example/Bitmap_cubical_complex)
add_subdirectory(example/Witness_complex)
add_subdirectory(example/Alpha_complex)
add_subdirectory(example/Rips_complex)
add_subdirectory(example/Spatial_searching)
add_subdirectory(example/Subsampling)
add_subdirectory(example/Tangential_complex)
add_subdirectory(example/Bottleneck_distance)
add_subdirectory(example/Nerve_GIC)

# data points generator
add_subdirectory(data/points/generator)

add_subdirectory(GudhUI)

# specific for cython module
add_subdirectory(${GUDHI_CYTHON_PATH})
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
#---------------------------------------------------------------------------------------
