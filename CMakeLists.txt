cmake_minimum_required(VERSION 2.6)
project(GUDHIdev)

include(CMakeGUDHIVersion.txt)

set(CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/src/cmake/modules/")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/src/cmake/modules/")
message("CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}")
message("CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")

enable_testing()

# For "make user_version"
include("${CMAKE_MODULE_PATH}/GUDHI_user_version_target.txt")
# For "make doxygen"
include(GUDHI_doxygen_target)
# This variable is used by Cython CMakeLists.txt to know its path
set(GUDHI_CYTHON_PATH "src/cython")
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
  # For programs to be more verbose
  message(STATUS "DEBUG_TRACES are activated")
  add_definitions(-DDEBUG_TRACES)
endif()

include_directories(src/common/include/)
include_directories(src/Alpha_complex/include/)
include_directories(src/Bitmap_cubical_complex/include/)
include_directories(src/Bottleneck_distance/include/)
include_directories(src/Contraction/include/)
include_directories(src/Hasse_complex/include/)
include_directories(src/Persistent_cohomology/include/)
include_directories(src/Rips_complex/include/)
include_directories(src/Simplex_tree/include/)
include_directories(src/Skeleton_blocker/include/)
include_directories(src/Spatial_searching/include/)
include_directories(src/Subsampling/include/)
include_directories(src/Tangential_complex/include/)
include_directories(src/Witness_complex/include/)
include_directories(src/Nerve_GIC/include/)

add_subdirectory(src/common/example)
add_subdirectory(src/common/test)
add_subdirectory(src/Simplex_tree/test)
add_subdirectory(src/Simplex_tree/example)
add_subdirectory(src/Persistent_cohomology/test)
add_subdirectory(src/Persistent_cohomology/example)
add_subdirectory(src/Persistent_cohomology/benchmark)
add_subdirectory(src/Skeleton_blocker/test)
add_subdirectory(src/Skeleton_blocker/example)
add_subdirectory(src/Contraction/example)
add_subdirectory(src/Witness_complex/test)
add_subdirectory(src/Witness_complex/example)
add_subdirectory(src/Bitmap_cubical_complex/test)
add_subdirectory(src/Bitmap_cubical_complex/example)
add_subdirectory(src/Alpha_complex/example)
add_subdirectory(src/Alpha_complex/test)
add_subdirectory(src/Spatial_searching/example)
add_subdirectory(src/Spatial_searching/test)
add_subdirectory(src/Subsampling/example)
add_subdirectory(src/Subsampling/test)
add_subdirectory(src/Tangential_complex/example)
add_subdirectory(src/Tangential_complex/test)
add_subdirectory(src/Tangential_complex/benchmark)
add_subdirectory(src/Bottleneck_distance/example)
add_subdirectory(src/Bottleneck_distance/test)
add_subdirectory(src/Bottleneck_distance/benchmark)
add_subdirectory(src/Rips_complex/example)
add_subdirectory(src/Rips_complex/test)
add_subdirectory(src/Nerve_GIC/example)
add_subdirectory(src/Nerve_GIC/test)

# data points generator
add_subdirectory(data/points/generator)

add_subdirectory(src/GudhUI)

# specific for cython module
add_subdirectory(${GUDHI_CYTHON_PATH})
