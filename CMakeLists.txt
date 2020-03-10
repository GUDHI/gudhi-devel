cmake_minimum_required(VERSION 3.1)

project(GUDHIdev)

include(CMakeGUDHIVersion.txt)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/src/cmake/modules/")

# Reset cache
set(GUDHI_MODULES "" CACHE INTERNAL "GUDHI_MODULES")
set(GUDHI_MISSING_MODULES "" CACHE INTERNAL "GUDHI_MISSING_MODULES")

# This variable is used by Cython CMakeLists.txt and by GUDHI_third_party_libraries to know its path
set(GUDHI_PYTHON_PATH "src/python")

# For third parties libraries management - To be done last as CGAL updates CMAKE_MODULE_PATH
include(GUDHI_third_party_libraries NO_POLICY_SCOPE)

include(GUDHI_compilation_flags)

# Add your new module in the list, order is not important
include(GUDHI_modules)

add_gudhi_module(common)
add_gudhi_module(Alpha_complex)
add_gudhi_module(Bitmap_cubical_complex)
add_gudhi_module(Bottleneck_distance)
add_gudhi_module(Contraction)
add_gudhi_module(Cech_complex)
add_gudhi_module(Hasse_complex)
add_gudhi_module(Persistence_representations)
add_gudhi_module(Persistent_cohomology)
add_gudhi_module(Rips_complex)
add_gudhi_module(Simplex_tree)
add_gudhi_module(Skeleton_blocker)
add_gudhi_module(Spatial_searching)
add_gudhi_module(Subsampling)
add_gudhi_module(Tangential_complex)
add_gudhi_module(Toplex_map)
add_gudhi_module(Witness_complex)
add_gudhi_module(Nerve_GIC)

# Include module CMake subdirectories
# GUDHI_SUB_DIRECTORIES is managed in CMAKE_MODULE_PATH/GUDHI_modules.cmake
foreach(GUDHI_MODULE ${GUDHI_MODULES})
  foreach(GUDHI_SUB_DIRECTORY ${GUDHI_SUB_DIRECTORIES})
    if(EXISTS ${CMAKE_SOURCE_DIR}/src/${GUDHI_MODULE}/${GUDHI_SUB_DIRECTORY}/CMakeLists.txt)
      add_subdirectory(src/${GUDHI_MODULE}/${GUDHI_SUB_DIRECTORY}/)
    endif()
  endforeach()
endforeach()

add_subdirectory(src/GudhUI)

if (WITH_GUDHI_PYTHON)
  # specific for cython module
  add_subdirectory(${GUDHI_PYTHON_PATH})
else()
  message("++ Python module will not be compiled because WITH_GUDHI_PYTHON is set to OFF")
  set(GUDHI_MISSING_MODULES ${GUDHI_MISSING_MODULES} "python")
endif()

# For "make user_version" - Requires GUDHI_modules to be performed
include(GUDHI_user_version_target)
# For "make doxygen" - Requires GUDHI_USER_VERSION_DIR to be set - Done in GUDHI_user_version_target for dev version
include(GUDHI_doxygen_target)

configure_file(${CMAKE_SOURCE_DIR}/.github/for_maintainers/new_gudhi_version_creation.md "${CMAKE_CURRENT_BINARY_DIR}/" @ONLY)

message("++ GUDHI_MODULES list is:\"${GUDHI_MODULES}\"")
message("++ GUDHI_MISSING_MODULES list is:\"${GUDHI_MISSING_MODULES}\"")
