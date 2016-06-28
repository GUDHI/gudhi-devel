cmake_minimum_required(VERSION 2.6)
project(Spatial_searching_examples)

if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.8.0)
    message(STATUS "CGAL version: ${CGAL_VERSION}.")

    message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
    include( ${EIGEN3_USE_FILE} )
    include_directories (BEFORE "../../include")

    add_executable( Spatial_searching_example_spatial_searching example_spatial_searching.cpp )
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Spatial_searching examples. Version 4.8.0 is required.")
  endif ()
else()
  message(WARNING "CGAL not found. It is required for the Spatial_searching examples.")
endif()
