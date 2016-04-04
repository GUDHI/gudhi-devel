cmake_minimum_required(VERSION 2.6)
project(GUDHIBottleneckExample)

# need CGAL 4.6
# cmake -DCGAL_DIR=~/workspace/CGAL-4.6-beta1 ../../..
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.6.0)
    message(STATUS "CGAL version: ${CGAL_VERSION}.")

    include( ${CGAL_USE_FILE} )
 
    find_package(Eigen3 3.1.0)
    if (EIGEN3_FOUND)
      message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
      include( ${EIGEN3_USE_FILE} )
      include_directories (BEFORE "../../include")

      add_executable (basic_example basic.cpp)
      #add_test(dtoffrw_tore3D ${CMAKE_CURRENT_BINARY_DIR}/dtoffrw ${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off 3)

    else()
      message(WARNING "Eigen3 not found. Version 3.1.0 is required.")
    endif()
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Alpha shapes feature. Version 4.6.0 is required.")
  endif ()
endif()
