cmake_minimum_required(VERSION 2.6)
project(GUDHIAlphaShapesExample)

# need CGAL 4.7
# cmake -DCGAL_DIR=~/workspace/CGAL-4.7-Ic-41 ../../..
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.7.0)
    message(STATUS "CGAL version: ${CGAL_VERSION}.")

    include( ${CGAL_USE_FILE} )
 
    find_package(Eigen3 3.1.0)
    if (EIGEN3_FOUND)
      message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
      include( ${EIGEN3_USE_FILE} )
      
      add_definitions(-DDEBUG_TRACES)

      add_executable ( dtoffrw Delaunay_triangulation_off_rw.cpp )
      target_link_libraries(dtoffrw  ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
      add_test(dtoffrw_tore3D ${CMAKE_CURRENT_BINARY_DIR}/dtoffrw ${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off 3)

      add_executable ( stfromdt Simplex_tree_from_delaunay_triangulation.cpp )
      target_link_libraries(stfromdt  ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
    else()
      message(WARNING "Eigen3 not found. Version 3.1.0 is required for Alpha shapes feature.")
    endif()
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Alpha shapes feature. Version 4.6.0 is required.")
  endif ()
endif()
