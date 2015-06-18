cmake_minimum_required(VERSION 2.6)
project(GUDHIDelaunayTriangulationOffFileReadWrite)

# need CGAL 4.6
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.6.0)
    message(STATUS "CGAL version: ${CGAL_VERSION}.")

    include( ${CGAL_USE_FILE} )
 
    find_package(Eigen3 3.1.0)
    if (EIGEN3_FOUND)
      message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
      include( ${EIGEN3_USE_FILE} )
      
      add_executable ( dtoffrw Delaunay_triangulation_off_rw.cpp )
      target_link_libraries(dtoffrw  ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
      add_test(dtoffrw ${CMAKE_CURRENT_BINARY_DIR}/dtoffrw ${CMAKE_SOURCE_DIR}/data/points/alphashapedoc.off ${CMAKE_CURRENT_BINARY_DIR}/result.off)

    else()
      message(WARNING "Eigen3 not found. Version 3.1.0 is required for Alpha shapes feature.")
    endif()
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Alpha shapes feature. Version 4.6.0 is required.")
  endif ()
endif()
