cmake_minimum_required(VERSION 2.6)
project(GUDHIDelaunayTriangulationOffFileReadWrite)

# need CGAL 4.7
if(CGAL_FOUND)
  add_executable ( cgal3Doffreader CGAL_3D_points_off_reader.cpp )
  target_link_libraries(cgal3Doffreader  ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
  add_test(cgal3Doffreader ${CMAKE_CURRENT_BINARY_DIR}/cgal3Doffreader ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off)

  if (NOT CGAL_VERSION VERSION_LESS 4.7.0)
    find_package(Eigen3 3.1.0)
    if (EIGEN3_FOUND)
      message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
      include( ${EIGEN3_USE_FILE} )
      
      add_executable ( cgaloffreader CGAL_points_off_reader.cpp )
      target_link_libraries(cgaloffreader  ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
      add_test(cgaloffreader ${CMAKE_CURRENT_BINARY_DIR}/cgaloffreader ${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off)

    else()
      message(WARNING "Eigen3 not found. Version 3.1.0 is required for Alpha shapes feature.")
    endif()
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Alpha shapes feature. Version 4.6.0 is required.")
  endif ()
endif()
