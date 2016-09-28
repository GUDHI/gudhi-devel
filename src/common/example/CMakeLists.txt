cmake_minimum_required(VERSION 2.6)
project(Common_examples)

add_executable ( vector_double_off_reader example_vector_double_points_off_reader.cpp )
target_link_libraries(vector_double_off_reader  ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
add_test(vector_double_off_reader ${CMAKE_CURRENT_BINARY_DIR}/vector_double_off_reader ${CMAKE_SOURCE_DIR}/data/points/SO3_10000.off)

if(CGAL_FOUND)
  add_executable ( cgal_3D_off_reader example_CGAL_3D_points_off_reader.cpp )
  target_link_libraries(cgal_3D_off_reader  ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
  add_test(cgal_3D_off_reader ${CMAKE_CURRENT_BINARY_DIR}/cgal_3D_off_reader ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off)

  # need CGAL 4.7
  if (NOT CGAL_VERSION VERSION_LESS 4.7.0)
    if (EIGEN3_FOUND)
      add_executable ( cgal_off_reader example_CGAL_points_off_reader.cpp )
      target_link_libraries(cgal_off_reader  ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
      add_test(cgal_off_reader ${CMAKE_CURRENT_BINARY_DIR}/cgal_off_reader ${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off)
    endif(EIGEN3_FOUND)
  endif (NOT CGAL_VERSION VERSION_LESS 4.7.0)
endif()
