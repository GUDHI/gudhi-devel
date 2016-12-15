cmake_minimum_required(VERSION 2.6)
project(Spatial_searching_examples)

if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.8.1)
    if (EIGEN3_FOUND)
      add_executable( Spatial_searching_example_spatial_searching example_spatial_searching.cpp )
      target_link_libraries(Spatial_searching_example_spatial_searching ${CGAL_LIBRARY})
      add_test(Spatial_searching_example_spatial_searching
               ${CMAKE_CURRENT_BINARY_DIR}/Spatial_searching_example_spatial_searching)
    endif()
  endif()
endif()
