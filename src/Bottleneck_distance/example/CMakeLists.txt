cmake_minimum_required(VERSION 2.6)
project(Bottleneck_distance_examples)

# requires CGAL 4.8
# cmake -DCGAL_DIR=~/workspace/CGAL-4.8 ../../..
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.8.0)
    if (EIGEN3_FOUND)
      add_executable (bottleneck_read_file_example bottleneck_read_file_example.cpp)
      add_executable (bottleneck_basic_example bottleneck_basic_example.cpp)

      add_test(bottleneck_basic_example ${CMAKE_CURRENT_BINARY_DIR}/bottleneck_basic_example)
    endif()
  endif ()
endif()
