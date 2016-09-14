cmake_minimum_required(VERSION 2.6)
project(Bottleneck_distance_examples)

# requires CGAL 4.8
# cmake -DCGAL_DIR=~/workspace/CGAL-4.8 ../../..
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.8.0)
    if (EIGEN3_FOUND)
      add_executable (bottleneck_example bottleneck_example.cpp)
    endif()
  endif ()
endif()
