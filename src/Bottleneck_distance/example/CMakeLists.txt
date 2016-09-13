cmake_minimum_required(VERSION 2.6)
project(Bottleneck_distance_examples)

# need CGAL 4.6
# cmake -DCGAL_DIR=~/workspace/CGAL-4.6-beta1 ../../..
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.6.0)
    if (EIGEN3_FOUND)
      add_executable (bottleneck_example bottleneck_example.cpp)
    endif()
  endif ()
endif()
