cmake_minimum_required(VERSION 2.6)
project(Bottleneck_distance_benchmark)


# requires CGAL 4.8
# cmake -DCGAL_DIR=~/workspace/CGAL-4.8 ../../..
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.8.0)
    add_executable ( bottleneck_chrono bottleneck_chrono.cpp )
    if (TBB_FOUND)
      target_link_libraries(bottleneck_chrono ${TBB_LIBRARIES})
    endif(TBB_FOUND)
  endif ()
endif()
