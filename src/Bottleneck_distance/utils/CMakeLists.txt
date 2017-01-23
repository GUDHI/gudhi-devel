cmake_minimum_required(VERSION 2.6)
project(Bottleneck_distance_utils)

# requires CGAL 4.8
# cmake -DCGAL_DIR=~/workspace/CGAL-4.8 ../../..
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.8.0)
    if (EIGEN3_FOUND)
      add_executable (alpha_rips_persistence_bottleneck_distance alpha_rips_persistence_bottleneck_distance.cpp)
      target_link_libraries(alpha_rips_persistence_bottleneck_distance ${Boost_SYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})
      if (TBB_FOUND)
        target_link_libraries(alpha_rips_persistence_bottleneck_distance ${TBB_LIBRARIES})
      endif(TBB_FOUND)
    endif(EIGEN3_FOUND)
  endif(NOT CGAL_VERSION VERSION_LESS 4.8.0)
endif(CGAL_FOUND)
