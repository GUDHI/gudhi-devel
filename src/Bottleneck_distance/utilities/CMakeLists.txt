cmake_minimum_required(VERSION 2.6)
project(Bottleneck_distance_utilities)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
  add_executable (bottleneck_read_file bottleneck_read_file.cpp)
  if (TBB_FOUND)
    target_link_libraries(bottleneck_read_file ${TBB_LIBRARIES})
  endif(TBB_FOUND)

  add_test(NAME Bottleneck_distance_utilities_Bottleneck_read_file
      COMMAND $<TARGET_FILE:bottleneck_read_file>
      "${CMAKE_SOURCE_DIR}/data/persistence_diagram/first.pers" "${CMAKE_SOURCE_DIR}/data/persistence_diagram/second.pers")

  install(TARGETS bottleneck_read_file DESTINATION bin)

endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
