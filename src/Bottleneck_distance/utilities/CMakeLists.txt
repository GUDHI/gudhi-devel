cmake_minimum_required(VERSION 2.6)
project(Bottleneck_distance_utilities)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
  add_executable (bottleneck_read_file_example bottleneck_read_file_example.cpp)
  if (TBB_FOUND)
    target_link_libraries(alpha_rips_persistence_bottleneck_distance ${TBB_LIBRARIES})
  endif(TBB_FOUND)

  add_test(NAME Bottleneck_read_file_example
      COMMAND $<TARGET_FILE:bottleneck_read_file_example>
      "${CMAKE_SOURCE_DIR}/data/persistence_diagram/first.pers" "${CMAKE_SOURCE_DIR}/data/persistence_diagram/second.pers")

  install(TARGETS bottleneck_read_file_example DESTINATION bin)

endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
