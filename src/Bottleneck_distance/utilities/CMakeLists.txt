project(Bottleneck_distance_utilities)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
  add_executable (bottleneck_distance bottleneck_distance.cpp)
  if (TBB_FOUND)
    target_link_libraries(bottleneck_distance ${TBB_LIBRARIES})
  endif(TBB_FOUND)

  add_test(NAME Bottleneck_distance_utilities_Bottleneck_read_file
      COMMAND $<TARGET_FILE:bottleneck_distance>
      "${CMAKE_SOURCE_DIR}/data/persistence_diagram/first.pers" "${CMAKE_SOURCE_DIR}/data/persistence_diagram/second.pers")

  install(TARGETS bottleneck_distance DESTINATION bin)

endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
