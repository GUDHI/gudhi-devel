add_executable_with_targets(bottleneck_distance bottleneck_distance.cpp CGAL::CGAL TBB::tbb)

if(TARGET CGAL::CGAL)
  add_test(NAME Bottleneck_distance_utilities_Bottleneck_read_file
      COMMAND $<TARGET_FILE:bottleneck_distance>
      "${CMAKE_SOURCE_DIR}/data/persistence_diagram/first.pers" "${CMAKE_SOURCE_DIR}/data/persistence_diagram/second.pers")

  install(TARGETS bottleneck_distance DESTINATION bin)
endif()
