project(Bottleneck_distance_utilities)

if(CGAL_FOUND AND EIGEN3_FOUND)
  add_executable (bottleneck_distance bottleneck_distance.cpp)
  target_link_libraries(bottleneck_distance CGAL::CGAL Eigen3::Eigen)
  if(TARGET TBB::tbb)
    target_link_libraries(bottleneck_distance TBB::tbb)
  endif()

  add_test(NAME Bottleneck_distance_utilities_Bottleneck_read_file
      COMMAND $<TARGET_FILE:bottleneck_distance>
      "${CMAKE_SOURCE_DIR}/data/persistence_diagram/first.pers" "${CMAKE_SOURCE_DIR}/data/persistence_diagram/second.pers")

  install(TARGETS bottleneck_distance DESTINATION bin)

endif(CGAL_FOUND AND EIGEN3_FOUND)
