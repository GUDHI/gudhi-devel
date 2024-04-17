project(Bottleneck_distance_benchmark)

if (NOT CGAL_VERSION VERSION_LESS 4.11.0)
  add_executable ( bottleneck_chrono bottleneck_chrono.cpp )
  target_link_libraries(bottleneck_chrono CGAL::CGAL)
  if(TARGET TBB::tbb)
    target_link_libraries(bottleneck_chrono TBB::tbb)
  endif()
endif(NOT CGAL_VERSION VERSION_LESS 4.11.0)
