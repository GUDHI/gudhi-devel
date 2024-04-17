project(Bottleneck_distance_benchmark)

if (CGAL_FOUND)
  add_executable ( bottleneck_chrono bottleneck_chrono.cpp )
  target_link_libraries(bottleneck_chrono CGAL::CGAL)
  if(TARGET TBB::tbb)
    target_link_libraries(bottleneck_chrono TBB::tbb)
  endif()
endif(CGAL_FOUND)
