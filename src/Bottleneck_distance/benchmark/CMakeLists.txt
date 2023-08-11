project(Bottleneck_distance_benchmark)

if(TARGET CGAL::CGAL)
  add_executable ( bottleneck_chrono bottleneck_chrono.cpp )
  target_link_libraries(bottleneck_chrono CGAL::CGAL)
  if(TARGET TBB::tbb)
    target_link_libraries(bottleneck_chrono TBB::tbb)
  endif()
endif(TARGET CGAL::CGAL)
