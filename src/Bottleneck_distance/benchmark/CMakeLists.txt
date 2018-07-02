project(Bottleneck_distance_benchmark)

if (NOT CGAL_VERSION VERSION_LESS 4.8.1)
  add_executable ( bottleneck_chrono bottleneck_chrono.cpp )
  if (TBB_FOUND)
    target_link_libraries(bottleneck_chrono ${TBB_LIBRARIES})
  endif(TBB_FOUND)
endif(NOT CGAL_VERSION VERSION_LESS 4.8.1)
