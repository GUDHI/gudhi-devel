if(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
  project(Subsampling_benchmark)

  add_executable ( choose_n_farthest_points choose_n_farthest_points.cpp )
  target_link_libraries(choose_n_farthest_points CGAL::CGAL)
endif()
