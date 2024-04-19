if(CGAL_FOUND AND EIGEN3_FOUND)
  project(Subsampling_benchmark)

  add_executable ( choose_n_farthest_points choose_n_farthest_points.cpp )
  target_link_libraries(choose_n_farthest_points CGAL::CGAL Eigen3::Eigen)
endif()
