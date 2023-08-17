project(Subsampling_benchmark)

if (CGAL_FOUND)
  add_executable ( choose_n_farthest_points choose_n_farthest_points.cpp )
  target_link_libraries(choose_n_farthest_points CGAL::CGAL)
  add_test(NAME Compare_choose_n_farthest_points COMMAND $<TARGET_FILE:choose_n_farthest_points>)
endif()
