project(Subsampling_examples)

if(CGAL_FOUND AND EIGEN3_FOUND)
  add_executable(Subsampling_example_pick_n_random_points example_pick_n_random_points.cpp)
  target_link_libraries(Subsampling_example_pick_n_random_points CGAL::CGAL)
  add_executable(Subsampling_example_choose_n_farthest_points example_choose_n_farthest_points.cpp)
  target_link_libraries(Subsampling_example_choose_n_farthest_points CGAL::CGAL)
  add_executable(Subsampling_example_sparsify_point_set example_sparsify_point_set.cpp)
  target_link_libraries(Subsampling_example_sparsify_point_set CGAL::CGAL)

  add_test(NAME Subsampling_example_pick_n_random_points
      COMMAND $<TARGET_FILE:Subsampling_example_pick_n_random_points>)
  add_test(NAME Subsampling_example_choose_n_farthest_points
      COMMAND $<TARGET_FILE:Subsampling_example_choose_n_farthest_points>)
  add_test(NAME Subsampling_example_sparsify_point_set
      COMMAND $<TARGET_FILE:Subsampling_example_sparsify_point_set>)
endif(CGAL_FOUND AND EIGEN3_FOUND)

add_executable(Subsampling_example_custom_distance example_custom_distance.cpp)
