project(Subsampling_examples)

if(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
  add_executable(Subsampling_example_pick_n_random_points example_pick_n_random_points.cpp)
  add_executable(Subsampling_example_choose_n_farthest_points example_choose_n_farthest_points.cpp)
  add_executable(Subsampling_example_custom_kernel example_custom_kernel.cpp)
  add_executable(Subsampling_example_sparsify_point_set example_sparsify_point_set.cpp)
  target_link_libraries(Subsampling_example_sparsify_point_set ${CGAL_LIBRARY})

  add_test(NAME Subsampling_example_pick_n_random_points
      COMMAND $<TARGET_FILE:Subsampling_example_pick_n_random_points>)
  add_test(NAME Subsampling_example_choose_n_farthest_points
      COMMAND $<TARGET_FILE:Subsampling_example_choose_n_farthest_points>)
  add_test(NAME Subsampling_example_sparsify_point_set
      COMMAND $<TARGET_FILE:Subsampling_example_sparsify_point_set>)

endif(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
