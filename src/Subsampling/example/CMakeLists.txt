add_executable_with_targets(Subsampling_example_pick_n_random_points example_pick_n_random_points.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Subsampling_example_choose_n_farthest_points example_choose_n_farthest_points.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Subsampling_example_sparsify_point_set example_sparsify_point_set.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)

if(TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  add_test(NAME Subsampling_example_pick_n_random_points
      COMMAND $<TARGET_FILE:Subsampling_example_pick_n_random_points>)
  add_test(NAME Subsampling_example_choose_n_farthest_points
      COMMAND $<TARGET_FILE:Subsampling_example_choose_n_farthest_points>)
  add_test(NAME Subsampling_example_sparsify_point_set
      COMMAND $<TARGET_FILE:Subsampling_example_sparsify_point_set>)
endif()

add_executable(Subsampling_example_custom_distance example_custom_distance.cpp)
add_test(NAME Subsampling_example_custom_distance
      COMMAND $<TARGET_FILE:Subsampling_example_custom_distance>)
