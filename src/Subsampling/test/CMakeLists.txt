add_executable_with_targets(Subsampling_test_pick_n_random_points test_pick_n_random_points.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Subsampling_test_choose_n_farthest_points test_choose_n_farthest_points.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Subsampling_test_sparsify_point_set test_sparsify_point_set.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)

if(TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  include(GUDHI_boost_test)

  gudhi_add_boost_test(Subsampling_test_pick_n_random_points)
  gudhi_add_boost_test(Subsampling_test_choose_n_farthest_points)
  gudhi_add_boost_test(Subsampling_test_sparsify_point_set)
endif()
