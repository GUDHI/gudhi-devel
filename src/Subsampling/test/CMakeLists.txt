project(Subsampling_tests)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
  include(GUDHI_test_coverage)

  add_executable( Subsampling_test_pick_n_random_points test_pick_n_random_points.cpp )
  target_link_libraries(Subsampling_test_pick_n_random_points ${CGAL_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

  add_executable( Subsampling_test_choose_n_farthest_points test_choose_n_farthest_points.cpp )
  target_link_libraries(Subsampling_test_choose_n_farthest_points ${CGAL_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

  add_executable(Subsampling_test_sparsify_point_set test_sparsify_point_set.cpp)
  target_link_libraries(Subsampling_test_sparsify_point_set ${CGAL_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

  gudhi_add_coverage_test(Subsampling_test_pick_n_random_points)
  gudhi_add_coverage_test(Subsampling_test_choose_n_farthest_points)
  gudhi_add_coverage_test(Subsampling_test_sparsify_point_set)
endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
