cmake_minimum_required(VERSION 2.6)
project(Subsampling_tests)

if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()
if (GPROF_PATH)
  # for gprof to make coverage reports - Jenkins
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
  add_executable( Subsampling_test_pick_n_random_points test_pick_n_random_points.cpp )
  target_link_libraries(Subsampling_test_pick_n_random_points ${CGAL_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

  add_executable( Subsampling_test_choose_n_farthest_points test_choose_n_farthest_points.cpp )
  target_link_libraries(Subsampling_test_choose_n_farthest_points ${CGAL_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

  add_executable(Subsampling_test_sparsify_point_set test_sparsify_point_set.cpp)
  target_link_libraries(Subsampling_test_sparsify_point_set ${CGAL_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

  add_test(Subsampling_test_pick_n_random_points ${CMAKE_CURRENT_BINARY_DIR}/Subsampling_test_pick_n_random_points
  # XML format for Jenkins xUnit plugin
  --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/Subsampling_test_pick_n_random_points_UT.xml --log_level=test_suite --report_level=no)

  add_test(Subsampling_test_choose_n_farthest_points ${CMAKE_CURRENT_BINARY_DIR}/Subsampling_test_choose_n_farthest_points
  # XML format for Jenkins xUnit plugin
  --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/Subsampling_test_choose_n_farthest_points_UT.xml --log_level=test_suite --report_level=no)

  add_test(Subsampling_test_sparsify_point_set ${CMAKE_CURRENT_BINARY_DIR}/Subsampling_test_sparsify_point_set
  # XML format for Jenkins xUnit plugin
  --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/Subsampling_test_sparsify_point_set_UT.xml --log_level=test_suite --report_level=no)
endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
