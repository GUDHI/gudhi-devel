cmake_minimum_required(VERSION 2.6)
project(Subsampling_examples)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
  add_executable(Subsampling_example_pick_n_random_points example_pick_n_random_points.cpp)
  add_executable(Subsampling_example_choose_n_farthest_points example_choose_n_farthest_points.cpp)
  add_executable(Subsampling_example_custom_kernel example_custom_kernel.cpp)
  add_executable(Subsampling_example_sparsify_point_set example_sparsify_point_set.cpp)
  target_link_libraries(Subsampling_example_sparsify_point_set ${CGAL_LIBRARY})

  add_test(Subsampling_example_pick_n_random_points
           ${CMAKE_CURRENT_BINARY_DIR}/Subsampling_example_pick_n_random_points)
  add_test(Subsampling_example_choose_n_farthest_points
           ${CMAKE_CURRENT_BINARY_DIR}/Subsampling_example_choose_n_farthest_points)
  add_test(Subsampling_example_sparsify_point_set
           ${CMAKE_CURRENT_BINARY_DIR}/Subsampling_example_sparsify_point_set)
endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
