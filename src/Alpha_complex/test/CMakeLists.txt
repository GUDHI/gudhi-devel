add_executable_with_targets(Alpha_complex_test_unit Alpha_complex_unit_test.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Alpha_complex_dim3_Epeck_dynamic_test_unit Alpha_complex_dim3_Epeck_dynamic_unit_test.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Alpha_complex_dim3_Epeck_static_test_unit Alpha_complex_dim3_Epeck_static_unit_test.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Alpha_complex_dim3_Epick_dynamic_test_unit Alpha_complex_dim3_Epick_dynamic_unit_test.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Alpha_complex_dim3_Epick_static_test_unit Alpha_complex_dim3_Epick_static_unit_test.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Delaunay_complex_Epeck_dynamic_test_unit Delaunay_complex_Epeck_dynamic_unit_test.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Delaunay_complex_Epeck_static_test_unit Delaunay_complex_Epeck_static_unit_test.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Delaunay_complex_Epick_dynamic_test_unit Delaunay_complex_Epick_dynamic_unit_test.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Delaunay_complex_Epick_static_test_unit Delaunay_complex_Epick_static_unit_test.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Incremental_delaunay_test Incremental_delaunay_test.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)

include(GUDHI_boost_test)
if (TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  # Do not forget to copy test files in current binary dir
  file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

  gudhi_add_boost_test(Alpha_complex_test_unit)
  gudhi_add_boost_test(Alpha_complex_dim3_Epeck_dynamic_test_unit)
  gudhi_add_boost_test(Alpha_complex_dim3_Epeck_static_test_unit)
  gudhi_add_boost_test(Alpha_complex_dim3_Epick_dynamic_test_unit)
  gudhi_add_boost_test(Alpha_complex_dim3_Epick_static_test_unit)
  gudhi_add_boost_test(Delaunay_complex_Epeck_dynamic_test_unit)
  gudhi_add_boost_test(Delaunay_complex_Epeck_static_test_unit)
  gudhi_add_boost_test(Delaunay_complex_Epick_dynamic_test_unit)
  gudhi_add_boost_test(Delaunay_complex_Epick_static_test_unit)
  gudhi_add_boost_test(Incremental_delaunay_test)
endif ()

add_executable_with_targets(Alpha_complex_3d_test_unit Alpha_complex_3d_unit_test.cpp CGAL::CGAL TBB::tbb)
add_executable_with_targets(Weighted_alpha_complex_3d_test_unit Weighted_alpha_complex_3d_unit_test.cpp CGAL::CGAL TBB::tbb)
add_executable_with_targets(Periodic_alpha_complex_3d_test_unit Periodic_alpha_complex_3d_unit_test.cpp CGAL::CGAL TBB::tbb)
add_executable_with_targets(Weighted_periodic_alpha_complex_3d_test_unit Weighted_periodic_alpha_complex_3d_unit_test.cpp CGAL::CGAL TBB::tbb)

if (TARGET CGAL::CGAL)
  gudhi_add_boost_test(Alpha_complex_3d_test_unit)
  gudhi_add_boost_test(Weighted_alpha_complex_3d_test_unit)
  gudhi_add_boost_test(Periodic_alpha_complex_3d_test_unit)
  gudhi_add_boost_test(Weighted_periodic_alpha_complex_3d_test_unit)
endif ()

add_executable_with_targets(Alpha_kernel_d_Epeck_dynamic_test_unit Alpha_kernel_d_Epeck_dynamic_unit_test.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Alpha_kernel_d_Epeck_static_test_unit Alpha_kernel_d_Epeck_static_unit_test.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Alpha_kernel_d_Epick_dynamic_test_unit Alpha_kernel_d_Epick_dynamic_unit_test.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Alpha_kernel_d_Epick_static_test_unit Alpha_kernel_d_Epick_static_unit_test.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)

add_executable_with_targets(Weighted_alpha_complex_test_unit Weighted_alpha_complex_unit_test.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)

add_executable_with_targets(Weighted_alpha_complex_non_visible_points_test_unit Weighted_alpha_complex_non_visible_points_unit_test.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)

add_executable_with_targets(Zero_weighted_alpha_complex_dynamic_test_unit Zero_weighted_alpha_complex_dynamic_unit_test.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Zero_weighted_alpha_complex_static_test_unit Zero_weighted_alpha_complex_static_unit_test.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)

if (TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  gudhi_add_boost_test(Alpha_kernel_d_Epeck_dynamic_test_unit)
  gudhi_add_boost_test(Alpha_kernel_d_Epeck_static_test_unit)
  gudhi_add_boost_test(Alpha_kernel_d_Epick_dynamic_test_unit)
  gudhi_add_boost_test(Alpha_kernel_d_Epick_static_test_unit)
  
  gudhi_add_boost_test(Weighted_alpha_complex_test_unit)

  gudhi_add_boost_test(Weighted_alpha_complex_non_visible_points_test_unit)

  gudhi_add_boost_test(Zero_weighted_alpha_complex_dynamic_test_unit)
  gudhi_add_boost_test(Zero_weighted_alpha_complex_static_test_unit)
endif ()
