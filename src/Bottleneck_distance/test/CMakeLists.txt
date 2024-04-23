include(GUDHI_boost_test)

add_executable_with_targets(Bottleneck_distance_test_unit bottleneck_unit_test.cpp CGAL::CGAL TBB::tbb)

if (TARGET CGAL::CGAL)
  gudhi_add_boost_test(Bottleneck_distance_test_unit)
endif()