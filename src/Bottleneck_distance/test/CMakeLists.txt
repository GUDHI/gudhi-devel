project(Bottleneck_distance_tests)

if (NOT CGAL_VERSION VERSION_LESS 4.11.0)
  include(GUDHI_test_coverage)

  add_executable ( Bottleneck_distance_test_unit bottleneck_unit_test.cpp )
  target_link_libraries(Bottleneck_distance_test_unit ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(Bottleneck_distance_test_unit ${TBB_LIBRARIES})
  endif(TBB_FOUND)

  gudhi_add_coverage_test(Bottleneck_distance_test_unit)

endif (NOT CGAL_VERSION VERSION_LESS 4.11.0)
