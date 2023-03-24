project(Bottleneck_distance_tests)

if (NOT CGAL_VERSION VERSION_LESS 4.11.0)
  include(GUDHI_boost_test)

  add_executable ( Bottleneck_distance_test_unit bottleneck_unit_test.cpp )
  if(TARGET TBB::tbb)
    target_link_libraries(Bottleneck_distance_test_unit TBB::tbb)
  endif()

  gudhi_add_boost_test(Bottleneck_distance_test_unit)

endif (NOT CGAL_VERSION VERSION_LESS 4.11.0)
