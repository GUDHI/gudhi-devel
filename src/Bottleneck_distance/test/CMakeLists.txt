project(Bottleneck_distance_tests)

if(CGAL_FOUND)
  include(GUDHI_boost_test)

  add_executable ( Bottleneck_distance_test_unit bottleneck_unit_test.cpp )
  target_link_libraries(Bottleneck_distance_test_unit CGAL::CGAL)
  if(TARGET TBB::tbb)
    target_link_libraries(Bottleneck_distance_test_unit TBB::tbb)
  endif()

  gudhi_add_boost_test(Bottleneck_distance_test_unit)

endif(CGAL_FOUND)
