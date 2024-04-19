project(Collapse_tests)

if (EIGEN3_FOUND)

  include(GUDHI_boost_test)
  
  add_executable ( Collapse_test_unit collapse_unit_test.cpp )
  target_link_libraries(Collapse_test_unit Eigen3::Eigen)
  if(TARGET TBB::tbb)
    target_link_libraries(Collapse_test_unit TBB::tbb)
  endif()
  gudhi_add_boost_test(Collapse_test_unit)

endif()