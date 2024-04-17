project(Collapse_tests)

if (EIGEN3_FOUND)

  include(GUDHI_boost_test)
  
  add_executable ( Collapse_test_unit collapse_unit_test.cpp )
  if(TARGET TBB::tbb)
    target_link_libraries(Collapse_test_unit TBB::tbb)
  endif()
  gudhi_add_boost_test(Collapse_test_unit)

endif()