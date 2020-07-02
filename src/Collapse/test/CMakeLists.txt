project(Collapse_tests)

include(GUDHI_boost_test)

add_executable ( Collapse_test_unit collapse_unit_test.cpp )
if (TBB_FOUND)
  target_link_libraries(Collapse_test_unit ${TBB_LIBRARIES})
endif()
gudhi_add_boost_test(Collapse_test_unit)
