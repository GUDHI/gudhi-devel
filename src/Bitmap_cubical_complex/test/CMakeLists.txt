project(Bitmap_cubical_complex_tests)

include(GUDHI_test_coverage)

add_executable ( Bitmap_cubical_complex_test_unit Bitmap_test.cpp )
target_link_libraries(Bitmap_cubical_complex_test_unit ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(Bitmap_cubical_complex_test_unit ${TBB_LIBRARIES})
endif()

gudhi_add_coverage_test(Bitmap_cubical_complex_test_unit)
