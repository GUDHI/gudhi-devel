cmake_minimum_required(VERSION 2.6)
project(kernel_tests)

include(GUDHI_test_coverage)

add_executable ( kernel_test_unit test_kernel.cpp )
target_link_libraries(kernel_test_unit ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(kernel_test_unit ${TBB_LIBRARIES})
endif()

gudhi_add_coverage_test(kernel_test_unit)
