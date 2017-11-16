cmake_minimum_required(VERSION 2.6)
project(Graph_induced_complex_tests)

include(GUDHI_test_coverage)

add_executable ( Graph_induced_complex_test_unit test_GIC.cpp )
target_link_libraries(Graph_induced_complex_test_unit ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(Graph_induced_complex_test_unit ${TBB_LIBRARIES})
endif()

file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

gudhi_add_coverage_test(Graph_induced_complex_test_unit)
