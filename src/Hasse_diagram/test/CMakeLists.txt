cmake_minimum_required(VERSION 2.6)
project(Hasse_diagrram_tests)

include(GUDHI_test_coverage)

add_executable ( Hasse_diagram_test_unit Hasse_diagram_test.cpp )
target_link_libraries(Hasse_diagram_test_unit ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(Hasse_diagram_test_unit ${TBB_LIBRARIES})
endif()

file(COPY "${CMAKE_SOURCE_DIR}/data/hasse_diagram/cw_decomposition_of_torus.hasse" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

gudhi_add_coverage_test(Hasse_diagram_test_unit)
