project(Toplex_map_tests)

include(GUDHI_test_coverage)

add_executable( Toplex_map_unit_test toplex_map_unit_test.cpp )
target_link_libraries(Toplex_map_unit_test ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

gudhi_add_coverage_test(Toplex_map_unit_test)
