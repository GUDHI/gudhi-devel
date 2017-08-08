cmake_minimum_required(VERSION 2.6)
project(Simplex_tree_tests)

include(GUDHI_test_coverage)

add_executable ( Simplex_tree_test_unit simplex_tree_unit_test.cpp )
target_link_libraries(Simplex_tree_test_unit ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(Simplex_tree_test_unit ${TBB_LIBRARIES})
endif()

# Do not forget to copy test files in current binary dir
file(COPY "simplex_tree_for_unit_test.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

gudhi_add_coverage_test(Simplex_tree_test_unit)
