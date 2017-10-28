cmake_minimum_required(VERSION 2.6)
project(Rips_complex_tests)

include(GUDHI_test_coverage)

add_executable ( Rips_complex_test_unit test_rips_complex.cpp )
target_link_libraries(Rips_complex_test_unit ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(Rips_complex_test_unit ${TBB_LIBRARIES})
endif()

# Do not forget to copy test files in current binary dir
file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY "${CMAKE_SOURCE_DIR}/data/distance_matrix/full_square_distance_matrix.csv" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

gudhi_add_coverage_test(Rips_complex_test_unit)
