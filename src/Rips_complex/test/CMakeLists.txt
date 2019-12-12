project(Rips_complex_tests)

include(GUDHI_boost_test)

add_executable ( Rips_complex_test_unit test_rips_complex.cpp )
if (TBB_FOUND)
  target_link_libraries(Rips_complex_test_unit ${TBB_LIBRARIES})
endif()

# Do not forget to copy test files in current binary dir
file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY "${CMAKE_SOURCE_DIR}/data/distance_matrix/full_square_distance_matrix.csv" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

gudhi_add_boost_test(Rips_complex_test_unit)
