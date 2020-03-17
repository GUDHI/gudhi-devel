project(Common_tests)

include(GUDHI_boost_test)

add_executable ( Common_test_points_off_reader test_points_off_reader.cpp )
add_executable ( Common_test_distance_matrix_reader test_distance_matrix_reader.cpp )
add_executable ( Common_test_persistence_intervals_reader test_persistence_intervals_reader.cpp )
add_executable ( Common_test_compute_edge_graph test_compute_edge_graph.cpp )

# Do not forget to copy test files in current binary dir
file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY "${CMAKE_SOURCE_DIR}/data/distance_matrix/lower_triangular_distance_matrix.csv" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY "${CMAKE_SOURCE_DIR}/data/distance_matrix/full_square_distance_matrix.csv" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY "${CMAKE_SOURCE_DIR}/src/common/test/persistence_intervals_with_dimension.pers" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY "${CMAKE_SOURCE_DIR}/src/common/test/persistence_intervals_with_field.pers" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY "${CMAKE_SOURCE_DIR}/src/common/test/persistence_intervals_without_dimension.pers" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

gudhi_add_boost_test(Common_test_points_off_reader)
gudhi_add_boost_test(Common_test_distance_matrix_reader)
gudhi_add_boost_test(Common_test_persistence_intervals_reader)
gudhi_add_boost_test(Common_test_compute_edge_graph)
