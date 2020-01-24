project(Rips_complex_utilities)

add_executable(rips_distance_matrix_persistence rips_distance_matrix_persistence.cpp)
target_link_libraries(rips_distance_matrix_persistence Boost::program_options)

add_executable(rips_persistence rips_persistence.cpp)
target_link_libraries(rips_persistence Boost::program_options)

add_executable(rips_correlation_matrix_persistence rips_correlation_matrix_persistence.cpp)
target_link_libraries(rips_correlation_matrix_persistence Boost::system Boost::program_options)

add_executable(sparse_rips_persistence sparse_rips_persistence.cpp)
target_link_libraries(sparse_rips_persistence Boost::program_options)

if (TBB_FOUND)
  target_link_libraries(rips_distance_matrix_persistence ${TBB_LIBRARIES})
  target_link_libraries(rips_persistence ${TBB_LIBRARIES})
  target_link_libraries(rips_correlation_matrix_persistence ${TBB_LIBRARIES})
  target_link_libraries(sparse_rips_persistence ${TBB_LIBRARIES})
endif()

add_test(NAME Rips_complex_utility_from_rips_distance_matrix COMMAND $<TARGET_FILE:rips_distance_matrix_persistence>
    "${CMAKE_SOURCE_DIR}/data/distance_matrix/full_square_distance_matrix.csv" "-r" "1.0" "-d" "3" "-p" "3" "-m" "0")
add_test(NAME Rips_complex_utility_from_rips_on_tore_3D COMMAND $<TARGET_FILE:rips_persistence>
    "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-r" "0.25" "-m" "0.5" "-d" "3" "-p" "3")
add_test(NAME Rips_complex_utility_from_rips_correlation_matrix COMMAND $<TARGET_FILE:rips_correlation_matrix_persistence>
    "${CMAKE_SOURCE_DIR}/data/correlation_matrix/lower_triangular_correlation_matrix.csv" "-c" "0.3" "-d" "3" "-p" "3" "-m" "0")
add_test(NAME Sparse_rips_complex_utility_on_tore_3D COMMAND $<TARGET_FILE:sparse_rips_persistence>
    "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-e" "0.5" "-m" "0.2" "-d" "3" "-p" "2")

install(TARGETS rips_distance_matrix_persistence DESTINATION bin)
install(TARGETS rips_persistence DESTINATION bin)
install(TARGETS rips_correlation_matrix_persistence DESTINATION bin)
install(TARGETS sparse_rips_persistence DESTINATION bin)
