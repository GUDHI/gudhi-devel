add_executable_with_targets(rips_distance_matrix_persistence rips_distance_matrix_persistence.cpp Boost::program_options TBB::tbb)
add_executable_with_targets(rips_persistence rips_persistence.cpp Boost::program_options TBB::tbb)
add_executable_with_targets(rips_correlation_matrix_persistence rips_correlation_matrix_persistence.cpp Boost::program_options TBB::tbb)
add_executable_with_targets(sparse_rips_persistence sparse_rips_persistence.cpp Boost::program_options TBB::tbb)

if(TARGET Boost::program_options)
  add_test(NAME Rips_complex_utility_from_rips_distance_matrix COMMAND $<TARGET_FILE:rips_distance_matrix_persistence>
      "${CMAKE_SOURCE_DIR}/data/distance_matrix/full_square_distance_matrix.csv" "-r" "1.0" "-d" "3" "-p" "3" "-m" "0")
  install(TARGETS rips_distance_matrix_persistence DESTINATION bin)

  add_test(NAME Rips_complex_utility_from_rips_on_tore_3D COMMAND $<TARGET_FILE:rips_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-r" "0.25" "-m" "0.5" "-d" "3" "-p" "3")
  install(TARGETS rips_persistence DESTINATION bin)

  add_test(NAME Rips_complex_utility_from_rips_correlation_matrix COMMAND $<TARGET_FILE:rips_correlation_matrix_persistence>
      "${CMAKE_SOURCE_DIR}/data/correlation_matrix/lower_triangular_correlation_matrix.csv" "-c" "0.3" "-d" "3" "-p" "3" "-m" "0")
  install(TARGETS rips_correlation_matrix_persistence DESTINATION bin)

  add_test(NAME Sparse_rips_complex_utility_on_tore_3D COMMAND $<TARGET_FILE:sparse_rips_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-e" "0.5" "-m" "0.2" "-d" "3" "-p" "2")
  install(TARGETS sparse_rips_persistence DESTINATION bin)
endif()
