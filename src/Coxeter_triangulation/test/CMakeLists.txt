add_executable_with_targets(Coxeter_triangulation_permutahedral_representation_test perm_rep_test.cpp Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Coxeter_triangulation_freudenthal_triangulation_test freud_triang_test.cpp Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Coxeter_triangulation_functions_test function_test.cpp Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Coxeter_triangulation_oracle_test oracle_test.cpp Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Coxeter_triangulation_manifold_tracing_test manifold_tracing_test.cpp Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Coxeter_triangulation_cell_complex_test cell_complex_test.cpp Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Coxeter_triangulation_random_orthogonal_matrix_function_test random_orthogonal_matrix_function_test.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)

if (TARGET Eigen3::Eigen)
  include(GUDHI_boost_test)

  gudhi_add_boost_test(Coxeter_triangulation_permutahedral_representation_test)
  gudhi_add_boost_test(Coxeter_triangulation_freudenthal_triangulation_test)
  gudhi_add_boost_test(Coxeter_triangulation_functions_test)
  gudhi_add_boost_test(Coxeter_triangulation_oracle_test)
  gudhi_add_boost_test(Coxeter_triangulation_manifold_tracing_test)
  gudhi_add_boost_test(Coxeter_triangulation_cell_complex_test)
  
  # because of random_orthogonal_matrix inclusion
  if (TARGET CGAL::CGAL)
    gudhi_add_boost_test(Coxeter_triangulation_random_orthogonal_matrix_function_test)
  endif()
endif()