project(Coxeter_triangulation_test)

include(GUDHI_boost_test)

if (NOT EIGEN3_VERSION VERSION_LESS 3.1.0)
  add_executable ( Coxeter_triangulation_permutahedral_representation_test perm_rep_test.cpp )
  gudhi_add_boost_test(Coxeter_triangulation_permutahedral_representation_test)
  
  add_executable ( Coxeter_triangulation_freudenthal_triangulation_test freud_triang_test.cpp )
  gudhi_add_boost_test(Coxeter_triangulation_freudenthal_triangulation_test)
  
  add_executable ( Coxeter_triangulation_functions_test function_test.cpp )
  gudhi_add_boost_test(Coxeter_triangulation_functions_test)
  
  # because of random_orthogonal_matrix inclusion
  if (NOT CGAL_VERSION VERSION_LESS 4.11.0)
    add_executable ( Coxeter_triangulation_random_orthogonal_matrix_function_test random_orthogonal_matrix_function_test.cpp )
    target_link_libraries(Coxeter_triangulation_random_orthogonal_matrix_function_test ${CGAL_LIBRARY})
    gudhi_add_boost_test(Coxeter_triangulation_random_orthogonal_matrix_function_test)
  endif()
  
  add_executable ( Coxeter_triangulation_oracle_test oracle_test.cpp )
  gudhi_add_boost_test(Coxeter_triangulation_oracle_test)
  
  add_executable ( Coxeter_triangulation_manifold_tracing_test manifold_tracing_test.cpp )
  gudhi_add_boost_test(Coxeter_triangulation_manifold_tracing_test)
  
  add_executable ( Coxeter_triangulation_cell_complex_test cell_complex_test.cpp )
  gudhi_add_boost_test(Coxeter_triangulation_cell_complex_test)
endif()