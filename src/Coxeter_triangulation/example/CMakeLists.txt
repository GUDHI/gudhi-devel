project(Coxeter_triangulation_example)

if (EIGEN3_FOUND)
  # because of random_orthogonal_matrix inclusion
  if (TARGET CGAL::CGAL)
    add_executable ( Coxeter_triangulation_manifold_tracing_flat_torus_with_boundary_example manifold_tracing_flat_torus_with_boundary.cpp )
    target_link_libraries(Coxeter_triangulation_manifold_tracing_flat_torus_with_boundary_example CGAL::CGAL)
    add_test(NAME Coxeter_triangulation_manifold_tracing_flat_torus_with_boundary_example
             COMMAND $<TARGET_FILE:Coxeter_triangulation_manifold_tracing_flat_torus_with_boundary_example>)
  endif()
  
  add_executable ( Coxeter_triangulation_manifold_tracing_custom_function_example manifold_tracing_custom_function.cpp )
  target_link_libraries(Coxeter_triangulation_manifold_tracing_custom_function_example)
  add_test(NAME Coxeter_triangulation_manifold_tracing_custom_function_example
           COMMAND $<TARGET_FILE:Coxeter_triangulation_manifold_tracing_custom_function_example>)
  
  add_executable ( Coxeter_triangulation_cell_complex_from_basic_circle_manifold_example cell_complex_from_basic_circle_manifold.cpp )
  target_link_libraries(Coxeter_triangulation_cell_complex_from_basic_circle_manifold_example)
  add_test(NAME Coxeter_triangulation_cell_complex_from_basic_circle_manifold_example
           COMMAND $<TARGET_FILE:Coxeter_triangulation_cell_complex_from_basic_circle_manifold_example>)
endif()