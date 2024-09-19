add_executable_with_targets(Coxeter_triangulation_manifold_tracing_flat_torus_with_boundary_example manifold_tracing_flat_torus_with_boundary.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Coxeter_triangulation_manifold_tracing_custom_function_example manifold_tracing_custom_function.cpp Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Coxeter_triangulation_cell_complex_from_basic_circle_manifold_example cell_complex_from_basic_circle_manifold.cpp Eigen3::Eigen TBB::tbb)

if (TARGET Eigen3::Eigen)
  # because of random_orthogonal_matrix inclusion
  if (TARGET CGAL::CGAL)
    add_test(NAME Coxeter_triangulation_manifold_tracing_flat_torus_with_boundary_example
             COMMAND $<TARGET_FILE:Coxeter_triangulation_manifold_tracing_flat_torus_with_boundary_example>)
  endif()
  
  add_test(NAME Coxeter_triangulation_manifold_tracing_custom_function_example
           COMMAND $<TARGET_FILE:Coxeter_triangulation_manifold_tracing_custom_function_example>)
  add_test(NAME Coxeter_triangulation_cell_complex_from_basic_circle_manifold_example
           COMMAND $<TARGET_FILE:Coxeter_triangulation_cell_complex_from_basic_circle_manifold_example>)
endif()