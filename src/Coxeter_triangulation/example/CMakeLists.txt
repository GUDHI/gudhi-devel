project(Coxeter_triangulation_example)

add_executable ( Coxeter_triangulation_manifold_tracing_flat_torus_with_boundary_example manifold_tracing_flat_torus_with_boundary.cpp )
if (TBB_FOUND)
   target_link_libraries(Coxeter_triangulation_manifold_tracing_flat_torus_with_boundary_example ${TBB_LIBRARIES})
endif()
add_test(NAME Coxeter_triangulation_manifold_tracing_flat_torus_with_boundary_example
         COMMAND $<TARGET_FILE:Coxeter_triangulation_manifold_tracing_flat_torus_with_boundary_example>)

add_executable ( Coxeter_triangulation_manifold_tracing_custom_function_example manifold_tracing_custom_function.cpp )
if (TBB_FOUND)
   target_link_libraries(Coxeter_triangulation_manifold_tracing_custom_function_example ${TBB_LIBRARIES})
endif()
add_test(NAME Coxeter_triangulation_manifold_tracing_custom_function_example
         COMMAND $<TARGET_FILE:Coxeter_triangulation_manifold_tracing_custom_function_example>)
