cmake_minimum_required(VERSION 2.6)
project(Coxeter_triangulation_example)

add_executable ( Coxeter_triangulation_manifold_tracing_flat_torus_with_boundary_example manifold_tracing_flat_torus_with_boundary.cpp )

if (TBB_FOUND)
   target_link_libraries(Coxeter_triangulation_manifold_tracing_flat_torus_with_boundary_example ${TBB_LIBRARIES})
endif()

