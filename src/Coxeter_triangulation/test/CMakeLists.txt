project(Coxeter_triangulation_test)

include(GUDHI_test_coverage)


add_executable ( Coxeter_triangulation_permutahedral_representation_test perm_rep_test.cpp )
target_link_libraries(Coxeter_triangulation_permutahedral_representation_test ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_executable ( Coxeter_triangulation_freudenthal_triangulation_test freud_triang_test.cpp )
add_executable ( Coxeter_triangulation_functions_test function_test.cpp )
add_executable ( Coxeter_triangulation_oracle_test oracle_test.cpp )
add_executable ( Coxeter_triangulation_manifold_tracing_test manifold_tracing_test.cpp )
add_executable ( Coxeter_triangulation_cell_complex_test cell_complex_test.cpp )
add_executable ( Coxeter_triangulation_projective_plane_test projective_plane_test.cpp )
add_executable ( Coxeter_triangulation_fun_on_RP2_test fun_on_RP2_test.cpp )
add_executable ( Coxeter_triangulation_CP2_test CP2_test.cpp )
add_executable ( Coxeter_triangulation_fun_on_CP2_test fun_on_CP2_test.cpp )

if (TBB_FOUND)
   target_link_libraries(Coxeter_triangulation_cell_complex_test ${TBB_LIBRARIES})
   target_link_libraries(Coxeter_triangulation_projective_plane_test ${TBB_LIBRARIES})
   target_link_libraries(Coxeter_triangulation_fun_on_RP2_test ${TBB_LIBRARIES})
   target_link_libraries(Coxeter_triangulation_CP2_test ${TBB_LIBRARIES})
   target_link_libraries(Coxeter_triangulation_fun_on_CP2_test ${TBB_LIBRARIES})
endif()
