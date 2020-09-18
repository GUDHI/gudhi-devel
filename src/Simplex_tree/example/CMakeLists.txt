project(Simplex_tree_examples)

add_executable ( Simplex_tree_example_from_cliques_of_graph simplex_tree_from_cliques_of_graph.cpp )
if (TBB_FOUND)
  target_link_libraries(Simplex_tree_example_from_cliques_of_graph ${TBB_LIBRARIES})
endif()
add_test(NAME Simplex_tree_example_from_cliques_of_graph_2 COMMAND $<TARGET_FILE:Simplex_tree_example_from_cliques_of_graph>
    "${CMAKE_SOURCE_DIR}/data/filtered_simplicial_complex/Klein_bottle_complex.fsc" "2")
add_test(NAME Simplex_tree_example_from_cliques_of_graph_3 COMMAND $<TARGET_FILE:Simplex_tree_example_from_cliques_of_graph>
    "${CMAKE_SOURCE_DIR}/data/filtered_simplicial_complex/Klein_bottle_complex.fsc" "3")

add_executable ( Simplex_tree_example_simple_simplex_tree simple_simplex_tree.cpp )
if (TBB_FOUND)
  target_link_libraries(Simplex_tree_example_simple_simplex_tree ${TBB_LIBRARIES})
endif()
add_test(NAME Simplex_tree_example_simple_simplex_tree COMMAND $<TARGET_FILE:Simplex_tree_example_simple_simplex_tree>)

add_executable ( Simplex_tree_example_mini_simplex_tree mini_simplex_tree.cpp )
if (TBB_FOUND)
  target_link_libraries(Simplex_tree_example_mini_simplex_tree ${TBB_LIBRARIES})
endif()
add_test(NAME Simplex_tree_example_mini_simplex_tree COMMAND $<TARGET_FILE:Simplex_tree_example_mini_simplex_tree>)

# An example with Simplex-tree using CGAL alpha_shapes_3
if(GMP_FOUND AND NOT CGAL_VERSION VERSION_LESS 4.11.0)
  add_executable ( Simplex_tree_example_alpha_shapes_3_from_off example_alpha_shapes_3_simplex_tree_from_off_file.cpp )
  target_link_libraries(Simplex_tree_example_alpha_shapes_3_from_off ${GMP_LIBRARIES} ${CGAL_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(Simplex_tree_example_alpha_shapes_3_from_off ${TBB_LIBRARIES})
  endif()
  add_test(NAME Simplex_tree_example_alpha_shapes_3_from_off COMMAND $<TARGET_FILE:Simplex_tree_example_alpha_shapes_3_from_off>
      "${CMAKE_SOURCE_DIR}/data/points/bunny_5000.off")

endif()

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
  add_executable ( Simplex_tree_example_cech_complex_cgal_mini_sphere_3d cech_complex_cgal_mini_sphere_3d.cpp )
  target_link_libraries(Simplex_tree_example_cech_complex_cgal_mini_sphere_3d Boost::program_options ${CGAL_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(Simplex_tree_example_cech_complex_cgal_mini_sphere_3d ${TBB_LIBRARIES})
  endif()
  add_test(NAME Simplex_tree_example_cech_complex_cgal_mini_sphere_3d COMMAND $<TARGET_FILE:Simplex_tree_example_cech_complex_cgal_mini_sphere_3d>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" -r 0.3 -d 3)
endif ()

add_executable ( Simplex_tree_example_graph_expansion_with_blocker graph_expansion_with_blocker.cpp )
if (TBB_FOUND)
  target_link_libraries(Simplex_tree_example_graph_expansion_with_blocker ${TBB_LIBRARIES})
endif()

add_test(NAME Simplex_tree_example_graph_expansion_with_blocker COMMAND $<TARGET_FILE:Simplex_tree_example_graph_expansion_with_blocker>)

add_executable ( Simplex_tree_oscillating_rips_zigzag_filtration_from_points oscillating_rips_zigzag_filtration_from_points.cpp )
if (TBB_FOUND)
  target_link_libraries(Simplex_tree_oscillating_rips_zigzag_filtration_from_points ${TBB_LIBRARIES})
endif()
