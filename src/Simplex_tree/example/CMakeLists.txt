cmake_minimum_required(VERSION 2.6)
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
add_test(NAME Simplex_tree_example_mini_simplex_tree COMMAND $<TARGET_FILE:Simplex_tree_example_mini_simplex_tree>)

install(TARGETS Simplex_tree_example_from_cliques_of_graph DESTINATION bin)
install(TARGETS Simplex_tree_example_simple_simplex_tree DESTINATION bin)
install(TARGETS Simplex_tree_example_mini_simplex_tree DESTINATION bin)

# An example with Simplex-tree using CGAL alpha_shapes_3
if(GMP_FOUND AND CGAL_FOUND)
  add_executable ( Simplex_tree_example_alpha_shapes_3_from_off example_alpha_shapes_3_simplex_tree_from_off_file.cpp )
  target_link_libraries(Simplex_tree_example_alpha_shapes_3_from_off ${GMP_LIBRARIES} ${CGAL_LIBRARY} ${Boost_SYSTEM_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(Simplex_tree_example_alpha_shapes_3_from_off ${TBB_LIBRARIES})
  endif()
  add_test(NAME Simplex_tree_example_alpha_shapes_3_from_off COMMAND $<TARGET_FILE:Simplex_tree_example_alpha_shapes_3_from_off>
      "${CMAKE_SOURCE_DIR}/data/points/bunny_5000.off")

  install(TARGETS Simplex_tree_example_alpha_shapes_3_from_off DESTINATION bin)

endif()

add_executable ( Simplex_tree_example_cech_complex_step_by_step cech_complex_step_by_step.cpp)
target_link_libraries(Simplex_tree_example_cech_complex_step_by_step ${Boost_PROGRAM_OPTIONS_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(Simplex_tree_example_cech_complex_step_by_step ${TBB_LIBRARIES})
endif()
add_test(NAME Simplex_tree_example_cech_complex_step_by_step COMMAND $<TARGET_FILE:Simplex_tree_example_cech_complex_step_by_step>
    "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" "-r" "12." "-d" "3")
install(TARGETS Simplex_tree_example_cech_complex_step_by_step DESTINATION bin)


#
# TO BE REMOVED !!
#

#add_executable ( rips_step_by_step rips_step_by_step.cpp)
#target_link_libraries(rips_step_by_step ${Boost_PROGRAM_OPTIONS_LIBRARY})
#if (TBB_FOUND)
#  target_link_libraries(rips_step_by_step ${TBB_LIBRARIES})
#endif()


add_executable ( block block.cpp )
if (TBB_FOUND)
  target_link_libraries(block ${TBB_LIBRARIES})
endif()
