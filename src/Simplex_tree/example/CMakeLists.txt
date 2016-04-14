cmake_minimum_required(VERSION 2.6)
project(GUDHISimplexTreeFromFile)

add_executable ( simplex_tree_from_cliques_of_graph simplex_tree_from_cliques_of_graph.cpp )
if (TBB_FOUND)
  target_link_libraries(simplex_tree_from_cliques_of_graph ${TBB_RELEASE_LIBRARY})
endif()
add_test(simplex_tree_from_cliques_of_graph_2 ${CMAKE_CURRENT_BINARY_DIR}/simplex_tree_from_cliques_of_graph ${CMAKE_SOURCE_DIR}/data/points/Klein_bottle_complex.txt 2)
add_test(simplex_tree_from_cliques_of_graph_3 ${CMAKE_CURRENT_BINARY_DIR}/simplex_tree_from_cliques_of_graph ${CMAKE_SOURCE_DIR}/data/points/Klein_bottle_complex.txt 3)

add_executable ( simple_simplex_tree simple_simplex_tree.cpp )
if (TBB_FOUND)
  target_link_libraries(simple_simplex_tree ${TBB_RELEASE_LIBRARY})
endif()
add_test(simple_simplex_tree ${CMAKE_CURRENT_BINARY_DIR}/simple_simplex_tree)

add_executable ( mini_simplex_tree mini_simplex_tree.cpp )
add_test(mini_simplex_tree ${CMAKE_CURRENT_BINARY_DIR}/mini_simplex_tree)


# An example with Simplex-tree using CGAL alpha_shapes_3
if(GMP_FOUND AND CGAL_FOUND)
  add_executable ( simplex_tree_from_alpha_shapes_3 simplex_tree_from_alpha_shapes_3.cpp )
  target_link_libraries(simplex_tree_from_alpha_shapes_3 ${GMP_LIBRARIES} ${CGAL_LIBRARY} ${Boost_SYSTEM_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(simplex_tree_from_alpha_shapes_3 ${TBB_RELEASE_LIBRARY})
  endif()
  add_test(simplex_tree_from_alpha_shapes_3 ${CMAKE_CURRENT_BINARY_DIR}/simplex_tree_from_alpha_shapes_3 ${CMAKE_SOURCE_DIR}/data/points/bunny_5000)
endif()
