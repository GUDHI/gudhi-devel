cmake_minimum_required(VERSION 2.6)
project(GUDHISimplexTreeFromFile)

add_executable ( simplex_tree_from_cliques_of_graph simplex_tree_from_cliques_of_graph.cpp )
add_test(simplex_tree_from_cliques_of_graph_2 ${CMAKE_CURRENT_BINARY_DIR}/simplex_tree_from_cliques_of_graph ${CMAKE_SOURCE_DIR}/data/points/Klein_bottle_complex.txt 2)
add_test(simplex_tree_from_cliques_of_graph_3 ${CMAKE_CURRENT_BINARY_DIR}/simplex_tree_from_cliques_of_graph ${CMAKE_SOURCE_DIR}/data/points/Klein_bottle_complex.txt 3)

add_executable ( simple_simplex_tree simple_simplex_tree.cpp )
add_test(simple_simplex_tree ${CMAKE_CURRENT_BINARY_DIR}/simple_simplex_tree)

add_executable ( mini_simplex_tree mini_simplex_tree.cpp )
add_test(mini_simplex_tree ${CMAKE_CURRENT_BINARY_DIR}/mini_simplex_tree)


# An example with Simplex-tree using CGAL alpha_shapes_3
if(GMP_FOUND AND CGAL_FOUND)
  include( ${CGAL_USE_FILE} )
  # In CMakeLists.txt, when include(${CGAL_USE_FILE}), CXX_FLAGS are overwritten.
  # cf. http://doc.cgal.org/latest/Manual/installation.html#title40
  # A workaround is to add "-std=c++11" again.
  # A fix would be to use https://cmake.org/cmake/help/v3.1/prop_gbl/CMAKE_CXX_KNOWN_FEATURES.html 
  # or even better https://cmake.org/cmake/help/v3.1/variable/CMAKE_CXX_STANDARD.html
  # but it implies to use cmake version 3.1 at least.
  if(NOT MSVC)
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG(-std=c++11 COMPILER_SUPPORTS_CXX11)
    if(COMPILER_SUPPORTS_CXX11)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    endif()
  endif()
  # - End of workaround

  INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR})
  add_executable ( simplex_tree_from_alpha_shapes_3 simplex_tree_from_alpha_shapes_3.cpp )
  target_link_libraries(simplex_tree_from_alpha_shapes_3 ${GMP_LIBRARIES} ${CGAL_LIBRARY})
  add_test(simplex_tree_from_alpha_shapes_3 ${CMAKE_CURRENT_BINARY_DIR}/simplex_tree_from_alpha_shapes_3 ${CMAKE_SOURCE_DIR}/data/points/bunny_5000)
endif()
