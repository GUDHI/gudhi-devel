cmake_minimum_required(VERSION 2.6)
project(GUDHIWitnessComplex)

add_executable ( simple_witness_complex simple_witness_complex.cpp )
add_test(simple_witness_complex ${CMAKE_CURRENT_BINARY_DIR}/simple_witness_complex)

add_executable( witness_complex_from_file witness_complex_from_file.cpp )
add_test( witness_complex_from_bunny &{CMAKE_CURRENT_BINARY_DIR}/witness_complex_from_file ${CMAKE_SOURCE_DIR}/data/points/bunny_5000 )   

add_executable( witness_complex_from_off witness_complex_from_off.cpp )

add_executable( witness_complex_from_wl_matrix witness_complex_from_wl_matrix.cpp )

# An example with Simplex-tree using CGAL alpha_shapes_3
#if(GMP_FOUND AND CGAL_FOUND)
#   message("CGAL_lib = ${CGAL_LIBRARIES_DIR}")
#   message("GMP_LIBRARIES = ${GMP_LIBRARIES}")   
#   INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR})
#   INCLUDE_DIRECTORIES(${CGAL_INCLUDE_DIRS})
#   add_executable ( simplex_tree_from_alpha_shapes_3 simplex_tree_from_alpha_shapes_3.cpp )
#   target_link_libraries(simplex_tree_from_alpha_shapes_3 ${GMP_LIBRARIES} ${CGAL_LIBRARY})
#   add_test(simplex_tree_from_alpha_shapes_3 ${CMAKE_CURRENT_BINARY_DIR}/simplex_tree_from_alpha_shapes_3 ${CMAKE_SOURCE_DIR}/data/points/bunny_5000)
#endif()  
