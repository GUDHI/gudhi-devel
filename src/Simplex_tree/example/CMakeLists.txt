cmake_minimum_required(VERSION 2.6)
project(GUDHISimplexTreeFromFile)

add_executable ( simplex_tree_from_file simplex_tree_from_file.cpp )
target_link_libraries(simplex_tree_from_file ${Boost_SYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY})
add_test(simplex_tree_from_file_2 ${CMAKE_CURRENT_BINARY_DIR}/simplex_tree_from_file ${CMAKE_SOURCE_DIR}/data/points/Klein_bottle_complex.txt 2)
add_test(simplex_tree_from_file_3 ${CMAKE_CURRENT_BINARY_DIR}/simplex_tree_from_file ${CMAKE_SOURCE_DIR}/data/points/Klein_bottle_complex.txt 3)

add_executable ( simple_simplex_tree simple_simplex_tree.cpp )
target_link_libraries(simple_simplex_tree ${Boost_SYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY})
add_test(simple_simplex_tree ${CMAKE_CURRENT_BINARY_DIR}/simple_simplex_tree)
   
# An example with Simplex-tree using CGAL alpha_shapes_3
if(GMPXX_FOUND AND GMP_FOUND AND CGAL_FOUND)
   message("CGAL_lib = ${CGAL_LIBRARIES_DIR}")
   message("GMPXX_LIBRARIES = ${GMPXX_LIBRARIES}")
   message("GMP_LIBRARIES = ${GMP_LIBRARIES}")

   add_executable ( simplex_tree_from_alpha_shapes_3 simplex_tree_from_alpha_shapes_3.cpp )
   target_link_libraries(simplex_tree_from_alpha_shapes_3 ${Boost_SYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} ${CGAL_LIBRARY})
   add_test(simplex_tree_from_alpha_shapes_3 ${CMAKE_CURRENT_BINARY_DIR}/simplex_tree_from_alpha_shapes_3 ${CMAKE_SOURCE_DIR}/data/points/bunny_5000)
endif()  
