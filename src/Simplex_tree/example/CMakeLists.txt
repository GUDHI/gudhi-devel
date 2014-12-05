cmake_minimum_required(VERSION 2.6)
project(GUDHISimplexTreeFromFile)

# NEED TO FIND BOOST NEEDED COMPONENTS TO LINK WITH
find_package(Boost 1.45.0 COMPONENTS system thread)
message("Boost_lib = ${Boost_LIBRARIES}")

add_executable ( simplex_tree_from_file simplex_tree_from_file.cpp )
target_link_libraries(simplex_tree_from_file ${Boost_LIBRARIES})
   
add_executable ( simple_simplex_tree simple_simplex_tree.cpp )
target_link_libraries(simple_simplex_tree ${Boost_LIBRARIES})
   
# An example with Simplex-tree using CGAL alpha_shapes_3
if(GMPXX_FOUND AND GMP_FOUND AND CGAL_FOUND)
   message("CGAL_lib = ${CGAL_LIBRARIES_DIR}")
   message("GMPXX_LIBRARIES = ${GMPXX_LIBRARIES}")
   message("GMP_LIBRARIES = ${GMP_LIBRARIES}")

   add_executable ( simplex_tree_from_alpha_shapes_3 simplex_tree_from_alpha_shapes_3.cpp )
   target_link_libraries(simplex_tree_from_alpha_shapes_3 ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} ${CGAL_LIBRARY})
endif()  
