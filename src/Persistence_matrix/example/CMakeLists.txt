project(Matrix_examples)

add_executable(Matrix_examples_representative_cycles representative_cycles_from_matrix.cpp)
if(TARGET TBB::tbb)
  target_link_libraries(Matrix_examples_representative_cycles TBB::tbb)
endif()
add_test(NAME Matrix_examples_representative_cycles COMMAND $<TARGET_FILE:Matrix_examples_representative_cycles>)

add_executable(Matrix_examples_simplex_tree_matrix simplex_tree_to_matrix.cpp)
if(TARGET TBB::tbb)
  target_link_libraries(Matrix_examples_simplex_tree_matrix TBB::tbb)
endif()
add_test(NAME Matrix_examples_simplex_tree_matrix COMMAND $<TARGET_FILE:Matrix_examples_simplex_tree_matrix>)

