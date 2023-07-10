project(Simplex_tree_benchmark)

add_executable(simplex_tree_cofaces_benchmark simplex_tree_cofaces_benchmark.cpp)
if(TARGET TBB::tbb)
  target_link_libraries(simplex_tree_cofaces_benchmark TBB::tbb)
endif()
