project(Simplex_tree_benchmark)

add_executable(simplex_tree_stars_benchmark simplex_tree_stars_benchmark.cpp)
if(TARGET TBB::tbb)
  target_link_libraries(simplex_tree_stars_benchmark TBB::tbb)
endif()
