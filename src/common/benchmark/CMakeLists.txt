project(common_benchmark)

add_executable(Graph_simplicial_complex_benchmark Graph_simplicial_complex_benchmark.cpp)

if (TARGET TBB::tbb)
  target_link_libraries(Graph_simplicial_complex_benchmark TBB::tbb)
endif()
