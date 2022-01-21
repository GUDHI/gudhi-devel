project(common_benchmark)

add_executable(Graph_simplicial_complex_benchmark Graph_simplicial_complex_benchmark.cpp)

if (TBB_FOUND)
  target_link_libraries(Graph_simplicial_complex_benchmark ${TBB_LIBRARIES})
endif()
