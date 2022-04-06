project(Toplex_map_benchmark)

add_executable(Toplex_map_benchmark benchmark_tm.cpp)

if (TBB_FOUND)
  target_link_libraries(Toplex_map_benchmark ${TBB_LIBRARIES})
endif()
