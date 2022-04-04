project(Bitmap_cubical_complex_benchmark)

add_executable(bitmap_cubical_complex_benchmark bitmap_cubical_complex_benchmark.cpp)
if (TBB_FOUND)
  target_link_libraries(bitmap_cubical_complex_benchmark ${TBB_LIBRARIES})
endif()
