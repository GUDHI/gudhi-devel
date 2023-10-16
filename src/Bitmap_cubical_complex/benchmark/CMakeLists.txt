project(Bitmap_cubical_complex_benchmark)

add_executable(bitmap_cubical_complex_benchmark bitmap_cubical_complex_benchmark.cpp)
if(TARGET TBB::tbb)
  target_link_libraries(bitmap_cubical_complex_benchmark TBB::tbb)
endif()
