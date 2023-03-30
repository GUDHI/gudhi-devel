project(Bitmap_cubical_complex_benchmark)

add_executable(bitmap_cubical_complex_benchmark bitmap_cubical_complex_benchmark.cpp)
if(TARGET TBB::tbb)
  target_link_libraries(bitmap_cubical_complex_benchmark TBB::tbb)
endif()

add_executable(persistence_2d persistence_2d.cpp)
if (TBB_FOUND)
  target_link_libraries(persistence_2d TBB::tbb)
endif()
add_test(NAME Compare_persistence_2d COMMAND $<TARGET_FILE:persistence_2d>)
