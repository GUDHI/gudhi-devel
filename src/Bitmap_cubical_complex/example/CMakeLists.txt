project(Bitmap_cubical_complex_examples)

add_executable ( Random_bitmap_cubical_complex Random_bitmap_cubical_complex.cpp )
if(TARGET TBB::tbb)
  target_link_libraries(Random_bitmap_cubical_complex TBB::tbb)
endif()
add_test(NAME Bitmap_cubical_complex_example_random COMMAND $<TARGET_FILE:Random_bitmap_cubical_complex>
    "2" "100" "100")
