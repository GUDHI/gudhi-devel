add_executable_with_targets(Random_bitmap_cubical_complex Random_bitmap_cubical_complex.cpp TBB::tbb)

add_test(NAME Bitmap_cubical_complex_example_random COMMAND $<TARGET_FILE:Random_bitmap_cubical_complex>
    "2" "100" "100")
