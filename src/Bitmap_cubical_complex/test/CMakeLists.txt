include(GUDHI_boost_test)

# Do not forget to copy test files in current binary dir
file(COPY "${CMAKE_SOURCE_DIR}/data/bitmap/sinusoid.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

add_executable_with_targets(Bitmap_cubical_complex_test_unit Bitmap_test.cpp TBB::tbb)

gudhi_add_boost_test(Bitmap_cubical_complex_test_unit)
