project(Bitmap_cubical_complex_tests)

include(GUDHI_boost_test)

# Do not forget to copy test files in current binary dir
file(COPY "${CMAKE_SOURCE_DIR}/data/bitmap/sinusoid.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

add_executable ( Bitmap_cubical_complex_test_unit Bitmap_test.cpp )
if (TBB_FOUND)
  target_link_libraries(Bitmap_cubical_complex_test_unit ${TBB_LIBRARIES})
endif()

gudhi_add_boost_test(Bitmap_cubical_complex_test_unit)
