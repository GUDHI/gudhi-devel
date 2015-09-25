cmake_minimum_required(VERSION 2.6)
project(GUDHISimplexTreeFromFile)

add_executable ( Bitmap_cubical_complex Bitmap_cubical_complex.cpp )
target_link_libraries(Bitmap_cubical_complex ${Boost_SYSTEM_LIBRARY})	
add_test(Bitmap_cubical_complex ${CMAKE_CURRENT_BINARY_DIR}/Bitmap_cubical_complex ${CMAKE_SOURCE_DIR}/data/bitmap/CubicalOneSphere.txt)
add_test(Bitmap_cubical_complex ${CMAKE_CURRENT_BINARY_DIR}/Bitmap_cubical_complex ${CMAKE_SOURCE_DIR}/data/bitmap/CubicalTwoSphere.txt)

add_executable ( Random_bitmap_cubical_complex Random_bitmap_cubical_complex.cpp )
target_link_libraries(Random_bitmap_cubical_complex ${Boost_SYSTEM_LIBRARY})	
add_test(Random_bitmap_cubical_complex ${CMAKE_CURRENT_BINARY_DIR}/Random_bitmap_cubical_complex 2 100 100)

