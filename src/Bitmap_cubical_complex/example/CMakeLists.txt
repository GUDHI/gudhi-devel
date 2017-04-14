cmake_minimum_required(VERSION 2.6)
project(Bitmap_cubical_complex_examples)

add_executable ( Bitmap_cubical_complex Bitmap_cubical_complex.cpp )
target_link_libraries(Bitmap_cubical_complex ${Boost_SYSTEM_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(Bitmap_cubical_complex ${TBB_LIBRARIES})
endif()

# Do not forget to copy test files in current binary dir
file(COPY "${CMAKE_SOURCE_DIR}/data/bitmap/CubicalOneSphere.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY "${CMAKE_SOURCE_DIR}/data/bitmap/CubicalTwoSphere.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

add_test(Bitmap_cubical_complex_one_sphere
    ${CMAKE_CURRENT_BINARY_DIR}/Bitmap_cubical_complex
    ${CMAKE_CURRENT_BINARY_DIR}/CubicalOneSphere.txt)

add_test(Bitmap_cubical_complex_two_sphere
    ${CMAKE_CURRENT_BINARY_DIR}/Bitmap_cubical_complex
    ${CMAKE_CURRENT_BINARY_DIR}/CubicalTwoSphere.txt)

add_executable ( Random_bitmap_cubical_complex Random_bitmap_cubical_complex.cpp )
target_link_libraries(Random_bitmap_cubical_complex ${Boost_SYSTEM_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(Random_bitmap_cubical_complex ${TBB_LIBRARIES})
endif()
add_test(Random_bitmap_cubical_complex ${CMAKE_CURRENT_BINARY_DIR}/Random_bitmap_cubical_complex 2 100 100)

add_executable ( Bitmap_cubical_complex_periodic_boundary_conditions Bitmap_cubical_complex_periodic_boundary_conditions.cpp )
target_link_libraries(Bitmap_cubical_complex_periodic_boundary_conditions ${Boost_SYSTEM_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(Bitmap_cubical_complex_periodic_boundary_conditions ${TBB_LIBRARIES})
endif()

# Do not forget to copy test files in current binary dir
file(COPY "${CMAKE_SOURCE_DIR}/data/bitmap/2d_torus.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY "${CMAKE_SOURCE_DIR}/data/bitmap/3d_torus.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

add_test(Bitmap_cubical_complex_periodic_2d_torus
    ${CMAKE_CURRENT_BINARY_DIR}/Bitmap_cubical_complex_periodic_boundary_conditions
    ${CMAKE_CURRENT_BINARY_DIR}/2d_torus.txt)
add_test(Bitmap_cubical_complex_periodic_3d_torus
    ${CMAKE_CURRENT_BINARY_DIR}/Bitmap_cubical_complex_periodic_boundary_conditions
    ${CMAKE_CURRENT_BINARY_DIR}/3d_torus.txt)

