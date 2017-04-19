cmake_minimum_required(VERSION 2.6)
project(Bitmap_cubical_complex_examples)

add_executable ( Bitmap_cubical_complex Bitmap_cubical_complex.cpp )
target_link_libraries(Bitmap_cubical_complex ${Boost_SYSTEM_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(Bitmap_cubical_complex ${TBB_LIBRARIES})
endif()

add_test(NAME Bitmap_cubical_complex_example_persistence_one_sphere COMMAND $<TARGET_FILE:Bitmap_cubical_complex>
    "${CMAKE_SOURCE_DIR}/data/bitmap/CubicalOneSphere.txt")

add_test(NAME Bitmap_cubical_complex_example_persistence_two_sphere COMMAND $<TARGET_FILE:Bitmap_cubical_complex>
    "${CMAKE_SOURCE_DIR}/data/bitmap/CubicalTwoSphere.txt")

add_executable ( Random_bitmap_cubical_complex Random_bitmap_cubical_complex.cpp )
target_link_libraries(Random_bitmap_cubical_complex ${Boost_SYSTEM_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(Random_bitmap_cubical_complex ${TBB_LIBRARIES})
endif()
add_test(NAME Bitmap_cubical_complex_example_random COMMAND $<TARGET_FILE:Random_bitmap_cubical_complex>
    "2" "100" "100")

add_executable ( Bitmap_cubical_complex_periodic_boundary_conditions Bitmap_cubical_complex_periodic_boundary_conditions.cpp )
target_link_libraries(Bitmap_cubical_complex_periodic_boundary_conditions ${Boost_SYSTEM_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(Bitmap_cubical_complex_periodic_boundary_conditions ${TBB_LIBRARIES})
endif()

add_test(NAME Bitmap_cubical_complex_example_periodic_boundary_conditions_2d_torus
    COMMAND $<TARGET_FILE:Bitmap_cubical_complex_periodic_boundary_conditions>
    "${CMAKE_SOURCE_DIR}/data/bitmap/2d_torus.txt")

add_test(NAME Bitmap_cubical_complex_example_periodic_boundary_conditions_3d_torus
    COMMAND $<TARGET_FILE:Bitmap_cubical_complex_periodic_boundary_conditions>
    "${CMAKE_SOURCE_DIR}/data/bitmap/3d_torus.txt")
