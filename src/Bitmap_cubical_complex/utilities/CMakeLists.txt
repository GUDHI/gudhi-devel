add_executable_with_targets(cubical_complex_persistence cubical_complex_persistence.cpp TBB::tbb)

add_test(NAME Bitmap_cubical_complex_utility_persistence_one_sphere COMMAND $<TARGET_FILE:cubical_complex_persistence>
    "${CMAKE_SOURCE_DIR}/data/bitmap/CubicalOneSphere.txt")

add_test(NAME Bitmap_cubical_complex_utility_persistence_two_sphere COMMAND $<TARGET_FILE:cubical_complex_persistence>
    "${CMAKE_SOURCE_DIR}/data/bitmap/CubicalTwoSphere.txt")

add_executable_with_targets(periodic_cubical_complex_persistence periodic_cubical_complex_persistence.cpp TBB::tbb)

add_test(NAME Bitmap_cubical_complex_utility_periodic_boundary_conditions_2d_torus
    COMMAND $<TARGET_FILE:periodic_cubical_complex_persistence>
    "${CMAKE_SOURCE_DIR}/data/bitmap/2d_torus.txt")

add_test(NAME Bitmap_cubical_complex_utility_periodic_boundary_conditions_3d_torus
    COMMAND $<TARGET_FILE:periodic_cubical_complex_persistence>
    "${CMAKE_SOURCE_DIR}/data/bitmap/3d_torus.txt")

install(TARGETS cubical_complex_persistence DESTINATION bin)
install(TARGETS periodic_cubical_complex_persistence DESTINATION bin)
