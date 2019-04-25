cmake_minimum_required(VERSION 2.6)
project(Coxeter_triangulation_examples)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
  add_executable ( Coxeter_triangulation_isosurface isosurface.cpp )
  add_executable ( Coxeter_triangulation_point_sample_reconstruction point_sample_reconstruction_fuzzy_sphere.cpp )
  add_executable ( Coxeter_triangulation_point_sample_reconstruction_knn point_sample_reconstruction_knn.cpp )
  # add_executable ( Coxeter_triangulation_point_sample_reconstruction point_sample_reconstruction.cpp )
  if (TBB_FOUND)
    target_link_libraries(Coxeter_triangulation_isosurface ${TBB_LIBRARIES})
    target_link_libraries(Coxeter_triangulation_point_sample_reconstruction ${TBB_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
    target_link_libraries(Coxeter_triangulation_point_sample_reconstruction_knn ${TBB_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
  endif()
  add_executable ( Coxeter_triangulation_coxeter_cech coxeter_cech.cpp )
  add_executable ( Coxeter_triangulation_test_suite test_suite.cpp )
  if (TBB_FOUND)
    target_link_libraries(Coxeter_triangulation_coxeter_cech ${TBB_LIBRARIES})
    target_link_libraries(Coxeter_triangulation_test_suite ${TBB_LIBRARIES})
    target_link_libraries(Coxeter_triangulation_test_suite ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    target_link_libraries(Coxeter_triangulation_test_suite ${Boost_PROGRAM_OPTIONS_LIBRARY})
  endif()
endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)

add_executable ( generate_sphere_product generate_sphere_product.cpp )
target_link_libraries(generate_sphere_product ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_link_libraries(generate_sphere_product ${Boost_PROGRAM_OPTIONS_LIBRARY})
target_link_libraries(generate_sphere_product ${TBB_LIBRARIES})

add_executable ( generate_projective_plane generate_projective_plane.cpp )
target_link_libraries(generate_projective_plane ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_link_libraries(generate_projective_plane ${Boost_PROGRAM_OPTIONS_LIBRARY})
target_link_libraries(generate_projective_plane ${TBB_LIBRARIES})

add_executable ( generate_torus_3 generate_torus_3.cpp )
target_link_libraries(generate_torus_3 ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_link_libraries(generate_torus_3 ${Boost_PROGRAM_OPTIONS_LIBRARY})
target_link_libraries(generate_torus_3 ${TBB_LIBRARIES})

add_executable ( generate_swiss_roll generate_swiss_roll.cpp )
target_link_libraries(generate_swiss_roll ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_link_libraries(generate_swiss_roll ${Boost_PROGRAM_OPTIONS_LIBRARY})
target_link_libraries(generate_swiss_roll ${TBB_LIBRARIES})

add_executable ( Hasse_diagram_weird hasse_weird.cpp )
target_link_libraries(Hasse_diagram_weird ${TBB_LIBRARIES})


# add_executable ( Witness_complex_example_nearest_landmark_table example_nearest_landmark_table.cpp )
# if (TBB_FOUND)
#   target_link_libraries(Witness_complex_example_nearest_landmark_table ${TBB_LIBRARIES})
# endif()
# add_test(NAME Witness_complex_example_nearest_landmark_table
#     COMMAND $<TARGET_FILE:Witness_complex_example_nearest_landmark_table>)

# # CGAL and Eigen3 are required for Euclidean version of Witness
# if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
#   add_executable( Witness_complex_example_off example_witness_complex_off.cpp )
#   add_executable ( Witness_complex_example_sphere example_witness_complex_sphere.cpp )
  
#   add_executable( Witness_complex_example_strong_off example_strong_witness_complex_off.cpp )
#   target_link_libraries(Witness_complex_example_strong_off)

#   add_test(NAME Witness_complex_example_off_test_torus
#       COMMAND $<TARGET_FILE:Witness_complex_example_off>
#       "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "20" "1.0" "3")
#   add_test(NAME Witness_complex_example_test_sphere_10
#       COMMAND $<TARGET_FILE:Witness_complex_example_sphere> "10")
#   add_test(NAME Witness_complex_example_strong_off_test_torus
#       COMMAND $<TARGET_FILE:Witness_complex_example_strong_off>
#       "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "20" "1.0" "3")
      
#   install(TARGETS Witness_complex_example_off DESTINATION bin)
#   install(TARGETS Witness_complex_example_sphere DESTINATION bin)
#   install(TARGETS Witness_complex_example_strong_off DESTINATION bin)

# endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
