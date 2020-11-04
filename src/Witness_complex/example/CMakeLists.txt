project(Witness_complex_examples)

add_executable ( Witness_complex_example_nearest_landmark_table example_nearest_landmark_table.cpp )
if (TBB_FOUND)
  target_link_libraries(Witness_complex_example_nearest_landmark_table ${TBB_LIBRARIES})
endif()
add_test(NAME Witness_complex_example_nearest_landmark_table
    COMMAND $<TARGET_FILE:Witness_complex_example_nearest_landmark_table>)

# CGAL and Eigen3 are required for Euclidean version of Witness
if(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
  add_executable( Witness_complex_example_off example_witness_complex_off.cpp )
  if (TBB_FOUND)
    target_link_libraries(Witness_complex_example_off ${TBB_LIBRARIES})
  endif()
  add_executable ( Witness_complex_example_sphere example_witness_complex_sphere.cpp )
  if (TBB_FOUND)
    target_link_libraries(Witness_complex_example_sphere ${TBB_LIBRARIES})
  endif()
  
  add_executable( Witness_complex_example_strong_off example_strong_witness_complex_off.cpp )
  if (TBB_FOUND)
    target_link_libraries(Witness_complex_example_strong_off ${TBB_LIBRARIES})
  endif()

  add_test(NAME Witness_complex_example_off_test_torus
      COMMAND $<TARGET_FILE:Witness_complex_example_off>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "20" "1.0" "3")
  add_test(NAME Witness_complex_example_test_sphere_10
      COMMAND $<TARGET_FILE:Witness_complex_example_sphere> "10")
  add_test(NAME Witness_complex_example_strong_off_test_torus
      COMMAND $<TARGET_FILE:Witness_complex_example_strong_off>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "20" "1.0" "3")

endif(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
