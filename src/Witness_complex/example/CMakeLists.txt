project(Witness_complex_examples)

add_executable ( Witness_complex_example_nearest_landmark_table example_nearest_landmark_table.cpp )
if(TARGET TBB::tbb)
  target_link_libraries(Witness_complex_example_nearest_landmark_table TBB::tbb)
endif()
add_test(NAME Witness_complex_example_nearest_landmark_table
    COMMAND $<TARGET_FILE:Witness_complex_example_nearest_landmark_table>)

# CGAL and Eigen3 are required for Euclidean version of Witness
if(CGAL_FOUND AND EIGEN3_FOUND)
  add_executable( Witness_complex_example_off example_witness_complex_off.cpp )
  target_link_libraries(Witness_complex_example_off CGAL::CGAL)
  if(TARGET TBB::tbb)
    target_link_libraries(Witness_complex_example_off TBB::tbb)
  endif()
  add_executable ( Witness_complex_example_sphere example_witness_complex_sphere.cpp )
  target_link_libraries(Witness_complex_example_sphere CGAL::CGAL)
  if(TARGET TBB::tbb)
    target_link_libraries(Witness_complex_example_sphere TBB::tbb)
  endif()
  
  add_executable( Witness_complex_example_strong_off example_strong_witness_complex_off.cpp )
  target_link_libraries(Witness_complex_example_strong_off CGAL::CGAL)
  if(TARGET TBB::tbb)
    target_link_libraries(Witness_complex_example_strong_off TBB::tbb)
  endif()

  add_test(NAME Witness_complex_example_off_test_torus
      COMMAND $<TARGET_FILE:Witness_complex_example_off>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "20" "1.0" "3")
  add_test(NAME Witness_complex_example_test_sphere_10
      COMMAND $<TARGET_FILE:Witness_complex_example_sphere> "10")
  add_test(NAME Witness_complex_example_strong_off_test_torus
      COMMAND $<TARGET_FILE:Witness_complex_example_strong_off>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "20" "1.0" "3")

endif(CGAL_FOUND AND EIGEN3_FOUND)
