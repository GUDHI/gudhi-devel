add_executable_with_targets(Witness_complex_example_nearest_landmark_table example_nearest_landmark_table.cpp TBB::tbb)
add_test(NAME Witness_complex_example_nearest_landmark_table
    COMMAND $<TARGET_FILE:Witness_complex_example_nearest_landmark_table>)

add_executable_with_targets(Witness_complex_example_off example_witness_complex_off.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Witness_complex_example_sphere example_witness_complex_sphere.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Witness_complex_example_strong_off example_strong_witness_complex_off.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)

# CGAL and Eigen3 are required for Euclidean version of Witness
if(TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  add_test(NAME Witness_complex_example_off_test_torus
      COMMAND $<TARGET_FILE:Witness_complex_example_off>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "20" "1.0" "3")
  add_test(NAME Witness_complex_example_test_sphere_10
      COMMAND $<TARGET_FILE:Witness_complex_example_sphere> "10")
  add_test(NAME Witness_complex_example_strong_off_test_torus
      COMMAND $<TARGET_FILE:Witness_complex_example_strong_off>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "20" "1.0" "3")
endif()
