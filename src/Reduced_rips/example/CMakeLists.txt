add_executable_with_targets(Reduced_rips_example_from_off example_reduced_rips_from_off.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Reduced_rips_example_from_points example_reduced_rips_from_points.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Reduced_rips_example_from_distance_matrix example_reduced_rips_from_distance_matrix.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)

if(TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  add_test(NAME Reduced_rips_example_from_off COMMAND $<TARGET_FILE:Reduced_rips_example_from_off>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off")
  add_test(NAME Reduced_rips_example_from_points COMMAND $<TARGET_FILE:Reduced_rips_example_from_points>)
  add_test(NAME Reduced_rips_example_from_distance_matrix
      COMMAND $<TARGET_FILE:Reduced_rips_example_from_distance_matrix>)
endif()
