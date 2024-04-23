add_executable_with_targets(Cech_complex_example_from_points cech_complex_example_from_points.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)

if (TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  add_test(NAME Cech_complex_example_from_points COMMAND $<TARGET_FILE:Cech_complex_example_from_points>)
endif()
