add_executable_with_targets(Tangential_complex_example_basic example_basic.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Tangential_complex_example_with_perturb example_with_perturb.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)

if(TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  add_test(NAME Tangential_complex_example_basic
      COMMAND $<TARGET_FILE:Tangential_complex_example_basic>)
  add_test(NAME Tangential_complex_example_with_perturb
      COMMAND $<TARGET_FILE:Tangential_complex_example_with_perturb>)
endif()
