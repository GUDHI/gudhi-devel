add_executable_with_targets(Tangential_complex_test_TC test_tangential_complex.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)

if(TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  include(GUDHI_boost_test)

  gudhi_add_boost_test(Tangential_complex_test_TC)
endif()
