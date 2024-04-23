include(GUDHI_boost_test)

add_executable_with_targets(Witness_complex_test_simple_witness_complex test_simple_witness_complex.cpp TBB::tbb)
gudhi_add_boost_test(Witness_complex_test_simple_witness_complex)

add_executable_with_targets(Witness_complex_test_euclidean_simple_witness_complex test_euclidean_simple_witness_complex.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
# CGAL and Eigen3 are required for Euclidean version of Witness
if(TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  gudhi_add_boost_test(Witness_complex_test_euclidean_simple_witness_complex)
endif()
