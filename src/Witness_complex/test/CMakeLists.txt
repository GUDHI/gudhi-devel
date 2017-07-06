cmake_minimum_required(VERSION 2.6)
project(Witness_complex_tests)

include(GUDHI_test_coverage)

add_executable ( Witness_complex_test_simple_witness_complex test_simple_witness_complex.cpp )
target_link_libraries(Witness_complex_test_simple_witness_complex ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(Witness_complex_test_simple_witness_complex ${TBB_LIBRARIES})
endif(TBB_FOUND)

gudhi_add_coverage_test(Witness_complex_test_simple_witness_complex)

# CGAL and Eigen3 are required for Euclidean version of Witness
if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
  add_executable ( Witness_complex_test_euclidean_simple_witness_complex test_euclidean_simple_witness_complex.cpp )
  target_link_libraries(Witness_complex_test_euclidean_simple_witness_complex ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(Witness_complex_test_euclidean_simple_witness_complex ${TBB_LIBRARIES})
  endif(TBB_FOUND)

  gudhi_add_coverage_test(Witness_complex_test_euclidean_simple_witness_complex)
endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
