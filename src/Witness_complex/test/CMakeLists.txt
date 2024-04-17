project(Witness_complex_tests)

include(GUDHI_boost_test)

add_executable ( Witness_complex_test_simple_witness_complex test_simple_witness_complex.cpp )
if(TARGET TBB::tbb)
  target_link_libraries(Witness_complex_test_simple_witness_complex TBB::tbb)
endif()

gudhi_add_boost_test(Witness_complex_test_simple_witness_complex)

# CGAL and Eigen3 are required for Euclidean version of Witness
if(CGAL_FOUND AND EIGEN3_FOUND)
  add_executable ( Witness_complex_test_euclidean_simple_witness_complex test_euclidean_simple_witness_complex.cpp )
  target_link_libraries(Witness_complex_test_euclidean_simple_witness_complex CGAL::CGAL)
  if(TARGET TBB::tbb)
    target_link_libraries(Witness_complex_test_euclidean_simple_witness_complex TBB::tbb)
  endif()

  gudhi_add_boost_test(Witness_complex_test_euclidean_simple_witness_complex)
endif(CGAL_FOUND AND EIGEN3_FOUND)
