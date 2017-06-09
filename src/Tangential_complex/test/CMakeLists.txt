cmake_minimum_required(VERSION 2.6)
project(Tangential_complex_tests)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
  include(GUDHI_test_coverage)

  add_executable( Tangential_complex_test_TC test_tangential_complex.cpp )
  target_link_libraries(Tangential_complex_test_TC ${CGAL_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(Tangential_complex_test_TC ${TBB_LIBRARIES})
  endif()

  gudhi_add_coverage_test(Tangential_complex_test_TC)
endif(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
