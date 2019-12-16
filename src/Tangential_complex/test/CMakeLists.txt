project(Tangential_complex_tests)

if(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
  include(GUDHI_boost_test)

  add_executable( Tangential_complex_test_TC test_tangential_complex.cpp )
  target_link_libraries(Tangential_complex_test_TC ${CGAL_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(Tangential_complex_test_TC ${TBB_LIBRARIES})
  endif()

  gudhi_add_boost_test(Tangential_complex_test_TC)
endif(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
