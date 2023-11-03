project(Tangential_complex_tests)

if(CGAL_FOUND AND TARGET Eigen3::Eigen)
  include(GUDHI_boost_test)

  add_executable( Tangential_complex_test_TC test_tangential_complex.cpp )
  target_link_libraries(Tangential_complex_test_TC CGAL::CGAL Eigen3::Eigen)
  if(TARGET TBB::tbb)
    target_link_libraries(Tangential_complex_test_TC TBB::tbb)
  endif()

  gudhi_add_boost_test(Tangential_complex_test_TC)
endif()
