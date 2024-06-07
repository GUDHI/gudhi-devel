add_executable_with_targets(Cech_complex_test_unit test_cech_complex.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)

if (TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  include(GUDHI_boost_test)

  # Do not forget to copy test files in current binary dir
  file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

  gudhi_add_boost_test(Cech_complex_test_unit)
endif()
