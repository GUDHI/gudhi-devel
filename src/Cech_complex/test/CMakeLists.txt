if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 5.0.1)
  include(GUDHI_boost_test)

  add_executable ( Cech_complex_test_unit test_cech_complex.cpp )
  if(TARGET TBB::tbb)
    target_link_libraries(Cech_complex_test_unit TBB::tbb)
  endif()

  # Do not forget to copy test files in current binary dir
  file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

  gudhi_add_boost_test(Cech_complex_test_unit)

endif()
