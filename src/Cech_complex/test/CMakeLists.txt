if (CGAL_FOUND AND EIGEN3_FOUND)
  include(GUDHI_boost_test)

  add_executable ( Cech_complex_test_unit test_cech_complex.cpp )
  target_link_libraries(Cech_complex_test_unit CGAL::CGAL)
  if(TARGET TBB::tbb)
    target_link_libraries(Cech_complex_test_unit TBB::tbb)
  endif()

  # Do not forget to copy test files in current binary dir
  file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

  gudhi_add_boost_test(Cech_complex_test_unit)

endif()
