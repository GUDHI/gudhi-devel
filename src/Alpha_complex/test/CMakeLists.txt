project(Alpha_complex_tests)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.7.0)
  include(GUDHI_test_coverage)

  add_executable ( Alpha_complex_test_unit Alpha_complex_unit_test.cpp )
  target_link_libraries(Alpha_complex_test_unit ${CGAL_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(Alpha_complex_test_unit ${TBB_LIBRARIES})
  endif()

  add_executable ( Alpha_complex_3d_test_unit Alpha_complex_3d_unit_test.cpp )
  target_link_libraries(Alpha_complex_3d_test_unit ${CGAL_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(Alpha_complex_3d_test_unit ${TBB_LIBRARIES})
  endif()

  # Do not forget to copy test files in current binary dir
  file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  file(COPY "${CMAKE_SOURCE_DIR}/data/points/grid_5_5_5_in_0_1.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

  gudhi_add_coverage_test(Alpha_complex_test_unit)
  gudhi_add_coverage_test(Alpha_complex_3d_test_unit)

endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.7.0)

