cmake_minimum_required(VERSION 2.6)
project(Alpha_complex_tests)

include(GUDHI_test_coverage)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.7.0)
  add_executable ( Alpha_complex_unit_test Alpha_complex_unit_test.cpp )
  target_link_libraries(Alpha_complex_unit_test  ${Boost_SYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY} ${CGAL_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(Alpha_complex_unit_test ${TBB_LIBRARIES})
  endif()

  # Do not forget to copy test files in current binary dir
  file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

  gudhi_add_coverage_test(Alpha_complex_unit_test)

endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.7.0)

