project(Graph_induced_complex_tests)

if (NOT CGAL_VERSION VERSION_LESS 4.8.1)
  include(GUDHI_test_coverage)

  add_executable ( Nerve_GIC_test_unit test_GIC.cpp )
  target_link_libraries(Nerve_GIC_test_unit ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(Nerve_GIC_test_unit ${TBB_LIBRARIES})
  endif()

  file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

  gudhi_add_coverage_test(Nerve_GIC_test_unit)

endif (NOT CGAL_VERSION VERSION_LESS 4.8.1)
