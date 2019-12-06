project(Graph_induced_complex_tests)

if (NOT CGAL_VERSION VERSION_LESS 4.11.0)
  include(GUDHI_boost_test)

  add_executable ( Nerve_GIC_test_unit test_GIC.cpp )
  if (TBB_FOUND)
    target_link_libraries(Nerve_GIC_test_unit ${TBB_LIBRARIES})
  endif()

  file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

  gudhi_add_boost_test(Nerve_GIC_test_unit)

endif (NOT CGAL_VERSION VERSION_LESS 4.11.0)
