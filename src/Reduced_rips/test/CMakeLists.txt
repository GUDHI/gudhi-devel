set(REDUCED_RIPS_TEST_UNITS correctness inputs scalar_types edge_cases)

foreach(unit ${REDUCED_RIPS_TEST_UNITS})
  add_executable_with_targets(Reduced_rips_test_${unit} test_${unit}.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
endforeach()

if(TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  include(GUDHI_boost_test)
  # Copy the data file used by the OFF-file integration test into the test working directory.
  file(COPY "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  foreach(unit ${REDUCED_RIPS_TEST_UNITS})
    gudhi_add_boost_test(Reduced_rips_test_${unit})
  endforeach()
endif()
