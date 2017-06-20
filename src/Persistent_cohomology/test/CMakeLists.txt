cmake_minimum_required(VERSION 2.6)
project(Persistent_cohomology_tests)

include(GUDHI_test_coverage)

add_executable ( Persistent_cohomology_test_unit persistent_cohomology_unit_test.cpp )
target_link_libraries(Persistent_cohomology_test_unit ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_executable ( Persistent_cohomology_test_betti_numbers betti_numbers_unit_test.cpp )
target_link_libraries(Persistent_cohomology_test_betti_numbers ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(Persistent_cohomology_test_unit ${TBB_LIBRARIES})
  target_link_libraries(Persistent_cohomology_test_betti_numbers ${TBB_LIBRARIES})
endif(TBB_FOUND)

# Do not forget to copy test results files in current binary dir
file(COPY "${CMAKE_SOURCE_DIR}/src/Persistent_cohomology/test/simplex_tree_file_for_unit_test.txt"
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

# Unitary tests
gudhi_add_coverage_test(Persistent_cohomology_test_unit)
gudhi_add_coverage_test(Persistent_cohomology_test_betti_numbers)

if(GMPXX_FOUND AND GMP_FOUND)
  add_executable ( Persistent_cohomology_test_unit_multi_field persistent_cohomology_unit_test_multi_field.cpp )
  target_link_libraries(Persistent_cohomology_test_unit_multi_field ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
  if (TBB_FOUND)
    target_link_libraries(Persistent_cohomology_test_unit_multi_field ${TBB_LIBRARIES})
  endif(TBB_FOUND)

  # Do not forget to copy test results files in current binary dir
  file(COPY "${CMAKE_SOURCE_DIR}/src/Persistent_cohomology/test/simplex_tree_file_for_multi_field_unit_test.txt"
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  # Unitary tests
  gudhi_add_coverage_test(Persistent_cohomology_test_unit_multi_field)

endif(GMPXX_FOUND AND GMP_FOUND)

