cmake_minimum_required(VERSION 2.6)
project(GUDHITestSimplexTree)

if(NOT MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} --coverage")	
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} --coverage")
endif()

add_executable ( persistent_cohomology_unit_test persistent_cohomology_unit_test.cpp )
target_link_libraries(persistent_cohomology_unit_test ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# Unitary tests
add_test(persistent_cohomology_unit_test ${CMAKE_CURRENT_BINARY_DIR}/persistent_cohomology_unit_test ${CMAKE_SOURCE_DIR}/src/Persistent_cohomology/test/simplex_tree_file_for_unit_test.txt)

if(GMPXX_FOUND AND GMP_FOUND)
   add_executable ( persistent_cohomology_unit_test_multi_field persistent_cohomology_unit_test_multi_field.cpp )
   target_link_libraries(persistent_cohomology_unit_test_multi_field ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})

   # Unitary tests
   add_test(persistent_cohomology_unit_test_multi_field ${CMAKE_CURRENT_BINARY_DIR}/persistent_cohomology_unit_test_multi_field ${CMAKE_SOURCE_DIR}/src/Persistent_cohomology/test/simplex_tree_file_for_unit_test_multi_field.txt)
endif()

if (LCOV_PATH)
  # Lcov code coverage of unitary test
  add_test(src/Persistent_cohomology/lcov/coverage.log ${CMAKE_SOURCE_DIR}/scripts/check_code_coverage.sh ${CMAKE_SOURCE_DIR}/src/Persistent_cohomology)
endif()

cpplint_add_tests("${CMAKE_SOURCE_DIR}/src/Persistent_cohomology/include/gudhi")
