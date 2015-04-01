cmake_minimum_required(VERSION 2.6)
project(GUDHIPersistentCohomologyUT)

if(NOT MSVC)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")	
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fprofile-arcs -ftest-coverage")
endif()

add_executable ( PersistentCohomologyUT persistent_cohomology_unit_test.cpp )
target_link_libraries(PersistentCohomologyUT ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# Unitary tests
add_test(NAME PersistentCohomologyUT 
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/PersistentCohomologyUT
    ${CMAKE_SOURCE_DIR}/src/Persistent_cohomology/test/simplex_tree_file_for_unit_test.txt
    # XML format for Jenkins xUnit plugin 
    --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/PersistentCohomologyUT.xml --log_level=test_suite --report_level=no)

if(GMPXX_FOUND AND GMP_FOUND)
  add_executable ( PersistentCohomologyMultiFieldUT persistent_cohomology_unit_test_multi_field.cpp )
  target_link_libraries(PersistentCohomologyMultiFieldUT ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})

  # Unitary tests
  add_test(NAME PersistentCohomologyMultiFieldUT
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/PersistentCohomologyMultiFieldUT
      ${CMAKE_SOURCE_DIR}/src/Persistent_cohomology/test/simplex_tree_file_for_multi_field_unit_test.txt
      # XML format for Jenkins xUnit plugin 
      --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/PersistentCohomologyMultiFieldUT.xml --log_level=test_suite --report_level=no)

endif()

cpplint_add_tests("${CMAKE_SOURCE_DIR}/src/Persistent_cohomology/include/gudhi")
