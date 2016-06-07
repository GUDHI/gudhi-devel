cmake_minimum_required(VERSION 2.6)
project(Persistent_cohomology_tests)

if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()
if (GPROF_PATH)
  # for gprof to make coverage reports - Jenkins
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

add_executable ( PersistentCohomologyUT persistent_cohomology_unit_test.cpp )
target_link_libraries(PersistentCohomologyUT ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_executable ( BettiNumbersUT betti_numbers_unit_test.cpp )
target_link_libraries(BettiNumbersUT ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(PersistentCohomologyUT ${TBB_RELEASE_LIBRARY})
  target_link_libraries(BettiNumbersUT ${TBB_RELEASE_LIBRARY})
endif()

# Unitary tests
add_test(NAME PersistentCohomologyUT 
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/PersistentCohomologyUT
    ${CMAKE_SOURCE_DIR}/src/Persistent_cohomology/test/simplex_tree_file_for_unit_test.txt
    # XML format for Jenkins xUnit plugin 
    --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/PersistentCohomologyUT.xml --log_level=test_suite --report_level=no)

add_test(NAME BettiNumbersUT 
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/BettiNumbersUT
    # XML format for Jenkins xUnit plugin 
    --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/BettiNumbersUT.xml --log_level=test_suite --report_level=no)

if(GMPXX_FOUND AND GMP_FOUND)
  add_executable ( PersistentCohomologyMultiFieldUT persistent_cohomology_unit_test_multi_field.cpp )
target_link_libraries(PersistentCohomologyMultiFieldUT ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
if (TBB_FOUND)
  target_link_libraries(PersistentCohomologyMultiFieldUT ${TBB_RELEASE_LIBRARY})
endif()

  # Unitary tests
  add_test(NAME PersistentCohomologyMultiFieldUT
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/PersistentCohomologyMultiFieldUT
      ${CMAKE_SOURCE_DIR}/src/Persistent_cohomology/test/simplex_tree_file_for_multi_field_unit_test.txt
      # XML format for Jenkins xUnit plugin 
      --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/PersistentCohomologyMultiFieldUT.xml --log_level=test_suite --report_level=no)

endif()

