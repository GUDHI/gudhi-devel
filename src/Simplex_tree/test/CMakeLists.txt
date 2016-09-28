cmake_minimum_required(VERSION 2.6)
project(Simplex_tree_tests)

if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()
if (GPROF_PATH)
  # for gprof to make coverage reports - Jenkins
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

add_executable ( SimplexTreeUT simplex_tree_unit_test.cpp )
target_link_libraries(SimplexTreeUT ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(SimplexTreeUT ${TBB_LIBRARIES})
endif()

# Do not forget to copy test files in current binary dir
file(COPY "simplex_tree_for_unit_test.txt" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

# Unitary tests
add_test(NAME SimplexTreeUT
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/SimplexTreeUT 
    # XML format for Jenkins xUnit plugin 
    --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/SimplexTreeUT.xml --log_level=test_suite --report_level=no)

