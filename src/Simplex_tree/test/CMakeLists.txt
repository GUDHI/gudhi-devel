cmake_minimum_required(VERSION 2.6)
project(GUDHISimplexTreeUT)

if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")	
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fprofile-arcs -ftest-coverage")
endif()

add_executable ( SimplexTreeUT simplex_tree_unit_test.cpp )
target_link_libraries(SimplexTreeUT ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# Unitary tests
add_test(NAME SimplexTreeUT
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/SimplexTreeUT 
    # XML format for Jenkins xUnit plugin 
    --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/SimplexTreeUT.xml --log_level=test_suite --report_level=no)

cpplint_add_tests("${CMAKE_SOURCE_DIR}/src/Simplex_tree/include/gudhi")
