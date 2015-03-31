cmake_minimum_required(VERSION 2.6)
project(GUDHISimplexTreeUT)

if(NOT MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} --coverage")	
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} --coverage")
endif()

add_executable ( SimplexTreeUT simplex_tree_unit_test.cpp )
target_link_libraries(SimplexTreeUT ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# Unitary tests
add_test(NAME SimplexTreeUT
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/SimplexTreeUT 
    # XML format for Jenkins xUnit plugin 
    --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/SimplexTreeUT.xml --log_level=test_suite --report_level=no)

if (LCOV_PATH)
  # Lcov code coverage of unitary test
  add_test(src/Simplex_tree/lcov/coverage.log ${CMAKE_SOURCE_DIR}/scripts/check_code_coverage.sh ${CMAKE_SOURCE_DIR}/src/Simplex_tree)
endif()

cpplint_add_tests("${CMAKE_SOURCE_DIR}/src/Simplex_tree/include/gudhi")
