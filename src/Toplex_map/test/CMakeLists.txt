cmake_minimum_required(VERSION 2.6)
project(Toplex_map_tests)

add_executable(chrono chrono.cpp)
add_executable ( ToplexMapUT toplex_map_unit_test.cpp )
target_link_libraries(ToplexMapUT ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})


# Unitary tests
add_test(NAME ToplexMapUT 
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/ToplexMapUT
    ${CMAKE_SOURCE_DIR}/src/Toplex_map/test/test.txt
    # XML format for Jenkins xUnit plugin 
    --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/ToplexMapUT.xml --log_level=test_suite --report_level=no)

