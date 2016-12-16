cmake_minimum_required(VERSION 2.6)
project(Witness_complex_tests)

if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()
if (GPROF_PATH)
  # for gprof to make coverage reports - Jenkins
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

add_executable ( Witness_complex_test_simple_witness_complex test_simple_witness_complex.cpp )
target_link_libraries(Witness_complex_test_simple_witness_complex ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_link_libraries(Witness_complex_test_simple_witness_complex ${TBB_LIBRARIES})

# Unitary tests definition and xml result file generation
add_test(NAME simple_witness_complex
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/Witness_complex_test_simple_witness_complex 
    # XML format for Jenkins xUnit plugin 
    --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/Witness_complex_test_simple_witness_complexUT.xml --log_level=test_suite --report_level=no)
