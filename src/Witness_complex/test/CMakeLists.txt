cmake_minimum_required(VERSION 2.6)
project(GUDHIWitnessComplexUT)

if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")	
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fprofile-arcs -ftest-coverage")
endif()
if (GPROF_PATH)
  # for gprof to make coverage reports - Jenkins
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pg")
endif()

add_executable ( simple_witness_complexUT simple_witness_complex.cpp )
target_link_libraries(simple_witness_complexUT ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# Unitary tests definition and xml result file generation
add_test(NAME simple_witness_complexUT
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/simple_witness_complexUT 
    # XML format for Jenkins xUnit plugin 
    --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/simple_witness_complexUT.xml --log_level=test_suite --report_level=no)

add_executable ( witness_complex_pointsUT witness_complex_points.cpp )
target_link_libraries(witness_complex_pointsUT ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# Unitary tests definition and xml result file generation
add_test(NAME witness_complex_pointsUT
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/witness_complex_pointsUT 
    # XML format for Jenkins xUnit plugin 
    --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/witness_complex_pointsUT.xml --log_level=test_suite --report_level=no)

