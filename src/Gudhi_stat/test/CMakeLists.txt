cmake_minimum_required(VERSION 2.6)
project(GUDHI_STAT_CCUT)

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







add_executable ( persistence_intervals_test persistence_intervals_test.cpp )
target_link_libraries(persistence_intervals_test ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# Unitary tests
add_test(NAME persistence_intervals_test
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/persistence_intervals_test	 
    # XML format for Jenkins xUnit plugin 
    --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/persistence_intervals_test.xml --log_level=test_suite --report_level=no)










add_executable ( vector_representation_test vector_representation_test.cpp )
target_link_libraries(vector_representation_test ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# Unitary tests
add_test(NAME vector_representation_test
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/vector_representation_test
    # XML format for Jenkins xUnit plugin 
    --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/vector_representation_test.xml --log_level=test_suite --report_level=no)




add_executable ( persistence_lanscapes_test persistence_lanscapes_test.cpp )
target_link_libraries(persistence_lanscapes_test ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# Unitary tests
add_test(NAME persistence_lanscapes_test
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/persistence_lanscapes_test
    # XML format for Jenkins xUnit plugin 
    --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/persistence_lanscapes_test.xml --log_level=test_suite --report_level=no)
