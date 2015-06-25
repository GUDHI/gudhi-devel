cmake_minimum_required(VERSION 2.6)
project(GUDHIBottleneckUT)

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
  message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
  message("CMAKE_CXX_FLAGS_DEBUG = ${CMAKE_CXX_FLAGS_DEBUG}")	
  message("CMAKE_CXX_FLAGS_RELEASE = ${CMAKE_CXX_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -std=c++11 -Wall -Wpedantic -Wsign-compare")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb -O0")	
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

add_executable ( BottleneckUT bottleneck_unit_test.cpp )
target_link_libraries(BottleneckUT ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# Unitary tests
add_test(NAME BottleneckUT 
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/BottleneckUT
  # XML format for Jenkins xUnit plugin 
  --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/BottleneckUT.xml --log_level=test_suite --report_level=no)

