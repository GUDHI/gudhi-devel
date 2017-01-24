cmake_minimum_required(VERSION 2.6)
project(Bottleneck_distance_tests)


if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()
if (GPROF_PATH)
  # for gprof to make coverage reports - Jenkins
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

# requires CGAL 4.8
# cmake -DCGAL_DIR=~/workspace/CGAL-4.8 ../../..
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.8.0)
    add_executable ( bottleneckUT bottleneck_unit_test.cpp )
    target_link_libraries(bottleneckUT ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

    # Unitary tests
    add_test(NAME bottleneckUT COMMAND ${CMAKE_CURRENT_BINARY_DIR}/bottleneckUT
      # XML format for Jenkins xUnit plugin
      --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/bottleneckUT.xml --log_level=test_suite --report_level=no)
  endif ()
endif()
