cmake_minimum_required(VERSION 2.6)
project(Common_tests)

if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()
if (GPROF_PATH)
  # for gprof to make coverage reports - Jenkins
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

add_executable ( poffreader_UT test_points_off_reader.cpp )
target_link_libraries(poffreader_UT ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# Do not forget to copy test files in current binary dir
file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

# Unitary tests
add_test(poffreader_UT ${CMAKE_CURRENT_BINARY_DIR}/poffreader_UT
         # XML format for Jenkins xUnit plugin 
         --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/poffreader_UT.xml --log_level=test_suite --report_level=no)


