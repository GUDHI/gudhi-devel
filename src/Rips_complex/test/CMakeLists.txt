cmake_minimum_required(VERSION 2.6)
project(Rips_complex_tests)

if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()
if (GPROF_PATH)
  # for gprof to make coverage reports - Jenkins
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

add_executable ( rips_complex_UT test_rips_complex.cpp )
target_link_libraries(rips_complex_UT  ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(rips_complex_UT ${TBB_LIBRARIES})
endif()

# Do not forget to copy test files in current binary dir
file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

add_test(rips_complex_UT ${CMAKE_CURRENT_BINARY_DIR}/rips_complex_UT
         # XML format for Jenkins xUnit plugin 
         --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/rips_complex_UT.xml --log_level=test_suite --report_level=no)
