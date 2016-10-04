cmake_minimum_required(VERSION 2.6)
project(Skeleton_blocker_tests)

if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()
if (GPROF_PATH)
  # for gprof to make coverage reports - Jenkins
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

add_executable ( test_skeleton_blocker_complex test_skeleton_blocker_complex.cpp )
target_link_libraries(test_skeleton_blocker_complex ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_executable ( test_skeleton_blocker_geometric_complex test_skeleton_blocker_geometric_complex.cpp )
target_link_libraries(test_skeleton_blocker_geometric_complex ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_executable ( test_skeleton_blocker_simplifiable test_skeleton_blocker_simplifiable.cpp )
target_link_libraries(test_skeleton_blocker_simplifiable ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# Do not forget to copy test files in current binary dir
file(COPY "test2.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

# Unitary tests
add_test(NAME test_skeleton_blocker_complex
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_skeleton_blocker_complex
    # XML format for Jenkins xUnit plugin 
    --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/skbl_complexUT.xml --log_level=test_suite --report_level=no)
add_test(NAME test_skeleton_blocker_geometric_complex
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_skeleton_blocker_geometric_complex
    # XML format for Jenkins xUnit plugin 
    --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/skbl_geometric_complexUT.xml --log_level=test_suite --report_level=no)
add_test(NAME test_skeleton_blocker_simplifiable
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_skeleton_blocker_simplifiable
    # XML format for Jenkins xUnit plugin 
    --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/skbl_simplifiableUT.xml --log_level=test_suite --report_level=no)
