
if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()
if (GPROF_PATH)
  # for gprof to make coverage reports - Jenkins
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

function(gudhi_add_coverage_test unitary_test)
  add_test(NAME ${unitary_test} COMMAND $<TARGET_FILE:${unitary_test}>
      "--log_format=XML" "--log_sink=${CMAKE_BINARY_DIR}/${unitary_test}_UT.xml" "--log_level=test_suite" "--report_level=no")
endfunction()
