
if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()
if (GPROF_PATH)
  # for gprof to make coverage reports - Jenkins
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

if (DEBUG_TRACES)
  # Make CTest more verbose with DEBUG_TRACES - no XML output
  set(GUDHI_UT_LOG_LEVEL "--log_level=all")
  set(GUDHI_UT_REPORT_LEVEL "--report_level=detailed")
else()
  set(GUDHI_UT_LOG_FORMAT "--log_format=XML")
  set(GUDHI_UT_LOG_SINK "--log_sink=${CMAKE_BINARY_DIR}/${unitary_test}_UT.xml")
  set(GUDHI_UT_LOG_LEVEL "--log_level=test_suite")
  set(GUDHI_UT_REPORT_LEVEL "--report_level=no")
endif()

function(gudhi_add_coverage_test unitary_test)
  add_test(NAME ${unitary_test} COMMAND $<TARGET_FILE:${unitary_test}>
      ${GUDHI_UT_LOG_FORMAT} ${GUDHI_UT_LOG_SINK}
      ${GUDHI_UT_LOG_LEVEL} ${GUDHI_UT_REPORT_LEVEL})
endfunction()
