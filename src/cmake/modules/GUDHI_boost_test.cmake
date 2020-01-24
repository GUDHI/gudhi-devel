if (WITH_GUDHI_BOOST_TEST_COVERAGE)
  # Make CTest output XML coverage report - WITH_GUDHI_BOOST_TEST_COVERAGE must be set - default is OFF
  if (GCOVR_PATH)
    # for gcovr to make coverage reports - Corbera Jenkins plugin
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
  endif()
  if (GPROF_PATH)
    # for gprof to make coverage reports - Jenkins
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
  endif()
  set(GUDHI_UT_LOG_FORMAT "--log_format=XML")
  set(GUDHI_UT_LOG_SINK "--log_sink=${CMAKE_BINARY_DIR}/${unitary_test}_UT.xml")
  set(GUDHI_UT_LOG_LEVEL "--log_level=test_suite")
  set(GUDHI_UT_REPORT_LEVEL "--report_level=no")
else (WITH_GUDHI_BOOST_TEST_COVERAGE)
  # Make CTest more verbose and color output
  set(GUDHI_UT_LOG_LEVEL "--color_output")
  set(GUDHI_UT_REPORT_LEVEL "--report_level=detailed")
endif(WITH_GUDHI_BOOST_TEST_COVERAGE)

function(gudhi_add_boost_test unitary_test)
  target_link_libraries(${unitary_test} Boost::unit_test_framework)
  add_test(NAME ${unitary_test} COMMAND $<TARGET_FILE:${unitary_test}>
      ${GUDHI_UT_LOG_FORMAT} ${GUDHI_UT_LOG_SINK}
      ${GUDHI_UT_LOG_LEVEL} ${GUDHI_UT_REPORT_LEVEL})
endfunction()
