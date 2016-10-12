cmake_minimum_required(VERSION 2.6)
project(Tangential_complex_tests)

if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()
if (GPROF_PATH)
  # for gprof to make coverage reports - Jenkins
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.8.0)
    if (EIGEN3_FOUND)
      message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
      include( ${EIGEN3_USE_FILE} )
      include_directories (BEFORE "../../include")

      add_executable( Tangential_complex_test_TC test_tangential_complex.cpp )
      target_link_libraries(Tangential_complex_test_TC ${CGAL_LIBRARY} ${Boost_DATE_TIME_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
      if (TBB_FOUND)
        target_link_libraries(Tangential_complex_test_TC ${TBB_LIBRARIES})
      endif()
      add_test(Tangential_complex_test_TC ${CMAKE_CURRENT_BINARY_DIR}/Tangential_complex_test_TC
      # XML format for Jenkins xUnit plugin 
      --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/Tangential_complex_UT.xml --log_level=test_suite --report_level=no)

    endif()
  endif ()
endif()