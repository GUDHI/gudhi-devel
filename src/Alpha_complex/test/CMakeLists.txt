cmake_minimum_required(VERSION 2.6)
project(Alpha_complex_tests)

if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()
if (GPROF_PATH)
  # for gprof to make coverage reports - Jenkins
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.7.0)
    if (EIGEN3_FOUND)
      add_executable ( AlphaComplexUT Alpha_complex_unit_test.cpp )
      target_link_libraries(AlphaComplexUT  ${Boost_SYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY} ${CGAL_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
      if (TBB_FOUND)
        target_link_libraries(AlphaComplexUT ${TBB_LIBRARIES})
      endif()

      # Do not forget to copy test files in current binary dir
      file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

      add_test(AlphaComplexUT ${CMAKE_CURRENT_BINARY_DIR}/AlphaComplexUT
      # XML format for Jenkins xUnit plugin 
      --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/AlphaComplexUT.xml --log_level=test_suite --report_level=no)

    endif(EIGEN3_FOUND)
  endif (NOT CGAL_VERSION VERSION_LESS 4.7.0)
endif(CGAL_FOUND)

