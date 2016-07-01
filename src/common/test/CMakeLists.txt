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

# need CGAL 4.7
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.7.0)
    find_package(Eigen3 3.1.0)
    if (EIGEN3_FOUND)
      message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
      include( ${EIGEN3_USE_FILE} )
      
      add_executable ( poffreader_UT points_off_reader_unit_test.cpp )
      target_link_libraries(poffreader_UT ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

      # Do not forget to copy test files in current binary dir
      file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

      # Unitary tests
      add_test(poffreader_UT ${CMAKE_CURRENT_BINARY_DIR}/poffreader_UT
      # XML format for Jenkins xUnit plugin 
      --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/poffreader_UT.xml --log_level=test_suite --report_level=no)

    else()
      message(WARNING "Eigen3 not found. Version 3.1.0 is required for Alpha shapes feature.")
    endif()
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Alpha shapes feature. Version 4.6.0 is required.")
  endif ()
endif()

