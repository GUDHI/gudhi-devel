cmake_minimum_required(VERSION 2.6)
project(Spatial_searching_tests)

if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()
if (GPROF_PATH)
  # for gprof to make coverage reports - Jenkins
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.8.1)
    if (EIGEN3_FOUND)
      add_executable( Spatial_searching_test_Kd_tree_search test_Kd_tree_search.cpp )
      target_link_libraries(Spatial_searching_test_Kd_tree_search 
        ${CGAL_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

      add_test(Spatial_searching_test_Kd_tree_search ${CMAKE_CURRENT_BINARY_DIR}/Spatial_searching_test_Kd_tree_search
      # XML format for Jenkins xUnit plugin 
      --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/Spatial_searching_UT.xml --log_level=test_suite --report_level=no)
    endif()
  endif ()
endif()
