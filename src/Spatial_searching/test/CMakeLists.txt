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
    message(STATUS "CGAL version: ${CGAL_VERSION}.")

    if (EIGEN3_FOUND)

      add_executable( Spatial_searching_test_Kd_tree_search test_Kd_tree_search.cpp )
      target_link_libraries(Spatial_searching_test_Kd_tree_search 
        ${CGAL_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    else()
      message(WARNING "Eigen3 not found. Version 3.1.0 is required for the Spatial_searching tests.")
    endif()
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Subsampling tests. Version 4.8.1 is required.")
  endif ()
else()
  message(WARNING "CGAL not found. It is required for the Spatial_searching tests.")
endif()
