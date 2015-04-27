cmake_minimum_required(VERSION 2.6)
project(GUDHIAlphaShapesUT)

# need CGAL 4.6
# cmake -DCGAL_DIR=~/workspace/CGAL-4.6-beta1 ../../..
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.6.0)
    message(STATUS "CGAL version: ${CGAL_VERSION}.")

    include( ${CGAL_USE_FILE} )
 
    find_package(Eigen3 3.1.0)
    if (EIGEN3_FOUND)
      message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
      include( ${EIGEN3_USE_FILE} )
      include_directories (BEFORE "../../include")

      if (GCOVR_PATH)
        # for gcovr to make coverage reports - Corbera Jenkins plugin
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fprofile-arcs -ftest-coverage")
      endif()
      if (GPROF_PATH)
        # for gprof to make coverage reports - Jenkins
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pg")
      endif()

      # uncomment to display debug traces
      # add_definitions(-DDEBUG_TRACES)
      add_executable ( AlphaShapesUT Alpha_shapes_unit_test.cpp )
      target_link_libraries(AlphaShapesUT  ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
      add_test(NAME AlphaShapesUT 
        COMMAND ${CMAKE_CURRENT_BINARY_DIR}/AlphaShapesUT
        # XML format for Jenkins xUnit plugin 
        --log_format=XML --log_sink=${CMAKE_SOURCE_DIR}/AlphaShapesUT.xml --log_level=test_suite --report_level=no)


    else()
      message(WARNING "Eigen3 not found. Version 3.1.0 is required for Alpha shapes feature.")
    endif()
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Alpha shapes feature. Version 4.6.0 is required.")
  endif ()
endif()

