cmake_minimum_required(VERSION 2.6)
project(Tangential_complex_benchmark)

if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()
if (GPROF_PATH)
  # for gprof to make coverage reports - Jenkins
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

# need CGAL 4.8
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.8.0)
    message(STATUS "CGAL version: ${CGAL_VERSION}.")

    find_package(Eigen3 3.1.0)
    if (EIGEN3_FOUND)
      message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
      include( ${EIGEN3_USE_FILE} )

      add_executable(Tangential_complex_benchmark benchmark_tc.cpp)
      target_link_libraries(Tangential_complex_benchmark 
		${Boost_DATE_TIME_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
      if (TBB_FOUND)
        target_link_libraries(Tangential_complex_benchmark ${TBB_LIBRARIES})
      endif()

      # Do not forget to copy test files in current binary dir
      #file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

    else()
      message(WARNING "Eigen3 not found. Version 3.1.0 is required for Tangential complex feature.")
    endif()
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Tangential complex feature. Version 4.8.0 is required.")
  endif ()
endif()

