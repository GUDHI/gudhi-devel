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
    if (EIGEN3_FOUND)
      add_executable(Tangential_complex_benchmark benchmark_tc.cpp)
      target_link_libraries(Tangential_complex_benchmark 
        ${Boost_DATE_TIME_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
      if (TBB_FOUND)
        target_link_libraries(Tangential_complex_benchmark ${TBB_LIBRARIES})
      endif(TBB_FOUND)
    endif(EIGEN3_FOUND)
  endif (NOT CGAL_VERSION VERSION_LESS 4.8.0)
endif(CGAL_FOUND)

