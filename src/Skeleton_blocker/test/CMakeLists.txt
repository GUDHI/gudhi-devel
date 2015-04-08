cmake_minimum_required(VERSION 2.6)
project(GUDHIskbl)

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

add_executable(TestSkeletonBlockerComplex TestSkeletonBlockerComplex.cpp)
add_executable(TestSimplifiable TestSimplifiable.cpp)
add_executable(TestGeometricComplex TestGeometricComplex.cpp)

add_test(TestSkeletonBlockerComplex ${CMAKE_CURRENT_BINARY_DIR}/TestSkeletonBlockerComplex)
add_test(TestSimplifiable ${CMAKE_CURRENT_BINARY_DIR}/TestSimplifiable)
add_test(TestGeometricComplex ${CMAKE_CURRENT_BINARY_DIR}/TestGeometricComplex)

