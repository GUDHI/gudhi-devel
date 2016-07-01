cmake_minimum_required(VERSION 2.6)
project(Skeleton_blocker_tests)

if (GCOVR_PATH)
  # for gcovr to make coverage reports - Corbera Jenkins plugin
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()
if (GPROF_PATH)
  # for gprof to make coverage reports - Jenkins
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

add_executable(TestSkeletonBlockerComplex TestSkeletonBlockerComplex.cpp)
add_executable(TestSimplifiable TestSimplifiable.cpp)
add_executable(TestGeometricComplex TestGeometricComplex.cpp)

# Do not forget to copy test files in current binary dir
file(COPY "test.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY "test2.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

add_test(TestSkeletonBlockerComplex ${CMAKE_CURRENT_BINARY_DIR}/TestSkeletonBlockerComplex)
add_test(TestSimplifiable ${CMAKE_CURRENT_BINARY_DIR}/TestSimplifiable)
add_test(TestGeometricComplex ${CMAKE_CURRENT_BINARY_DIR}/TestGeometricComplex)

