cmake_minimum_required(VERSION 2.6)
project(GUDHIskbl)

if(NOT MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} --coverage")	
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} --coverage")
endif()

add_executable(TestSkeletonBlockerComplex TestSkeletonBlockerComplex.cpp)
add_executable(TestSimplifiable TestSimplifiable.cpp)

add_test(TestSkeletonBlockerComplex ${CMAKE_CURRENT_BINARY_DIR}/TestSkeletonBlockerComplex)
add_test(TestSimplifiable ${CMAKE_CURRENT_BINARY_DIR}/TestSimplifiable)

if (LCOV_PATH)
  # Lcov code coverage of unitary test
  add_test(src/Skeleton_blocker/lcov/coverage.log ${CMAKE_SOURCE_DIR}/scripts/check_code_coverage.sh ${CMAKE_SOURCE_DIR}/src/Skeleton_blocker)
endif()

cpplint_add_tests("${CMAKE_SOURCE_DIR}/src/Skeleton_blocker/include/gudhi")
