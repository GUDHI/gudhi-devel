cmake_minimum_required(VERSION 2.6)
project(GUDHIskbl)

if(NOT MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} --coverage")	
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} --coverage")
endif()

add_executable(TestContraction TestContraction.cpp)

# Unitary tests
add_test(TestContraction ${CMAKE_CURRENT_BINARY_DIR}/TestContraction ${CMAKE_SOURCE_DIR}/data/meshes/SO3_50000.off 0.2)

if (LCOV_PATH)
  # Lcov code coverage of unitary test
  add_test(src/Contraction/lcov/coverage.log ${CMAKE_SOURCE_DIR}/scripts/check_code_coverage.sh ${CMAKE_SOURCE_DIR}/src/Contraction)
endif()

cpplint_add_tests("${CMAKE_SOURCE_DIR}/src/Contraction/include/gudhi")
