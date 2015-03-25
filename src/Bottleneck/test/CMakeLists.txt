cmake_minimum_required(VERSION 2.6)
project(GUDHIBottleneckUnitTest)

if(NOT MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} --coverage")	
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} --coverage")
endif()

add_executable ( BottleneckUnitTest bottleneck_unit_test.cpp )
target_link_libraries(BottleneckUnitTest ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

# Unitary tests
add_test(BottleneckUnitTest ${CMAKE_CURRENT_BINARY_DIR}/BottleneckUnitTest)

if (LCOV_PATH)
  # Lcov code coverage of unitary test
  add_test(src/Bottleneck/lcov/coverage.log ${CMAKE_SOURCE_DIR}/scripts/check_code_coverage.sh ${CMAKE_SOURCE_DIR}/src/Bottleneck)
endif()

cpplint_add_tests("${CMAKE_SOURCE_DIR}/src/Bottleneck/include/gudhi")
