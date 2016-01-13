cmake_minimum_required(VERSION 2.6)
project(GUDHITestWitnessComplex)

#if(NOT MSVC)
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
#  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} --coverage")	
#  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} --coverage")
#endif()

add_executable ( simple_witness_complex simple_witness_complex.cpp )
add_test(test ${CMAKE_CURRENT_BINARY_DIR}/simple_witness_complex)

add_executable ( witness_complex_points witness_complex_points.cpp )
add_test(test ${CMAKE_CURRENT_BINARY_DIR}/witness_complex_points)

#cpplint_add_tests("${CMAKE_SOURCE_DIR}/src/Simplex_tree/include/gudhi")
