cmake_minimum_required(VERSION 2.6)
project(Bottleneck_examples)

add_executable ( RandomDiagrams random_diagrams.cpp )
add_test(RandomDiagrams ${CMAKE_CURRENT_BINARY_DIR}/RandomDiagrams)
