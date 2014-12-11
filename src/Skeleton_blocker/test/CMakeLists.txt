cmake_minimum_required(VERSION 2.6)
project(GUDHIskbl)

add_executable(TestSkeletonBlockerComplex TestSkeletonBlockerComplex.cpp)
add_executable(TestSimplifiable TestSimplifiable.cpp)

add_test(TestSkeletonBlockerComplex ${CMAKE_CURRENT_BINARY_DIR}/TestSkeletonBlockerComplex)
add_test(TestSimplifiable ${CMAKE_CURRENT_BINARY_DIR}/TestSimplifiable)
