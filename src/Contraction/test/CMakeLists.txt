cmake_minimum_required(VERSION 2.6)
project(GUDHIskbl)

add_executable(TestContraction TestContraction.cpp)

add_test(TestContraction ${CMAKE_CURRENT_BINARY_DIR}/TestContraction ${CMAKE_SOURCE_DIR}/data/meshes/SO3_50000.off 0.2)
