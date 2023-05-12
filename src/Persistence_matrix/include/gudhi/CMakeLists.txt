cmake_minimum_required(VERSION 3.5)

project(Matrix LANGUAGES CXX)

set(CMAKE_CXX_FLAGS_RELEASE "-Wall -g -O3 -march=native")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Boost 1.6 REQUIRED COMPONENTS system filesystem chrono timer REQUIRED)
set(Boost_CXX_FLAGS "-DBOOST_SYSTEM_NO_DEPRECATED")
set(CMAKE_REQUIRED_FLAGS ${Boost_CXX_FLAGS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Boost_CXX_FLAGS}")
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

include_directories(. chain_matrix boundary_matrix column_types)
add_executable(matrix main.cpp)

