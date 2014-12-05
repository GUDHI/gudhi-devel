cmake_minimum_required(VERSION 2.6)
project(GUDHITestSimplexTree)

# NEED TO FIND BOOST NEEDED COMPONENTS TO LINK WITH
find_package(Boost 1.45.0 COMPONENTS system unit_test_framework)
message("Boost_lib = ${Boost_LIBRARIES}")

include_directories(${Boost_INCLUDE_DIRS})
add_executable ( TestSimplexTree UnitTestSimplexTree.cpp )
target_link_libraries(TestSimplexTree ${Boost_LIBRARIES})

