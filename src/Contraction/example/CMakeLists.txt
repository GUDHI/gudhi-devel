cmake_minimum_required(VERSION 2.6)
project(GUDHIskbl)

add_executable(RipsContraction Rips_contraction.cpp)
target_link_libraries(RipsContraction ${Boost_TIMER_LIBRARY} ${Boost_SYSTEM_LIBRARY})	

add_test(RipsContraction ${CMAKE_CURRENT_BINARY_DIR}/RipsContraction ${CMAKE_SOURCE_DIR}/data/meshes/sphere3D_2646.off 0.2)
