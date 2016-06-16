cmake_minimum_required(VERSION 2.6)
project(GUDHIskbl)


add_executable(RipsContraction Rips_contraction.cpp)
add_executable(GarlandHeckbert Garland_heckbert.cpp)

target_link_libraries(RipsContraction ${Boost_TIMER_LIBRARY} ${Boost_SYSTEM_LIBRARY})	
target_link_libraries(GarlandHeckbert ${Boost_TIMER_LIBRARY} ${Boost_SYSTEM_LIBRARY})	


add_test(RipsContraction.tore3D.0.2 ${CMAKE_CURRENT_BINARY_DIR}/RipsContraction ${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off 0.2)
# TODO(DS) : These tests are too long under Windows
#add_test(RipsContraction.sphere.0.2 ${CMAKE_CURRENT_BINARY_DIR}/RipsContraction ${CMAKE_SOURCE_DIR}/data/points/sphere3D_2646.off 0.2)
#add_test(RipsContraction.S0310000 ${CMAKE_CURRENT_BINARY_DIR}/RipsContraction ${CMAKE_SOURCE_DIR}/data/points/SO3_10000.off 0.3)
