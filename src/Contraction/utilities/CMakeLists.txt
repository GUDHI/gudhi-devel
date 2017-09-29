cmake_minimum_required(VERSION 2.6)
project(Contraction_utilities)

add_executable(GarlandHeckbert Garland_heckbert.cpp)

install(TARGETS GarlandHeckbert DESTINATION bin)
