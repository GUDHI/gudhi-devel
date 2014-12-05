cmake_minimum_required(VERSION 2.6)
project(GUDHIskbl)
add_executable ( RipsContraction Rips_contraction.cpp)
target_link_libraries( RipsContraction ${Boost_TIMER_LIBRARY})	


