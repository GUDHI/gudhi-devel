cmake_minimum_required(VERSION 2.6)
project(GUDHI)


find_package(Boost REQUIRED COMPONENTS system filesystem unit_test_framework chrono timer REQUIRED)

if(MSVC)
	# Turn off some VC++ warnings
	SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4668 /wd4311 /wd4800 /wd4820 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -std=c++11 -Wall -Wpedantic -lboost_system -lgmpxx -lgmp")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")	
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
endif()


find_package(Boost) 
if(NOT Boost_FOUND)
  message(FATAL_ERROR "NOTICE: This demo requires Boost and will not be compiled.")
else()
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})    
    LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})                      
	include_directories(include/)
	
	add_subdirectory(example/Simplex_tree)                   
	add_subdirectory(example/Persistent_cohomology)                   
	add_subdirectory(example/Skeleton_blocker)                   
	add_subdirectory(example/Contraction)                   
endif()  
