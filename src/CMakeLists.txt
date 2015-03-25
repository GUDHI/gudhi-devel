cmake_minimum_required(VERSION 2.6)
project(GUDHI)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules/")

find_package(Boost REQUIRED COMPONENTS system filesystem program_options chrono timer REQUIRED)


if (NOT CMAKE_BUILD_TYPE)
 set(CMAKE_BUILD_TYPE "Release")
endif()
if(MSVC)
	SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4668 /wd4311 /wd4800 /wd4820 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")
else()
  list(APPEND CMAKE_CXX_FLAGS "-std=c++11")
endif() 

# BOOST ISSUE result_of vs C++11
add_definitions(-DBOOST_RESULT_OF_USE_DECLTYPE)

find_package(Boost) 
find_package(GMP)
if(GMP_FOUND)
        find_package(GMPXX)
endif()
find_package(CGAL)

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
  add_subdirectory(example/Hasse_complex)
  add_subdirectort(example/Witness_complex)
  add_subdirectory(example/Alpha_shapes)
  add_subdirectory(example/Bottleneck)

endif()  
