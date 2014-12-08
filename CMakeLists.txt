cmake_minimum_required(VERSION 2.6)
project(GUDHI)

#find_package(Boost 1.45.0 COMPONENTS system filesystem unit_test_framework chrono timer)
find_package(Boost REQUIRED COMPONENTS system filesystem unit_test_framework chrono timer REQUIRED)

set(CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/src/cmake/modules/")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/src/cmake/modules/")
message("CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}")
message("CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")

enable_testing()

if(MSVC)
  # Turn off some VC++ warnings
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4668 /wd4311 /wd4800 /wd4820 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -std=c++11 -Wall -Wpedantic -lboost_system -lboost_unit_test_framework -lgmpxx -lgmp")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")	
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
endif()

# BOOST ISSUE result_of vs C++11
add_definitions(-DBOOST_RESULT_OF_USE_DECLTYPE)

find_package(Boost)
find_package(GMP)
find_package(GMPXX)
find_package(CGAL)

# Required programs for unitary tests purpose
FIND_PROGRAM( LCOV_PATH lcov )
if (LCOV_PATH)
  message("lcov found in ${LCOV_PATH}")
endif()

FIND_PROGRAM( PYTHON_PATH python )
if (PYTHON_PATH)
  message("python found in ${PYTHON_PATH}")
endif()


if(NOT Boost_FOUND)
  message(FATAL_ERROR "NOTICE: This demo requires Boost and will not be compiled.")
else()
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
  include_directories(src/common/include/)
  include_directories(src/Contraction/include/)
  include_directories(src/Hasse_complex/include/)
  include_directories(src/Persistent_cohomology/include/)
  include_directories(src/Simplex_tree/include/)
  include_directories(src/Skeleton_blocker/include/)

  add_subdirectory(src/Simplex_tree/test)
  add_subdirectory(src/Simplex_tree/example)
  add_subdirectory(src/Persistent_cohomology/example)
  add_subdirectory(src/Skeleton_blocker/example) 
  add_subdirectory(src/Contraction/example)
endif()  
