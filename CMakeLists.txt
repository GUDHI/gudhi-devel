cmake_minimum_required(VERSION 2.6)
project(GUDHIdev)

include(CMakeGUDHIVersion.txt)

if (NOT CMAKE_BUILD_TYPE)
 # Set default build type to Release
 set(CMAKE_BUILD_TYPE "Release")
endif()

enable_testing()

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/src/cmake/modules/")
message("CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")

# Generate GUDHI official version file
configure_file(GUDHIVersion.cmake.in "${PROJECT_BINARY_DIR}/GUDHIVersion.cmake" @ONLY)

find_package(Boost REQUIRED COMPONENTS system filesystem unit_test_framework chrono timer program_options thread REQUIRED)

if(MSVC)
  # Turn off some VC++ warnings
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4668 /wd4311 /wd4800 /wd4820 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -std=c++11 -Wall -Wpedantic -Wsign-compare")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb -O1")	
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
endif()

set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)

find_package(GMP)
if(GMP_FOUND)
	find_package(GMPXX)
endif()

find_package(CGAL)

# Find TBB package for parallel sort - not mandatory, just optional.
set(TBB_FIND_QUIETLY ON)
find_package(TBB)

# Required programs for unitary tests purpose
FIND_PROGRAM( GCOVR_PATH gcovr )
if (GCOVR_PATH)
  message("gcovr found in ${GCOVR_PATH}")
endif()
# Required programs for unitary tests purpose
FIND_PROGRAM( GPROF_PATH gprof )
if (GPROF_PATH)
  message("gprof found in ${GPROF_PATH}")
endif()

FIND_PROGRAM( DIFF_PATH diff )
if (DIFF_PATH)
  message("diff found in ${DIFF_PATH}")
endif()


if(NOT Boost_FOUND)
  message(FATAL_ERROR "NOTICE: This demo requires Boost and will not be compiled.")
else()
  # BOOST ISSUE result_of vs C++11
  add_definitions(-DBOOST_RESULT_OF_USE_DECLTYPE)
  # BOOST ISSUE with Libraries name resolution under Windows
  add_definitions(-DBOOST_ALL_NO_LIB)

  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})                      

  message(STATUS "boost include dirs:" ${Boost_INCLUDE_DIRS})
  message(STATUS "boost library dirs:" ${Boost_LIBRARY_DIRS})

  if (DEBUG_TRACES)
    # For programs to be more verbose
    message(STATUS "DEBUG_TRACES are activated")
    add_definitions(-DDEBUG_TRACES)
  endif()

  include_directories(src/common/include/)
  include_directories(src/Alpha_complex/include/)
  include_directories(src/Bottleneck/include/)
  include_directories(src/Contraction/include/)
  include_directories(src/Hasse_complex/include/)
  include_directories(src/Persistent_cohomology/include/)
  include_directories(src/Simplex_tree/include/)
  include_directories(src/Skeleton_blocker/include/)

  add_subdirectory(src/common/example)
  add_subdirectory(src/common/test)
  add_subdirectory(src/Simplex_tree/test)
  add_subdirectory(src/Simplex_tree/example)
  add_subdirectory(src/Persistent_cohomology/test)
  add_subdirectory(src/Persistent_cohomology/example)
  add_subdirectory(src/Skeleton_blocker/test) 
  add_subdirectory(src/Skeleton_blocker/example) 
  add_subdirectory(src/Contraction/example)
  add_subdirectory(src/Alpha_complex/example)
  add_subdirectory(src/Alpha_complex/test)

  # data points generator
  add_subdirectory(data/points/generator)

  # Please let GudhUI in last compilation position as QT is known to modify CMAKE_CXX_FLAGS
  # GudhUI
  add_subdirectory(src/GudhUI)

endif()  


