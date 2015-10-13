cmake_minimum_required(VERSION 2.6)
project(GUDHI)

include("CMakeGUDHIVersion.txt")

enable_testing()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules/")

find_package(Boost REQUIRED COMPONENTS system filesystem program_options chrono timer REQUIRED)

if (NOT CMAKE_BUILD_TYPE)
 set(CMAKE_BUILD_TYPE "Release")
endif()
if(MSVC)
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4668 /wd4311 /wd4800 /wd4820 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)

find_package(GMP)
if(GMP_FOUND)
  find_package(GMPXX)
endif()

find_package(CGAL)

if(NOT Boost_FOUND)
  message(FATAL_ERROR "NOTICE: This demo requires Boost and will not be compiled.")
else()
  # BOOST ISSUE result_of vs C++11
  add_definitions(-DBOOST_RESULT_OF_USE_DECLTYPE)
  # BOOST ISSUE with Libraries name resolution under Windows
  add_definitions(-DBOOST_ALL_NO_LIB)

  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

  #---------------------------------------------------------------------------------------
  # Gudhi compilation part
  include_directories(include)

  add_subdirectory(example/common)
  add_subdirectory(example/Simplex_tree)                   
  add_subdirectory(example/Persistent_cohomology)                   
  add_subdirectory(example/Skeleton_blocker)                   
  add_subdirectory(example/Contraction)                   
  add_subdirectory(example/Alpha_complex)

  # data points generator
  add_subdirectory(data/points/generator)

  # Please let GudhUI in last compilation position as QT is known to modify CMAKE_CXX_FLAGS
  # GudhUI
  add_subdirectory(GudhUI)
  #---------------------------------------------------------------------------------------

  #---------------------------------------------------------------------------------------
  # GUDHIConfig.cmake
  # Export the package for use from the build-tree
  # (this registers the build-tree with a global CMake-registry)
  export(PACKAGE GUDHI)

  # Create the GUDHIConfig.cmake and GUDHIConfigVersion files
  set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include")
  configure_file(GUDHIConfig.cmake.in "${PROJECT_BINARY_DIR}/GUDHIConfig.cmake" @ONLY)
  configure_file(GUDHIConfigVersion.cmake.in "${PROJECT_BINARY_DIR}/GUDHIConfigVersion.cmake" @ONLY)

  #---------------------------------------------------------------------------------------

  #---------------------------------------------------------------------------------------
  # Gudhi installation part

  # Install the GUDHIConfig.cmake and GUDHIConfigVersion.cmake
  install(FILES
    "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/GUDHIConfig.cmake"
    "${PROJECT_BINARY_DIR}/GUDHIConfigVersion.cmake"
    DESTINATION share/gudhi)

  # install the include file on "make install"
  install(DIRECTORY include/gudhi DESTINATION include)
  #---------------------------------------------------------------------------------------

endif()  
