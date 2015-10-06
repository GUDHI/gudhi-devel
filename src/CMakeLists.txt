cmake_minimum_required(VERSION 2.6)
project(GUDHI)

include("CMakeGUDHIVersion.txt")

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
  # GUDHIConfig.cmake 1st part

  # Offer the user the choice of overriding the installation directories - i.e. "cmake -DGUDHI_DIR=~/gudhi/GUDHI_1.2.0 ."
  set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
  set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
  set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
  if(WIN32 AND NOT CYGWIN)
    set(DEF_INSTALL_CMAKE_DIR CMake)
  else()
    set(DEF_INSTALL_CMAKE_DIR lib/CMake/GUDHI)
  endif()
  set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Installation directory for CMake files")

  # Make relative paths absolute (needed later on)
  foreach(p LIB BIN INCLUDE CMAKE)
    set(var INSTALL_${p}_DIR)
    if(NOT IS_ABSOLUTE "${${var}}")
      set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
    endif()
  endforeach()
  #---------------------------------------------------------------------------------------

  #---------------------------------------------------------------------------------------
  # Gudhi compilation part
  include_directories(include)

  add_subdirectory(example/Simplex_tree)                   
  add_subdirectory(example/Persistent_cohomology)                   
  add_subdirectory(example/Skeleton_blocker)                   
  add_subdirectory(example/Contraction)                   
  add_subdirectory(example/Hasse_complex)
  add_subdirectory(example/Bottleneck)

  # GudhUI
  add_subdirectory(GudhUI)

  # data points generator
  add_subdirectory(data/points/generator)
  #---------------------------------------------------------------------------------------

  #---------------------------------------------------------------------------------------
  # GUDHIConfig.cmake 2nd part
  # Export the package for use from the build-tree
  # (this registers the build-tree with a global CMake-registry)
  export(PACKAGE GUDHI)

  # Create the GUDHIConfig.cmake and GUDHIConfigVersion files
  file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")
  # ... for the build tree
  set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include" "${PROJECT_BINARY_DIR}/include")
  configure_file(GUDHIConfig.cmake.in "${PROJECT_BINARY_DIR}/GUDHIConfig.cmake" @ONLY)
  # ... for the install tree
  set(CONF_INCLUDE_DIRS "\${GUDHI_CMAKE_DIR}/${REL_INCLUDE_DIR}")
  configure_file(GUDHIConfig.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/GUDHIConfig.cmake" @ONLY)
  # ... for both
  configure_file(GUDHIConfigVersion.cmake.in "${PROJECT_BINARY_DIR}/GUDHIConfigVersion.cmake" @ONLY)

  configure_file(GUDHIVersion.cmake.in "${PROJECT_BINARY_DIR}/GUDHIVersion.cmake" @ONLY)

  # Install the GUDHIConfig.cmake and GUDHIConfigVersion.cmake
  install(FILES
    "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/GUDHIConfig.cmake"
    "${PROJECT_BINARY_DIR}/GUDHIConfigVersion.cmake"
    DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
  #---------------------------------------------------------------------------------------

  #---------------------------------------------------------------------------------------
  # Gudhi installation part

  # install the include file on "make install"
  install(DIRECTORY include/gudhi DESTINATION include)
  #---------------------------------------------------------------------------------------

endif()  
