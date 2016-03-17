cmake_minimum_required(VERSION 2.6)
project(GUDHI)

include("CMakeGUDHIVersion.txt")

enable_testing()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules/")

find_package(Boost REQUIRED COMPONENTS system filesystem program_options chrono timer REQUIRED)

if(NOT Boost_FOUND)
  message(FATAL_ERROR "NOTICE: This demo requires Boost and will not be compiled.")
else()

  find_package(GMP)
  if(GMP_FOUND)
    INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR})
    find_package(GMPXX)
    if(GMPXX_FOUND)
      INCLUDE_DIRECTORIES(${GMPXX_INCLUDE_DIR})
    endif()
  endif()

  # In CMakeLists.txt, when include(${CGAL_USE_FILE}), CMAKE_CXX_FLAGS are overwritten.
  # cf. http://doc.cgal.org/latest/Manual/installation.html#title40
  # A workaround is to include(${CGAL_USE_FILE}) before adding "-std=c++11".
  # A fix would be to use https://cmake.org/cmake/help/v3.1/prop_gbl/CMAKE_CXX_KNOWN_FEATURES.html 
  # or even better https://cmake.org/cmake/help/v3.1/variable/CMAKE_CXX_STANDARD.html
  # but it implies to use cmake version 3.1 at least.
  find_package(CGAL)
  if(CGAL_FOUND)
    include( ${CGAL_USE_FILE} )
  endif()
  
  if(MSVC)
    # Turn off some VC++ warnings
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4668 /wd4311 /wd4800 /wd4820 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wpedantic -Wsign-compare")
  endif()

  if(CMAKE_BUILD_TYPE MATCHES Debug)
    message("++ Debug compilation flags are: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
  else()
    message("++ Release compilation flags are: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
  endif()

  set(Boost_USE_STATIC_LIBS        ON)
  set(Boost_USE_MULTITHREADED      ON)
  set(Boost_USE_STATIC_RUNTIME    OFF)

  # Find TBB package for parallel sort - not mandatory, just optional.
  set(TBB_FIND_QUIETLY ON)
  find_package(TBB)

  # BOOST ISSUE result_of vs C++11
  add_definitions(-DBOOST_RESULT_OF_USE_DECLTYPE)
  # BOOST ISSUE with Libraries name resolution under Windows
  add_definitions(-DBOOST_ALL_NO_LIB)

  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

  #---------------------------------------------------------------------------------------
  # Gudhi compilation part
  include_directories(include)

  #add_subdirectory(example/Simplex_tree)                   
  #add_subdirectory(example/Persistent_cohomology)                   
  #add_subdirectory(example/Skeleton_blocker)                   
  #add_subdirectory(example/Contraction)                   
  add_subdirectory(example/Bitmap_cubical_complex)

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
    "${PROJECT_BINARY_DIR}/GUDHIConfig.cmake"
    "${PROJECT_BINARY_DIR}/GUDHIConfigVersion.cmake"
    DESTINATION share/gudhi)

  # install the include file on "make install"
  install(DIRECTORY include/gudhi DESTINATION include)
  #---------------------------------------------------------------------------------------

endif()  
