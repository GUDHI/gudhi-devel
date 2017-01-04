cmake_minimum_required(VERSION 2.6)
project(GUDHIdev)

include(CMakeGUDHIVersion.txt)
# Generate GUDHI official version file
configure_file(GUDHIVersion.cmake.in "${CMAKE_SOURCE_DIR}/GUDHIVersion.cmake" @ONLY)

find_package(Boost REQUIRED COMPONENTS system filesystem unit_test_framework chrono timer date_time program_options thread REQUIRED)

if(NOT Boost_FOUND)
  message(FATAL_ERROR "NOTICE: This program requires Boost and will not be compiled.")
else()

  set(CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/src/cmake/modules/")
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/src/cmake/modules/")
  message("CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}")
  message("CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")

  enable_testing()

  # For "make user_version"
  include(${CMAKE_MODULE_PATH}/GUDHI_user_version_target.txt)
  # For "make doxygen"
  include(${CMAKE_MODULE_PATH}/GUDHI_doxygen_target.txt)

  find_package(GMP)
  if(GMP_FOUND)
    message(STATUS "GMP_LIBRARIES = ${GMP_LIBRARIES}")
    INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR})
    find_package(GMPXX)
    if(GMPXX_FOUND)
      message(STATUS "GMPXX_LIBRARIES = ${GMPXX_LIBRARIES}")
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
  # Only CGAL versions > 4.4 supports what Gudhi uses from CGAL
  if (CGAL_VERSION VERSION_LESS 4.4.0)
    message("CGAL version ${CGAL_VERSION} is considered too old to be used by Gudhi.")
    unset(CGAL_FOUND)
  endif()
  if(CGAL_FOUND)
    message(STATUS "CGAL version: ${CGAL_VERSION}.")
    include( ${CGAL_USE_FILE} )

    if (NOT CGAL_VERSION VERSION_LESS 4.8.0)
      include_directories(BEFORE "src/common/include/gudhi_patches")
    endif()
  endif()
  
  if(MSVC)
    # Turn off some VC++ warnings
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4668 /wd4311 /wd4800 /wd4820 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -pedantic")
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
  if (TBB_FOUND)
    include(${TBB_USE_FILE})
    message("TBB found in ${TBB_LIBRARY_DIRS}")
    add_definitions(-DGUDHI_USE_TBB)
  endif()

  find_package(Eigen3 3.1.0)
  if (EIGEN3_FOUND)
    message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
    include( ${EIGEN3_USE_FILE} )
    #include_directories (BEFORE "../../include")
  endif (EIGEN3_FOUND)

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
  include_directories(src/Bitmap_cubical_complex/include/)
  include_directories(src/Bottleneck_distance/include/)
  include_directories(src/Contraction/include/)
  include_directories(src/Hasse_complex/include/)
  include_directories(src/Persistent_cohomology/include/)
  include_directories(src/Simplex_tree/include/)
  include_directories(src/Skeleton_blocker/include/)
  include_directories(src/Spatial_searching/include/)
  include_directories(src/Subsampling/include/)
  include_directories(src/Tangential_complex/include/)
  include_directories(src/Witness_complex/include/)

  add_subdirectory(src/common/example)
  add_subdirectory(src/common/test)
  add_subdirectory(src/Simplex_tree/test)
  add_subdirectory(src/Simplex_tree/example)
  add_subdirectory(src/Persistent_cohomology/test)
  add_subdirectory(src/Persistent_cohomology/example)
  add_subdirectory(src/Skeleton_blocker/test) 
  add_subdirectory(src/Skeleton_blocker/example) 
  add_subdirectory(src/Contraction/example)
  add_subdirectory(src/Witness_complex/test)
  add_subdirectory(src/Witness_complex/example)
  add_subdirectory(src/Bitmap_cubical_complex/test)
  add_subdirectory(src/Bitmap_cubical_complex/example)
  add_subdirectory(src/Alpha_complex/example)
  add_subdirectory(src/Alpha_complex/test)
  add_subdirectory(src/Spatial_searching/example)
  add_subdirectory(src/Spatial_searching/test)
  add_subdirectory(src/Subsampling/example)
  add_subdirectory(src/Subsampling/test)
  add_subdirectory(src/Tangential_complex/example)
  add_subdirectory(src/Tangential_complex/test)
  add_subdirectory(src/Tangential_complex/benchmark)
  add_subdirectory(src/Bottleneck_distance/example)
  add_subdirectory(src/Bottleneck_distance/test)
  add_subdirectory(src/Bottleneck_distance/benchmark)

  # data points generator
  add_subdirectory(data/points/generator)

  # Please let GudhUI in last compilation position as QT is known to modify CMAKE_CXX_FLAGS
  # GudhUI
  add_subdirectory(src/GudhUI)

endif()  


