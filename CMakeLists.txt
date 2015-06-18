cmake_minimum_required(VERSION 2.6)
project(GUDHIdev)

find_package(Boost REQUIRED COMPONENTS system filesystem unit_test_framework chrono timer program_options thread REQUIRED)

set(CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/src/cmake/modules/")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/src/cmake/modules/")
message("CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}")
message("CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")

enable_testing()

if(MSVC)
  # Turn off some VC++ warnings
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4668 /wd4311 /wd4800 /wd4820 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -std=c++11 -Wall -Wpedantic -Wsign-compare")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb -O0")	
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
endif()

set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)

find_package(Boost)
find_package(GMP)
if(GMP_FOUND)
	find_package(GMPXX)
endif()

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

FIND_PROGRAM( DIFF_PATH diff )
if (DIFF_PATH)
  message("diff found in ${DIFF_PATH}")
endif()

# Function to add_test cpplint on each header file of the Gudhi module 
function(cpplint_add_tests the_directory)
  if (PYTHON_PATH)
    # Cpplint tests on coding style
    file(GLOB files "${the_directory}/*.h" "${the_directory}/*/*.h")
    foreach(filename ${files})
      message(${filename})
      add_test("${filename}.cpplint" ${CMAKE_SOURCE_DIR}/scripts/check_google_style.sh ${filename}  ${CMAKE_SOURCE_DIR}/scripts/cpplint.py)
    endforeach()
  endif()
endfunction(cpplint_add_tests)


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
  add_subdirectory(src/Hasse_complex/example)
  add_subdirectory(src/Alpha_complex/example)
  add_subdirectory(src/Alpha_complex/test)
  add_subdirectory(src/Bottleneck/example)
  add_subdirectory(src/Bottleneck/test)

  # data points generator
  add_subdirectory(data/points/generator)

endif()  


