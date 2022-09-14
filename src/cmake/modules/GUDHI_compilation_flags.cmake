# This files manage compilation flags required by GUDHI

include(TestCXXAcceptsFlag)

# add a compiler flag only if it is accepted
macro(add_cxx_compiler_flag _flag)
  string(REPLACE "-" "_" "/" _flag_var ${_flag})
  check_cxx_accepts_flag("${_flag}" CXX_COMPILER_${_flag_var}_OK)
  if(CXX_COMPILER_${_flag_var}_OK)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_flag}")
  endif()
endmacro()

set (CMAKE_CXX_STANDARD 17)

enable_testing()

if(MSVC)
  add_cxx_compiler_flag("/W3")
else()
  add_cxx_compiler_flag("-Wall")
  # Only for dev version
  if(PROJECT_NAME STREQUAL "GUDHIdev")
    add_cxx_compiler_flag("-pedantic")
  endif()
endif()

if (DEBUG_TRACES)
  # For programs to be more verbose
  message(STATUS "DEBUG_TRACES are activated")
  add_definitions(-DDEBUG_TRACES)
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug)
  message("++ Debug compilation flags are: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
else()
  message("++ Release compilation flags are: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

option(WITH_GUDHI_BOOST_TEST_COVERAGE "Report xml coverage files on boost tests" OFF)
