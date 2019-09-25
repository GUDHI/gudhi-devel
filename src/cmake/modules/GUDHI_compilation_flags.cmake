# This files manage compilation flags required by GUDHI

include(TestCXXAcceptsFlag)
include(CheckCXXSourceCompiles)

# add a compiler flag only if it is accepted
macro(add_cxx_compiler_flag _flag)
  string(REPLACE "-" "_" "/" _flag_var ${_flag})
  check_cxx_accepts_flag("${_flag}" CXX_COMPILER_${_flag_var}_OK)
  if(CXX_COMPILER_${_flag_var}_OK)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_flag}")
  endif()
endmacro()

function(can_cgal_use_cxx11_thread_local)
  # This is because of https://github.com/CGAL/cgal/blob/master/Installation/include/CGAL/tss.h
  # CGAL is using boost thread if thread_local is not ready (requires XCode 8 for Mac).
  # The test in https://github.com/CGAL/cgal/blob/master/Installation/include/CGAL/config.h
  #   #if __has_feature(cxx_thread_local) || \
  #       ( (__GNUC__ * 100 + __GNUC_MINOR__) >= 408 && __cplusplus >= 201103L ) || \
  #       ( _MSC_VER >= 1900 )
  #   #define CGAL_CAN_USE_CXX11_THREAD_LOCAL
  #   #endif
  set(CGAL_CAN_USE_CXX11_THREAD_LOCAL "
      int main() {
      #ifndef __has_feature
        #define __has_feature(x) 0  // Compatibility with non-clang compilers.
      #endif
      #if __has_feature(cxx_thread_local) || \
          ( (__GNUC__ * 100 + __GNUC_MINOR__) >= 408 && __cplusplus >= 201103L ) || \
          ( _MSC_VER >= 1900 )
        bool has_feature_thread_local = true;
      #else
        // Explicit error of compilation for CMake test purpose - has_feature_thread_local is not defined
      #endif
        bool result = has_feature_thread_local;
      }  ")
  check_cxx_source_compiles("${CGAL_CAN_USE_CXX11_THREAD_LOCAL}" CGAL_CAN_USE_CXX11_THREAD_LOCAL_RESULT)
endfunction()

set (CMAKE_CXX_STANDARD 14)

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

set(GUDHI_CAN_USE_CXX11_THREAD_LOCAL "
    int main() {
      thread_local int result = 0;
      return result;
    }  ")
check_cxx_source_compiles("${GUDHI_CAN_USE_CXX11_THREAD_LOCAL}" GUDHI_CAN_USE_CXX11_THREAD_LOCAL_RESULT)
if (GUDHI_CAN_USE_CXX11_THREAD_LOCAL_RESULT)
  add_definitions(-DGUDHI_CAN_USE_CXX11_THREAD_LOCAL)
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug)
  message("++ Debug compilation flags are: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
else()
  message("++ Release compilation flags are: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
endif()
