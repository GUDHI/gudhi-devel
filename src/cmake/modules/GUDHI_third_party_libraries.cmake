# This files manage third party libraries required by GUDHI

find_package(Boost REQUIRED COMPONENTS system filesystem unit_test_framework program_options thread)

if(NOT Boost_FOUND)
  message(FATAL_ERROR "NOTICE: This program requires Boost and will not be compiled.")
endif(NOT Boost_FOUND)

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
    # HACK to detect CGAL version 4.8.0
    # CGAL version 4.8, 4.8.1 and 4.8.2 are identified as version 4.8.1000)
    # cf. https://github.com/CGAL/cgal/issues/1559
    # Limit the HACK between CGAL versions 4.8 and 4.9 because of file read
    if (NOT CGAL_VERSION VERSION_GREATER 4.9.0)
      foreach(CGAL_INCLUDE_DIR ${CGAL_INCLUDE_DIRS})
        if (EXISTS "${CGAL_INCLUDE_DIR}/CGAL/version.h")
          FILE(READ "${CGAL_INCLUDE_DIR}/CGAL/version.h" contents)
          STRING(REGEX REPLACE "\n" ";" contents "${contents}")
          foreach(Line ${contents})
            if("${Line}" STREQUAL "#define CGAL_VERSION 4.8")
              set(CGAL_VERSION 4.8.0)
              message (">>>>> HACK CGAL version to ${CGAL_VERSION}")
            endif("${Line}" STREQUAL "#define CGAL_VERSION 4.8")
          endforeach(Line ${contents})
        endif (EXISTS "${CGAL_INCLUDE_DIR}/CGAL/version.h")
      endforeach(CGAL_INCLUDE_DIR ${CGAL_INCLUDE_DIRS})
    endif(NOT CGAL_VERSION VERSION_GREATER 4.9.0)

    # For dev version
    include_directories(BEFORE "src/common/include/gudhi_patches")
    # For user version
    include_directories(BEFORE "include/gudhi_patches")
  endif()
endif()

# Find TBB package for parallel sort - not mandatory, just optional.
set(TBB_FIND_QUIETLY ON)
find_package(TBB)
if (TBB_FOUND)
  include(${TBB_USE_FILE})
  message("TBB found in ${TBB_LIBRARY_DIRS}")
  add_definitions(-DGUDHI_USE_TBB)
endif()

set(CGAL_WITH_EIGEN3_VERSION 0.0.0)
find_package(Eigen3 3.1.0)
if (EIGEN3_FOUND)
  message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
  include( ${EIGEN3_USE_FILE} )
  set(CGAL_WITH_EIGEN3_VERSION ${CGAL_VERSION})
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
# problem with Visual Studio link on Boost program_options
add_definitions( -DBOOST_ALL_DYN_LINK )
# problem on Mac with boost_system and boost_thread
add_definitions( -DBOOST_SYSTEM_NO_DEPRECATED )

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

message(STATUS "boost include dirs:" ${Boost_INCLUDE_DIRS})
message(STATUS "boost library dirs:" ${Boost_LIBRARY_DIRS})

# Find the correct Python interpreter.
# Can be set with -DPYTHON_EXECUTABLE=/usr/bin/python3 or -DPython_ADDITIONAL_VERSIONS=3 for instance.
find_package(Cython)

if(NOT GUDHI_CYTHON_PATH)
  message(FATAL_ERROR "ERROR: GUDHI_CYTHON_PATH is not valid.")
endif(NOT GUDHI_CYTHON_PATH)

if(PYTHONINTERP_FOUND AND CYTHON_FOUND)
  # Unitary tests are available through py.test
  find_program( PYTEST_PATH py.test )
  # Default found version 2
  if(PYTHON_VERSION_MAJOR EQUAL 2)
    # Documentation generation is available through sphinx
    find_program( SPHINX_PATH sphinx-build )
  elseif(PYTHON_VERSION_MAJOR EQUAL 3)
    # No sphinx-build in Pyton3, just hack it
    set(SPHINX_PATH "${CMAKE_SOURCE_DIR}/${GUDHI_CYTHON_PATH}/doc/python3-sphinx-build")
  else()
    message(FATAL_ERROR "ERROR: Try to compile the Cython interface. Python version ${PYTHON_VERSION_STRING} is not valid.")
  endif(PYTHON_VERSION_MAJOR EQUAL 2)
  # get PYTHON_SITE_PACKAGES relative path from a python command line
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c "from distutils.sysconfig import get_python_lib; print (get_python_lib(prefix='', plat_specific=True))"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE)
endif(PYTHONINTERP_FOUND AND CYTHON_FOUND)

