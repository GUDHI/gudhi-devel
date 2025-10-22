# This files manage third party libraries required by GUDHI

find_package(Boost 1.71.0 QUIET OPTIONAL_COMPONENTS filesystem unit_test_framework program_options)

# Boost_FOUND is not reliable
if(NOT Boost_VERSION)
    message(FATAL_ERROR "NOTICE: This program requires Boost and will not be compiled.")
endif(NOT Boost_VERSION)
message("++ BOOST version ${Boost_VERSION}. Includes found in ${Boost_INCLUDE_DIRS}, libraries found in ${Boost_LIBRARY_DIRS}")
include_directories(${Boost_INCLUDE_DIRS})

find_package(GMP)
if(GMP_FOUND)
    INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR})
    find_package(GMPXX)
    if(GMPXX_FOUND)
        INCLUDE_DIRECTORIES(${GMPXX_INCLUDE_DIR})
    endif()
endif()

# from windows vcpkg eigen 3.4.0 #2 : build fails with
# error C2440: '<function-style-cast>': cannot convert from 'Eigen::EigenBase<Derived>::Index' to '__gmp_expr<mpq_t,mpq_t>'
# cf. https://gitlab.com/libeigen/eigen/-/issues/2476
# Workaround is to compile with '-DEIGEN_DEFAULT_DENSE_INDEX_TYPE=int'
if (FORCE_EIGEN_DEFAULT_DENSE_INDEX_TYPE_TO_INT)
    message("++ User explicit demand to force EIGEN_DEFAULT_DENSE_INDEX_TYPE to int")
    add_definitions(-DEIGEN_DEFAULT_DENSE_INDEX_TYPE=int)
endif()

find_package(CGAL 5.1.0)

if (TARGET CGAL::CGAL)
    message("++ CGAL version: ${CGAL_VERSION}. Includes found in ${CGAL_INCLUDE_DIRS}")
endif ()

find_package(Eigen3 3.3...5 NO_MODULE) # Any version >=3.3.0 but <6.0.0
if(TARGET Eigen3::Eigen)
    # Not mandatory as it is set by Eigen3Config.cmake
    get_target_property(EIGEN3_INCLUDE_DIRS Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
    message("++ Eigen 3 version ${EIGEN3_VERSION_STRING}. Includes found in ${EIGEN3_INCLUDE_DIRS}")
endif()

option(WITH_GUDHI_USE_TBB "Build with Intel TBB parallelization" ON)

# Find TBB package for parallel sort - not mandatory, just optional.
if(WITH_GUDHI_USE_TBB)
    find_package(TBB CONFIG)
    if(TARGET TBB::tbb)
        # Specific windows case with its debug/release management
        if(CMAKE_BUILD_TYPE MATCHES Debug)
            get_target_property(TBB_LIBRARY TBB::tbb IMPORTED_LOCATION_DEBUG)
            get_target_property(TBB_MALLOC_LIBRARY TBB::tbbmalloc IMPORTED_LOCATION_DEBUG)
        else()
            get_target_property(TBB_LIBRARY TBB::tbb IMPORTED_LOCATION_RELEASE)
            get_target_property(TBB_MALLOC_LIBRARY TBB::tbbmalloc IMPORTED_LOCATION_RELEASE)
        endif()
        # Generic case
        if (NOT TBB_LIBRARY)
            get_target_property(TBB_LIBRARY TBB::tbb LOCATION)
            get_target_property(TBB_MALLOC_LIBRARY TBB::tbbmalloc LOCATION)
        endif()
        # TBB Error management
        if (TBB_VERSION VERSION_LESS 2019.0.11007)
            # TBBTargets.cmake was introduced in 2019.7, so this case should not happen
            # cf. https://github.com/oneapi-src/oneTBB/blob/2019_U7/CHANGES
            message(WARNING "++ TBB found but version ${TBB_VERSION} is too old - GUDHI cannot compile with TBB")
        else()
            if (NOT TBB_LIBRARY)
                message(WARNING "++ TBB found but not ${TBB_LIBRARY} - GUDHI cannot compile with TBB")
            else()
                # A correct version of TBB was found
                get_target_property(TBB_INCLUDE_DIRS TBB::tbb INTERFACE_INCLUDE_DIRECTORIES)
                get_filename_component(TBB_LIBRARY_DIRS ${TBB_LIBRARY} DIRECTORY)
                message("++ TBB version ${TBB_VERSION}. Includes found in ${TBB_INCLUDE_DIRS}, libraries found in ${TBB_LIBRARY_DIRS}")
                add_definitions(-DGUDHI_USE_TBB)
                if(MSVC)
                    # cf. https://github.com/oneapi-src/oneTBB/issues/573
                    add_definitions(-DNOMINMAX)
                endif()
            endif()
        endif()
    endif()
endif()

function(add_executable_with_targets)
    if (ARGC LESS_EQUAL 2)
        message (FATAL_ERROR "add_executable_with_targets requires at least 2 arguments.")
    endif()

    list(POP_FRONT ARGN EXECUTABLE_NAME EXECUTABLE_SOURCE)
    message(DEBUG "${EXECUTABLE_NAME} - ${EXECUTABLE_SOURCE}")
    # Do not add_executable if one of the target is not here, except for TBB that is optional
    foreach(USER_TARGET IN LISTS ARGN)
        if(NOT TARGET ${USER_TARGET})
            if(NOT ${USER_TARGET} STREQUAL "TBB::tbb")
                return()
            endif()
        endif()
    endforeach()
    add_executable(${EXECUTABLE_NAME} ${EXECUTABLE_SOURCE})
    foreach(USER_TARGET IN LISTS ARGN)
        message(DEBUG "target_link_libraries(${EXECUTABLE_NAME} ${USER_TARGET})")
        # TARGET_NAME_IF_EXISTS is specific to TBB case (optional)
        target_link_libraries(${EXECUTABLE_NAME} $<TARGET_NAME_IF_EXISTS:${USER_TARGET}>)
    endforeach()
endfunction()

# Required programs for unitary tests purpose
FIND_PROGRAM( GCOVR_PATH gcovr )
if (GCOVR_PATH)
    message("gcovr found in ${GCOVR_PATH}")
endif()
FIND_PROGRAM( GPROF_PATH gprof )
if (GPROF_PATH)
    message("gprof found in ${GPROF_PATH}")
endif()
FIND_PROGRAM( DIFF_PATH diff )
if (DIFF_PATH)
    message("diff found in ${DIFF_PATH}")
endif()
FIND_PROGRAM( GNUPLOT_PATH gnuplot )
if (GNUPLOT_PATH)
    message("gnuplot found in ${GNUPLOT_PATH}")
endif()

# BOOST ISSUE result_of vs C++11
add_definitions(-DBOOST_RESULT_OF_USE_DECLTYPE)
# BOOST ISSUE with Libraries name resolution under Windows
add_definitions(-DBOOST_ALL_NO_LIB)
# problem with Visual Studio link on Boost program_options
add_definitions( -DBOOST_ALL_DYN_LINK )
# problem on Mac with boost_system and boost_thread
add_definitions( -DBOOST_SYSTEM_NO_DEPRECATED )

if (WITH_GUDHI_PYTHON)
    # Find the correct Python interpreter.
    # Can be set with -DPython_EXECUTABLE=/usr/bin/python3 for instance.
    # Default Python_FIND_STRATEGY to LOCATION: Stops lookup as soon as a version satisfying version constraints is founded
    # (as opposed to VERSION: Try to find the most recent version in all specified locations.)
    cmake_policy(SET CMP0094 NEW)

    # cf. https://nanobind.readthedocs.io/en/latest/building.html #preliminaries
    if (CMAKE_VERSION VERSION_LESS 3.18)
        set(DEV_MODULE Development)
    else()
        set(DEV_MODULE Development.Module)
    endif()

    find_package( Python COMPONENTS Interpreter ${DEV_MODULE} NumPy)

    # find_python_module tries to import module in Python interpreter and to retrieve its version number
    # returns ${PYTHON_MODULE_NAME_UP}_VERSION and ${PYTHON_MODULE_NAME_UP}_FOUND
    function( find_python_module PYTHON_MODULE_NAME )
        string(TOUPPER ${PYTHON_MODULE_NAME} PYTHON_MODULE_NAME_UP)
        # Modify tracebacklimit as the exception is quite verbose when module is not found
        execute_process(
            COMMAND ${Python_EXECUTABLE}  -c "import sys; sys.tracebacklimit = 0; from importlib.metadata import version; print(version('${PYTHON_MODULE_NAME}'))"
            RESULT_VARIABLE PYTHON_MODULE_RESULT
            OUTPUT_VARIABLE PYTHON_MODULE_VERSION
            ERROR_VARIABLE PYTHON_MODULE_ERROR)
        if(PYTHON_MODULE_RESULT EQUAL 0)
            # Remove all carriage returns as it can be multiline
            string(REGEX REPLACE "\n" " " PYTHON_MODULE_VERSION "${PYTHON_MODULE_VERSION}")
            message ("++ Python module ${PYTHON_MODULE_NAME} - Version ${PYTHON_MODULE_VERSION} found")

            set(${PYTHON_MODULE_NAME_UP}_VERSION ${PYTHON_MODULE_VERSION} PARENT_SCOPE)
            set(${PYTHON_MODULE_NAME_UP}_FOUND TRUE PARENT_SCOPE)
        else()
            # Remove all carriage returns as it can be multiline
            string(REGEX REPLACE "\n" " " PYTHON_MODULE_ERROR "${PYTHON_MODULE_ERROR}")
            message ("## Python module ${PYTHON_MODULE_NAME} NOT found [${PYTHON_MODULE_RESULT}]: '${PYTHON_MODULE_ERROR}'")
            unset(${PYTHON_MODULE_NAME_UP}_VERSION PARENT_SCOPE)
            set(${PYTHON_MODULE_NAME_UP}_FOUND FALSE PARENT_SCOPE)
        endif()
    endfunction( find_python_module )

    if( TARGET Python::Interpreter )
        # Mandatory
        find_python_module("nanobind")
        find_python_module("scikit_build_core")
        find_python_module("numpy")
        # Tests purpose
        find_python_module("pytest")
        # Documentation purpose
        find_python_module("sphinx")
        find_python_module("sphinx_paramlinks")
        find_python_module("pydata_sphinx_theme")
        find_python_module("sphinxcontrib.bibtex")
        find_python_module("sphinx-autodoc-typehints")
        # Optional third parties libraries (also required for documentation)
        find_python_module("matplotlib")
        find_python_module("scipy")
        find_python_module("scikit-learn")
        find_python_module("POT")
        find_python_module("torch")
        find_python_module("eagerpy")
        find_python_module("hnswlib")
        find_python_module("tensorflow")
        find_python_module("networkx")
        # Specific case for PyKeops on Windows where import fails because it uses fcntl (not available on Windows)
        if (NOT WIN32)
            find_python_module("pykeops")
        endif()
    endif()


    if(NOT GUDHI_PYTHON_PATH)
        message(FATAL_ERROR "ERROR: GUDHI_PYTHON_PATH is not valid.")
    endif(NOT GUDHI_PYTHON_PATH)

endif (WITH_GUDHI_PYTHON)
