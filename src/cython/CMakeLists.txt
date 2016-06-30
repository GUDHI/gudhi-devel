cmake_minimum_required(VERSION 2.8)
project(Cython)

FIND_PROGRAM( PYTHON_PATH python )
FIND_PROGRAM( CYTHON_PATH cython )

if(PYTHON_PATH AND CYTHON_PATH)
  # Gudhi and CGAL compilation option
  if(MSVC)
    set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'/fp:strict', ")
  else(MSVC)
    set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-std=c++11', ")
  endif(MSVC)
  if(CMAKE_COMPILER_IS_GNUCXX)
    set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-frounding-math', ")
  endif(CMAKE_COMPILER_IS_GNUCXX)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-fp-model strict', ")
  endif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")

  find_package(Eigen3 3.1.0)
  
  if (EIGEN3_FOUND)
    # No problem, even if no CGAL found
    set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-DCGAL_EIGEN3_ENABLED', ")
    set(GUDHI_CYTHON_INCLUDE_DIRS "${GUDHI_CYTHON_INCLUDE_DIRS}'${EIGEN3_INCLUDE_DIR}', ")
  endif (EIGEN3_FOUND)
  
  # Copy recursively src and test repositories before packages finding
  # Some tests files are removed in case some packages are not found
  file(COPY src DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY test DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  
  if (CGAL_FOUND)
    if (NOT CGAL_VERSION VERSION_LESS 4.7.0)
      # If CGAL_VERSION >= 4.7.0, include alpha complex
      set(GUDHI_CYTHON_LIBRARIES "${GUDHI_CYTHON_LIBRARIES}'CGAL', ")
      set(GUDHI_CYTHON_LIBRARY_DIRS "${GUDHI_CYTHON_LIBRARY_DIRS}'${CGAL_LIBRARIES_DIR}', ")
      # GMP and GMPXX are not required, but if present, CGAL will link with them.
      if(GMP_FOUND)
        set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-DCGAL_USE_GMP', ")
        set(GUDHI_CYTHON_LIBRARIES "${GUDHI_CYTHON_LIBRARIES}'gmp', ")
        set(GUDHI_CYTHON_LIBRARY_DIRS "${GUDHI_CYTHON_LIBRARY_DIRS}'${GMP_LIBRARIES_DIR}', ")
        if(GMPXX_FOUND)
          set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-DCGAL_USE_GMPXX', ")
          set(GUDHI_CYTHON_LIBRARIES "${GUDHI_CYTHON_LIBRARIES}'gmpxx', ")
          set(GUDHI_CYTHON_LIBRARY_DIRS "${GUDHI_CYTHON_LIBRARY_DIRS}'${GMPXX_LIBRARIES_DIR}', ")
        endif(GMPXX_FOUND)
      endif(GMP_FOUND)
  
      set(GUDHI_CYTHON_ALPHA_COMPLEX "include 'src/cython/alpha_complex.pyx'")
    else (NOT CGAL_VERSION VERSION_LESS 4.7.0)
      # Remove alpha complex unitary tests
      file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/test/test_alpha_complex.py)
    endif (NOT CGAL_VERSION VERSION_LESS 4.7.0)
  endif (CGAL_FOUND)
  
  # Loop on INCLUDE_DIRECTORIES PROPERTY
  get_property(GUDHI_INCLUDE_DIRECTORIES DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
  foreach(GUDHI_INCLUDE_DIRECTORY ${GUDHI_INCLUDE_DIRECTORIES})
    set(GUDHI_CYTHON_INCLUDE_DIRS "${GUDHI_CYTHON_INCLUDE_DIRS}'${GUDHI_INCLUDE_DIRECTORY}', ")
  endforeach()
  set(GUDHI_CYTHON_INCLUDE_DIRS "${GUDHI_CYTHON_INCLUDE_DIRS}'${CMAKE_SOURCE_DIR}/${GUDHI_CYTHON_PATH}/src/cpp', ")
  
  # Generate cythonize_gudhi.py file to cythonize Gudhi
  configure_file(cythonize_gudhi.py.in "${CMAKE_CURRENT_BINARY_DIR}/cythonize_gudhi.py" @ONLY)
  # Generate gudhi.pyx - Gudhi cython file 
  configure_file(gudhi.pyx.in "${CMAKE_CURRENT_BINARY_DIR}/gudhi.pyx" @ONLY)
  
  add_custom_command(
      OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/gudhi.so"
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMAND python "${CMAKE_CURRENT_BINARY_DIR}/cythonize_gudhi.py" "build_ext" "--inplace")
  
  add_custom_target(cythonize_gudhi ALL DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/gudhi.so")

  find_program( PYTEST_PATH py.test )
  if(PYTEST_PATH)
    add_test(
        NAME gudhi_cython_py_test
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND ${PYTEST_PATH})
    set_tests_properties(gudhi_cython_py_test PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")
  endif(PYTEST_PATH)

endif(PYTHON_PATH AND CYTHON_PATH)
