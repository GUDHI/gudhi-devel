if (NOT TARGET Python::Interpreter)
  message("++ Python module will not be compiled because no Python interpreter was found")
  set(GUDHI_MISSING_MODULES ${GUDHI_MISSING_MODULES} "python" CACHE INTERNAL "GUDHI_MISSING_MODULES")
  return()
endif()
if (NOT NUMPY_FOUND)
  message("++ Python module will not be compiled because numpy was not found")
  set(GUDHI_MISSING_MODULES ${GUDHI_MISSING_MODULES} "python" CACHE INTERNAL "GUDHI_MISSING_MODULES")
  return()
endif()
if(NOT CYTHON_FOUND)
  message("++ Python module will not be compiled because cython was not found")
  set(GUDHI_MISSING_MODULES ${GUDHI_MISSING_MODULES} "python" CACHE INTERNAL "GUDHI_MISSING_MODULES")
  return()
endif()
if(NOT PYBIND11_FOUND)
  message("++ Python module will not be compiled because pybind11 was not found")
  set(GUDHI_MISSING_MODULES ${GUDHI_MISSING_MODULES} "python" CACHE INTERNAL "GUDHI_MISSING_MODULES")
  return()
endif()

function( add_GUDHI_PYTHON_lib THE_LIB )
  if(EXISTS ${THE_LIB})
    get_filename_component(THE_LIB_FILE_NAME ${THE_LIB} NAME_WE)
    if(WIN32)
      message("++ ${THE_LIB} => THE_LIB_FILE_NAME = ${THE_LIB_FILE_NAME}")
      set(GUDHI_PYTHON_LIBRARIES "${GUDHI_PYTHON_LIBRARIES}'${THE_LIB_FILE_NAME}', " PARENT_SCOPE)
    else(WIN32)
      STRING(REGEX REPLACE "lib" "" UNIX_LIB_FILE_NAME ${THE_LIB_FILE_NAME})
      message("++ ${THE_LIB} => UNIX_LIB_FILE_NAME = ${UNIX_LIB_FILE_NAME}")
      set(GUDHI_PYTHON_LIBRARIES "${GUDHI_PYTHON_LIBRARIES}'${UNIX_LIB_FILE_NAME}', " PARENT_SCOPE)
    endif(WIN32)
  endif(EXISTS ${THE_LIB})
endfunction( add_GUDHI_PYTHON_lib )

function( add_GUDHI_PYTHON_lib_dir)
  # Argument may be a list (specifically on windows with release/debug paths)
  foreach(THE_LIB_DIR IN LISTS ARGN)
    # deals when it is not set - error on windows
    if(EXISTS ${THE_LIB_DIR})
      set(GUDHI_PYTHON_LIBRARY_DIRS "${GUDHI_PYTHON_LIBRARY_DIRS}'${THE_LIB_DIR}', " PARENT_SCOPE)
    else()
      message("add_GUDHI_PYTHON_lib_dir - '${THE_LIB_DIR}' does not exist")
    endif()
  endforeach()
endfunction( add_GUDHI_PYTHON_lib_dir )

# THE_TEST is the python test file name (without .py extension) containing tests functions
function( add_gudhi_py_test THE_TEST )
  if(PYTEST_FOUND)
    # use ${Python_EXECUTABLE} -B, otherwise a __pycache__ directory is created in sources by python
    # use py.test no cache provider, otherwise a .cache file is created in sources by py.test
    add_test(NAME ${THE_TEST}_py_test
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
             COMMAND ${Python_EXECUTABLE} -B -m pytest -p no:cacheprovider ${CMAKE_CURRENT_SOURCE_DIR}/test/${THE_TEST}.py)
  endif()
endfunction( add_gudhi_py_test )

function( add_gudhi_symbolic_links GLOBBING_EXPRESSION )
  # No symbolic links if build in sources (aka. 'cmake .' in sources)
  if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
    return()
  endif()
  file(GLOB GUDHI_GLOB_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${GLOBBING_EXPRESSION})
  foreach(GUDHI_GLOB_FILENAME ${GUDHI_GLOB_FILES})
    # COPY_ON_ERROR is activated because on windows 10, user needs to be root or in dev mode to be able to make symbolic links
    # This case can be problematic, because if you modify the sources, build is still done on the copy
    file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/${GUDHI_GLOB_FILENAME}" "${CMAKE_CURRENT_BINARY_DIR}/${GUDHI_GLOB_FILENAME}" COPY_ON_ERROR SYMBOLIC)
  endforeach()
endfunction( add_gudhi_symbolic_links )

# Set gudhi.__debug_info__
# WARNING : to be done before setup.py.in configure_file
function( add_gudhi_debug_info DEBUG_INFO )
  set(GUDHI_PYTHON_DEBUG_INFO "${GUDHI_PYTHON_DEBUG_INFO}    \"${DEBUG_INFO}\\n\" \\\n" PARENT_SCOPE)
endfunction( add_gudhi_debug_info )

function(disable_python_documentation MISSING_MODULE)
  set(GUDHI_PYTHON_DOCUMENTATION_MISSING_MODULES "${GUDHI_PYTHON_DOCUMENTATION_MISSING_MODULES}'${MISSING_MODULE}', " PARENT_SCOPE)
  set(GUDHI_MISSING_MODULES ${GUDHI_MISSING_MODULES} "python-documentation" CACHE INTERNAL "GUDHI_MISSING_MODULES")
endfunction( disable_python_documentation )

add_gudhi_debug_info("Python version ${Python_VERSION}")
add_gudhi_debug_info("Pybind11 version ${PYBIND11_VERSION}")
# PyBind11 modules
set(GUDHI_PYTHON_MODULES "${GUDHI_PYTHON_MODULES}'bottleneck', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'hera', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'clustering', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'datasets', ")

# Cython modules
set(GUDHI_PYTHON_MODULES "${GUDHI_PYTHON_MODULES}'off_utils', ")
set(GUDHI_PYTHON_MODULES "${GUDHI_PYTHON_MODULES}'simplex_tree', ")
set(GUDHI_PYTHON_MODULES "${GUDHI_PYTHON_MODULES}'rips_complex', ")
set(GUDHI_PYTHON_MODULES "${GUDHI_PYTHON_MODULES}'cubical_complex', ")
set(GUDHI_PYTHON_MODULES "${GUDHI_PYTHON_MODULES}'periodic_cubical_complex', ")
set(GUDHI_PYTHON_MODULES "${GUDHI_PYTHON_MODULES}'persistence_graphical_tools', ")
set(GUDHI_PYTHON_MODULES "${GUDHI_PYTHON_MODULES}'reader_utils', ")
set(GUDHI_PYTHON_MODULES "${GUDHI_PYTHON_MODULES}'witness_complex', ")
set(GUDHI_PYTHON_MODULES "${GUDHI_PYTHON_MODULES}'strong_witness_complex', ")
set(GUDHI_PYTHON_MODULES "${GUDHI_PYTHON_MODULES}'nerve_gic', ")
set(GUDHI_PYTHON_MODULES "${GUDHI_PYTHON_MODULES}'subsampling', ")
set(GUDHI_PYTHON_MODULES "${GUDHI_PYTHON_MODULES}'tangential_complex', ")
set(GUDHI_PYTHON_MODULES "${GUDHI_PYTHON_MODULES}'delaunay_complex', ")
set(GUDHI_PYTHON_MODULES "${GUDHI_PYTHON_MODULES}'euclidean_witness_complex', ")
set(GUDHI_PYTHON_MODULES "${GUDHI_PYTHON_MODULES}'euclidean_strong_witness_complex', ")
# Modules that should not be auto-imported in __init__.py
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'representations', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'sklearn', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'tensorflow', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'wasserstein', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'point_cloud', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'weighted_rips_complex', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'dtm_rips_complex', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'cover_complex', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'flag_filtration', ")

add_gudhi_debug_info("Cython version ${CYTHON_VERSION}")
add_gudhi_debug_info("Numpy version ${NUMPY_VERSION}")

if(PYTEST_FOUND)
  add_gudhi_debug_info("Pytest version ${PYTEST_VERSION}")
endif()
if(MATPLOTLIB_FOUND)
  add_gudhi_debug_info("Matplotlib version ${MATPLOTLIB_VERSION}")
else()
  disable_python_documentation("matplotlib")
endif()
if(SCIPY_FOUND)
  add_gudhi_debug_info("Scipy version ${SCIPY_VERSION}")
else()
  disable_python_documentation("scipy")
endif()
if(SCIKIT-LEARN_FOUND)
  add_gudhi_debug_info("Scikit-learn version ${SCIKIT-LEARN_VERSION}")
else()
  disable_python_documentation("scikit-learn")
endif()
if(POT_FOUND)
  add_gudhi_debug_info("POT version ${POT_VERSION}")
else()
  disable_python_documentation("ot")
endif()
if(HNSWLIB_FOUND)
  add_gudhi_debug_info("HNSWlib version ${HNSWLIB_VERSION}")
endif()
if(TORCH_FOUND)
  add_gudhi_debug_info("PyTorch version ${TORCH_VERSION}")
endif()
if(PYKEOPS_FOUND)
  add_gudhi_debug_info("PyKeOps version ${PYKEOPS_VERSION}")
endif()
if(EAGERPY_FOUND)
  add_gudhi_debug_info("EagerPy version ${EAGERPY_VERSION}")
endif()
if(TENSORFLOW_FOUND)
  add_gudhi_debug_info("TensorFlow version ${TENSORFLOW_VERSION}")
endif()
if(SPHINX_FOUND)
  add_gudhi_debug_info("Sphinx version ${SPHINX_VERSION}")
  # Calling sphinx-build executable (that could be found find_program) may use a different version of python and fail
  set(SPHINX_PATH "${Python_EXECUTABLE}" "-m" "sphinx.cmd.build")
else()
  disable_python_documentation("sphinx")
endif()
if(SPHINX_PARAMLINKS_FOUND)
  add_gudhi_debug_info("Sphinx-paramlinks version ${SPHINX_PARAMLINKS_VERSION}")
else()
  disable_python_documentation("sphinx-paramlinks")
endif()
if(PYDATA_SPHINX_THEME_FOUND)
  add_gudhi_debug_info("pydata_sphinx_theme version ${PYDATA_SPHINX_THEME_VERSION}")
else()
  disable_python_documentation("pydata_sphinx_theme")
endif()
if(SPHINXCONTRIB.BIBTEX_FOUND)
  add_gudhi_debug_info("sphinxcontrib-bibtex version ${SPHINXCONTRIB.BIBTEX_VERSION}")
else()
  disable_python_documentation("sphinxcontrib-bibtex")
endif()
if(SPHINX-AUTODOC-TYPEHINTS_FOUND)
  add_gudhi_debug_info("sphinx-autodoc-typehints version ${SPHINX-AUTODOC-TYPEHINTS_VERSION}")
else()
  disable_python_documentation("sphinx-autodoc-typehints")
endif()
if(NETWORKX_FOUND)
  add_gudhi_debug_info("NetworkX version ${NETWORKX_VERSION}")
endif()
if(NOT TARGET CGAL::CGAL)
  disable_python_documentation("CGAL")
endif()
if (NOT TARGET Eigen3::Eigen)
  disable_python_documentation("Eigen3")
endif()

set(GUDHI_PYTHON_EXTRA_COMPILE_ARGS "${GUDHI_PYTHON_EXTRA_COMPILE_ARGS}'-DBOOST_RESULT_OF_USE_DECLTYPE', ")
set(GUDHI_PYTHON_EXTRA_COMPILE_ARGS "${GUDHI_PYTHON_EXTRA_COMPILE_ARGS}'-DBOOST_ALL_NO_LIB', ")
set(GUDHI_PYTHON_EXTRA_COMPILE_ARGS "${GUDHI_PYTHON_EXTRA_COMPILE_ARGS}'-DBOOST_SYSTEM_NO_DEPRECATED', ")

# Gudhi and CGAL compilation option
if(MSVC)
  set(GUDHI_PYTHON_EXTRA_COMPILE_ARGS "${GUDHI_PYTHON_EXTRA_COMPILE_ARGS}'/std:c++17', ")
  set(GUDHI_PYTHON_EXTRA_COMPILE_ARGS "${GUDHI_PYTHON_EXTRA_COMPILE_ARGS}'/fp:strict', ")
  # cf. https://github.com/oneapi-src/oneTBB/issues/573
  set(GUDHI_PYTHON_EXTRA_COMPILE_ARGS "${GUDHI_PYTHON_EXTRA_COMPILE_ARGS}'-DNOMINMAX', ")
else(MSVC)
  set(GUDHI_PYTHON_EXTRA_COMPILE_ARGS "${GUDHI_PYTHON_EXTRA_COMPILE_ARGS}'-std=c++17', ")
endif(MSVC)
if(CMAKE_COMPILER_IS_GNUCXX)
  set(GUDHI_PYTHON_EXTRA_COMPILE_ARGS "${GUDHI_PYTHON_EXTRA_COMPILE_ARGS}'-frounding-math', ")
  set(GUDHI_PYBIND11_EXTRA_COMPILE_ARGS "${GUDHI_PYBIND11_EXTRA_COMPILE_ARGS}'-fvisibility=hidden', ")
endif(CMAKE_COMPILER_IS_GNUCXX)
if (CMAKE_CXX_COMPILER_ID MATCHES Intel)
  set(GUDHI_PYTHON_EXTRA_COMPILE_ARGS "${GUDHI_PYTHON_EXTRA_COMPILE_ARGS}'-fp-model strict', ")
endif(CMAKE_CXX_COMPILER_ID MATCHES Intel)
if (DEBUG_TRACES)
  # For programs to be more verbose
  set(GUDHI_PYTHON_EXTRA_COMPILE_ARGS "${GUDHI_PYTHON_EXTRA_COMPILE_ARGS}'-DDEBUG_TRACES', ")
endif()

if (TARGET Eigen3::Eigen)
  add_gudhi_debug_info("Eigen3 version ${EIGEN3_VERSION_STRING}")
  # No problem, even if no CGAL found
  set(GUDHI_PYTHON_EXTRA_COMPILE_ARGS "${GUDHI_PYTHON_EXTRA_COMPILE_ARGS}'-DCGAL_EIGEN3_ENABLED', ")
  set(GUDHI_PYTHON_INCLUDE_DIRS "${GUDHI_PYTHON_INCLUDE_DIRS}'${EIGEN3_INCLUDE_DIRS}', ")
endif (TARGET Eigen3::Eigen)

set(GUDHI_CYTHON_MODULES "${GUDHI_CYTHON_MODULES}'off_utils', ")
set(GUDHI_CYTHON_MODULES "${GUDHI_CYTHON_MODULES}'simplex_tree', ")
set(GUDHI_CYTHON_MODULES "${GUDHI_CYTHON_MODULES}'rips_complex', ")
set(GUDHI_CYTHON_MODULES "${GUDHI_CYTHON_MODULES}'cubical_complex', ")
set(GUDHI_CYTHON_MODULES "${GUDHI_CYTHON_MODULES}'periodic_cubical_complex', ")
set(GUDHI_CYTHON_MODULES "${GUDHI_CYTHON_MODULES}'reader_utils', ")
set(GUDHI_CYTHON_MODULES "${GUDHI_CYTHON_MODULES}'witness_complex', ")
set(GUDHI_CYTHON_MODULES "${GUDHI_CYTHON_MODULES}'strong_witness_complex', ")
set(GUDHI_PYBIND11_MODULES "${GUDHI_PYBIND11_MODULES}'clustering/_tomato', ")
set(GUDHI_PYBIND11_MODULES "${GUDHI_PYBIND11_MODULES}'hera/wasserstein', ")
set(GUDHI_PYBIND11_MODULES "${GUDHI_PYBIND11_MODULES}'hera/bottleneck', ")
set(GUDHI_CYTHON_MODULES "${GUDHI_CYTHON_MODULES}'nerve_gic', ")
set(GUDHI_CYTHON_MODULES "${GUDHI_CYTHON_MODULES}'subsampling', ")
if (TARGET CGAL::CGAL)
  set(GUDHI_PYBIND11_MODULES "${GUDHI_PYBIND11_MODULES}'bottleneck', ")
endif ()
if (TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  set(GUDHI_CYTHON_MODULES "${GUDHI_CYTHON_MODULES}'tangential_complex', ")
  set(GUDHI_CYTHON_MODULES "${GUDHI_CYTHON_MODULES}'euclidean_witness_complex', ")
  set(GUDHI_CYTHON_MODULES "${GUDHI_CYTHON_MODULES}'euclidean_strong_witness_complex', ")
  set(GUDHI_CYTHON_MODULES "${GUDHI_CYTHON_MODULES}'delaunay_complex', ")
  set(GUDHI_PYBIND11_MODULES "${GUDHI_PYBIND11_MODULES}'datasets/generators/_points', ")
endif ()
set(GUDHI_PYBIND11_MODULES "${GUDHI_PYBIND11_MODULES}'_pers_cub_low_dim', ")
set(GUDHI_PYBIND11_MODULES "${GUDHI_PYBIND11_MODULES}'_edge_collapse', ")
set(GUDHI_PYBIND11_MODULES "${GUDHI_PYBIND11_MODULES}'_ripser', ")

# from windows vcpkg eigen 3.4.0#2 : build fails with
# error C2440: '<function-style-cast>': cannot convert from 'Eigen::EigenBase<Derived>::Index' to '__gmp_expr<mpq_t,mpq_t>'
# cf. https://gitlab.com/libeigen/eigen/-/issues/2476
# Workaround is to compile with '-DEIGEN_DEFAULT_DENSE_INDEX_TYPE=int'
if (FORCE_EIGEN_DEFAULT_DENSE_INDEX_TYPE_TO_INT)
  set(GUDHI_PYTHON_EXTRA_COMPILE_ARGS "${GUDHI_PYTHON_EXTRA_COMPILE_ARGS}'-DEIGEN_DEFAULT_DENSE_INDEX_TYPE=int', ")
endif()

add_gudhi_debug_info("Boost version ${Boost_VERSION}")
set(GUDHI_PYTHON_INCLUDE_DIRS "${GUDHI_PYTHON_INCLUDE_DIRS}'${Boost_INCLUDE_DIRS}', ")
if(TARGET CGAL::CGAL)
  foreach(CGAL_INCLUDE_DIR ${CGAL_INCLUDE_DIRS})
    set(GUDHI_PYTHON_INCLUDE_DIRS "${GUDHI_PYTHON_INCLUDE_DIRS}'${CGAL_INCLUDE_DIR}', ")
  endforeach()
  if(NOT CGAL_VERSION VERSION_LESS 5.3.0)
    # CGAL_HEADER_ONLY has been dropped for CGAL >= 5.3. Only the header-only version is supported.
    set(CGAL_HEADER_ONLY True)
  endif(NOT CGAL_VERSION VERSION_LESS 5.3.0)
  # Add CGAL compilation args
  if(CGAL_HEADER_ONLY)
    add_gudhi_debug_info("CGAL header only version ${CGAL_VERSION}")
    set(GUDHI_PYTHON_EXTRA_COMPILE_ARGS "${GUDHI_PYTHON_EXTRA_COMPILE_ARGS}'-DCGAL_HEADER_ONLY', ")
  else(CGAL_HEADER_ONLY)
    add_gudhi_debug_info("CGAL version ${CGAL_VERSION}")
    add_GUDHI_PYTHON_lib("${CGAL_LIBRARY}")
    add_GUDHI_PYTHON_lib_dir(${CGAL_LIBRARIES_DIR})
    message("** Add CGAL ${CGAL_LIBRARIES_DIR}")
    # If CGAL is not header only, CGAL library may link with boost system,
    if(CMAKE_BUILD_TYPE MATCHES Debug)
      add_GUDHI_PYTHON_lib("${Boost_SYSTEM_LIBRARY_DEBUG}")
    else()
      add_GUDHI_PYTHON_lib("${Boost_SYSTEM_LIBRARY_RELEASE}")
    endif()
    add_GUDHI_PYTHON_lib_dir(${Boost_LIBRARY_DIRS})
    message("** Add Boost ${Boost_LIBRARY_DIRS}")
  endif(CGAL_HEADER_ONLY)
  # GMP and GMPXX are not required, but if present, CGAL will link with them.
  if(GMP_FOUND)
    add_gudhi_debug_info("GMP_LIBRARIES = ${GMP_LIBRARIES}")
    set(GUDHI_PYTHON_EXTRA_COMPILE_ARGS "${GUDHI_PYTHON_EXTRA_COMPILE_ARGS}'-DCGAL_USE_GMP', ")
    add_GUDHI_PYTHON_lib("${GMP_LIBRARIES}")
    if(NOT GMP_LIBRARIES_DIR)
      get_filename_component(GMP_LIBRARIES_DIR ${GMP_LIBRARIES} PATH)
      message("GMP_LIBRARIES_DIR from GMP_LIBRARIES set to ${GMP_LIBRARIES_DIR}")
    endif(NOT GMP_LIBRARIES_DIR)
    add_GUDHI_PYTHON_lib_dir(${GMP_LIBRARIES_DIR})
    message("** Add gmp ${GMP_LIBRARIES_DIR}")
    # When FORCE_CGAL_NOT_TO_BUILD_WITH_GMPXX is set, not defining CGAL_USE_GMPXX is sufficient enough
    if(GMPXX_FOUND)
      add_gudhi_debug_info("GMPXX_LIBRARIES = ${GMPXX_LIBRARIES}")
      set(GUDHI_PYTHON_EXTRA_COMPILE_ARGS "${GUDHI_PYTHON_EXTRA_COMPILE_ARGS}'-DCGAL_USE_GMPXX', ")
      add_GUDHI_PYTHON_lib("${GMPXX_LIBRARIES}")
      add_GUDHI_PYTHON_lib_dir(${GMPXX_LIBRARIES_DIR})
      message("** Add gmpxx ${GMPXX_LIBRARIES_DIR}")
    endif()
  endif(GMP_FOUND)
  if(MPFR_FOUND)
    add_gudhi_debug_info("MPFR_LIBRARIES = ${MPFR_LIBRARIES}")
    set(GUDHI_PYTHON_EXTRA_COMPILE_ARGS "${GUDHI_PYTHON_EXTRA_COMPILE_ARGS}'-DCGAL_USE_MPFR', ")
    add_GUDHI_PYTHON_lib("${MPFR_LIBRARIES}")
    # In case CGAL is not header only, all MPFR variables are set except MPFR_LIBRARIES_DIR - Just set it
    if(NOT MPFR_LIBRARIES_DIR)
      get_filename_component(MPFR_LIBRARIES_DIR ${MPFR_LIBRARIES} PATH)
      message("MPFR_LIBRARIES_DIR from MPFR_LIBRARIES set to ${MPFR_LIBRARIES_DIR}")
    endif(NOT MPFR_LIBRARIES_DIR)
    add_GUDHI_PYTHON_lib_dir(${MPFR_LIBRARIES_DIR})
    message("** Add mpfr ${MPFR_LIBRARIES_DIR}")
  endif(MPFR_FOUND)
endif(TARGET CGAL::CGAL)

# Specific for Mac
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  # CMAKE_OSX_DEPLOYMENT_TARGET is always set on OSx
  # but it can be modified by the user by setting MACOSX_DEPLOYMENT_TARGET environment variable
  if(CMAKE_OSX_DEPLOYMENT_TARGET)
    set(GUDHI_PYTHON_EXTRA_COMPILE_ARGS "${GUDHI_PYTHON_EXTRA_COMPILE_ARGS}'-mmacosx-version-min=${CMAKE_OSX_DEPLOYMENT_TARGET}', ")
    set(GUDHI_PYTHON_EXTRA_LINK_ARGS "${GUDHI_PYTHON_EXTRA_LINK_ARGS}'-mmacosx-version-min=${CMAKE_OSX_DEPLOYMENT_TARGET}', ")
  endif()
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Strip dynamic libraries in release mode for smaller pip packages under linux
if(CMAKE_COMPILER_IS_GNUCXX)
  if(CMAKE_BUILD_TYPE MATCHES Release)
    set(GUDHI_PYTHON_EXTRA_LINK_ARGS "${GUDHI_PYTHON_EXTRA_LINK_ARGS}'-s', ")
  endif(CMAKE_BUILD_TYPE MATCHES Release)
endif(CMAKE_COMPILER_IS_GNUCXX)

# Loop on INCLUDE_DIRECTORIES PROPERTY
get_property(GUDHI_INCLUDE_DIRECTORIES DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(GUDHI_INCLUDE_DIRECTORY ${GUDHI_INCLUDE_DIRECTORIES})
  set(GUDHI_PYTHON_INCLUDE_DIRS "${GUDHI_PYTHON_INCLUDE_DIRS}'${GUDHI_INCLUDE_DIRECTORY}', ")
endforeach()
set(GUDHI_PYTHON_INCLUDE_DIRS "${GUDHI_PYTHON_INCLUDE_DIRS}'${CMAKE_SOURCE_DIR}/${GUDHI_PYTHON_PATH}/include', ")

if (TBB_FOUND AND WITH_GUDHI_USE_TBB)
  add_gudhi_debug_info("TBB version ${TBB_VERSION} found and used")
  add_gudhi_debug_info("TBB_LIBRARY = ${TBB_LIBRARY}")
  add_gudhi_debug_info("TBB_MALLOC_LIBRARY = ${TBB_MALLOC_LIBRARY}")

  set(GUDHI_PYTHON_EXTRA_COMPILE_ARGS "${GUDHI_PYTHON_EXTRA_COMPILE_ARGS}'-DGUDHI_USE_TBB', ")
  add_GUDHI_PYTHON_lib("${TBB_LIBRARY}")
  add_GUDHI_PYTHON_lib("${TBB_MALLOC_LIBRARY}")
  add_GUDHI_PYTHON_lib_dir(${TBB_LIBRARY_DIRS})
  message("** Add tbb ${TBB_LIBRARY_DIRS}")
  set(GUDHI_PYTHON_INCLUDE_DIRS "${GUDHI_PYTHON_INCLUDE_DIRS}'${TBB_INCLUDE_DIRS}', ")
endif()

if(DEBUG_TRACES)
  set(GUDHI_PYTHON_EXTRA_COMPILE_ARGS "${GUDHI_PYTHON_EXTRA_COMPILE_ARGS}'-DDEBUG_TRACES', ")
endif(DEBUG_TRACES)

if(UNIX AND WITH_GUDHI_PYTHON_RUNTIME_LIBRARY_DIRS)
  set( GUDHI_PYTHON_RUNTIME_LIBRARY_DIRS "${GUDHI_PYTHON_LIBRARY_DIRS}")
endif(UNIX AND WITH_GUDHI_PYTHON_RUNTIME_LIBRARY_DIRS)

# doc/pip_introduction.inc is a restructured text file (extension changed because of sphinx warning) that describes
# the GUDHI project on https://pypi.org/project/gudhi/
file(READ "doc/pypi_introduction.inc" GUDHI_LONG_DESCRIPTION_WITH_UNPROPER_CR)
string(REPLACE "\n" "\\n' \\\n        '" GUDHI_LONG_DESCRIPTION ${GUDHI_LONG_DESCRIPTION_WITH_UNPROPER_CR})

# Generate setup.py file to cythonize Gudhi - This file must be named setup.py by convention
configure_file(setup.py.in "${CMAKE_CURRENT_BINARY_DIR}/setup.py" @ONLY)

# Generate gudhi/__init__.py
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/gudhi")
configure_file("gudhi/__init__.py.in" "${CMAKE_CURRENT_BINARY_DIR}/gudhi/__init__.py" @ONLY)

# Other .py files
file(COPY "gudhi/persistence_graphical_tools.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gudhi")
file(COPY "gudhi/representations" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gudhi/")
file(COPY "gudhi/tensorflow" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gudhi/")
file(COPY "gudhi/wasserstein" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gudhi")
file(COPY "gudhi/point_cloud" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gudhi")
file(COPY "gudhi/clustering" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gudhi" FILES_MATCHING PATTERN "*.py")
file(COPY "gudhi/weighted_rips_complex.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gudhi")
file(COPY "gudhi/dtm_rips_complex.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gudhi")
file(COPY "gudhi/cover_complex.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gudhi")
file(COPY "gudhi/_kepler_mapper.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gudhi")
file(COPY "gudhi/_kepler_mapper.py.COPYRIGHT" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gudhi")
file(COPY "gudhi/hera/__init__.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gudhi/hera")
file(COPY "gudhi/datasets" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gudhi" FILES_MATCHING PATTERN "*.py")
file(COPY "gudhi/sklearn" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gudhi/")
file(COPY "gudhi/flag_filtration" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/gudhi/")

# Some files for pip package
file(COPY "pyproject.toml" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

# Symbolic links to cython and Pybind11 sources to be built. "MAKE_DIRECTORY gudhi" must be done before
add_gudhi_symbolic_links("gudhi/*.pxd")
add_gudhi_symbolic_links("gudhi/*.pyx")
add_gudhi_symbolic_links("gudhi/*.cc")
add_gudhi_symbolic_links("gudhi/clustering/_tomato.cc")
add_gudhi_symbolic_links("gudhi/hera/*.cc")
add_gudhi_symbolic_links("gudhi/datasets/generators/_points.cc")

add_custom_command(
    OUTPUT gudhi.so
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/setup.py" "build_ext" "--inplace")

add_custom_target(python ALL DEPENDS gudhi.so
                  COMMENT "Do not forget to add ${CMAKE_CURRENT_BINARY_DIR}/ to your PYTHONPATH before using examples or tests")

# Path separator management for windows
if (WIN32)
  set(GUDHI_PYTHON_PATH_ENV "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR};$ENV{PYTHONPATH}")
else(WIN32)
  set(GUDHI_PYTHON_PATH_ENV "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{PYTHONPATH}")
endif(WIN32)

# Documentation generation is available through sphinx - requires all modules
# Make it first as sphinx test is by far the longest test which is nice when testing in parallel
if(GUDHI_PYTHON_DOCUMENTATION_MISSING_MODULES)
  message("++ Python documentation will not be compilable because of the following missing modules: ${GUDHI_PYTHON_DOCUMENTATION_MISSING_MODULES}")
else()
  set (GUDHI_SPHINX_MESSAGE "Generating API documentation with Sphinx in ${CMAKE_CURRENT_BINARY_DIR}/sphinx/")
  # User warning - Sphinx is a static pages generator, and configured to work fine with user_version
  # Images and biblio warnings because not found on developer version
  if (GUDHI_PYTHON_PATH STREQUAL "src/python")
    set (GUDHI_SPHINX_MESSAGE "${GUDHI_SPHINX_MESSAGE} \n WARNING : Sphinx is configured for user version, you run it on developer version. Images and biblio will miss")
  endif()
  # sphinx target requires gudhi.so, because conf.py reads gudhi version from it
  add_custom_target(sphinx
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
      COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
      ${SPHINX_PATH} -b html ${CMAKE_CURRENT_SOURCE_DIR}/doc ${CMAKE_CURRENT_BINARY_DIR}/sphinx
      DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/gudhi.so"
      COMMENT "${GUDHI_SPHINX_MESSAGE}" VERBATIM)
  add_test(NAME sphinx_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
           ${SPHINX_PATH} -b doctest ${CMAKE_CURRENT_SOURCE_DIR}/doc ${CMAKE_CURRENT_BINARY_DIR}/doctest)
  # Set missing or not modules
  set(GUDHI_MODULES ${GUDHI_MODULES} "python-documentation" CACHE INTERNAL "GUDHI_MODULES")
endif()

# Cubical
add_test(NAME periodic_cubical_complex_barcode_persistence_from_perseus_file_example_py_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
         COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
         ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/periodic_cubical_complex_barcode_persistence_from_perseus_file_example.py"
         --no-barcode -f ${CMAKE_SOURCE_DIR}/data/bitmap/CubicalTwoSphere.txt)

add_test(NAME random_cubical_complex_persistence_example_py_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
         COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
         ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/random_cubical_complex_persistence_example.py"
         10 10 10)

add_gudhi_py_test(test_cubical_complex)

# Datasets are fetched for these tests
if(SCIKIT-LEARN_FOUND AND WITH_GUDHI_REMOTE_TEST)
  add_gudhi_py_test(test_sklearn_cubical_persistence)

  add_test(NAME cubical_complex_sklearn_itf_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
           ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/cubical_complex_sklearn_itf.py")
endif()

# Test examples
if (TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  # Bottleneck and Alpha
  add_test(NAME alpha_rips_persistence_bottleneck_distance_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
           ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/alpha_rips_persistence_bottleneck_distance.py"
           -f ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off -t 0.15 -d 3)
endif (TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
if (TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  # Tangential
  add_test(NAME tangential_complex_plain_homology_from_off_file_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
           ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/tangential_complex_plain_homology_from_off_file_example.py"
           --no-diagram -i 2 -f ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off)

  add_gudhi_py_test(test_tangential_complex)

  # Witness complex
  add_test(NAME euclidean_strong_witness_complex_diagram_persistence_from_off_file_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
           ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/euclidean_strong_witness_complex_diagram_persistence_from_off_file_example.py"
           --no-diagram -f ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off -a 1.0 -n 20 -d 2)

  add_test(NAME euclidean_witness_complex_diagram_persistence_from_off_file_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
           ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/euclidean_witness_complex_diagram_persistence_from_off_file_example.py"
           --no-diagram -f ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off -a 1.0 -n 20 -d 2)

endif (TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
if (TARGET CGAL::CGAL)
  # Bottleneck
  add_test(NAME bottleneck_basic_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
           ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/bottleneck_basic_example.py")

  add_gudhi_py_test(test_bottleneck_distance)
endif (TARGET CGAL::CGAL)

# Cover complex
file(COPY ${CMAKE_SOURCE_DIR}/data/points/human.off DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY ${CMAKE_SOURCE_DIR}/data/points/COIL_database/lucky_cat.off DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
file(COPY ${CMAKE_SOURCE_DIR}/data/points/COIL_database/lucky_cat_PCA1 DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
add_test(NAME cover_complex_nerve_example_py_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
         COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
         ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/nerve_of_a_covering.py"
         -f human.off -c 2 -r 10 -g 0.3)

add_test(NAME cover_complex_coordinate_gic_example_py_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
         COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
         ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/coordinate_graph_induced_complex.py"
         -f human.off -c 0 -v)

add_test(NAME cover_complex_functional_gic_example_py_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
         COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
         ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/functional_graph_induced_complex.py"
         -o lucky_cat.off
         -f lucky_cat_PCA1 -v)

add_test(NAME cover_complex_voronoi_gic_example_py_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
         COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
         ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/voronoi_graph_induced_complex.py"
         -f human.off -n 700 -v)

if (TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  # Alpha
  add_test(NAME alpha_complex_from_points_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
           ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/alpha_complex_from_points_example.py")
  add_test(NAME alpha_complex_from_generated_points_on_sphere_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
           ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/alpha_complex_from_generated_points_on_sphere_example.py")
  add_test(NAME alpha_complex_diagram_persistence_from_off_file_example_py_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
         COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
         ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/alpha_complex_diagram_persistence_from_off_file_example.py"
         --no-diagram -f ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off)
  add_gudhi_py_test(test_alpha_complex)
  add_gudhi_py_test(test_delaunay_complex)

  # Cech scikit-learn like interface
  if(SCIKIT-LEARN_FOUND)
    add_test(NAME cech_complex_sklearn_itf_py_test
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
             COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
             ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/cech_complex_sklearn_itf.py" "--no-plot")

    add_gudhi_py_test(test_sklearn_cech_persistence)
  endif()
endif (TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)

if (TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  # Euclidean witness
  add_gudhi_py_test(test_euclidean_witness_complex)

  # Datasets generators
  add_gudhi_py_test(test_datasets_generators) # TODO separate full python datasets generators in another test file independent from CGAL ?

endif (TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)

# Rips
add_test(NAME rips_complex_diagram_persistence_from_distance_matrix_file_example_py_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
         COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
         ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/rips_complex_diagram_persistence_from_distance_matrix_file_example.py"
         --no-diagram -f ${CMAKE_SOURCE_DIR}/data/distance_matrix/lower_triangular_distance_matrix.csv -s , -e 12.0 -d 3)

add_test(NAME rips_complex_diagram_persistence_from_off_file_example_py_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
         COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
         ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/example/rips_complex_diagram_persistence_from_off_file_example.py
         --no-diagram -f ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off  -e 0.25 -d 3)

add_test(NAME rips_complex_from_points_example_py_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
         COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
         ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/example/rips_complex_from_points_example.py)

add_gudhi_py_test(test_rips_complex)

if(SCIKIT-LEARN_FOUND)
  # test_sklearn_rips_persistence is using gudhi.datasets.generators
  if (TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
    add_gudhi_py_test(test_sklearn_rips_persistence)
  endif()

  # Datasets are fetched for these tests
  if(WITH_GUDHI_REMOTE_TEST)
    add_test(NAME rips_complex_sklearn_itf_py_test
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
             COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
             ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/rips_complex_sklearn_itf.py" "--no-plot")
  endif()
endif()

# Simplex tree
add_test(NAME simplex_tree_example_py_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
         COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
         ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/example/simplex_tree_example.py)

add_gudhi_py_test(test_simplex_tree)
add_gudhi_py_test(test_simplex_generators)
add_gudhi_py_test(test_simplex_tree_serialization)

# Edge collapse
if(SCIPY_FOUND)
  add_gudhi_py_test(test_collapse_edges)
endif()

# Subsampling
add_gudhi_py_test(test_subsampling)

# Witness
add_test(NAME witness_complex_from_nearest_landmark_table_py_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
         COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
         ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/example/witness_complex_from_nearest_landmark_table.py)

add_gudhi_py_test(test_witness_complex)

# Reader utils
add_gudhi_py_test(test_reader_utils)
add_gudhi_py_test(test_off)

# Wasserstein
if(POT_FOUND)
  # EagerPy dependency because of enable_autodiff=True
  if(EAGERPY_FOUND)
    add_gudhi_py_test(test_wasserstein_distance)
  endif()
  add_gudhi_py_test(test_wasserstein_barycenter)
  if(TORCH_FOUND AND TENSORFLOW_FOUND AND EAGERPY_FOUND)
    add_gudhi_py_test(test_wasserstein_with_tensors)
  endif()
endif()

# Representations
if(SCIKIT-LEARN_FOUND AND MATPLOTLIB_FOUND AND POT_FOUND AND TARGET CGAL::CGAL)
  add_gudhi_py_test(test_representations)
endif()

# Differentiation
if(TENSORFLOW_FOUND)
  add_gudhi_py_test(test_diff)
endif()

# Perslay
if(TENSORFLOW_FOUND AND SCIKIT-LEARN_FOUND)
  add_gudhi_py_test(test_perslay)
endif()

# Betti curves
if(SCIKIT-LEARN_FOUND AND SCIPY_FOUND)
  add_gudhi_py_test(test_betti_curve_representations)
endif()

# Representations preprocessing
if(SCIKIT-LEARN_FOUND)
  add_gudhi_py_test(test_representations_preprocessing)
endif()

# Time Delay
add_gudhi_py_test(test_time_delay)

# DTM
if(SCIPY_FOUND AND SCIKIT-LEARN_FOUND AND TORCH_FOUND AND HNSWLIB_FOUND AND PYKEOPS_FOUND AND EAGERPY_FOUND)
  add_gudhi_py_test(test_knn)
  add_gudhi_py_test(test_dtm)
endif()

# Tomato
if(SCIPY_FOUND AND SCIKIT-LEARN_FOUND)
  add_gudhi_py_test(test_tomato)
endif()

# Weighted Rips
if(SCIPY_FOUND)
  add_gudhi_py_test(test_weighted_rips_complex)
endif()

# DTM Rips
if(SCIPY_FOUND)
  add_gudhi_py_test(test_dtm_rips_complex)
endif()

# Cover complex
if(SCIKIT-LEARN_FOUND)
  add_gudhi_py_test(test_cover_complex)
endif()

# Fetch remote datasets
if(WITH_GUDHI_REMOTE_TEST)
  add_gudhi_py_test(test_remote_datasets)
endif()

# persistence graphical tools
if(MATPLOTLIB_FOUND)
  add_gudhi_py_test(test_persistence_graphical_tools)
endif()

# Set missing or not modules
set(GUDHI_MODULES ${GUDHI_MODULES} "python" CACHE INTERNAL "GUDHI_MODULES")
