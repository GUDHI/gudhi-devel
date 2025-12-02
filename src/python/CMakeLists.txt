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
if(NOT NANOBIND_FOUND)
  message("++ Python module will not be compiled because Nanobind was not found")
  set(GUDHI_MISSING_MODULES ${GUDHI_MISSING_MODULES} "python" CACHE INTERNAL "GUDHI_MISSING_MODULES")
  return()
endif()

function(nanobind_add_module_with_targets)
    if (ARGC LESS 2)
        message (FATAL_ERROR "add_executable_with_targets requires at least 2 arguments.")
    endif()

    list(POP_FRONT ARGN EXECUTABLE_NAME EXECUTABLE_SOURCE)
    # while loop also ends if EXECUTABLE_SOURCE is empty
    while(NOT ${EXECUTABLE_SOURCE} STREQUAL "TARGETS")
      list(APPEND EXE_SOURCES "${EXECUTABLE_SOURCE}")
      list(POP_FRONT ARGN EXECUTABLE_SOURCE)
    endwhile()

    # Do not add_executable if one of the target is not here, except for TBB that is optional
    foreach(USER_TARGET IN LISTS ARGN)
        if(NOT TARGET ${USER_TARGET})
            if(NOT ${USER_TARGET} STREQUAL "TBB::tbb")
                return()
            endif()
        endif()
    endforeach()
    message(DEBUG "nanobind_add_module(${EXECUTABLE_NAME} ${EXE_SOURCES}) + target_link_libraries(${ARGN})")
    # Maintain a list of nanobind modules - Mainly for sphinx purposes to know when the documentation is "buildable"
    set(GUDHI_NANOBIND_MODULES ${GUDHI_NANOBIND_MODULES} ${EXECUTABLE_NAME} CACHE INTERNAL "GUDHI_NANOBIND_MODULES")

    # cf. https://nanobind.readthedocs.io/en/latest/api_cmake.html#high-level-interface
    # could be:
    # nanobind_add_module(${EXECUTABLE_NAME} NB_STATIC NOMINSIZE ${EXE_SOURCES})
    # but
    # cf. https://nanobind.readthedocs.io/en/latest/api_cmake.html#low-level-interface
    # To avoid nanobind_set_visibility(my_ext) for the "random" package
    nanobind_build_library(nanobind STATIC)
    add_library(${EXECUTABLE_NAME} MODULE ${EXE_SOURCES})
    target_link_libraries(${EXECUTABLE_NAME} PRIVATE nanobind)
    # .. strip unneeded symbols and debug info from the binary (only active in release builds)
    nanobind_strip(${EXECUTABLE_NAME})
    # .. disable the stack protector
    nanobind_disable_stack_protector(${EXECUTABLE_NAME})
    # .. set the Python extension suffix
    nanobind_extension(${EXECUTABLE_NAME})
    # .. set important compilation flags
    nanobind_compile_options(${EXECUTABLE_NAME})
    # .. set important linker flags
    nanobind_link_options(${EXECUTABLE_NAME})
    # Statically link against libstdc++/libgcc when targeting musllinux
    nanobind_musl_static_libcpp(${EXECUTABLE_NAME})

    if(MSVC)
      target_compile_options(${EXECUTABLE_NAME} PRIVATE "-DNOMINMAX")
    endif()
    foreach(USER_TARGET IN LISTS ARGN)
        # TARGET_NAME_IF_EXISTS is specific to TBB case (optional)
        target_link_libraries(${EXECUTABLE_NAME} PRIVATE $<TARGET_NAME_IF_EXISTS:${USER_TARGET}>)
    endforeach()
    if(MSVC)
        # cf. https://nanobind.readthedocs.io/en/latest/faq.html#importing-fails-due-to-missing-lib-nanobind-dylib-so-dll
        message(DEBUG "${EXECUTABLE_NAME} set LIBRARY_OUTPUT_DIRECTORIES to ${CMAKE_CURRENT_BINARY_DIR}")
        set_target_properties(${EXECUTABLE_NAME} PROPERTIES
                              LIBRARY_OUTPUT_DIRECTORY                ${CMAKE_CURRENT_BINARY_DIR}
                              LIBRARY_OUTPUT_DIRECTORY_RELEASE        ${CMAKE_CURRENT_BINARY_DIR}
                              LIBRARY_OUTPUT_DIRECTORY_DEBUG          ${CMAKE_CURRENT_BINARY_DIR}
                              LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_CURRENT_BINARY_DIR}
                              LIBRARY_OUTPUT_DIRECTORY_MINSIZEREL     ${CMAKE_CURRENT_BINARY_DIR})
    endif(MSVC)
endfunction()

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
add_gudhi_debug_info("Nanobind version ${NANOBIND_VERSION}")
add_gudhi_debug_info("Numpy version ${NUMPY_VERSION}")

# Python modules for __all__ list in __init__.py
set(GUDHI_PYTHON_MODULES "${GUDHI_PYTHON_MODULES}'bottleneck', ")
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
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'hera', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'clustering', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'datasets', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'representations', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'sklearn', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'tensorflow', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'wasserstein', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'point_cloud', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'weighted_rips_complex', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'dtm_rips_complex', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'cover_complex', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'flag_filtration', ")
set(GUDHI_PYTHON_MODULES_EXTRA "${GUDHI_PYTHON_MODULES_EXTRA}'random', ")


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
else()
  disable_python_documentation("tensorflow")
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

if (TARGET Eigen3::Eigen)
  add_gudhi_debug_info("Eigen3 version ${EIGEN3_VERSION_STRING}")
endif (TARGET Eigen3::Eigen)

add_gudhi_debug_info("Boost version ${Boost_VERSION}")
if(TARGET CGAL::CGAL)
  add_gudhi_debug_info("CGAL version ${CGAL_VERSION}")
  if(GMP_FOUND)
    add_gudhi_debug_info("GMP_LIBRARIES = ${GMP_LIBRARIES}")
    if(GMPXX_FOUND)
      add_gudhi_debug_info("GMPXX_LIBRARIES = ${GMPXX_LIBRARIES}")
    endif(GMPXX_FOUND)
  endif(GMP_FOUND)
  if(MPFR_FOUND)
    add_gudhi_debug_info("MPFR_LIBRARIES = ${MPFR_LIBRARIES}")
  endif(MPFR_FOUND)
endif(TARGET CGAL::CGAL)

if (TBB_FOUND AND WITH_GUDHI_USE_TBB)
  add_gudhi_debug_info("TBB version ${TBB_VERSION} found and used")
  add_gudhi_debug_info("TBB_LIBRARY = ${TBB_LIBRARY}")
  add_gudhi_debug_info("TBB_MALLOC_LIBRARY = ${TBB_MALLOC_LIBRARY}")
endif()

# configure all python files from sources to binary dir forces CMake to copy them when modified
file(GLOB_RECURSE GUDHI_PYTHON_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "gudhi/*.py")
foreach(GUDHI_PYTHON_FILE ${GUDHI_PYTHON_FILES})
  if(NOT GUDHI_PYTHON_FILE MATCHES "gudhi/__init__.py")
    configure_file("${GUDHI_PYTHON_FILE}" "${CMAKE_CURRENT_BINARY_DIR}/${GUDHI_PYTHON_FILE}" COPYONLY)
    message(DEBUG "Python configure_file ${GUDHI_PYTHON_FILE}")
  endif()
endforeach()

# Generate gudhi/__init__.py
configure_file("gudhi/__init__.py.in" "${CMAKE_CURRENT_BINARY_DIR}/gudhi/__init__.py" @ONLY)

# Detect the installed nanobind package and import it into CMake
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

# Add src/python/include directory
include_directories(include)
add_subdirectory("gudhi")
add_subdirectory("gudhi/clustering")
add_subdirectory("gudhi/datasets/generators")

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
  # sphinx target requires python modules to be built, because conf.py reads gudhi version from it
  add_custom_target(sphinx
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
      COMMAND ${CMAKE_COMMAND} -E env "${GUDHI_PYTHON_PATH_ENV}"
      ${SPHINX_PATH} -b html ${CMAKE_CURRENT_SOURCE_DIR}/doc ${CMAKE_CURRENT_BINARY_DIR}/sphinx
      DEPENDS "${GUDHI_NANOBIND_MODULES}"
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
add_gudhi_py_test(test_cubical_persistence_low_dim)

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
if (TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  add_gudhi_py_test(test_subsampling)
endif()

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

# random
if (TARGET CGAL::CGAL)
  add_gudhi_py_test(test_random)
endif()

# Set missing or not modules
set(GUDHI_MODULES ${GUDHI_MODULES} "python" CACHE INTERNAL "GUDHI_MODULES")

install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/gudhi" DESTINATION ".")
