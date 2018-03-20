cmake_minimum_required(VERSION 2.8)
project(Cython)

function( add_gudhi_cython_lib THE_LIB )
  if(EXISTS ${THE_LIB})
    get_filename_component(THE_LIB_FILE_NAME ${THE_LIB} NAME_WE)
    if(WIN32)
      message("++ ${THE_LIB} => THE_LIB_FILE_NAME = ${THE_LIB_FILE_NAME}")
      set(GUDHI_CYTHON_LIBRARIES "${GUDHI_CYTHON_LIBRARIES}'${THE_LIB_FILE_NAME}', " PARENT_SCOPE)
    else(WIN32)
      STRING(REGEX REPLACE "lib" "" UNIX_LIB_FILE_NAME ${THE_LIB_FILE_NAME})
      message("++ ${THE_LIB} => UNIX_LIB_FILE_NAME = ${UNIX_LIB_FILE_NAME}")
      set(GUDHI_CYTHON_LIBRARIES "${GUDHI_CYTHON_LIBRARIES}'${UNIX_LIB_FILE_NAME}', " PARENT_SCOPE)
    endif(WIN32)
  endif(EXISTS ${THE_LIB})
endfunction( add_gudhi_cython_lib )

# THE_TEST is the python test file name (without .py extension) containing tests functions
function( add_gudhi_py_test THE_TEST )
  # use ${PYTHON_EXECUTABLE} -B, otherwise a __pycache__ directory is created in sources by python
  # use py.test no cache provider, otherwise a .cache file is created in sources by py.test
  add_test(NAME ${THE_TEST}_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${PYTHON_EXECUTABLE} -B -m pytest -p no:cacheprovider ${CMAKE_CURRENT_SOURCE_DIR}/test/${THE_TEST}.py)
endfunction( add_gudhi_py_test )


if(CYTHON_FOUND)
  message("++ ${PYTHON_EXECUTABLE} v.${PYTHON_VERSION_STRING} - Cython is ${CYTHON_EXECUTABLE} - Sphinx is ${SPHINX_PATH}")

  set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-DBOOST_RESULT_OF_USE_DECLTYPE', ")
  set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-DBOOST_ALL_NO_LIB', ")
  set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-DBOOST_SYSTEM_NO_DEPRECATED', ")

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

  if (NOT CGAL_CAN_USE_CXX11_THREAD_LOCAL_RESULT)
    add_gudhi_cython_lib(${Boost_THREAD_LIBRARY})
    set(GUDHI_CYTHON_LIBRARY_DIRS "${GUDHI_CYTHON_LIBRARY_DIRS}'${Boost_LIBRARY_DIRS}', ")
  endif()

  # Gudhi and CGAL compilation option
  if(MSVC)
    set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'/fp:strict', ")
  else(MSVC)
    set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-std=c++11', ")
  endif(MSVC)
  if(CMAKE_COMPILER_IS_GNUCXX)
    set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-frounding-math', ")
  endif(CMAKE_COMPILER_IS_GNUCXX)
  if (CMAKE_CXX_COMPILER_ID MATCHES Intel)
    set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-fp-model strict', ")
  endif(CMAKE_CXX_COMPILER_ID MATCHES Intel)
  if (DEBUG_TRACES)
    # For programs to be more verbose
    set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-DDEBUG_TRACES', ")
  endif()

  if (EIGEN3_FOUND)
    # No problem, even if no CGAL found
    set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-DCGAL_EIGEN3_ENABLED', ")
  endif (EIGEN3_FOUND)

  if (NOT CGAL_VERSION VERSION_LESS 4.8.1)
    set(GUDHI_CYTHON_BOTTLENECK_DISTANCE "include '${CMAKE_CURRENT_SOURCE_DIR}/cython/bottleneck_distance.pyx'")
    set(GUDHI_CYTHON_NERVE_GIC "include '${CMAKE_CURRENT_SOURCE_DIR}/cython/nerve_gic.pyx'")
  endif (NOT CGAL_VERSION VERSION_LESS 4.8.1)
  if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
    set(GUDHI_CYTHON_SUBSAMPLING "include '${CMAKE_CURRENT_SOURCE_DIR}/cython/subsampling.pyx'")
    set(GUDHI_CYTHON_TANGENTIAL_COMPLEX "include '${CMAKE_CURRENT_SOURCE_DIR}/cython/tangential_complex.pyx'")
  endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
  if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.7.0)
    set(GUDHI_CYTHON_ALPHA_COMPLEX "include '${CMAKE_CURRENT_SOURCE_DIR}/cython/alpha_complex.pyx'")
  endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.7.0)
  if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
    set(GUDHI_CYTHON_EUCLIDEAN_WITNESS_COMPLEX
        "include '${CMAKE_CURRENT_SOURCE_DIR}/cython/euclidean_witness_complex.pyx'\ninclude '${CMAKE_CURRENT_SOURCE_DIR}/cython/euclidean_strong_witness_complex.pyx'\n")
  endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)

  if(CGAL_FOUND)
    # Add CGAL compilation args
    if(CGAL_HEADER_ONLY)
      set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-DCGAL_HEADER_ONLY', ")
    else(CGAL_HEADER_ONLY)
      add_gudhi_cython_lib(${CGAL_LIBRARIES})
      set(GUDHI_CYTHON_LIBRARY_DIRS "${GUDHI_CYTHON_LIBRARY_DIRS}'${CGAL_LIBRARIES_DIR}', ")
      # If CGAL is not header only, CGAL library may link with boost system,
      add_gudhi_cython_lib(${Boost_SYSTEM_LIBRARY})
      set(GUDHI_CYTHON_LIBRARY_DIRS "${GUDHI_CYTHON_LIBRARY_DIRS}'${Boost_LIBRARY_DIRS}', ")
    endif(CGAL_HEADER_ONLY)
    # GMP and GMPXX are not required, but if present, CGAL will link with them.
    if(GMP_FOUND)
      set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-DCGAL_USE_GMP', ")
      add_gudhi_cython_lib(${GMP_LIBRARIES})
      set(GUDHI_CYTHON_LIBRARY_DIRS "${GUDHI_CYTHON_LIBRARY_DIRS}'${GMP_LIBRARIES_DIR}', ")
      if(GMPXX_FOUND)
        set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-DCGAL_USE_GMPXX', ")
        add_gudhi_cython_lib(${GMPXX_LIBRARIES})
        set(GUDHI_CYTHON_LIBRARY_DIRS "${GUDHI_CYTHON_LIBRARY_DIRS}'${GMPXX_LIBRARIES_DIR}', ")
      endif(GMPXX_FOUND)
    endif(GMP_FOUND)
  endif(CGAL_FOUND)

  # Specific for Mac
  if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-mmacosx-version-min=10.12', ")
      set(GUDHI_CYTHON_EXTRA_LINK_ARGS "${GUDHI_CYTHON_EXTRA_LINK_ARGS}'-mmacosx-version-min=10.12', ")
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

  # Loop on INCLUDE_DIRECTORIES PROPERTY
  get_property(GUDHI_INCLUDE_DIRECTORIES DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
  foreach(GUDHI_INCLUDE_DIRECTORY ${GUDHI_INCLUDE_DIRECTORIES})
    set(GUDHI_CYTHON_INCLUDE_DIRS "${GUDHI_CYTHON_INCLUDE_DIRS}'${GUDHI_INCLUDE_DIRECTORY}', ")
  endforeach()
  set(GUDHI_CYTHON_INCLUDE_DIRS "${GUDHI_CYTHON_INCLUDE_DIRS}'${CMAKE_SOURCE_DIR}/${GUDHI_CYTHON_PATH}/include', ")

  if (TBB_FOUND)
    set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-DGUDHI_USE_TBB', ")
    add_gudhi_cython_lib(${TBB_RELEASE_LIBRARY})
    add_gudhi_cython_lib(${TBB_MALLOC_RELEASE_LIBRARY})
    set(GUDHI_CYTHON_LIBRARY_DIRS "${GUDHI_CYTHON_LIBRARY_DIRS}'${TBB_LIBRARY_DIRS}', ")
    set(GUDHI_CYTHON_INCLUDE_DIRS "${GUDHI_CYTHON_INCLUDE_DIRS}'${TBB_INCLUDE_DIRS}', ")
  endif()

  if(UNIX)
    set( GUDHI_CYTHON_RUNTIME_LIBRARY_DIRS "${GUDHI_CYTHON_LIBRARY_DIRS}")
  endif(UNIX)

  # Generate setup.py file to cythonize Gudhi - This file must be named setup.py by convention
  configure_file(setup.py.in "${CMAKE_CURRENT_BINARY_DIR}/setup.py" @ONLY)
  # Generate gudhi.pyx - Gudhi cython file
  configure_file(gudhi.pyx.in "${CMAKE_CURRENT_BINARY_DIR}/gudhi.pyx" @ONLY)

  add_custom_command(
      OUTPUT gudhi.so
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/setup.py" "build_ext" "--inplace")

  add_custom_target(cython ALL DEPENDS gudhi.so
                    COMMENT "Do not forget to add ${CMAKE_CURRENT_BINARY_DIR}/ to your PYTHONPATH before using examples or tests")

  # For installation purpose
  # TODO(VR) : files matching pattern mechanism is copying all cython directory
  install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" DESTINATION "${PYTHON_SITE_PACKAGES}/" FILES_MATCHING
      PATTERN "*.so"
      PATTERN "*.dylib"
      PATTERN "*.pyd")

  # Test examples
  if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
    # Bottleneck and Alpha
    add_test(NAME alpha_rips_persistence_bottleneck_distance_py_test
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
             COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
             ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/alpha_rips_persistence_bottleneck_distance.py"
             -f ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off -t 0.15 -d 3)

    # Tangential
    add_test(NAME tangential_complex_plain_homology_from_off_file_example_py_test
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
             COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
             ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/tangential_complex_plain_homology_from_off_file_example.py"
             --no-diagram -f ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off)

    add_gudhi_py_test(test_tangential_complex)

    # Witness complex AND Subsampling
    add_test(NAME euclidean_strong_witness_complex_diagram_persistence_from_off_file_example_py_test
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
             COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
             ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/euclidean_strong_witness_complex_diagram_persistence_from_off_file_example.py"
             --no-diagram -f ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off -a 1.0 -n 20 -d 2)

    add_test(NAME euclidean_witness_complex_diagram_persistence_from_off_file_example_py_test
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
             COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
             ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/euclidean_witness_complex_diagram_persistence_from_off_file_example.py"
             --no-diagram -f ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off -a 1.0 -n 20 -d 2)

    # Subsampling
    add_gudhi_py_test(test_subsampling)

  endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
  if (NOT CGAL_VERSION VERSION_LESS 4.8.1)
    # Bottleneck
    add_test(NAME bottleneck_basic_example_py_test
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
             COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
             ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/bottleneck_basic_example.py")

    add_gudhi_py_test(test_bottleneck_distance)

    # Cover complex
    add_test(NAME cover_complex_nerve_example_py_test
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
             COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
             ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/nerve_of_a_covering.py"
             -f ${CMAKE_SOURCE_DIR}/data/points/human.off -c 2 -r 10 -g 0.3)

    add_test(NAME cover_complex_coordinate_gic_example_py_test
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
             COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
             ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/coordinate_graph_induced_complex.py"
             -f ${CMAKE_SOURCE_DIR}/data/points/human.off -c 0 -v)

    add_test(NAME cover_complex_functional_gic_example_py_test
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
             COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
             ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/functional_graph_induced_complex.py"
             -o ${CMAKE_SOURCE_DIR}/data/points/COIL_database/lucky_cat.off
             -f ${CMAKE_SOURCE_DIR}/data/points/COIL_database/lucky_cat_PCA1 -v)

    add_test(NAME cover_complex_voronoi_gic_example_py_test
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
             COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
             ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/voronoi_graph_induced_complex.py"
             -f ${CMAKE_SOURCE_DIR}/data/points/human.off -n 700 -v)

    add_gudhi_py_test(test_cover_complex)
  endif (NOT CGAL_VERSION VERSION_LESS 4.8.1)

  if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.7.0)
    # Alpha
    add_test(NAME alpha_complex_from_points_example_py_test
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
             COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
             ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/alpha_complex_from_points_example.py")

    add_test(NAME alpha_complex_diagram_persistence_from_off_file_example_py_test
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
             COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
             ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/alpha_complex_diagram_persistence_from_off_file_example.py"
             --no-diagram -f ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off -a 0.6)

    add_gudhi_py_test(test_alpha_complex)

  endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.7.0)

  if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
    # Euclidean witness
    add_gudhi_py_test(test_euclidean_witness_complex)

  endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)

  # Cubical
  add_test(NAME periodic_cubical_complex_barcode_persistence_from_perseus_file_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
           ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/periodic_cubical_complex_barcode_persistence_from_perseus_file_example.py"
           --no-barcode -f ${CMAKE_SOURCE_DIR}/data/bitmap/CubicalTwoSphere.txt)

  add_test(NAME random_cubical_complex_persistence_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
           ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/random_cubical_complex_persistence_example.py"
           10 10 10)

  add_gudhi_py_test(test_cubical_complex)

  # Rips
  add_test(NAME rips_complex_diagram_persistence_from_distance_matrix_file_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
           ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/example/rips_complex_diagram_persistence_from_distance_matrix_file_example.py"
           --no-diagram -f ${CMAKE_SOURCE_DIR}/data/distance_matrix/lower_triangular_distance_matrix.csv -e 12.0 -d 3)

  add_test(NAME rips_complex_diagram_persistence_from_off_file_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
           ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/example/rips_complex_diagram_persistence_from_off_file_example.py
           --no-diagram -f ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off  -e 0.25 -d 3)

  add_test(NAME rips_complex_from_points_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
           ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/example/rips_complex_from_points_example.py)

  add_gudhi_py_test(test_rips_complex)

  # Simplex tree
  add_test(NAME simplex_tree_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
           ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/example/simplex_tree_example.py)

  add_gudhi_py_test(test_simplex_tree)

  # Witness
  add_test(NAME witness_complex_from_nearest_landmark_table_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
           ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/example/witness_complex_from_nearest_landmark_table.py)

  add_gudhi_py_test(test_witness_complex)

  # Reader utils
  add_gudhi_py_test(test_reader_utils)

  # Documentation generation is available through sphinx - requires all modules
  if(SPHINX_PATH AND NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
    set (GUDHI_SPHINX_MESSAGE "Generating API documentation with Sphinx in ${CMAKE_CURRENT_BINARY_DIR}/sphinx/")
    # User warning - Sphinx is a static pages generator, and configured to work fine with user_version
    # Images and biblio warnings because not found on developper version
    if (GUDHI_CYTHON_PATH STREQUAL "src/cython")
      set (GUDHI_SPHINX_MESSAGE "${GUDHI_SPHINX_MESSAGE} \n WARNING : Sphinx is configured for user version, you run it on developper version. Images and biblio will miss")
    endif()
    # sphinx target requires gudhi.so, because conf.py reads gudhi version from it
    add_custom_target(sphinx
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc
        COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
        ${SPHINX_PATH} -b html ${CMAKE_CURRENT_SOURCE_DIR}/doc ${CMAKE_CURRENT_BINARY_DIR}/sphinx
        DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/gudhi.so"
        COMMENT "${GUDHI_SPHINX_MESSAGE}" VERBATIM)

    add_test(NAME sphinx_py_test
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
             COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
             ${SPHINX_PATH} -b doctest ${CMAKE_CURRENT_SOURCE_DIR}/doc ${CMAKE_CURRENT_BINARY_DIR}/doctest)

  endif(SPHINX_PATH AND NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.8.1)
endif(CYTHON_FOUND)
