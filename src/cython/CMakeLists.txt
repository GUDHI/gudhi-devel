cmake_minimum_required(VERSION 2.8)
project(Cython)

if(CYTHON_PATH)
  message("${PYTHON_EXECUTABLE} v.${PYTHON_VERSION_STRING} - Cython is ${CYTHON_PATH} - py.test is ${PYTEST_PATH} - Sphinx is ${SPHINX_PATH}")

  set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-DBOOST_RESULT_OF_USE_DECLTYPE', ")
  set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-DBOOST_ALL_NO_LIB', ")
  set(GUDHI_CYTHON_LIBRARY_DIRS "${GUDHI_CYTHON_LIBRARY_DIRS}'${Boost_LIBRARY_DIRS}', ")
  if(WIN32)
    set( returnValue "" )
    find_the_lib (${returnValue} ${Boost_SYSTEM_LIBRARY})
    set(BOOST_SYSTEM_LIB_NAME ${returnValue})
  else()
      set(BOOST_SYSTEM_LIB_NAME "boost_system")
      endif()
  set(GUDHI_CYTHON_LIBRARIES "${GUDHI_CYTHON_LIBRARIES}'${BOOST_SYSTEM_LIB_NAME}', ")

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
  if (DEBUG_TRACES)
    # For programs to be more verbose
    set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-DDEBUG_TRACES', ")
  endif()

  if (EIGEN3_FOUND)
    # No problem, even if no CGAL found
    set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-DCGAL_EIGEN3_ENABLED', ")
  endif (EIGEN3_FOUND)

  # Copy recursively include, cython, example, doc and test repositories before packages finding
  # Some tests and doc files are removed in case some packages are not found
  file(COPY include DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY cython DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY example DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY test DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY doc DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  # Developper version for doc images
  file(GLOB GUDHI_DEV_DOC_IMAGES "${CMAKE_SOURCE_DIR}/src/*/doc/*.png")
  file(COPY ${GUDHI_DEV_DOC_IMAGES} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/doc/img")
  file(GLOB GUDHI_DEV_DOC_IMAGES "${CMAKE_SOURCE_DIR}/src/*/doc/*.svg")
  file(COPY ${GUDHI_DEV_DOC_IMAGES} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/doc/img")
  # User version for doc images
  file(GLOB GUDHI_USER_DOC_IMAGES "${CMAKE_SOURCE_DIR}/doc/*/*.png")
  file(COPY ${GUDHI_USER_DOC_IMAGES} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/doc/img")
  file(GLOB GUDHI_USER_DOC_IMAGES "${CMAKE_SOURCE_DIR}/doc/*/*.svg")
  file(COPY ${GUDHI_USER_DOC_IMAGES} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/doc/img")
  # Biblio
  file(GLOB GUDHI_BIB_FILES "${CMAKE_SOURCE_DIR}/biblio/*.bib")
  file(COPY ${GUDHI_BIB_FILES} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/doc/")
  # Cubical complex perseus doc example
  file(GLOB GUDHI_CUBICAL_PERSEUS_FILES "${CMAKE_SOURCE_DIR}/data/bitmap/*cubicalcomplexdoc.txt")
  file(COPY ${GUDHI_CUBICAL_PERSEUS_FILES} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/doc/")
  file(COPY "${CMAKE_SOURCE_DIR}/data/points/alphacomplexdoc.off" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/doc/")
  file(COPY "${CMAKE_SOURCE_DIR}/data/distance_matrix/full_square_distance_matrix.csv" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/doc/")
  # Persistence graphical tools examples
  file(COPY "${CMAKE_SOURCE_DIR}/data/bitmap/3d_torus.txt" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/doc/")
  file(COPY "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/doc/")

  if (NOT CGAL_VERSION VERSION_LESS 4.8.1)
    # If CGAL_VERSION >= 4.8.1, include subsampling
    set(GUDHI_CYTHON_SUBSAMPLING "include 'cython/subsampling.pyx'")
    set(GUDHI_CYTHON_TANGENTIAL_COMPLEX "include 'cython/tangential_complex.pyx'")
    set(GUDHI_CYTHON_BOTTLENECK_DISTANCE "include 'cython/bottleneck_distance.pyx'")
  else (NOT CGAL_VERSION VERSION_LESS 4.8.1)
    # Remove subsampling unitary tests
    file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/test/test_subsampling.py)
    file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/doc/subsampling_ref.rst")
    file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/doc/subsampling_sum.rst")
    file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/doc/subsampling_user.rst")
    # Remove tangential complex and bottleneck unitary tests
    file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/test/test_tangential_complex.py)
    file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/test/test_bottleneck_distance.py)
    file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/doc/bottleneck_distance_ref.rst")
    file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/doc/bottleneck_distance_sum.rst")
    file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/doc/bottleneck_distance_user.rst")
    file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/doc/tangential_complex_ref.rst")
    file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/doc/tangential_complex_sum.rst")
    file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/doc/tangential_complex_user.rst")
  endif (NOT CGAL_VERSION VERSION_LESS 4.8.1)
  if (NOT CGAL_VERSION VERSION_LESS 4.7.0)
    # If CGAL_VERSION >= 4.7.0, include alpha
    set(GUDHI_CYTHON_ALPHA_COMPLEX "include 'cython/alpha_complex.pyx'")
  else (NOT CGAL_VERSION VERSION_LESS 4.7.0)
    # Remove alpha complex unitary tests
    file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/test/test_alpha_complex.py)
    file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/doc/alpha_complex_ref.rst")
    file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/doc/alpha_complex_sum.rst")
    file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/doc/alpha_complex_user.rst")
  endif (NOT CGAL_VERSION VERSION_LESS 4.7.0)
  if (NOT CGAL_VERSION VERSION_LESS 4.6.0)
    # If CGAL_VERSION >= 4.6.0, include euclidean versions of witness complex
    set(GUDHI_CYTHON_EUCLIDEAN_WITNESS_COMPLEX
        "include 'cython/euclidean_witness_complex.pyx'\ninclude 'cython/euclidean_strong_witness_complex.pyx'\n")
  else (NOT CGAL_VERSION VERSION_LESS 4.6.0)
    # Remove alpha complex unitary tests
    file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/test/test_euclidean_witness_complex.py)
    file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/doc/euclidean_witness_complex_ref.rst")
    file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/doc/euclidean_strong_witness_complex_ref.rst")
  endif (NOT CGAL_VERSION VERSION_LESS 4.6.0)

  if(CGAL_FOUND)
    # Add CGAL compilation args
    if(CGAL_HEADER_ONLY)
      set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-DCGAL_HEADER_ONLY', ")
    else(CGAL_HEADER_ONLY)
      if(WIN32)
        set(GUDHI_CYTHON_LIBRARIES "${GUDHI_CYTHON_LIBRARIES}'CGAL-vc140-mt-4.7', ")
      else(WIN32)
        set(GUDHI_CYTHON_LIBRARIES "${GUDHI_CYTHON_LIBRARIES}'CGAL', ")
      endif(WIN32)
      set(GUDHI_CYTHON_LIBRARY_DIRS "${GUDHI_CYTHON_LIBRARY_DIRS}'${CGAL_LIBRARIES_DIR}', ")
    endif(CGAL_HEADER_ONLY)
    # GMP and GMPXX are not required, but if present, CGAL will link with them.
    if(GMP_FOUND)
      set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-DCGAL_USE_GMP', ")
      if(WIN32)
        set(GUDHI_CYTHON_LIBRARIES "${GUDHI_CYTHON_LIBRARIES}'libgmp-10', ")
      else(WIN32)
        set(GUDHI_CYTHON_LIBRARIES "${GUDHI_CYTHON_LIBRARIES}'gmp', ")
      endif(WIN32)
        set(GUDHI_CYTHON_LIBRARY_DIRS "${GUDHI_CYTHON_LIBRARY_DIRS}'${GMP_LIBRARIES_DIR}', ")
      if(GMPXX_FOUND)
        set(GUDHI_CYTHON_EXTRA_COMPILE_ARGS "${GUDHI_CYTHON_EXTRA_COMPILE_ARGS}'-DCGAL_USE_GMPXX', ")
        set(GUDHI_CYTHON_LIBRARIES "${GUDHI_CYTHON_LIBRARIES}'gmpxx', ")
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
    set(GUDHI_CYTHON_LIBRARIES "${GUDHI_CYTHON_LIBRARIES}'tbb', 'tbbmalloc', ")
    set(GUDHI_CYTHON_LIBRARY_DIRS "${GUDHI_CYTHON_LIBRARY_DIRS}'${TBB_LIBRARY_DIRS}', ")
    set(GUDHI_CYTHON_INCLUDE_DIRS "${GUDHI_CYTHON_INCLUDE_DIRS}'${TBB_INCLUDE_DIRS}', ")
  endif()

  # set sphinx-build in make files
  configure_file(doc/Makefile.in "${CMAKE_CURRENT_BINARY_DIR}/doc/Makefile" @ONLY)
  configure_file(doc/make.bat.in "${CMAKE_CURRENT_BINARY_DIR}/doc/make.bat" @ONLY)

  # Generate cythonize_gudhi.py file to cythonize Gudhi
  configure_file(cythonize_gudhi.py.in "${CMAKE_CURRENT_BINARY_DIR}/cythonize_gudhi.py" @ONLY)
  # Generate gudhi.pyx - Gudhi cython file
  configure_file(gudhi.pyx.in "${CMAKE_CURRENT_BINARY_DIR}/gudhi.pyx" @ONLY)

  add_custom_command(
      OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/gudhi.so"
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/cythonize_gudhi.py" "build_ext" "--inplace")

  add_custom_target(cython ALL DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/gudhi.so"
                    COMMENT "Do not forget to add ${CMAKE_CURRENT_BINARY_DIR}/ to your PYTHONPATH before using examples or tests")

  if(UNIX)
    set( ENV{PYTHONPATH} $ENV{PYTHONPATH}:${CMAKE_CURRENT_BINARY_DIR}/ )
  endif(UNIX)

  # Test examples
  add_test(NAME alpha_complex_from_points_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/example/alpha_complex_from_points_example.py")
  set_tests_properties(alpha_complex_from_points_example_py_test PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

  add_test(NAME alpha_complex_diagram_persistence_from_off_file_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/example/alpha_complex_diagram_persistence_from_off_file_example.py"
           --no-diagram -f ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off -a 0.6)
  set_tests_properties(alpha_complex_diagram_persistence_from_off_file_example_py_test PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

  add_test(NAME alpha_rips_persistence_bottleneck_distance_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/example/alpha_rips_persistence_bottleneck_distance.py"
           -f ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off -t 0.15 -d 3)
  set_tests_properties(alpha_rips_persistence_bottleneck_distance_py_test PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

  add_test(NAME bottleneck_basic_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/example/bottleneck_basic_example.py")
  set_tests_properties(bottleneck_basic_example_py_test PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

  add_test(NAME euclidean_strong_witness_complex_diagram_persistence_from_off_file_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/example/euclidean_strong_witness_complex_diagram_persistence_from_off_file_example.py"
           --no-diagram -f ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off -a 1.0 -n 20 -d 2)
  set_tests_properties(euclidean_strong_witness_complex_diagram_persistence_from_off_file_example_py_test PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

  add_test(NAME euclidean_witness_complex_diagram_persistence_from_off_file_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/example/euclidean_witness_complex_diagram_persistence_from_off_file_example.py"
           --no-diagram -f ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off -a 1.0 -n 20 -d 2)
  set_tests_properties(euclidean_witness_complex_diagram_persistence_from_off_file_example_py_test PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

  add_test(NAME periodic_cubical_complex_barcode_persistence_from_perseus_file_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/example/periodic_cubical_complex_barcode_persistence_from_perseus_file_example.py"
           --no-barcode -f ${CMAKE_SOURCE_DIR}/data/bitmap/CubicalTwoSphere.txt)
  set_tests_properties(periodic_cubical_complex_barcode_persistence_from_perseus_file_example_py_test PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

  add_test(NAME random_cubical_complex_persistence_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/example/random_cubical_complex_persistence_example.py"
           10 10 10)
  set_tests_properties(random_cubical_complex_persistence_example_py_test PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

  add_test(NAME rips_complex_diagram_persistence_from_distance_matrix_file_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/example/rips_complex_diagram_persistence_from_distance_matrix_file_example.py"
           --no-diagram -f ${CMAKE_SOURCE_DIR}/data/distance_matrix/lower_triangular_distance_matrix.csv -e 12.0 -d 3)
  set_tests_properties(rips_complex_diagram_persistence_from_distance_matrix_file_example_py_test PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

  add_test(NAME rips_complex_diagram_persistence_from_off_file_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/example/rips_complex_diagram_persistence_from_off_file_example.py
           --no-diagram -f ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off  -e 0.25 -d 3)
  set_tests_properties(rips_complex_diagram_persistence_from_off_file_example_py_test PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

  add_test(NAME rips_complex_from_points_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/example/rips_complex_from_points_example.py)
  set_tests_properties(rips_complex_from_points_example_py_test PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

  add_test(NAME simplex_tree_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/example/simplex_tree_example.py)
  set_tests_properties(simplex_tree_example_py_test PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

  add_test(NAME tangential_complex_plain_homology_from_off_file_example_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/example/tangential_complex_plain_homology_from_off_file_example.py"
           --no-diagram -f ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off)
  set_tests_properties(tangential_complex_plain_homology_from_off_file_example_py_test PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

  add_test(NAME witness_complex_from_nearest_landmark_table_py_test
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
           COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/example/witness_complex_from_nearest_landmark_table.py)
  set_tests_properties(witness_complex_from_nearest_landmark_table_py_test PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

  # Unitary tests are available through py.test
  if(PYTEST_PATH)
    add_test(
        NAME gudhi_cython_py_test
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND ${PYTHON_EXECUTABLE} "${PYTEST_PATH}")
    set_tests_properties(gudhi_cython_py_test PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")
  endif(PYTEST_PATH)

  # Documentation generation is available through sphinx
  if(SPHINX_PATH)
    if (UNIX)
      add_custom_target(sphinx
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
          DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/gudhi.so"
          COMMAND make html doctest)
    else (UNIX)
      add_custom_target(sphinx
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
          COMMAND make.bat html doctest)
    endif (UNIX)
  endif(SPHINX_PATH)
endif(CYTHON_PATH)
