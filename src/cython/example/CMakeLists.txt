cmake_minimum_required(VERSION 2.8)
project(Cython_examples)

if(CYTHON_PATH)
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
endif(CYTHON_PATH)
