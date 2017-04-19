cmake_minimum_required(VERSION 2.6)
project(Witness_complex_examples)

add_executable ( Witness_complex_example_nearest_landmark_table example_nearest_landmark_table.cpp )
target_link_libraries(Witness_complex_example_nearest_landmark_table ${Boost_SYSTEM_LIBRARY})
if (TBB_FOUND)
  target_link_libraries(Witness_complex_example_nearest_landmark_table ${TBB_LIBRARIES})
endif()
add_test(NAME Witness_complex_example_nearest_landmark_table
    COMMAND $<TARGET_FILE:Witness_complex_example_nearest_landmark_table>)

# CGAL and Eigen3 are required for Euclidean version of Witness
if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
  add_executable( Witness_complex_example_off example_witness_complex_off.cpp )
  target_link_libraries(Witness_complex_example_off ${Boost_SYSTEM_LIBRARY})
  add_executable( Witness_complex_example_strong_off example_strong_witness_complex_off.cpp )
  target_link_libraries(Witness_complex_example_strong_off ${Boost_SYSTEM_LIBRARY})
  add_executable ( Witness_complex_example_sphere example_witness_complex_sphere.cpp )
  target_link_libraries(Witness_complex_example_sphere ${Boost_SYSTEM_LIBRARY})

  add_executable ( Witness_complex_example_witness_persistence example_witness_complex_persistence.cpp )
  target_link_libraries(Witness_complex_example_witness_persistence ${Boost_SYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})

  add_executable ( Witness_complex_example_strong_witness_persistence example_strong_witness_persistence.cpp )
  target_link_libraries(Witness_complex_example_strong_witness_persistence ${Boost_SYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})

  if (TBB_FOUND)
    target_link_libraries(Witness_complex_example_witness_persistence ${TBB_LIBRARIES})
    target_link_libraries(Witness_complex_example_strong_witness_persistence ${TBB_LIBRARIES})
  endif()

  add_test(NAME Witness_complex_example_off_test_torus
      COMMAND $<TARGET_FILE:Witness_complex_example_off>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "20" "1.0" "3")
  add_test(NAME Witness_complex_example_strong_off_test_torus
      COMMAND $<TARGET_FILE:Witness_complex_example_strong_off>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "20" "1.0" "3")
  add_test(NAME Witness_complex_example_test_sphere_10
      COMMAND $<TARGET_FILE:Witness_complex_example_sphere> "10")
  add_test(NAME Witness_complex_example_test_torus_persistence
      COMMAND $<TARGET_FILE:Witness_complex_example_witness_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-l" "20" "-a" "0.5")
  add_test(NAME Witness_complex_example_strong_test_torus_persistence
      COMMAND $<TARGET_FILE:Witness_complex_example_strong_witness_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-l" "20" "-a" "0.5")
endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
