cmake_minimum_required(VERSION 2.6)
project(Witness_complex_examples)


if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.6.0)
    if (EIGEN3_FOUND)
      add_executable( Witness_complex_example_off example_witness_complex_off.cpp )
      add_executable ( Witness_complex_example_sphere example_witness_complex_sphere.cpp )

      add_executable ( Witness_complex_example_witness_persistence example_witness_complex_persistence.cpp )
      target_link_libraries(Witness_complex_example_witness_persistence ${Boost_PROGRAM_OPTIONS_LIBRARY})

      add_executable ( Witness_complex_example_strong_witness_persistence example_strong_witness_persistence.cpp )
      target_link_libraries(Witness_complex_example_strong_witness_persistence ${Boost_PROGRAM_OPTIONS_LIBRARY})
      if (TBB_FOUND)
        target_link_libraries(Witness_complex_example_witness_persistence ${TBB_LIBRARIES})
        target_link_libraries(Witness_complex_example_strong_witness_persistence ${TBB_LIBRARIES})
      endif()

      add_test(Witness_complex_test_torus_persistence
               ${CMAKE_CURRENT_BINARY_DIR}/Witness_complex_example_witness_persistence
               ${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off -l 20 -a 0.5)
      add_test(Witness_complex_test_torus_persistence
               ${CMAKE_CURRENT_BINARY_DIR}/Witness_complex_example_strong_witness_persistence
               ${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off -l 20 -a 0.5)
    endif(EIGEN3_FOUND)
  endif (NOT CGAL_VERSION VERSION_LESS 4.6.0)
endif()
