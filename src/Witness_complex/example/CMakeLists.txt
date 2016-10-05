cmake_minimum_required(VERSION 2.6)
project(Witness_complex_examples)


if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.6.0)
    if (EIGEN3_FOUND)
      add_executable( Witness_complex_example_off example_witness_complex_off.cpp )
      add_executable ( Witness_complex_example_sphere example_witness_complex_sphere.cpp )
      #add_executable ( relaxed_witness_persistence relaxed_witness_persistence.cpp )
    endif(EIGEN3_FOUND)
  endif (NOT CGAL_VERSION VERSION_LESS 4.6.0)
endif()
