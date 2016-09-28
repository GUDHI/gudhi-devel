cmake_minimum_required(VERSION 2.6)
project(Witness_complex_examples)

# A simple example
   add_executable( witness_complex_from_file witness_complex_from_file.cpp )
   add_test( witness_complex_from_bunny ${CMAKE_CURRENT_BINARY_DIR}/witness_complex_from_file ${CMAKE_SOURCE_DIR}/data/points/bunny_5000 100)   

if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.6.0)
    if (EIGEN3_FOUND)
      add_executable ( witness_complex_sphere witness_complex_sphere.cpp )
      target_link_libraries(witness_complex_sphere ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
      add_test( witness_complex_sphere_10 ${CMAKE_CURRENT_BINARY_DIR}/witness_complex_sphere 10)
      add_executable ( relaxed_witness_persistence relaxed_witness_persistence.cpp )
    endif(EIGEN3_FOUND)
  endif (NOT CGAL_VERSION VERSION_LESS 4.6.0)
endif()
