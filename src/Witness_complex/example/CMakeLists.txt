cmake_minimum_required(VERSION 2.6)
project(GUDHIWitnessComplex)

# A simple example
   add_executable ( simple_witness_complex simple_witness_complex.cpp )
   add_test(simple_witness_complex ${CMAKE_CURRENT_BINARY_DIR}/simple_witness_complex)

   add_executable( witness_complex_from_file witness_complex_from_file.cpp )
   #target_link_libraries(witness_complex_from_file ${EIGEN3_LIBRARIES} ${CGAL_LIBRARY})
   add_test( witness_complex_from_bunny &{CMAKE_CURRENT_BINARY_DIR}/witness_complex_from_file ${CMAKE_SOURCE_DIR}/data/points/bunny_5000 100)   
   
