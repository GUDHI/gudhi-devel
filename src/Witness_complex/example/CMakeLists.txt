cmake_minimum_required(VERSION 2.6)
project(GUDHIWitnessComplex)
#cmake -DCGAL_DIR=~/GitDrive/CGAL/ ../../..

#if (CGAL_FOUND)
   #message(STATUS "CGAL version: ${CGAL_VERSION}.")
   #include( ${CGAL_USE_FILE} )

   #find_package(Eigen3 3.1.0)
   #include( ${EIGEN3_USE_FILE} )

   #INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
   #INCLUDE_DIRECTORIES(${CGAL_INCLUDE_DIRS})
   add_executable ( simple_witness_complex simple_witness_complex.cpp )
   add_test(simple_witness_complex ${CMAKE_CURRENT_BINARY_DIR}/simple_witness_complex)

   #add_executable( witness_complex_from_file witness_complex_from_file.cpp )
   #target_link_libraries(witness_complex_from_file ${EIGEN3_LIBRARIES} ${CGAL_LIBRARY})
   #add_test( witness_complex_from_bunny &{CMAKE_CURRENT_BINARY_DIR}/witness_complex_from_file ${CMAKE_SOURCE_DIR}/data/points/bunny_5000 100)   
   
   add_executable( witness_complex_from_off witness_complex_from_off.cpp )
   
   add_executable( witness_complex_from_wl_matrix witness_complex_from_wl_matrix.cpp )
#endif()

# An example with Simplex-tree using CGAL alpha_shapes_3

include( ${CGAL_USE_FILE} )
find_package(Eigen3 3.1.0)
if(GMP_FOUND AND CGAL_FOUND)
   message("CGAL_lib = ${CGAL_LIBRARIES_DIR}")
   message("GMP_LIBRARIES = ${GMP_LIBRARIES}")
   message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
   #message("EIGEN3_LIBRARIES = ${EIGEN3_LIBRARIES}")
   INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIRS})
   INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR})
   INCLUDE_DIRECTORIES(${CGAL_INCLUDE_DIRS})
   add_executable (witness_complex_from_file witness_complex_from_file.cpp )
   target_link_libraries(witness_complex_from_file ${GMP_LIBRARIES} ${EIGEN3_LIBRARIES} ${CGAL_LIBRARY})
   add_test(witness_complex_from_file ${CMAKE_CURRENT_BINARY_DIR}/witness_complex_from_file ${CMAKE_SOURCE_DIR}/data/points/bunny_5000 100)
endif()
