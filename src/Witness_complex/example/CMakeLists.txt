cmake_minimum_required(VERSION 2.6)
project(GUDHIWitnessComplex)

# A simple example
   add_executable( witness_complex_from_file witness_complex_from_file.cpp )
   add_test( witness_complex_from_bunny &{CMAKE_CURRENT_BINARY_DIR}/witness_complex_from_file ${CMAKE_SOURCE_DIR}/data/points/bunny_5000 100)   
   
   add_executable( witness_complex_bench witness_complex_bench.cpp )
   add_test( witness_complex_bench_from_bunny &{CMAKE_CURRENT_BINARY_DIR}/witness_complex_bench ${CMAKE_SOURCE_DIR}/data/points/bunny_5000 100)   

if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.6.0)
    message(STATUS "CGAL version: ${CGAL_VERSION}.")

    include( ${CGAL_USE_FILE} )

    find_package(Eigen3 3.1.0)
    if (EIGEN3_FOUND)
      message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
      include( ${EIGEN3_USE_FILE} )
      message(STATUS "Eigen3 use file: ${EIGEN3_USE_FILE}.")
      include_directories (BEFORE "../../include")

      add_executable ( witness_complex_bench2 witness_complex_bench2.cpp )
      target_link_libraries(witness_complex_bench2 ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
    else()
      message(WARNING "Eigen3 not found. Version 3.1.0 is required for Alpha shapes feature.")      
    endif()
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Alpha shapes feature. Version 4.6.0 is required.")
  endif ()
endif()
