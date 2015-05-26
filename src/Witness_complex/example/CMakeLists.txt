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

   add_executable( witness_complex_from_file witness_complex_from_file.cpp )
   #target_link_libraries(witness_complex_from_file ${EIGEN3_LIBRARIES} ${CGAL_LIBRARY})
   add_test( witness_complex_from_bunny &{CMAKE_CURRENT_BINARY_DIR}/witness_complex_from_file ${CMAKE_SOURCE_DIR}/data/points/bunny_5000 100)   
   
   add_executable( witness_complex_from_off witness_complex_from_off.cpp )
   
   add_executable( witness_complex_from_wl_matrix witness_complex_from_wl_matrix.cpp )
#endif()

# An example with Simplex-tree using CGAL alpha_shapes_3

#find_package(Eigen3 3.1.0)
#if(GMP_FOUND AND CGAL_FOUND)
#   message("CGAL_lib = ${CGAL_LIBRARIES_DIR}")
#   message("GMP_LIBRARIES = ${GMP_LIBRARIES}")
#   message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
#   #message("EIGEN3_LIBRARIES = ${EIGEN3_LIBRARIES}")
#   INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIRS})
#   INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR})
#   INCLUDE_DIRECTORIES(${CGAL_INCLUDE_DIRS})
#   add_executable (witness_complex_knn_landmarks witness_complex_knn_landmarks.cpp )
#   target_link_libraries(witness_complex_knn_landmarks ${EIGEN3_LIBRARIES} ${CGAL_LIBRARY})
#   add_test(witness_complex_knn_landmarks ${CMAKE_CURRENT_BINARY_DIR}/witness_complex_knn_landmarks ${CMAKE_SOURCE_DIR}/data/points/bunny_5000 100)
#endif()

# need CGAL 4.6
# cmake -DCGAL_DIR=~/workspace/CGAL-4.6-beta1 ../../..
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.6.0)
    message(STATUS "CGAL version: ${CGAL_VERSION}.")

    include( ${CGAL_USE_FILE} )

    find_package(Eigen3 3.1.0)
    if (EIGEN3_FOUND)
      message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
      include( ${EIGEN3_USE_FILE} )
      include_directories (BEFORE "../../include")

      add_executable ( witness_complex_knn_landmarks witness_complex_knn_landmarks.cpp )
      target_link_libraries(witness_complex_knn_landmarks ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
      add_test(witness_complex_knn_landmarks ${CMAKE_CURRENT_BINARY_DIR}/witness_complex_knn_landmarks ${CMAKE_SOURCE_DIR}/data/points/bunny_5000 100)
      #add_executable ( witness_complex_perturbations witness_complex_perturbations.cpp )
      #target_link_libraries(witness_complex_perturbations ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
      #add_test(witness_complex_perturbations ${CMAKE_CURRENT_BINARY_DIR}/witness_complex_perturbations ${CMAKE_SOURCE_DIR}/data/points/bunny_5000 100)
    else()
      message(WARNING "Eigen3 not found. Version 3.1.0 is required for Alpha shapes feature.")      
    endif()
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Alpha shapes feature. Version 4.6.0 is required.")
  endif ()
endif()
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.6.0)
    message(STATUS "CGAL version: ${CGAL_VERSION}.")

    include( ${CGAL_USE_FILE} )

    find_package(Eigen3 3.1.0)
    if (EIGEN3_FOUND)
      message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
      include( ${EIGEN3_USE_FILE} )
      include_directories (BEFORE "../../include")
      add_executable ( witness_complex_perturbations witness_complex_perturbations.cpp )
      target_link_libraries(witness_complex_perturbations ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
      add_test(witness_complex_perturbations ${CMAKE_CURRENT_BINARY_DIR}/witness_complex_perturbations ${CMAKE_SOURCE_DIR}/data/points/bunny_5000 100)
      add_executable ( witness_complex_flat_torus witness_complex_flat_torus.cpp )
      target_link_libraries(witness_complex_flat_torus ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
      add_test(witness_complex_flat_torus ${CMAKE_CURRENT_BINARY_DIR}/witness_complex_flat_torus ${CMAKE_SOURCE_DIR}/data/points/bunny_5000 100)
    else()
      message(WARNING "Eigen3 not found. Version 3.1.0 is required for Alpha shapes feature.")      
    endif()
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Alpha shapes feature. Version 4.6.0 is required.")
  endif ()
endif()
