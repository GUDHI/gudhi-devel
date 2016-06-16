cmake_minimum_required(VERSION 2.6)
project(GUDHIExPersCohom)

# problem with Visual Studio link on Boost program_options
add_definitions( -DBOOST_ALL_NO_LIB )
add_definitions( -DBOOST_ALL_DYN_LINK )

add_executable(plain_homology plain_homology.cpp)
target_link_libraries(plain_homology ${Boost_SYSTEM_LIBRARY})

add_executable(persistence_from_simple_simplex_tree persistence_from_simple_simplex_tree.cpp)
target_link_libraries(persistence_from_simple_simplex_tree ${Boost_SYSTEM_LIBRARY})

add_executable(rips_persistence rips_persistence.cpp)
target_link_libraries(rips_persistence ${Boost_SYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})

add_executable(rips_persistence_via_boundary_matrix rips_persistence_via_boundary_matrix.cpp)
target_link_libraries(rips_persistence_via_boundary_matrix ${Boost_SYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})

add_executable(persistence_from_file persistence_from_file.cpp)
target_link_libraries(persistence_from_file ${Boost_SYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})

if (TBB_FOUND)
  target_link_libraries(plain_homology ${TBB_LIBRARIES})
  target_link_libraries(persistence_from_simple_simplex_tree ${TBB_LIBRARIES})
  target_link_libraries(rips_persistence ${TBB_LIBRARIES})
  target_link_libraries(rips_persistence_via_boundary_matrix ${TBB_LIBRARIES})
  target_link_libraries(persistence_from_file ${TBB_LIBRARIES})
endif()

add_test(plain_homology ${CMAKE_CURRENT_BINARY_DIR}/plain_homology)
add_test(persistence_from_simple_simplex_tree ${CMAKE_CURRENT_BINARY_DIR}/persistence_from_simple_simplex_tree 1 0)
add_test(rips_persistence_3 ${CMAKE_CURRENT_BINARY_DIR}/rips_persistence ${CMAKE_SOURCE_DIR}/data/points/Kl.txt -r 0.2 -d 3 -p 3 -m 100)
add_test(rips_persistence_via_boundary_matrix_3 ${CMAKE_CURRENT_BINARY_DIR}/rips_persistence_via_boundary_matrix ${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.txt -r 0.3 -d 3 -p 3 -m 100)
add_test(persistence_from_file_3_2_0 ${CMAKE_CURRENT_BINARY_DIR}/persistence_from_file ${CMAKE_SOURCE_DIR}/data/points/bunny_5000.st -p 2 -m 0)
add_test(persistence_from_file_3_3_100 ${CMAKE_CURRENT_BINARY_DIR}/persistence_from_file ${CMAKE_SOURCE_DIR}/data/points/bunny_5000.st -p 3 -m 100)
   
if(GMP_FOUND)
   message("GMP_LIBRARIES = ${GMP_LIBRARIES}")

   if(GMPXX_FOUND)
      message("GMPXX_LIBRARIES = ${GMPXX_LIBRARIES}")

      add_executable(rips_multifield_persistence rips_multifield_persistence.cpp )
      target_link_libraries(rips_multifield_persistence ${Boost_SYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
      add_executable ( performance_rips_persistence performance_rips_persistence.cpp )
      target_link_libraries(performance_rips_persistence ${Boost_SYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
      if (TBB_FOUND)
        target_link_libraries(rips_multifield_persistence ${TBB_LIBRARIES})
        target_link_libraries(performance_rips_persistence ${TBB_LIBRARIES})
      endif(TBB_FOUND)

   add_test(rips_multifield_persistence_2_71 ${CMAKE_CURRENT_BINARY_DIR}/rips_multifield_persistence ${CMAKE_SOURCE_DIR}/data/points/Kl.txt -r 0.2 -d 3 -p 2 -q 71 -m 100)
   endif(GMPXX_FOUND)

  if(CGAL_FOUND)
    add_executable(alpha_complex_3d_persistence alpha_complex_3d_persistence.cpp)
    target_link_libraries(alpha_complex_3d_persistence ${Boost_SYSTEM_LIBRARY} ${GMP_LIBRARIES} ${CGAL_LIBRARY})

    if (TBB_FOUND)
      target_link_libraries(alpha_complex_3d_persistence ${TBB_LIBRARIES})
    endif(TBB_FOUND)
    add_test(alpha_complex_3d_persistence_2_0_5 ${CMAKE_CURRENT_BINARY_DIR}/alpha_complex_3d_persistence ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off 2 0.45)


    if (NOT CGAL_VERSION VERSION_LESS 4.7.0)
      message(STATUS "CGAL version: ${CGAL_VERSION}.")

      find_package(Eigen3 3.1.0)
      if (EIGEN3_FOUND)
        message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
        include( ${EIGEN3_USE_FILE} )

        add_executable (alpha_complex_persistence alpha_complex_persistence.cpp)
        target_link_libraries(alpha_complex_persistence  ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})

        add_executable(periodic_alpha_complex_3d_persistence periodic_alpha_complex_3d_persistence.cpp)
        target_link_libraries(periodic_alpha_complex_3d_persistence ${Boost_SYSTEM_LIBRARY} ${GMP_LIBRARIES} ${CGAL_LIBRARY})

        add_executable(custom_persistence_sort custom_persistence_sort.cpp)
        target_link_libraries(custom_persistence_sort ${Boost_SYSTEM_LIBRARY} ${GMP_LIBRARIES} ${CGAL_LIBRARY})

        if (TBB_FOUND)
          target_link_libraries(alpha_complex_persistence ${TBB_LIBRARIES})
          target_link_libraries(periodic_alpha_complex_3d_persistence ${TBB_LIBRARIES})
          target_link_libraries(custom_persistence_sort ${TBB_LIBRARIES})
        endif(TBB_FOUND)
        add_test(alpha_complex_persistence_2_0_45 ${CMAKE_CURRENT_BINARY_DIR}/alpha_complex_persistence ${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off -m 0.45 -p 2)
        add_test(periodic_alpha_complex_3d_persistence_2_0 ${CMAKE_CURRENT_BINARY_DIR}/periodic_alpha_complex_3d_persistence ${CMAKE_SOURCE_DIR}/data/points/grid_10_10_10_in_0_1.off ${CMAKE_SOURCE_DIR}/data/points/iso_cuboid_3_in_0_1.txt 2 0)
        add_test(custom_persistence_sort ${CMAKE_CURRENT_BINARY_DIR}/custom_persistence_sort)

      else()
        message(WARNING "Eigen3 not found. Version 3.1.0 is required for Alpha shapes feature.")
      endif(EIGEN3_FOUND)
    else()
      message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Alpha shapes feature. Version 4.6.0 is required.")
    endif ()
  else()
    # message(WARNING "CGAL not found.")
  endif(CGAL_FOUND)
   
else()
  # message(WARNING "GMP not found.")
endif(GMP_FOUND)
