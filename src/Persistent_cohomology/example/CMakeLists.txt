cmake_minimum_required(VERSION 2.6)
project(Persistent_cohomology_examples)

add_executable(plain_homology plain_homology.cpp)

add_executable(persistence_from_simple_simplex_tree persistence_from_simple_simplex_tree.cpp)

add_executable(rips_distance_matrix_persistence rips_distance_matrix_persistence.cpp)
target_link_libraries(rips_distance_matrix_persistence ${Boost_PROGRAM_OPTIONS_LIBRARY})

add_executable(rips_persistence rips_persistence.cpp)
target_link_libraries(rips_persistence ${Boost_PROGRAM_OPTIONS_LIBRARY})

add_executable(rips_persistence_step_by_step rips_persistence_step_by_step.cpp)
target_link_libraries(rips_persistence_step_by_step ${Boost_PROGRAM_OPTIONS_LIBRARY})

add_executable(rips_persistence_via_boundary_matrix rips_persistence_via_boundary_matrix.cpp)
target_link_libraries(rips_persistence_via_boundary_matrix ${Boost_PROGRAM_OPTIONS_LIBRARY})

add_executable(persistence_from_file persistence_from_file.cpp)
target_link_libraries(persistence_from_file ${Boost_PROGRAM_OPTIONS_LIBRARY})

if (TBB_FOUND)
  target_link_libraries(plain_homology ${TBB_LIBRARIES})
  target_link_libraries(persistence_from_simple_simplex_tree ${TBB_LIBRARIES})
  target_link_libraries(rips_distance_matrix_persistence ${TBB_LIBRARIES})
  target_link_libraries(rips_persistence ${TBB_LIBRARIES})
  target_link_libraries(rips_persistence_step_by_step ${TBB_LIBRARIES})
  target_link_libraries(rips_persistence_via_boundary_matrix ${TBB_LIBRARIES})
  target_link_libraries(persistence_from_file ${TBB_LIBRARIES})
endif()

add_test(NAME Persistent_cohomology_example_plain_homology COMMAND $<TARGET_FILE:plain_homology>)
add_test(NAME Persistent_cohomology_example_from_simple_simplex_tree COMMAND $<TARGET_FILE:persistence_from_simple_simplex_tree>
    "1" "0")
add_test(NAME Persistent_cohomology_example_from_rips_distance_matrix COMMAND $<TARGET_FILE:rips_distance_matrix_persistence>
    "${CMAKE_SOURCE_DIR}/data/distance_matrix/full_square_distance_matrix.csv" "-r" "1.0" "-d" "3" "-p" "3" "-m" "0")
add_test(NAME Persistent_cohomology_example_from_rips_on_tore_3D COMMAND $<TARGET_FILE:rips_persistence>
    "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-r" "0.25" "-m" "0.5" "-d" "3" "-p" "3")
add_test(NAME Persistent_cohomology_example_from_rips_step_by_step_on_tore_3D COMMAND $<TARGET_FILE:rips_persistence_step_by_step>
    "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-r" "0.25" "-m" "0.5" "-d" "3" "-p" "3")
add_test(NAME Persistent_cohomology_example_via_boundary_matrix COMMAND $<TARGET_FILE:rips_persistence_via_boundary_matrix>
    "${CMAKE_SOURCE_DIR}/data/points/Kl.off" "-r" "0.16" "-d" "3" "-p" "3" "-m" "100")
add_test(NAME Persistent_cohomology_example_from_file_3_2_0 COMMAND $<TARGET_FILE:persistence_from_file>
    "${CMAKE_SOURCE_DIR}/data/filtered_simplicial_complex/bunny_5000_complex.fsc" "-p" "2" "-m" "0")
add_test(NAME Persistent_cohomology_example_from_file_3_3_100 COMMAND $<TARGET_FILE:persistence_from_file>
    "${CMAKE_SOURCE_DIR}/data/filtered_simplicial_complex/bunny_5000_complex.fsc" "-p" "3" "-m" "100")

install(TARGETS plain_homology DESTINATION bin)
install(TARGETS persistence_from_simple_simplex_tree DESTINATION bin)
install(TARGETS rips_distance_matrix_persistence DESTINATION bin)
install(TARGETS rips_persistence DESTINATION bin)
install(TARGETS rips_persistence_step_by_step DESTINATION bin)
install(TARGETS rips_persistence_via_boundary_matrix DESTINATION bin)
install(TARGETS persistence_from_file DESTINATION bin)

if(GMP_FOUND)
   if(GMPXX_FOUND)
      add_executable(rips_multifield_persistence rips_multifield_persistence.cpp )
      target_link_libraries(rips_multifield_persistence
          ${Boost_PROGRAM_OPTIONS_LIBRARY} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
      if (TBB_FOUND)
        target_link_libraries(rips_multifield_persistence ${TBB_LIBRARIES})
      endif(TBB_FOUND)
      add_test(NAME Persistent_cohomology_example_multifield_2_71 COMMAND $<TARGET_FILE:rips_multifield_persistence>
          "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-r" "0.25" "-m" "0.5" "-d" "3" "-p" "2" "-q" "71")
      install(TARGETS rips_multifield_persistence DESTINATION bin)
   endif(GMPXX_FOUND)
endif(GMP_FOUND)

if(CGAL_FOUND)
  add_executable(alpha_complex_3d_persistence alpha_complex_3d_persistence.cpp)
  target_link_libraries(alpha_complex_3d_persistence ${CGAL_LIBRARY})
  add_executable(exact_alpha_complex_3d_persistence exact_alpha_complex_3d_persistence.cpp)
  target_link_libraries(exact_alpha_complex_3d_persistence ${CGAL_LIBRARY})

  if (TBB_FOUND)
    target_link_libraries(alpha_complex_3d_persistence ${TBB_LIBRARIES})
    target_link_libraries(exact_alpha_complex_3d_persistence ${TBB_LIBRARIES})
  endif(TBB_FOUND)
  add_test(NAME Persistent_cohomology_example_alpha_complex_3d COMMAND $<TARGET_FILE:alpha_complex_3d_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "2" "0.45")
  add_test(NAME Persistent_cohomology_example_exact_alpha_complex_3d COMMAND $<TARGET_FILE:exact_alpha_complex_3d_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "2" "0.45")

  install(TARGETS alpha_complex_3d_persistence DESTINATION bin)
  install(TARGETS exact_alpha_complex_3d_persistence DESTINATION bin)

  if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.7.0)
    add_executable (alpha_complex_persistence alpha_complex_persistence.cpp)
    target_link_libraries(alpha_complex_persistence
        ${CGAL_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})

    add_executable(periodic_alpha_complex_3d_persistence periodic_alpha_complex_3d_persistence.cpp)
    target_link_libraries(periodic_alpha_complex_3d_persistence ${CGAL_LIBRARY})

    add_executable(custom_persistence_sort custom_persistence_sort.cpp)
    target_link_libraries(custom_persistence_sort ${CGAL_LIBRARY})

    if (TBB_FOUND)
      target_link_libraries(alpha_complex_persistence ${TBB_LIBRARIES})
      target_link_libraries(periodic_alpha_complex_3d_persistence ${TBB_LIBRARIES})
      target_link_libraries(custom_persistence_sort ${TBB_LIBRARIES})
    endif(TBB_FOUND)
    add_test(NAME Persistent_cohomology_example_alpha_complex COMMAND $<TARGET_FILE:alpha_complex_persistence>
        "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-p" "2" "-m" "0.45")
    add_test(NAME Persistent_cohomology_example_periodic_alpha_complex_3d COMMAND $<TARGET_FILE:periodic_alpha_complex_3d_persistence>
        "${CMAKE_SOURCE_DIR}/data/points/grid_10_10_10_in_0_1.off" "${CMAKE_SOURCE_DIR}/data/points/iso_cuboid_3_in_0_1.txt" "2" "0")
    add_test(NAME Persistent_cohomology_example_custom_persistence_sort COMMAND $<TARGET_FILE:custom_persistence_sort>)

    install(TARGETS alpha_complex_persistence DESTINATION bin)
    install(TARGETS periodic_alpha_complex_3d_persistence DESTINATION bin)
    install(TARGETS custom_persistence_sort DESTINATION bin)

  endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.7.0)

  if (NOT CGAL_VERSION VERSION_LESS 4.11.0)
    add_executable(weighted_periodic_alpha_complex_3d_persistence weighted_periodic_alpha_complex_3d_persistence.cpp)
    target_link_libraries(weighted_periodic_alpha_complex_3d_persistence ${CGAL_LIBRARY})
    if (TBB_FOUND)
      target_link_libraries(weighted_periodic_alpha_complex_3d_persistence ${TBB_LIBRARIES})
    endif(TBB_FOUND)

    install(TARGETS weighted_periodic_alpha_complex_3d_persistence DESTINATION bin)

    add_executable(test_cgal test_cgal.cpp)
    target_link_libraries(test_cgal ${CGAL_LIBRARY})
    if (TBB_FOUND)
      target_link_libraries(test_cgal ${TBB_LIBRARIES})
    endif(TBB_FOUND)


    # On CGAL 4.11, the weighted alpha shape 3 interface was changed - the new interface
    add_executable(weighted_alpha_complex_3d_persistence weighted_alpha_complex_3d_persistence.cpp)
  else(NOT CGAL_VERSION VERSION_LESS 4.11.0)
    # On CGAL 4.11, the weighted alpha shape 3 interface was changed - the old interface
    add_executable(weighted_alpha_complex_3d_persistence weighted_alpha_complex_3d_persistence_old_cgal_interface.cpp)
  endif (NOT CGAL_VERSION VERSION_LESS 4.11.0)

  # cf. if (NOT CGAL_VERSION VERSION_LESS 4.11.0) and else(NOT CGAL_VERSION VERSION_LESS 4.11.0)
  # same target name on purpose
  target_link_libraries(weighted_alpha_complex_3d_persistence ${CGAL_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(weighted_alpha_complex_3d_persistence ${TBB_LIBRARIES})
  endif(TBB_FOUND)

  install(TARGETS weighted_alpha_complex_3d_persistence DESTINATION bin)

  add_test(NAME Persistent_cohomology_example_weighted_alpha_complex_3d COMMAND $<TARGET_FILE:weighted_alpha_complex_3d_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.weights" "2" "0.45")

endif(CGAL_FOUND)
