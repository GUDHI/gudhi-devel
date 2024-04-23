add_executable_with_targets(plain_homology plain_homology.cpp TBB::tbb)
add_test(NAME Persistent_cohomology_example_plain_homology COMMAND $<TARGET_FILE:plain_homology>)

add_executable_with_targets(persistence_from_simple_simplex_tree persistence_from_simple_simplex_tree.cpp TBB::tbb)
add_test(NAME Persistent_cohomology_example_from_simple_simplex_tree COMMAND $<TARGET_FILE:persistence_from_simple_simplex_tree>
    "2" "0")

add_executable_with_targets(rips_persistence_step_by_step rips_persistence_step_by_step.cpp Boost::program_options TBB::tbb)
add_executable_with_targets(rips_persistence_via_boundary_matrix rips_persistence_via_boundary_matrix.cpp Boost::program_options TBB::tbb)
add_executable_with_targets(persistence_from_file persistence_from_file.cpp Boost::program_options TBB::tbb)
if(TARGET Boost::program_options)
  add_test(NAME Persistent_cohomology_example_from_rips_step_by_step_on_tore_3D COMMAND $<TARGET_FILE:rips_persistence_step_by_step>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-r" "0.25" "-m" "0.5" "-d" "3" "-p" "3")

  add_test(NAME Persistent_cohomology_example_via_boundary_matrix COMMAND $<TARGET_FILE:rips_persistence_via_boundary_matrix>
      "${CMAKE_SOURCE_DIR}/data/points/Kl.off" "-r" "0.16" "-d" "3" "-p" "3" "-m" "100")

  add_test(NAME Persistent_cohomology_example_from_file_3_2_0 COMMAND $<TARGET_FILE:persistence_from_file>
      "${CMAKE_SOURCE_DIR}/data/filtered_simplicial_complex/Klein_bottle_complex.fsc" "-p" "2" "-m" "0")
  add_test(NAME Persistent_cohomology_example_from_file_3_3_100 COMMAND $<TARGET_FILE:persistence_from_file>
      "${CMAKE_SOURCE_DIR}/data/filtered_simplicial_complex/Klein_bottle_complex.fsc" "-p" "3" "-m" "100")
endif()

if(GMP_FOUND)
  if(GMPXX_FOUND)
    if(TARGET Boost::program_options)
      add_executable(rips_multifield_persistence rips_multifield_persistence.cpp )
      target_link_libraries(rips_multifield_persistence
          Boost::program_options ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
      if(TARGET TBB::tbb)
        target_link_libraries(rips_multifield_persistence TBB::tbb)
      endif()
      add_test(NAME Persistent_cohomology_example_multifield_2_71 COMMAND $<TARGET_FILE:rips_multifield_persistence>
          "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-r" "0.25" "-m" "0.5" "-d" "3" "-p" "2" "-q" "71")
    endif()
  endif(GMPXX_FOUND)
endif(GMP_FOUND)

add_executable_with_targets(custom_persistence_sort custom_persistence_sort.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
if(TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  add_test(NAME Persistent_cohomology_example_custom_persistence_sort COMMAND $<TARGET_FILE:custom_persistence_sort>)
endif()

add_executable(persistence_on_a_line persistence_on_a_line.cpp)
add_test(NAME persistence_on_a_line COMMAND $<TARGET_FILE:persistence_on_a_line>)
