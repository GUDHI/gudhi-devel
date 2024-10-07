add_executable_with_targets(cech_persistence cech_persistence.cpp CGAL::CGAL Eigen3::Eigen Boost::program_options TBB::tbb)
add_executable_with_targets(delaunay_cech_persistence delaunay_cech_persistence.cpp CGAL::CGAL Eigen3::Eigen Boost::program_options TBB::tbb)

if (TARGET CGAL::CGAL AND TARGET Eigen3::Eigen AND TARGET Boost::program_options)
  add_test(NAME Cech_complex_utility_from_rips_on_tore_3D_safe COMMAND $<TARGET_FILE:cech_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-r" "0.25" "-m" "0.5" "-d" "3" "-p" "3" "-o" "safe.pers")
  add_test(NAME Cech_complex_utility_from_rips_on_tore_3D_fast COMMAND $<TARGET_FILE:cech_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-r" "0.25" "-m" "0.5" "-d" "3" "-p" "3" "-o" "fast.pers" "-f")
  add_test(NAME Cech_complex_utility_from_rips_on_tore_3D_exact COMMAND $<TARGET_FILE:cech_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-r" "0.25" "-m" "0.5" "-d" "3" "-p" "3" "-o" "exact.pers" "-e")

  if (DIFF_PATH)
    add_test(Cech_complex_utilities_diff_exact ${DIFF_PATH}
             "exact.pers" "safe.pers")
    set_tests_properties(Cech_complex_utilities_diff_exact PROPERTIES DEPENDS
             "Cech_complex_utility_from_rips_on_tore_3D_safe;Cech_complex_utility_from_rips_on_tore_3D_exact")

    add_test(Cech_complex_utilities_diff_fast ${DIFF_PATH}
             "fast.pers" "safe.pers")
    set_tests_properties(Cech_complex_utilities_diff_fast PROPERTIES DEPENDS
             "Cech_complex_utility_from_rips_on_tore_3D_safe;Cech_complex_utility_from_rips_on_tore_3D_fast")
  endif()

  install(TARGETS cech_persistence DESTINATION bin)

  add_test(NAME Delaunay_Cech_complex_utility_from_rips_on_tore_3D_safe COMMAND $<TARGET_FILE:delaunay_cech_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-r" "0.25" "-m" "0.5" "-p" "3" "-o" "delaunay_safe.pers")
  add_test(NAME Delaunay_Cech_complex_utility_from_rips_on_tore_3D_fast COMMAND $<TARGET_FILE:delaunay_cech_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-r" "0.25" "-m" "0.5" "-p" "3" "-o" "delaunay_fast.pers" "-f")
  add_test(NAME Delaunay_Cech_complex_utility_from_rips_on_tore_3D_exact COMMAND $<TARGET_FILE:delaunay_cech_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-r" "0.25" "-m" "0.5" "-p" "3" "-o" "delaunay_exact.pers" "-e")
  # Same with sqrt filtration values - "-s"
  add_test(NAME Delaunay_Cech_complex_utility_from_rips_on_tore_3D_sqrt_safe COMMAND $<TARGET_FILE:delaunay_cech_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-s" "-r" "0.25" "-m" "0.5" "-p" "3" "-o" "delaunay_sqrt_safe.pers")
  add_test(NAME Delaunay_Cech_complex_utility_from_rips_on_tore_3D_sqrt_fast COMMAND $<TARGET_FILE:delaunay_cech_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-s" "-r" "0.25" "-m" "0.5" "-p" "3" "-o" "delaunay_sqrt_fast.pers" "-f")
  add_test(NAME Delaunay_Cech_complex_utility_from_rips_on_tore_3D_sqrt_exact COMMAND $<TARGET_FILE:delaunay_cech_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-s" "-r" "0.25" "-m" "0.5" "-p" "3" "-o" "delaunay_sqrt_exact.pers" "-e")

  if (DIFF_PATH)
    add_test(Delaunay_Cech_complex_utilities_diff_exact ${DIFF_PATH}
             "delaunay_exact.pers" "delaunay_safe.pers")
    set_tests_properties(Delaunay_Cech_complex_utilities_diff_exact PROPERTIES DEPENDS
             "Delaunay_Cech_complex_utility_from_rips_on_tore_3D_safe;Delaunay_Cech_complex_utility_from_rips_on_tore_3D_exact")

    add_test(Delaunay_Cech_complex_utilities_diff_fast ${DIFF_PATH}
             "delaunay_fast.pers" "delaunay_safe.pers")
    set_tests_properties(Delaunay_Cech_complex_utilities_diff_fast PROPERTIES DEPENDS
             "Delaunay_Cech_complex_utility_from_rips_on_tore_3D_safe;Delaunay_Cech_complex_utility_from_rips_on_tore_3D_fast")

    # Same with sqrt filtration values - "-s"
    add_test(Delaunay_Cech_complex_utilities_diff_sqrt_exact ${DIFF_PATH}
             "delaunay_sqrt_exact.pers" "delaunay_sqrt_safe.pers")
    set_tests_properties(Delaunay_Cech_complex_utilities_diff_sqrt_exact PROPERTIES DEPENDS
             "Delaunay_Cech_complex_utility_from_rips_on_tore_3D_sqrt_safe;Delaunay_Cech_complex_utility_from_rips_on_tore_3D_sqrt_exact")

    add_test(Delaunay_Cech_complex_utilities_diff_sqrt_fast ${DIFF_PATH}
             "delaunay_sqrt_fast.pers" "delaunay_sqrt_safe.pers")
    set_tests_properties(Delaunay_Cech_complex_utilities_diff_sqrt_fast PROPERTIES DEPENDS
             "Delaunay_Cech_complex_utility_from_rips_on_tore_3D_sqrt_safe;Delaunay_Cech_complex_utility_from_rips_on_tore_3D_sqrt_fast")
  endif()

  install(TARGETS delaunay_cech_persistence DESTINATION bin)

endif()
