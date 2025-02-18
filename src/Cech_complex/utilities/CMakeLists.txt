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

  # With square filtration values
  add_test(NAME Delaunay_Cech_complex_utility_on_tore_3D_square_safe COMMAND $<TARGET_FILE:delaunay_cech_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-r" "0.25" "-m" "0.5" "-p" "3" "-o" "delaunay_square_safe.pers")
  add_test(NAME Delaunay_Cech_complex_utility_on_tore_3D_square_fast COMMAND $<TARGET_FILE:delaunay_cech_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-r" "0.25" "-m" "0.5" "-p" "3" "-o" "delaunay_square_fast.pers" "-f")
  add_test(NAME Delaunay_Cech_complex_utility_on_tore_3D_square_exact COMMAND $<TARGET_FILE:delaunay_cech_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-r" "0.25" "-m" "0.5" "-p" "3" "-o" "delaunay_square_exact.pers" "-e")
  # Same without square filtration values - '-s off'
  add_test(NAME Delaunay_Cech_complex_utility_on_tore_3D_safe COMMAND $<TARGET_FILE:delaunay_cech_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-s" "off" "-r" "0.25" "-m" "0.5" "-p" "3" "-o" "delaunay_safe.pers")
  add_test(NAME Delaunay_Cech_complex_utility_on_tore_3D_fast COMMAND $<TARGET_FILE:delaunay_cech_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-s" "off" "-r" "0.25" "-m" "0.5" "-p" "3" "-o" "delaunay_fast.pers" "-f")
  add_test(NAME Delaunay_Cech_complex_utility_on_tore_3D_exact COMMAND $<TARGET_FILE:delaunay_cech_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-s" "off" "-r" "0.25" "-m" "0.5" "-p" "3" "-o" "delaunay_exact.pers" "-e")

  if (DIFF_PATH)
    # With square filtration values
    add_test(Delaunay_Cech_complex_utilities_diff_square_exact ${DIFF_PATH}
             "delaunay_square_exact.pers" "delaunay_square_safe.pers")
    set_tests_properties(Delaunay_Cech_complex_utilities_diff_square_exact PROPERTIES DEPENDS
             "Delaunay_Cech_complex_utility_on_tore_3D_square_safe;Delaunay_Cech_complex_utility_on_tore_3D_square_exact")

    add_test(Delaunay_Cech_complex_utilities_diff_square_fast ${DIFF_PATH}
             "delaunay_square_fast.pers" "delaunay_square_safe.pers")
    set_tests_properties(Delaunay_Cech_complex_utilities_diff_square_fast PROPERTIES DEPENDS
             "Delaunay_Cech_complex_utility_on_tore_3D_square_safe;Delaunay_Cech_complex_utility_on_tore_3D_square_fast")

    # Same without square filtration values - '-s off'
    add_test(Delaunay_Cech_complex_utilities_diff_exact ${DIFF_PATH}
             "delaunay_exact.pers" "delaunay_safe.pers")
    set_tests_properties(Delaunay_Cech_complex_utilities_diff_exact PROPERTIES DEPENDS
             "Delaunay_Cech_complex_utility_on_tore_3D_safe;Delaunay_Cech_complex_utility_on_tore_3D_exact")

    add_test(Delaunay_Cech_complex_utilities_diff_fast ${DIFF_PATH}
             "delaunay_fast.pers" "delaunay_safe.pers")
    set_tests_properties(Delaunay_Cech_complex_utilities_diff_fast PROPERTIES DEPENDS
             "Delaunay_Cech_complex_utility_on_tore_3D_safe;Delaunay_Cech_complex_utility_on_tore_3D_fast")
  endif()

  install(TARGETS delaunay_cech_persistence DESTINATION bin)

endif()
