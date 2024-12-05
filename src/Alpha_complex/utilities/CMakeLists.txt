add_executable_with_targets(alpha_complex_persistence alpha_complex_persistence.cpp CGAL::CGAL Eigen3::Eigen Boost::program_options TBB::tbb)

if (TARGET CGAL::CGAL AND TARGET Eigen3::Eigen AND TARGET Boost::program_options)
  add_test(NAME Alpha_complex_utilities_safe_alpha_complex_persistence COMMAND $<TARGET_FILE:alpha_complex_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-p" "2" "-m" "0.45" "-o" "safe.pers")
  add_test(NAME Alpha_complex_utilities_fast_alpha_complex_persistence COMMAND $<TARGET_FILE:alpha_complex_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-p" "2" "-m" "0.45" "-o" "fast.pers" "-f")
  add_test(NAME Alpha_complex_utilities_exact_alpha_complex_persistence COMMAND $<TARGET_FILE:alpha_complex_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-p" "2" "-m" "0.45" "-o" "exact.pers" "-e")
  # Same with sqrt filtration values - "-s"
  add_test(NAME Alpha_complex_utilities_safe_sqrt_alpha_complex_persistence COMMAND $<TARGET_FILE:alpha_complex_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-s" "-p" "2" "-m" "0.45" "-o" "safe_sqrt.pers")
  add_test(NAME Alpha_complex_utilities_fast_sqrt_alpha_complex_persistence COMMAND $<TARGET_FILE:alpha_complex_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-s" "-p" "2" "-m" "0.45" "-o" "fast_sqrt.pers" "-f")
  add_test(NAME Alpha_complex_utilities_exact_sqrt_alpha_complex_persistence COMMAND $<TARGET_FILE:alpha_complex_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-s" "-p" "2" "-m" "0.45" "-o" "exact_sqrt.pers" "-e")
  if (DIFF_PATH)
    add_test(Alpha_complex_utilities_diff_exact_alpha_complex ${DIFF_PATH}
             "exact.pers" "safe.pers")
    set_tests_properties(Alpha_complex_utilities_diff_exact_alpha_complex PROPERTIES DEPENDS
             "Alpha_complex_utilities_exact_alpha_complex_persistence;Alpha_complex_utilities_safe_alpha_complex_persistence")

    add_test(Alpha_complex_utilities_diff_fast_alpha_complex ${DIFF_PATH}
             "fast.pers" "safe.pers")
    set_tests_properties(Alpha_complex_utilities_diff_fast_alpha_complex PROPERTIES DEPENDS
             "Alpha_complex_utilities_fast_alpha_complex_persistence;Alpha_complex_utilities_safe_alpha_complex_persistence")

    # Same with sqrt filtration values - "-s"
    add_test(Alpha_complex_utilities_diff_exact_sqrt_alpha_complex ${DIFF_PATH}
             "exact_sqrt.pers" "safe_sqrt.pers")
    set_tests_properties(Alpha_complex_utilities_diff_exact_sqrt_alpha_complex PROPERTIES DEPENDS
             "Alpha_complex_utilities_exact_sqrt_alpha_complex_persistence;Alpha_complex_utilities_safe_sqrt_alpha_complex_persistence")

    add_test(Alpha_complex_utilities_diff_fast_sqrt_alpha_complex ${DIFF_PATH}
             "fast_sqrt.pers" "safe_sqrt.pers")
    set_tests_properties(Alpha_complex_utilities_diff_fast_sqrt_alpha_complex PROPERTIES DEPENDS
             "Alpha_complex_utilities_fast_sqrt_alpha_complex_persistence;Alpha_complex_utilities_safe_sqrt_alpha_complex_persistence")
  endif()

  install(TARGETS alpha_complex_persistence DESTINATION bin)
endif()

add_executable_with_targets(alpha_complex_3d_persistence alpha_complex_3d_persistence.cpp CGAL::CGAL Boost::program_options TBB::tbb)
if (TARGET CGAL::CGAL AND TARGET Boost::program_options)
  add_test(NAME Alpha_complex_utilities_alpha_complex_3d COMMAND $<TARGET_FILE:alpha_complex_3d_persistence>
          "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off"
          "-p" "2" "-m" "0.45" "-o" "safe_3d.pers")

  add_test(NAME Alpha_complex_utilities_exact_alpha_complex_3d COMMAND $<TARGET_FILE:alpha_complex_3d_persistence>
          "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off"
          "-p" "2" "-m" "0.45" "-o" "exact_3d.pers" "-e")

  add_test(NAME Alpha_complex_utilities_fast_alpha_complex_3d COMMAND $<TARGET_FILE:alpha_complex_3d_persistence>
          "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off"
          "-p" "2" "-m" "0.45" "-o" "fast_3d.pers" "-f")

  if (DIFF_PATH)
    add_test(Alpha_complex_utilities_diff_exact_alpha_complex_3d ${DIFF_PATH}
            "exact_3d.pers" "safe_3d.pers")
    set_tests_properties(Alpha_complex_utilities_diff_exact_alpha_complex_3d PROPERTIES DEPENDS
            "Alpha_complex_utilities_exact_alpha_complex_3d;Alpha_complex_utilities_alpha_complex_3d")
    add_test(Alpha_complex_utilities_diff_fast_alpha_complex_3d ${DIFF_PATH}
            "fast_3d.pers" "safe_3d.pers")
    set_tests_properties(Alpha_complex_utilities_diff_fast_alpha_complex_3d PROPERTIES DEPENDS
            "Alpha_complex_utilities_fast_alpha_complex_3d;Alpha_complex_utilities_alpha_complex_3d")
  endif()

  add_test(NAME Alpha_complex_utilities_periodic_alpha_complex_3d_persistence COMMAND $<TARGET_FILE:alpha_complex_3d_persistence>
          "${CMAKE_SOURCE_DIR}/data/points/grid_10_10_10_in_0_1.off"
          "-c" "${CMAKE_SOURCE_DIR}/data/points/iso_cuboid_3_in_0_1.txt"
          "-p" "2" "-m" "0")

  add_test(NAME Alpha_complex_utilities_weighted_alpha_complex_3d COMMAND $<TARGET_FILE:alpha_complex_3d_persistence>
          "${CMAKE_SOURCE_DIR}/data/points/grid_10_10_10_in_0_1.off"
          "-w" "${CMAKE_SOURCE_DIR}/data/points/grid_10_10_10_in_0_1.weights"
          "-p" "2" "-m" "0")

  add_test(NAME Alpha_complex_utilities_weighted_periodic_alpha_complex_3d COMMAND $<TARGET_FILE:alpha_complex_3d_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/grid_10_10_10_in_0_1.off"
      "-w" "${CMAKE_SOURCE_DIR}/data/points/grid_10_10_10_in_0_1.weights"
      "-c" "${CMAKE_SOURCE_DIR}/data/points/iso_cuboid_3_in_0_1.txt"
      "-p" "2" "-m" "0" "-e")

  install(TARGETS alpha_complex_3d_persistence DESTINATION bin)
endif()
