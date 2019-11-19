project(Alpha_complex_utilities)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
  add_executable (alpha_complex_persistence alpha_complex_persistence.cpp)
  target_link_libraries(alpha_complex_persistence ${CGAL_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})

  if (TBB_FOUND)
    target_link_libraries(alpha_complex_persistence ${TBB_LIBRARIES})
  endif(TBB_FOUND)
  add_test(NAME Alpha_complex_utilities_safe_alpha_complex_persistence COMMAND $<TARGET_FILE:alpha_complex_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-p" "2" "-m" "0.45" "-o" "safe.pers")
  add_test(NAME Alpha_complex_utilities_fast_alpha_complex_persistence COMMAND $<TARGET_FILE:alpha_complex_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-p" "2" "-m" "0.45" "-o" "fast.pers" "-f")
  add_test(NAME Alpha_complex_utilities_exact_alpha_complex_persistence COMMAND $<TARGET_FILE:alpha_complex_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-p" "2" "-m" "0.45" "-o" "exact.pers" "-e")
  if (DIFF_PATH)
    add_test(Alpha_complex_utilities_diff_exact_alpha_complex ${DIFF_PATH}
             "exact.pers" "safe.pers")
    add_test(Alpha_complex_utilities_diff_fast_alpha_complex ${DIFF_PATH}
             "fast.pers" "safe.pers")
  endif()
  
  install(TARGETS alpha_complex_persistence DESTINATION bin)

  add_executable(alpha_complex_3d_persistence alpha_complex_3d_persistence.cpp)
  target_link_libraries(alpha_complex_3d_persistence ${CGAL_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})
  if (TBB_FOUND)
    target_link_libraries(alpha_complex_3d_persistence ${TBB_LIBRARIES})
  endif(TBB_FOUND)

  add_test(NAME Alpha_complex_utilities_alpha_complex_3d COMMAND $<TARGET_FILE:alpha_complex_3d_persistence>
          "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off"
          "-p" "2" "-m" "0.45" "-o" "safe_3d.pers")

  add_test(NAME Alpha_complex_utilities_exact_alpha_complex_3d COMMAND $<TARGET_FILE:alpha_complex_3d_persistence>
          "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off"
          "-p" "2" "-m" "0.45" "-o" "exact_3d.pers" "-e")

  add_test(NAME Alpha_complex_utilities_safe_alpha_complex_3d COMMAND $<TARGET_FILE:alpha_complex_3d_persistence>
          "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off"
          "-p" "2" "-m" "0.45" "-o" "fast_3d.pers" "-f")

  if (DIFF_PATH)
    add_test(Alpha_complex_utilities_diff_exact_alpha_complex_3d ${DIFF_PATH}
            "exact_3d.pers" "safe_3d.pers")
    add_test(Alpha_complex_utilities_diff_fast_alpha_complex_3d ${DIFF_PATH}
            "fast_3d.pers" "safe_3d.pers")
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

endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
