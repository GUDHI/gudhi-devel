if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 5.0.1)
  project(Cech_complex_utilities)

  if (TARGET Boost::program_options)
    add_executable(cech_persistence cech_persistence.cpp)
    target_link_libraries(cech_persistence Boost::program_options)

    if(TARGET TBB::tbb)
      target_link_libraries(cech_persistence TBB::tbb)
    endif()

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
  endif()
endif()
