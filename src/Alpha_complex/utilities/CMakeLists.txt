project(Alpha_complex_utilities)

if(CGAL_FOUND)
  add_executable(alpha_complex_3d_persistence alpha_complex_3d_persistence.cpp)
  target_link_libraries(alpha_complex_3d_persistence ${CGAL_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})
  add_executable(exact_alpha_complex_3d_persistence exact_alpha_complex_3d_persistence.cpp)
  target_link_libraries(exact_alpha_complex_3d_persistence ${CGAL_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})
  add_executable(weighted_alpha_complex_3d_persistence weighted_alpha_complex_3d_persistence.cpp)
  target_link_libraries(weighted_alpha_complex_3d_persistence ${CGAL_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})

  if (TBB_FOUND)
    target_link_libraries(alpha_complex_3d_persistence ${TBB_LIBRARIES})
    target_link_libraries(exact_alpha_complex_3d_persistence ${TBB_LIBRARIES})
    target_link_libraries(weighted_alpha_complex_3d_persistence ${TBB_LIBRARIES})
  endif(TBB_FOUND)
  
  add_test(NAME Alpha_complex_utilities_alpha_complex_3d_persistence COMMAND $<TARGET_FILE:alpha_complex_3d_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-p" "2" "-m" "0.45")
  add_test(NAME Alpha_complex_utilities_exact_alpha_complex_3d COMMAND $<TARGET_FILE:exact_alpha_complex_3d_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off"  "-p" "2" "-m" "0.45")
  add_test(NAME Alpha_complex_utilities_weighted_alpha_complex_3d COMMAND $<TARGET_FILE:weighted_alpha_complex_3d_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.weights"  "-p" "2" "-m" "0.45")

  install(TARGETS alpha_complex_3d_persistence DESTINATION bin)
  install(TARGETS exact_alpha_complex_3d_persistence DESTINATION bin)
  install(TARGETS weighted_alpha_complex_3d_persistence DESTINATION bin)

  if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.7.0)
    add_executable (alpha_complex_persistence alpha_complex_persistence.cpp)
    target_link_libraries(alpha_complex_persistence ${CGAL_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})

    add_executable(periodic_alpha_complex_3d_persistence periodic_alpha_complex_3d_persistence.cpp)
    target_link_libraries(periodic_alpha_complex_3d_persistence ${CGAL_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})

    if (TBB_FOUND)
      target_link_libraries(alpha_complex_persistence ${TBB_LIBRARIES})
      target_link_libraries(periodic_alpha_complex_3d_persistence ${TBB_LIBRARIES})
    endif(TBB_FOUND)
    add_test(NAME Alpha_complex_utilities_alpha_complex_persistence COMMAND $<TARGET_FILE:alpha_complex_persistence>
        "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-p" "2" "-m" "0.45")
    add_test(NAME Alpha_complex_utilities_periodic_alpha_complex_3d_persistence COMMAND $<TARGET_FILE:periodic_alpha_complex_3d_persistence>
        "${CMAKE_SOURCE_DIR}/data/points/grid_10_10_10_in_0_1.off" "${CMAKE_SOURCE_DIR}/data/points/iso_cuboid_3_in_0_1.txt" "-p" "2" "-m" "0")

    install(TARGETS alpha_complex_persistence DESTINATION bin)
    install(TARGETS periodic_alpha_complex_3d_persistence DESTINATION bin)

  endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.7.0)

  if (NOT CGAL_VERSION VERSION_LESS 4.11.0)
    add_executable(weighted_periodic_alpha_complex_3d_persistence weighted_periodic_alpha_complex_3d_persistence.cpp)
    target_link_libraries(weighted_periodic_alpha_complex_3d_persistence ${CGAL_LIBRARY})
    if (TBB_FOUND)
      target_link_libraries(weighted_periodic_alpha_complex_3d_persistence ${TBB_LIBRARIES})
    endif(TBB_FOUND)

    add_test(NAME Alpha_complex_utilities_weigted_periodic_alpha_complex_3d COMMAND $<TARGET_FILE:weighted_periodic_alpha_complex_3d_persistence>
        "${CMAKE_SOURCE_DIR}/data/points/grid_10_10_10_in_0_1.off" "${CMAKE_SOURCE_DIR}/data/points/grid_10_10_10_in_0_1.weights"
        "${CMAKE_SOURCE_DIR}/data/points/iso_cuboid_3_in_0_1.txt" "3" "1.0")

    install(TARGETS weighted_periodic_alpha_complex_3d_persistence DESTINATION bin)

  endif (NOT CGAL_VERSION VERSION_LESS 4.11.0)

endif(CGAL_FOUND)
