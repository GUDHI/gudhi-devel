cmake_minimum_required(VERSION 2.6)
project(Persistent_cohomology_utilities)

add_executable(rips_distance_matrix_persistence rips_distance_matrix_persistence.cpp)
target_link_libraries(rips_distance_matrix_persistence ${Boost_PROGRAM_OPTIONS_LIBRARY})

add_executable(rips_persistence rips_persistence.cpp)
target_link_libraries(rips_persistence ${Boost_PROGRAM_OPTIONS_LIBRARY})

if (TBB_FOUND)
  target_link_libraries(rips_distance_matrix_persistence ${TBB_LIBRARIES})
  target_link_libraries(rips_persistence ${TBB_LIBRARIES})
endif()

add_test(NAME Persistent_cohomology_example_from_rips_distance_matrix COMMAND $<TARGET_FILE:rips_distance_matrix_persistence>
    "${CMAKE_SOURCE_DIR}/data/distance_matrix/full_square_distance_matrix.csv" "-r" "1.0" "-d" "3" "-p" "3" "-m" "0")
add_test(NAME Persistent_cohomology_example_from_rips_on_tore_3D COMMAND $<TARGET_FILE:rips_persistence>
    "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-r" "0.25" "-m" "0.5" "-d" "3" "-p" "3")

install(TARGETS rips_distance_matrix_persistence DESTINATION bin)
install(TARGETS rips_persistence DESTINATION bin)

if(CGAL_FOUND)
  add_executable(alpha_complex_3d_persistence alpha_complex_3d_persistence.cpp)
  target_link_libraries(alpha_complex_3d_persistence ${CGAL_LIBRARY})

  if (TBB_FOUND)
    target_link_libraries(alpha_complex_3d_persistence ${TBB_LIBRARIES})
  endif(TBB_FOUND)
  add_test(NAME Persistent_cohomology_example_alpha_complex_3d COMMAND $<TARGET_FILE:alpha_complex_3d_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "2" "0.45")

  install(TARGETS alpha_complex_3d_persistence DESTINATION bin)

  if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.7.0)
    add_executable (alpha_complex_persistence alpha_complex_persistence.cpp)
    target_link_libraries(alpha_complex_persistence
        ${CGAL_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})

    add_executable(periodic_alpha_complex_3d_persistence periodic_alpha_complex_3d_persistence.cpp)
    target_link_libraries(periodic_alpha_complex_3d_persistence ${CGAL_LIBRARY})

    if (TBB_FOUND)
      target_link_libraries(alpha_complex_persistence ${TBB_LIBRARIES})
      target_link_libraries(periodic_alpha_complex_3d_persistence ${TBB_LIBRARIES})
    endif(TBB_FOUND)
    add_test(NAME Persistent_cohomology_example_alpha_complex COMMAND $<TARGET_FILE:alpha_complex_persistence>
        "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-p" "2" "-m" "0.45")
    add_test(NAME Persistent_cohomology_example_periodic_alpha_complex_3d COMMAND $<TARGET_FILE:periodic_alpha_complex_3d_persistence>
        "${CMAKE_SOURCE_DIR}/data/points/grid_10_10_10_in_0_1.off" "${CMAKE_SOURCE_DIR}/data/points/iso_cuboid_3_in_0_1.txt" "2" "0")

    install(TARGETS alpha_complex_persistence DESTINATION bin)
    install(TARGETS periodic_alpha_complex_3d_persistence DESTINATION bin)

  endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.7.0)
endif(CGAL_FOUND)
