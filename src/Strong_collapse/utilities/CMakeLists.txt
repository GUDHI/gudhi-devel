project(Strong_collapse_utilities)

add_executable(strong_collapse_distance_matrix strong_collapse_distance_matrix.cpp)
target_link_libraries(strong_collapse_distance_matrix ${Boost_PROGRAM_OPTIONS_LIBRARY})

if (TBB_FOUND)
  target_link_libraries(strong_collapse_distance_matrix ${TBB_LIBRARIES})
endif()

add_test(NAME Strong_collapse_utility_distance_matrix_dry_run COMMAND $<TARGET_FILE:strong_collapse_distance_matrix>
    "-i" "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "--dry-run")
add_test(NAME Strong_collapse_utility_distance_matrix_torus COMMAND $<TARGET_FILE:strong_collapse_distance_matrix>
    "-i" "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-s" "0.0467216" "1.0467216" "2.0467216" "3.0467216" "4.0467216" "5.0467216" "6.0467216" "-o" "md.csv")

install(TARGETS strong_collapse_distance_matrix DESTINATION bin)

add_executable(strong_collapse_rips_persistence strong_collapse_rips_persistence.cpp)
target_link_libraries(strong_collapse_rips_persistence ${Boost_PROGRAM_OPTIONS_LIBRARY})

add_test(NAME Strong_collapse_utility_rips_persistence_dry_run COMMAND $<TARGET_FILE:strong_collapse_rips_persistence>
        "-i" "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "--dry-run")
add_test(NAME Strong_collapse_utility_rips_persistence_torus COMMAND $<TARGET_FILE:strong_collapse_rips_persistence>
        "-i" "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-s" "0.0467216" "1.0467216" "2.0467216" "3.0467216" "4.0467216" "5.0467216" "6.0467216" "-o" "pers.txt")

if (TBB_FOUND)
  target_link_libraries(strong_collapse_rips_persistence ${TBB_LIBRARIES})
endif()

install(TARGETS strong_collapse_rips_persistence DESTINATION bin)
