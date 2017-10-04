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
