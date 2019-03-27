project(Strong_collapse_utilities)

add_executable(strong_collapse_distance_matrix strong_collapse_distance_matrix.cpp)
target_link_libraries(strong_collapse_distance_matrix ${Boost_PROGRAM_OPTIONS_LIBRARY})

if (TBB_FOUND)
  target_link_libraries(strong_collapse_distance_matrix ${TBB_LIBRARIES})
endif()

add_test(NAME Strong_collapse_utility_distance_matrix COMMAND $<TARGET_FILE:strong_collapse_distance_matrix>
    "-i" "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-s" "1.0")

install(TARGETS strong_collapse_distance_matrix DESTINATION bin)
