project(Collapse_utilities)

# From a point cloud
add_executable ( point_cloud_edge_collapse_rips_persistence point_cloud_edge_collapse_rips_persistence.cpp )
target_link_libraries(point_cloud_edge_collapse_rips_persistence Boost::program_options)

if (TBB_FOUND)
  target_link_libraries(point_cloud_edge_collapse_rips_persistence ${TBB_LIBRARIES})
endif()
add_test(NAME point_cloud_edge_collapse_rips_persistence COMMAND $<TARGET_FILE:point_cloud_edge_collapse_rips_persistence>
         "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-r" "0.25" "-m" "0.5" "-d" "3" "-p" "3")

install(TARGETS point_cloud_edge_collapse_rips_persistence DESTINATION bin)

# From a distance matrix
add_executable ( distance_matrix_edge_collapse_rips_persistence distance_matrix_edge_collapse_rips_persistence.cpp )
target_link_libraries(distance_matrix_edge_collapse_rips_persistence Boost::program_options)

if (TBB_FOUND)
  target_link_libraries(distance_matrix_edge_collapse_rips_persistence ${TBB_LIBRARIES})
endif()
add_test(NAME distance_matrix_edge_collapse_rips_persistence COMMAND $<TARGET_FILE:distance_matrix_edge_collapse_rips_persistence>
         "-i" "${CMAKE_SOURCE_DIR}/data/distance_matrix/tore3D_300_distance_matrix.csv" "-e" "1e40" "-o" "csv_result.pers")

install(TARGETS distance_matrix_edge_collapse_rips_persistence DESTINATION bin)
