project(Collapse_utilities)

if (NOT EIGEN3_VERSION VERSION_LESS 3.1.0)
  if (TARGET Boost::program_options)
    # From a point cloud
    add_executable ( point_cloud_edge_collapse_rips_persistence point_cloud_edge_collapse_rips_persistence.cpp )
    target_link_libraries(point_cloud_edge_collapse_rips_persistence Boost::program_options)
  
    if (TBB_FOUND)
      target_link_libraries(point_cloud_edge_collapse_rips_persistence ${TBB_LIBRARIES})
    endif()
    add_test(NAME Edge_collapse_utilities_point_cloud_rips_persistence COMMAND $<TARGET_FILE:point_cloud_edge_collapse_rips_persistence>
             "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-r" "0.25" "-m" "0.5" "-d" "3" "-p" "3" "-o" "off_results.pers")
  
    install(TARGETS point_cloud_edge_collapse_rips_persistence DESTINATION bin)
  
    # From a distance matrix
    add_executable ( distance_matrix_edge_collapse_rips_persistence distance_matrix_edge_collapse_rips_persistence.cpp )
    target_link_libraries(distance_matrix_edge_collapse_rips_persistence Boost::program_options)
  
    if (TBB_FOUND)
      target_link_libraries(distance_matrix_edge_collapse_rips_persistence ${TBB_LIBRARIES})
    endif()
    add_test(NAME Edge_collapse_utilities_distance_matrix_rips_persistence COMMAND $<TARGET_FILE:distance_matrix_edge_collapse_rips_persistence>
             "${CMAKE_SOURCE_DIR}/data/distance_matrix/tore3D_1307_distance_matrix.csv" "-r" "0.25" "-m" "0.5" "-d" "3" "-p" "3" "-o" "csv_results.pers")
  
    install(TARGETS distance_matrix_edge_collapse_rips_persistence DESTINATION bin)
  
    if (DIFF_PATH)
      add_test(Edge_collapse_utilities_diff_persistence ${DIFF_PATH}
               "off_results.pers" "csv_results.pers")
      set_tests_properties(Edge_collapse_utilities_diff_persistence PROPERTIES DEPENDS
                           "Edge_collapse_utilities_point_cloud_rips_persistence;Edge_collapse_utilities_distance_matrix_rips_persistence")
    endif()
  endif()

endif()