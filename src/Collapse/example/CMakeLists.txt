project(Edge_collapse_examples)

# Point cloud
add_executable ( Edge_collapse_example_basic edge_collapse_basic_example.cpp )

if (TBB_FOUND)
  target_link_libraries(Edge_collapse_example_basic ${TBB_LIBRARIES})
endif()

add_test(NAME Edge_collapse_example_basic COMMAND $<TARGET_FILE:Edge_collapse_example_basic>)

# Point cloud
add_executable ( Edge_collapse_conserve_persistence edge_collapse_conserve_persistence.cpp )

if (TBB_FOUND)
  target_link_libraries(Edge_collapse_conserve_persistence ${TBB_LIBRARIES})
endif()

add_test(NAME Edge_collapse_conserve_persistence_1 COMMAND $<TARGET_FILE:Edge_collapse_conserve_persistence>
         "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "1.")

         add_test(NAME Edge_collapse_conserve_persistence_2 COMMAND $<TARGET_FILE:Edge_collapse_conserve_persistence>
         "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "2.")
