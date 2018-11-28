project(Strong_collapse_examples)

add_executable (proximity_graph proximity_graph.cpp)
target_link_libraries(proximity_graph ${Boost_PROGRAM_OPTIONS_LIBRARY})

if (TBB_FOUND)
  target_link_libraries(proximity_graph ${TBB_LIBRARIES})
endif()

add_test(NAME Strong_collapse_example_proximity_graph_on_tore_3D COMMAND $<TARGET_FILE:proximity_graph>
        "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-r" "12.")
