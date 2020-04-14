project(Edge_collapse_examples)

# Point cloud
add_executable ( Edge_collapse_example_basic edge_collapse_basic_example.cpp )

if (TBB_FOUND)
  target_link_libraries(Edge_collapse_example_basic ${TBB_LIBRARIES})
endif()

add_test(NAME Edge_collapse_example_basic COMMAND $<TARGET_FILE:Edge_collapse_example_basic>)
