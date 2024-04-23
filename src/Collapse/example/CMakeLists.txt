add_executable_with_targets(Edge_collapse_example_basic edge_collapse_basic_example.cpp Eigen3::Eigen TBB::tbb)
add_executable_with_targets(Edge_collapse_conserve_persistence edge_collapse_conserve_persistence.cpp Eigen3::Eigen TBB::tbb)

if (TARGET Eigen3::Eigen)
  add_test(NAME Edge_collapse_example_basic COMMAND $<TARGET_FILE:Edge_collapse_example_basic>)
  
  add_test(NAME Edge_collapse_conserve_persistence_1 COMMAND $<TARGET_FILE:Edge_collapse_conserve_persistence>
           "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "0.2")
  
  add_test(NAME Edge_collapse_conserve_persistence_2 COMMAND $<TARGET_FILE:Edge_collapse_conserve_persistence>
           "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "1.8")
endif()