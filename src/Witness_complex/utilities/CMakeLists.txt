add_executable_with_targets(Witness_complex_strong_witness_persistence strong_witness_persistence.cpp CGAL::CGAL Eigen3::Eigen Boost::program_options TBB::tbb)
add_executable_with_targets(Witness_complex_weak_witness_persistence weak_witness_persistence.cpp CGAL::CGAL Eigen3::Eigen Boost::program_options TBB::tbb)

# CGAL and Eigen3 are required for Euclidean version of Witness
if(TARGET CGAL::CGAL AND TARGET Eigen3::Eigen AND TARGET Boost::program_options)
  add_test(NAME Witness_complex_strong_test_torus_persistence
      COMMAND $<TARGET_FILE:Witness_complex_strong_witness_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-l" "20" "-a" "0.5")
  add_test(NAME Witness_complex_weak_test_torus_persistence
      COMMAND $<TARGET_FILE:Witness_complex_weak_witness_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-l" "20" "-a" "0.5")

  install(TARGETS Witness_complex_strong_witness_persistence DESTINATION bin)
  install(TARGETS Witness_complex_weak_witness_persistence DESTINATION bin)
endif()
