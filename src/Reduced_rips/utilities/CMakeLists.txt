add_executable_with_targets(reduced_rips_persistence reduced_rips_persistence.cpp CGAL::CGAL Eigen3::Eigen Boost::program_options TBB::tbb)

if(TARGET CGAL::CGAL AND TARGET Eigen3::Eigen AND TARGET Boost::program_options)
  add_test(NAME Reduced_rips_utility_from_off COMMAND $<TARGET_FILE:reduced_rips_persistence>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" "-k" "0")
  install(TARGETS reduced_rips_persistence DESTINATION bin)
endif()
