add_executable_with_targets(bottleneck_basic_example bottleneck_basic_example.cpp CGAL::CGAL TBB::tbb)

if (TARGET CGAL::CGAL)  
  add_test(NAME Bottleneck_distance_example_basic COMMAND $<TARGET_FILE:bottleneck_basic_example>)
endif ()

add_executable_with_targets(alpha_rips_persistence_bottleneck_distance alpha_rips_persistence_bottleneck_distance.cpp CGAL::CGAL Eigen3::Eigen Boost::program_options TBB::tbb)

if (TARGET CGAL::CGAL AND TARGET  Eigen3::Eigen AND TARGET Boost::program_options)
  add_test(NAME Bottleneck_distance_example_alpha_rips_persistence_bottleneck
      COMMAND $<TARGET_FILE:alpha_rips_persistence_bottleneck_distance>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-r" "0.15" "-m" "0.12" "-d" "3" "-p" "3")
endif ()
