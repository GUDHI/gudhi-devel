project(Bottleneck_distance_examples)

if (NOT CGAL_VERSION VERSION_LESS 4.11.0)
  add_executable (bottleneck_basic_example bottleneck_basic_example.cpp)

  if(TARGET TBB::tbb)
    target_link_libraries(bottleneck_basic_example TBB::tbb)
  endif()
  
  add_test(NAME Bottleneck_distance_example_basic COMMAND $<TARGET_FILE:bottleneck_basic_example>)

endif (NOT CGAL_VERSION VERSION_LESS 4.11.0)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
  if (TARGET Boost::program_options)
    add_executable (alpha_rips_persistence_bottleneck_distance alpha_rips_persistence_bottleneck_distance.cpp)
    target_link_libraries(alpha_rips_persistence_bottleneck_distance Boost::program_options)

    if(TARGET TBB::tbb)
      target_link_libraries(alpha_rips_persistence_bottleneck_distance TBB::tbb)
    endif()

    add_test(NAME Bottleneck_distance_example_alpha_rips_persistence_bottleneck
        COMMAND $<TARGET_FILE:alpha_rips_persistence_bottleneck_distance>
        "${CMAKE_SOURCE_DIR}/data/points/tore3D_1307.off" "-r" "0.15" "-m" "0.12" "-d" "3" "-p" "3")
  endif()
endif (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
