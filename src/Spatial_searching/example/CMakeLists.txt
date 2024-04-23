add_executable_with_targets(Spatial_searching_example_spatial_searching example_spatial_searching.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)
if(TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  add_test(NAME Spatial_searching_example_spatial_searching
      COMMAND $<TARGET_FILE:Spatial_searching_example_spatial_searching>)
endif()
