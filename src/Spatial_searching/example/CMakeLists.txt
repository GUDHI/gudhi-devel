project(Spatial_searching_examples)

if(CGAL_FOUND AND EIGEN3_FOUND)
  add_executable( Spatial_searching_example_spatial_searching example_spatial_searching.cpp )
  target_link_libraries(Spatial_searching_example_spatial_searching CGAL::CGAL Eigen3::Eigen)
  add_test(NAME Spatial_searching_example_spatial_searching
      COMMAND $<TARGET_FILE:Spatial_searching_example_spatial_searching>)
endif(CGAL_FOUND AND EIGEN3_FOUND)
