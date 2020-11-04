project(Spatial_searching_examples)

if(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
  add_executable( Spatial_searching_example_spatial_searching example_spatial_searching.cpp )
  target_link_libraries(Spatial_searching_example_spatial_searching ${CGAL_LIBRARY})
  add_test(NAME Spatial_searching_example_spatial_searching
      COMMAND $<TARGET_FILE:Spatial_searching_example_spatial_searching>)
endif(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
