project(Computational_geometry_examples)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.11.0)
  add_executable ( points_genericity points_genericity.cpp )
  
  add_test(NAME Computational_geometry_example_points_genericity COMMAND $<TARGET_FILE:points_genericity>
      "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off")
endif()