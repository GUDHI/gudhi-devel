cmake_minimum_required(VERSION 2.6)
project(data_points_generator)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
  add_executable ( data_points_generator hypergenerator.cpp )
  target_link_libraries(data_points_generator ${Boost_SYSTEM_LIBRARY})
  add_test(NAME data_points_generator_on_sphere_1000_3_15.2 COMMAND $<TARGET_FILE:data_points_generator>
      "on" "sphere" "onSphere.off" "1000" "3" "15.2")
  add_test(NAME data_points_generator_in_sphere_100_2 COMMAND $<TARGET_FILE:data_points_generator>
      "in" "sphere" "inSphere.off" "100" "2")

  # on cube is not available in CGAL
  add_test(NAME data_points_generator_in_cube_10000_3_5.8 COMMAND $<TARGET_FILE:data_points_generator>
      "in" "cube" "inCube.off" "10000" "3" "5.8")
endif(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
