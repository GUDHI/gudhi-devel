cmake_minimum_required(VERSION 2.6)
project(off_file_from_shape_generator)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
  add_executable ( off_file_from_shape_generator off_file_from_shape_generator.cpp )
  target_link_libraries(off_file_from_shape_generator ${Boost_SYSTEM_LIBRARY})
  add_test(NAME off_file_from_shape_generator_on_sphere_1000_3_15.2 COMMAND $<TARGET_FILE:off_file_from_shape_generator>
      "on" "sphere" "onSphere.off" "1000" "3" "15.2")
  add_test(NAME off_file_from_shape_generator_in_sphere_100_2 COMMAND $<TARGET_FILE:off_file_from_shape_generator>
      "in" "sphere" "inSphere.off" "100" "2")

  # on cube is not available in CGAL
  add_test(NAME off_file_from_shape_generator_in_cube_10000_3_5.8 COMMAND $<TARGET_FILE:off_file_from_shape_generator>
      "in" "cube" "inCube.off" "10000" "3" "5.8")

  install(TARGETS off_file_from_shape_generator DESTINATION bin)

endif(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
