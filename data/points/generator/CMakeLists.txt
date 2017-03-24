cmake_minimum_required(VERSION 2.6)
project(Data_points_generator)

if (NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
  add_executable ( hypergenerator hypergenerator.cpp )
  target_link_libraries(hypergenerator ${Boost_SYSTEM_LIBRARY})
  add_test(hypergenerator_on_sphere_3000_10_5.0 ${CMAKE_CURRENT_BINARY_DIR}/hypergenerator on sphere onSphere.off 3000 10 5.0)
  add_test(hypergenerator_on_sphere_10000_3 ${CMAKE_CURRENT_BINARY_DIR}/hypergenerator on sphere onSphere.off 10000 3)
  add_test(hypergenerator_in_sphere_7000_12_10.8 ${CMAKE_CURRENT_BINARY_DIR}/hypergenerator in sphere inSphere.off 7000 12 10.8)
  add_test(hypergenerator_in_sphere_50000_2 ${CMAKE_CURRENT_BINARY_DIR}/hypergenerator in sphere inSphere.off 50000 2)
  # on cube is not available in CGAL
  add_test(hypergenerator_in_cube_7000_12_10.8 ${CMAKE_CURRENT_BINARY_DIR}/hypergenerator in cube inCube.off 7000 12 10.8)
  add_test(hypergenerator_in_cube_50000_2 ${CMAKE_CURRENT_BINARY_DIR}/hypergenerator in cube inCube.off 50000 3)
endif(NOT CGAL_WITH_EIGEN3_VERSION VERSION_LESS 4.6.0)
