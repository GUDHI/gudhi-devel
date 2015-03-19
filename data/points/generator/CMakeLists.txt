cmake_minimum_required(VERSION 2.6)
project(GUDHIPointsGenerator)

if(CGAL_FOUND)
  include( ${CGAL_USE_FILE} )

  find_package(Eigen3 3.1.0)
  if (EIGEN3_FOUND)
    include( ${EIGEN3_USE_FILE} )
    include_directories (BEFORE "../../include")

    add_executable ( hypergenerator hypergenerator.cpp )
    add_test(hypergenerator_on_sphere_3000_10_5.0 ${CMAKE_CURRENT_BINARY_DIR}/hypergenerator on sphere onSphere.off 3000 10 5.0)
    add_test(hypergenerator_on_sphere_10000_3 ${CMAKE_CURRENT_BINARY_DIR}/hypergenerator on sphere onSphere.off 10000 3)
    add_test(hypergenerator_in_sphere_7000_12_10.8 ${CMAKE_CURRENT_BINARY_DIR}/hypergenerator in sphere inSphere.off 7000 12 10.8)
    add_test(hypergenerator_in_sphere_50000_2 ${CMAKE_CURRENT_BINARY_DIR}/hypergenerator in sphere inSphere.off 50000 2)
    # on cube is not available in CGAL
    add_test(hypergenerator_in_cube_7000_12_10.8 ${CMAKE_CURRENT_BINARY_DIR}/hypergenerator in cube inCube.off 7000 12 10.8)
    add_test(hypergenerator_in_cube_50000_2 ${CMAKE_CURRENT_BINARY_DIR}/hypergenerator in cube inCube.off 50000 3)
  endif()
endif()
