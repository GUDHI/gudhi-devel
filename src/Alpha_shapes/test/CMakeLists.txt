cmake_minimum_required(VERSION 2.6)
project(GUDHIAlphaShapesTest)

# need CGAL 4.6
# cmake -DCGAL_DIR=~/workspace/CGAL-4.6-beta1 ../../..
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.6.0)
    message(STATUS "CGAL version: ${CGAL_VERSION}.")

    include( ${CGAL_USE_FILE} )
 
    find_package(Eigen3 3.1.0)
    if (EIGEN3_FOUND)
      message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
      include( ${EIGEN3_USE_FILE} )
      include_directories (BEFORE "../../include")

      add_definitions(-DDEBUG_TRACES)
      add_executable ( AlphaShapesUnitTest Alpha_shapes_unit_test.cpp )
      target_link_libraries(AlphaShapesUnitTest  ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
      add_test(AlphaShapesUnitTest ${CMAKE_CURRENT_BINARY_DIR}/AlphaShapesUnitTest)

    else()
      message(WARNING "Eigen3 not found. Version 3.1.0 is required for Alpha shapes feature.")
    endif()
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Alpha shapes feature. Version 4.6.0 is required.")
  endif ()
endif()

cpplint_add_tests("${CMAKE_SOURCE_DIR}/src/Alpha_shapes/include/gudhi")
