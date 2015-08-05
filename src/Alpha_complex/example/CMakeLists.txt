cmake_minimum_required(VERSION 2.6)
project(GUDHIAlphaShapesExample)

# need CGAL 4.7
# cmake -DCGAL_DIR=~/workspace/CGAL-4.7-Ic-41 ../../..
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.7.0)
    message(STATUS "CGAL version: ${CGAL_VERSION}.")

    include( ${CGAL_USE_FILE} )
 
    find_package(Eigen3 3.1.0)
    if (EIGEN3_FOUND)
      message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
      include( ${EIGEN3_USE_FILE} )
      if (CMAKE_BUILD_TYPE MATCHES Debug)
        # For programs to be more verbose
        add_definitions(-DDEBUG_TRACES)
      endif()
      
      add_executable ( alphaoffreader Alpha_complex_from_off.cpp )
      target_link_libraries(alphaoffreader  ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})

      add_executable ( alphapoints Alpha_complex_from_points.cpp )
      target_link_libraries(alphapoints  ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
    else()
      message(WARNING "Eigen3 not found. Version 3.1.0 is required for Alpha shapes feature.")
    endif()
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Alpha shapes feature. Version 4.6.0 is required.")
  endif ()
endif()
