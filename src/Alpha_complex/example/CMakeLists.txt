cmake_minimum_required(VERSION 2.6)
project(GUDHIAlphaShapesExample)

add_executable ( flat flat.cpp )

# need CGAL 4.7
# cmake -DCGAL_DIR=~/workspace/CGAL-4.7-Ic-41 ../../..
if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.7.0)
    message(STATUS "CGAL version: ${CGAL_VERSION}.")

    include( ${CGAL_USE_FILE} )
    # In CMakeLists.txt, when include(${CGAL_USE_FILE}), CXX_FLAGS are overwritten.
    # cf. http://doc.cgal.org/latest/Manual/installation.html#title40
    # A workaround is to add "-std=c++11" again.
    # A fix would be to use https://cmake.org/cmake/help/v3.1/prop_gbl/CMAKE_CXX_KNOWN_FEATURES.html 
    # or even better https://cmake.org/cmake/help/v3.1/variable/CMAKE_CXX_STANDARD.html
    # but it implies to use cmake version 3.1 at least.
    if(NOT MSVC)
      include(CheckCXXCompilerFlag)
      CHECK_CXX_COMPILER_FLAG(-std=c++11 COMPILER_SUPPORTS_CXX11)
      if(COMPILER_SUPPORTS_CXX11)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
      endif()
    endif()
    # - End of workaround
 
    find_package(Eigen3 3.1.0)
    if (EIGEN3_FOUND)
      message(STATUS "Eigen3 version: ${EIGEN3_VERSION}.")
      include( ${EIGEN3_USE_FILE} )
      
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
