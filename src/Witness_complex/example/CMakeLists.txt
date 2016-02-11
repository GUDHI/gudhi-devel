cmake_minimum_required(VERSION 2.6)
project(GUDHIWitnessComplex)

# A simple example
   add_executable( witness_complex_from_file witness_complex_from_file.cpp )
   add_test( witness_complex_from_bunny ${CMAKE_CURRENT_BINARY_DIR}/witness_complex_from_file ${CMAKE_SOURCE_DIR}/data/points/bunny_5000 100)   

if(CGAL_FOUND)
  if (NOT CGAL_VERSION VERSION_LESS 4.6.0)
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
      message(STATUS "Eigen3 use file: ${EIGEN3_USE_FILE}.")
      include_directories (BEFORE "../../include")

      add_executable ( witness_complex_sphere witness_complex_sphere.cpp )
      target_link_libraries(witness_complex_sphere ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
      add_test( witness_complex_sphere_10 ${CMAKE_CURRENT_BINARY_DIR}/witness_complex_sphere 10)   
      add_executable ( relaxed_witness_persistence relaxed_witness_persistence.cpp )
      target_link_libraries(relaxed_witness_persistence ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})
      add_test( relaxed_witness_persistence_10 ${CMAKE_CURRENT_BINARY_DIR}/relaxed_witness_persistence 10)   
    else()
      message(WARNING "Eigen3 not found. Version 3.1.0 is required for witness_complex_sphere example.")      
    endif()
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile witness_complex_sphere example. Version 4.6.0 is required.")
  endif ()
endif()
