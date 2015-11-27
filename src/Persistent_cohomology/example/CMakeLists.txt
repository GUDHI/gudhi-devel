cmake_minimum_required(VERSION 2.6)
project(GUDHIExPersCohom)

# problem with Visual Studio link on Boost program_options
add_definitions( -DBOOST_ALL_NO_LIB )
add_definitions( -DBOOST_ALL_DYN_LINK )

add_executable(plain_homology plain_homology.cpp)
target_link_libraries(plain_homology ${Boost_SYSTEM_LIBRARY})
add_test(plain_homology ${CMAKE_CURRENT_BINARY_DIR}/plain_homology)

add_executable(persistence_from_simple_simplex_tree persistence_from_simple_simplex_tree.cpp)
target_link_libraries(persistence_from_simple_simplex_tree ${Boost_SYSTEM_LIBRARY})
add_test(persistence_from_simple_simplex_tree ${CMAKE_CURRENT_BINARY_DIR}/persistence_from_simple_simplex_tree 1 0)

add_executable(rips_persistence rips_persistence.cpp)
target_link_libraries(rips_persistence ${Boost_SYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})
   
add_test(rips_persistence_3 ${CMAKE_CURRENT_BINARY_DIR}/rips_persistence ${CMAKE_SOURCE_DIR}/data/points/Kl.txt -r 0.25 -d 3 -p 3 -m 100)
   
add_executable(persistence_from_file persistence_from_file.cpp)
target_link_libraries(persistence_from_file ${Boost_SYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY})
add_test(persistence_from_file_3_2_0 ${CMAKE_CURRENT_BINARY_DIR}/persistence_from_file ${CMAKE_SOURCE_DIR}/data/points/bunny_5000.st -p 2 -m 0)
add_test(persistence_from_file_3_3_100 ${CMAKE_CURRENT_BINARY_DIR}/persistence_from_file ${CMAKE_SOURCE_DIR}/data/points/bunny_5000.st -p 3 -m 100)
   
if(GMPXX_FOUND AND GMP_FOUND)
   message("GMPXX_LIBRARIES = ${GMPXX_LIBRARIES}")
   message("GMP_LIBRARIES = ${GMP_LIBRARIES}")

   add_executable(rips_multifield_persistence rips_multifield_persistence.cpp )
   target_link_libraries(rips_multifield_persistence ${Boost_SYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
   add_test(rips_multifield_persistence_2_71 ${CMAKE_CURRENT_BINARY_DIR}/rips_multifield_persistence ${CMAKE_SOURCE_DIR}/data/points/Kl.txt -r 0.25 -d 3 -p 2 -q 71 -m 100)

   add_executable ( performance_rips_persistence performance_rips_persistence.cpp )
   target_link_libraries(performance_rips_persistence ${Boost_SYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})

   if(CGAL_FOUND)
     add_executable(alpha_shapes_persistence alpha_shapes_persistence.cpp)
     target_link_libraries(alpha_shapes_persistence ${Boost_SYSTEM_LIBRARY} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} ${CGAL_LIBRARY})
     add_test(alpha_shapes_persistence_2_0_5 ${CMAKE_CURRENT_BINARY_DIR}/alpha_shapes_persistence ${CMAKE_SOURCE_DIR}/data/points/bunny_5000 2 0.5)
     #add_test(alpha_shapes_persistence_3_3_100 ${CMAKE_CURRENT_BINARY_DIR}/alpha_shapes_persistence ${CMAKE_SOURCE_DIR}/data/points/bunny_5000.st -p 3 -m 100)




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
      
      add_executable (alphacomplexpersistence alpha_complex_persistence.cpp)
      target_link_libraries(alphacomplexpersistence  ${Boost_SYSTEM_LIBRARY} ${CGAL_LIBRARY})

    else()
      message(WARNING "Eigen3 not found. Version 3.1.0 is required for Alpha shapes feature.")
    endif()
  else()
    message(WARNING "CGAL version: ${CGAL_VERSION} is too old to compile Alpha shapes feature. Version 4.6.0 is required.")
  endif ()



   endif()
   
endif()
