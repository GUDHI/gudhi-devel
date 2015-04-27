cmake_minimum_required(VERSION 2.6)
project(GUDHIExPersCohom)

# problem with Visual Studio link on Boost program_options
if (NOT MSVC)
   add_executable(persistence_from_simple_simplex_tree persistence_from_simple_simplex_tree.cpp)
   target_link_libraries(persistence_from_simple_simplex_tree ${Boost_SYSTEM_LIBRARY})

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
         # uncomment to display debug traces
         # add_definitions(-DDEBUG_TRACES)
         add_executable(alpha_shapes_persistence alpha_shapes_persistence.cpp)
         target_link_libraries(alpha_shapes_persistence ${Boost_SYSTEM_LIBRARY} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES} ${CGAL_LIBRARY})
         add_test(alpha_shapes_persistence_2_0_5 ${CMAKE_CURRENT_BINARY_DIR}/alpha_shapes_persistence ${CMAKE_SOURCE_DIR}/data/points/bunny_5000 2 0.5)
         #add_test(alpha_shapes_persistence_3_3_100 ${CMAKE_CURRENT_BINARY_DIR}/alpha_shapes_persistence ${CMAKE_SOURCE_DIR}/data/points/bunny_5000.st -p 3 -m 100)
      endif()
   
   endif()
endif()
