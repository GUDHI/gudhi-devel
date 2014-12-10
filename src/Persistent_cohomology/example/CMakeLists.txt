cmake_minimum_required(VERSION 2.6)
project(GUDHIExPersCohom)

# NEED TO FIND BOOST NEEDED COMPONENTS TO LINK WITH
find_package(Boost 1.45.0 COMPONENTS system program_options)
message("Boost_lib = ${Boost_LIBRARIES}")

add_executable(rips_persistence rips_persistence.cpp)
target_link_libraries(rips_persistence ${Boost_LIBRARIES})
add_test(rips_persistence_2 ${CMAKE_CURRENT_BINARY_DIR}/rips_persistence ${CMAKE_SOURCE_DIR}/data/points/Kl.txt -r 0.25 -d 3 -p 2 -m 100)
add_test(rips_persistence_3 ${CMAKE_CURRENT_BINARY_DIR}/rips_persistence ${CMAKE_SOURCE_DIR}/data/points/Kl.txt -r 0.25 -d 3 -p 3 -m 100)

add_executable(rips_persistence_from_alpha_shape_3 rips_persistence_from_alpha_shape_3.cpp )
target_link_libraries(rips_persistence_from_alpha_shape_3 ${Boost_LIBRARIES})

if(GMPXX_FOUND AND GMP_FOUND)
   message("GMPXX_LIBRARIES = ${GMPXX_LIBRARIES}")
   message("GMP_LIBRARIES = ${GMP_LIBRARIES}")

   add_executable(rips_multifield_persistence rips_multifield_persistence.cpp )
   target_link_libraries(rips_multifield_persistence  ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
   add_test(rips_multifield_persistence_2_3 ${CMAKE_CURRENT_BINARY_DIR}/rips_multifield_persistence ${CMAKE_SOURCE_DIR}/data/points/Kl.txt -r 0.25 -d 3 -p 2 -q 3 -m 100)
   add_test(rips_multifield_persistence_2_71 ${CMAKE_CURRENT_BINARY_DIR}/rips_multifield_persistence ${CMAKE_SOURCE_DIR}/data/points/Kl.txt -r 0.25 -d 3 -p 2 -q 71 -m 100)

   add_executable ( performance_rips_persistence performance_rips_persistence.cpp )
   target_link_libraries(performance_rips_persistence  ${Boost_LIBRARIES} ${GMPXX_LIBRARIES} ${GMP_LIBRARIES})
   add_test(rips_persistence_from_alpha_shape_3_2_0 ${CMAKE_CURRENT_BINARY_DIR}/rips_persistence_from_alpha_shape_3 ${CMAKE_SOURCE_DIR}/data/points/bunny_5000.st -p 2 -m 0)
   add_test(rips_persistence_from_alpha_shape_3_3_100 ${CMAKE_CURRENT_BINARY_DIR}/rips_persistence_from_alpha_shape_3 ${CMAKE_SOURCE_DIR}/data/points/bunny_5000.st -p 3 -m 100)
endif()
